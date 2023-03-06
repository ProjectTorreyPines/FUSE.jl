#= ========== =#
#  HFS sizing  #
#= ========== =#
Base.@kwdef mutable struct FUSEparameters__ActorHFSsizing{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    j_tolerance::Entry{T} = Entry(T, "-", "Tolerance on the conductor current limits"; default=0.4)
    stress_tolerance::Entry{T} = Entry(T, "-", "Tolerance on the structural stresses limits"; default=0.2)
    fixed_aspect_ratio::Entry{Bool} = Entry(Bool, "-", "Raise an error if aspect_ratio changes more than 10%"; default=true)
    unconstrained_flattop_duration::Entry{Bool} = Entry(Bool, "-", "Maximize flux_duration without targeting a specific value"; default=true)
    do_plot::Entry{Bool} = Entry(Bool, "-", "plot"; default=false)
    verbose::Entry{Bool} = Entry(Bool, "-", "verbose"; default=false)
end

mutable struct ActorHFSsizing <: ReactorAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorHFSsizing
    stresses_actor::ActorStresses
    fluxswing_actor::ActorFluxSwing
end

"""
    ActorHFSsizing(dd::IMAS.dd, act::ParametersAllActors; kw...)

Actor that resizes the High Field Side of the tokamak radial build
* takes into account the OH maximum allowed superconductor current/Field
* takes into account the stresses on the center stack
    
!!! note 
    Manipulates radial build information in `dd.build.layer`
"""
function ActorHFSsizing(dd::IMAS.dd, act::ParametersAllActors; kw_ActorFluxSwing=Dict(), kw_ActorStresses=Dict(), kw...)
    par = act.ActorHFSsizing(kw...)
    actor = ActorHFSsizing(dd, par, act; kw_ActorFluxSwing, kw_ActorStresses)
    if par.do_plot
        p = plot(dd.build)
    end
    step(actor)
    finalize(actor)
    if par.do_plot
        display(plot!(p, dd.build; cx=false))
    end
    return actor
end

function ActorHFSsizing(dd::IMAS.dd, par::FUSEparameters__ActorHFSsizing, act::ParametersAllActors; kw_ActorFluxSwing=Dict(), kw_ActorStresses=Dict(), kw...)
    par = act.ActorHFSsizing(kw...)
    fluxswing_actor = ActorFluxSwing(dd, act.ActorFluxSwing; kw_ActorFluxSwing...)
    stresses_actor = ActorStresses(dd, act.ActorStresses; kw_ActorStresses...)
    return ActorHFSsizing(dd, par, stresses_actor, fluxswing_actor)
end

function _step(
    actor::ActorHFSsizing;
    j_tolerance::Real=actor.par.j_tolerance,
    stress_tolerance::Real=actor.par.stress_tolerance,
    fixed_aspect_ratio::Bool=actor.par.fixed_aspect_ratio,
    unconstrained_flattop_duration::Bool=actor.par.unconstrained_flattop_duration,
    verbose::Bool=actor.par.verbose,
    debug_plot=false
)

    function target_value(value, target, tolerance) # relative error with tolerance
        return abs((value .* (1.0 .+ tolerance) .- target) ./ (abs(target) + 1.0))
    end

    function assign_PL_OH_TF_thicknesses(x0, what)
        x0 = map(abs, x0)
        c_extra = 0.0

        if what == :oh
            OH.thickness = x0[1]
            if length(x0) == 2
                fraction_stainless, c_extra = mirror_bound_w_cost(x0[2], 0.5, 1.0 - dd.build.oh.technology.fraction_void - 0.05)
                dd.build.oh.technology.fraction_stainless = fraction_stainless
            end

        elseif what == :tf
            TFhfs.thickness = x0[1]
            if length(x0) == 2
                fraction_stainless, c_extra = mirror_bound_w_cost(x0[2], 0.5, 1.0 - dd.build.oh.technology.fraction_void - 0.05)
                dd.build.tf.technology.fraction_stainless = fraction_stainless
            end

        else
            OH.thickness, TFhfs.thickness = x0
            if length(x0) == 4
                fraction_stainless, c_extra = mirror_bound_w_cost(x0[3], 0.5, 1.0 - dd.build.oh.technology.fraction_void - 0.05)
                dd.build.oh.technology.fraction_stainless = fraction_stainless
                fraction_stainless, c_extra = mirror_bound_w_cost(x0[4], 0.5, 1.0 - dd.build.oh.technology.fraction_void - 0.05)
                dd.build.tf.technology.fraction_stainless = fraction_stainless
            end
        end

        plug.thickness += old_plasma_start_radius - plasma.start_radius
        plug.thickness = max(OH.thickness / 4.0, plug.thickness)

        TFlfs.thickness = TFhfs.thickness
        return c_extra
    end

    function cost(x0, what)
        # assign optimization arguments and evaluate coils currents and stresses
        c_extra = assign_PL_OH_TF_thicknesses(x0, what)
        _step(actor.fluxswing_actor; operate_at_j_crit=unconstrained_flattop_duration, j_tolerance, only=what)
        _step(actor.stresses_actor)

        # OH and plug sizing based on stresses
        c_joh = c_soh = c_spl = 0.0
        if what ∈ [:oh, :all]
            c_joh1 = target_value(dd.build.oh.critical_j, dd.build.oh.max_j, -j_tolerance)
            c_joh2 = target_value(dd.build.oh.max_j, dd.build.oh.critical_j, j_tolerance)
            c_joh = norm([c_joh1, c_joh2])
            c_soh = target_value(maximum(dd.solid_mechanics.center_stack.stress.vonmises.oh), stainless_steel.yield_strength, stress_tolerance)
            if !ismissing(dd.solid_mechanics.center_stack.stress.vonmises, :pl)
                c_spl = target_value(maximum(dd.solid_mechanics.center_stack.stress.vonmises.pl), stainless_steel.yield_strength, stress_tolerance)
            end
        end

        # TF sizing based on stresses
        c_jtf = c_stf = 0.0
        if what ∈ [:tf, :all]
            c_jtf1 = target_value(dd.build.tf.critical_j, dd.build.tf.max_j, -j_tolerance)
            c_jtf2 = target_value(dd.build.tf.max_j, dd.build.tf.critical_j, j_tolerance)
            c_jtf = norm([c_jtf1, c_jtf2])
            c_stf = target_value(maximum(dd.solid_mechanics.center_stack.stress.vonmises.tf), stainless_steel.yield_strength, stress_tolerance)
        end

        if debug_plot
            push!(C_JOH, c_joh)
            push!(C_SOH, c_soh)
            push!(C_JTF, c_jtf)
            push!(C_STF, c_stf)
        end

        # total cost
        return norm(vcat([c_joh, c_jtf], [c_soh, c_stf, c_spl], [c_extra]))
    end

    @assert actor.stresses_actor.dd === actor.fluxswing_actor.dd
    dd = actor.stresses_actor.dd
    target_B0 = maximum(abs.(dd.equilibrium.vacuum_toroidal_field.b0))

    # init
    plug = dd.build.layer[1]
    OH = IMAS.get_build(dd.build, type=_oh_)
    TFhfs = IMAS.get_build(dd.build, type=_tf_, fs=_hfs_)
    TFlfs = IMAS.get_build(dd.build, type=_tf_, fs=_lfs_)
    iplasma = IMAS.get_build(dd.build, type=_plasma_, return_index=true)
    plasma = dd.build.layer[iplasma]

    old_R0 = (TFhfs.end_radius + TFlfs.start_radius) / 2.0
    old_plasma_start_radius = plasma.start_radius
    old_a = plasma.thickness / 2.0
    old_ϵ = old_R0 / old_a

    if debug_plot
        C_JOH = []
        C_SOH = []
        C_JTF = []
        C_STF = []
    end

    # initialize all dd fields
    step(actor.fluxswing_actor; operate_at_j_crit=unconstrained_flattop_duration, j_tolerance)
    step(actor.stresses_actor)

    dd.build.oh.technology.fraction_stainless = 0.5
    dd.build.tf.technology.fraction_stainless = 0.5

    # plug and OH optimization (w/fraction)
    old_thicknesses = [layer.thickness for layer in dd.build.layer]
    res = Optim.optimize(
        x0 -> cost(x0, :oh),
        [OH.thickness, dd.build.oh.technology.fraction_stainless],
        Optim.NelderMead(),
        Optim.Options(time_limit=60);
        autodiff=:forward
    )
    assign_PL_OH_TF_thicknesses(res.minimizer, :oh)
    _step(actor.fluxswing_actor; operate_at_j_crit=unconstrained_flattop_duration, j_tolerance, only=:oh)
    _step(actor.stresses_actor)
    if verbose
        display(res)
    end

    # TF optimization (w/fraction)
    old_thicknesses = [layer.thickness for layer in dd.build.layer]
    res = Optim.optimize(
        x0 -> cost(x0, :tf),
        [TFhfs.thickness, dd.build.tf.technology.fraction_stainless],
        Optim.NelderMead(),
        Optim.Options(time_limit=60);
        autodiff=:forward
    )
    assign_PL_OH_TF_thicknesses(res.minimizer, :tf)
    _step(actor.fluxswing_actor; operate_at_j_crit=unconstrained_flattop_duration, j_tolerance, only=:tf)
    _step(actor.stresses_actor)
    if verbose
        display(res)
    end

    # combined plug+OH+TF optimization
    res = nothing
    if (dd.solid_mechanics.center_stack.bucked == 1 || dd.solid_mechanics.center_stack.noslip == 1 || dd.solid_mechanics.center_stack.plug == 1)
        old_thicknesses = [layer.thickness for layer in dd.build.layer]
        res = Optim.optimize(
            x0 -> cost(x0, :all),
            [OH.thickness, TFhfs.thickness, dd.build.oh.technology.fraction_stainless, dd.build.tf.technology.fraction_stainless],
            Optim.NelderMead(),
            Optim.Options(time_limit=60, iterations=1000);
            autodiff=:forward
        )
        assign_PL_OH_TF_thicknesses(res.minimizer, :all)
        _step(actor.fluxswing_actor; operate_at_j_crit=unconstrained_flattop_duration, j_tolerance)
        _step(actor.stresses_actor)
        if verbose
            display(res)
        end
    end

    R0 = (TFhfs.end_radius + TFlfs.start_radius) / 2.0
    a = plasma.thickness / 2.0
    ϵ = R0 / a

    if debug_plot
        p = plot(yscale=:log10)
        plot!(p, C_JOH ./ (C_JOH .> 0.0), label="Jcrit OH")
        plot!(p, C_SOH ./ (C_SOH .> 0.0), label="Stresses OH")
        plot!(p, C_JTF ./ (C_JTF .> 0.0), label="Jcrit TF")
        plot!(p, C_STF ./ (C_STF .> 0.0), label="Stresses TF")
        display(p)
    end

    if verbose
        R0 = (TFhfs.end_radius + TFlfs.start_radius) / 2.0
        @show target_B0
        @show dd.build.tf.max_b_field * TFhfs.end_radius / R0

        @show dd.build.oh.flattop_duration
        @show dd.requirements.flattop_duration

        @show dd.build.oh.max_j
        @show dd.build.oh.critical_j

        @show dd.build.tf.max_j
        @show dd.build.tf.critical_j

        @show maximum(dd.solid_mechanics.center_stack.stress.vonmises.oh)
        @show stainless_steel.yield_strength

        @show maximum(dd.solid_mechanics.center_stack.stress.vonmises.tf)
        @show stainless_steel.yield_strength

        @show ϵ
        @show old_ϵ
    end

    function rel_error(value, target) # relative error with tolerance
        return abs((value .- target) ./ target)
    end

    max_B0 = dd.build.tf.max_b_field / TFhfs.end_radius * R0
    @assert target_B0 < max_B0 "TF cannot achieve requested B0 ($target_B0 --> $max_B0)"

    @assert dd.build.oh.max_j < dd.build.oh.critical_j
    @assert dd.build.tf.max_j < dd.build.tf.critical_j
    @assert maximum(dd.solid_mechanics.center_stack.stress.vonmises.oh) < stainless_steel.yield_strength
    @assert maximum(dd.solid_mechanics.center_stack.stress.vonmises.tf) < stainless_steel.yield_strength
    if !unconstrained_flattop_duration
        @assert rel_error(dd.build.oh.flattop_duration, dd.requirements.flattop_duration) < 0.1 "Relative error on flattop duration is more than 10% ($(dd.build.oh.flattop_duration) --> $(dd.requirements.flattop_duration))"
    end
    if fixed_aspect_ratio
        @assert rel_error(ϵ, old_ϵ) < 0.1 "ActorHFSsizing: plasma aspect ratio changed more than 10% ($old_ϵ --> $ϵ)"
    end

    return actor
end