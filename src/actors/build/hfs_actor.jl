#= ========== =#
#  HFS sizing  #
#= ========== =#
Base.@kwdef mutable struct FUSEparameters__ActorHFSsizing{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    j_tolerance::Entry{T} = Entry(T, "-", "Tolerance on the conductor current limits"; default=0.4)
    stress_tolerance::Entry{T} = Entry(T, "-", "Tolerance on the structural stresses limits"; default=0.2)
    aspect_ratio_tolerance::Entry{T} = Entry(T, "-", "Tolerance on the aspect_ratio change"; default=0.01)
    operate_at_j_crit::Entry{Bool} = Entry(Bool, "-", "Maximize flux_duration without targeting a specific value"; default=true)
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
    par = act.ActorHFSsizing
    actor = ActorHFSsizing(dd, par, act; kw_ActorFluxSwing, kw_ActorStresses, kw...)
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
    aspect_ratio_tolerance::Real=actor.par.aspect_ratio_tolerance,
    operate_at_j_crit::Bool=actor.par.operate_at_j_crit,
    verbose::Bool=actor.par.verbose
)

    #Relative error with tolerance
    #NOTE: we divide by (abs(target) + 1.0) because critical currents can drop to 0.0!
    function target_value(value, target, tolerance)
        return abs((value .* (1.0 .+ tolerance) .- target) ./ (abs(target) + 1.0))
    end

    function assign_PL_OH_TF(x0)
        # assign optimization arguments
        OH.thickness, c_extra1 = mirror_bound_w_cost(x0[1], 0.0, 100.0)
        TFhfs.thickness, c_extra2 = mirror_bound_w_cost(x0[2], 0.0, 100.0)
        dd.build.oh.technology.fraction_stainless, c_extra3 = mirror_bound_w_cost(x0[3], 0.45, 1.0 - dd.build.oh.technology.fraction_void - 0.05)
        dd.build.tf.technology.fraction_stainless, c_extra4 = mirror_bound_w_cost(x0[4], 0.45, 1.0 - dd.build.tf.technology.fraction_void - 0.05)

        # NOTE: the plug expands/contracts to keep original plasma radius constant
        #       but it does not contract more than 1/4 of the OH thickness
        plug.thickness += old_plasma_start_radius - plasma.start_radius
        plug.thickness = max(OH.thickness / 4.0, plug.thickness)

        return c_extra1 + c_extra2 + c_extra3 + c_extra4
    end

    function cost(x0)
        # assign optimization arguments
        c_extra = assign_PL_OH_TF(x0)

        # evaluate coils currents and stresses
        _step(actor.fluxswing_actor; operate_at_j_crit, j_tolerance)
        _step(actor.stresses_actor)

        # OH sizing
        c_joh = target_value(dd.build.oh.max_j, dd.build.oh.critical_j, j_tolerance) # we want max_j to be 20% below critical_j
        c_soh = target_value(maximum(dd.solid_mechanics.center_stack.stress.vonmises.oh), stainless_steel.yield_strength, stress_tolerance)  # we want stress to be 20% below yield_strength
       
        # plug sizing
        if !ismissing(dd.solid_mechanics.center_stack.stress.vonmises, :pl)
            c_spl = target_value(maximum(dd.solid_mechanics.center_stack.stress.vonmises.pl), stainless_steel.yield_strength, stress_tolerance)
        else
            c_spl = 0.0
        end

        # TF sizing
        c_jtf = c_stf = 0.0
        c_jtf = target_value(dd.build.tf.max_j, dd.build.tf.critical_j, j_tolerance) # we want max_j to be j_tolerance% below critical_j
        c_stf = target_value(maximum(dd.solid_mechanics.center_stack.stress.vonmises.tf), stainless_steel.yield_strength, stress_tolerance)  # we want stress to be stress_tolerance% below yield_strength

        if verbose
            push!(C_JOH, c_joh)
            push!(C_SOH, c_soh)
            push!(C_JTF, c_jtf)
            push!(C_STF, c_stf)
            push!(C_SPL, c_spl)
            push!(C_XTR, c_extra)
        end

        # total cost
        return norm([c_joh; c_jtf; c_spl; c_soh; c_stf; c_extra])
    end

    @assert actor.stresses_actor.dd === actor.fluxswing_actor.dd
    dd = actor.stresses_actor.dd

    # initialize
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

    dd.build.oh.technology.fraction_stainless = 0.5
    dd.build.tf.technology.fraction_stainless = 0.5

    if verbose
        C_JOH = Float64[]
        C_SOH = Float64[]
        C_JTF = Float64[]
        C_STF = Float64[]
        C_SPL = Float64[]
        C_XTR = Float64[]
    end

    # optimization
    res = Optim.optimize(
        x0 -> cost(x0),
        [OH.thickness, TFhfs.thickness, dd.build.oh.technology.fraction_stainless, dd.build.tf.technology.fraction_stainless],
        Optim.NelderMead(),
        Optim.Options(iterations=1000);
        autodiff=:forward
    )
    assign_PL_OH_TF(res.minimizer)
    step(actor.fluxswing_actor; operate_at_j_crit, j_tolerance)
    step(actor.stresses_actor)
    if verbose
        display(res)
    end

    R0 = (TFhfs.end_radius + TFlfs.start_radius) / 2.0
    a = plasma.thickness / 2.0
    ϵ = R0 / a

    if verbose
        p = plot(yscale=:log10, legend=:topright)
        plot!(p, C_JOH, label="cost Jcrit OH")
        plot!(p, C_JTF, label="cost Jcrit TF")
        if !ismissing(dd.solid_mechanics.center_stack.stress.vonmises, :pl)
            plot!(p, C_SPL, label="cost stresses PL")
        end
        plot!(p, C_SOH, label="cost stresses OH")
        plot!(p, C_STF, label="cost stresses TF")
        if sum(C_XTR) > 0.0
            plot!(p, C_XTR ./ (C_XTR .> 0.0), label="cost optimizer penalty")
        end
        display(p)
    end

    target_B0 = maximum(abs.(dd.equilibrium.vacuum_toroidal_field.b0))
    R0_of_B0 = dd.equilibrium.vacuum_toroidal_field.r0

    if verbose
        @show [OH.thickness, dd.build.oh.technology.fraction_stainless]
        @show [TFhfs.thickness, dd.build.tf.technology.fraction_stainless]
        println()
        @show target_B0
        @show dd.build.tf.max_b_field * TFhfs.end_radius / R0_of_B0
        println()
        @show dd.build.oh.flattop_duration
        @show dd.requirements.flattop_duration
        println()
        @show dd.build.oh.max_j
        @show dd.build.oh.critical_j
        println()
        @show dd.build.tf.max_j
        @show dd.build.tf.critical_j
        println()
        @show maximum(dd.solid_mechanics.center_stack.stress.vonmises.oh)
        @show stainless_steel.yield_strength
        println()
        @show maximum(dd.solid_mechanics.center_stack.stress.vonmises.tf)
        @show stainless_steel.yield_strength
        println()
        @show ϵ
        @show old_ϵ
    end

    function rel_error(value, target) # relative error with tolerance
        return abs((value .- target) ./ target)
    end

    max_B0 = dd.build.tf.max_b_field / TFhfs.end_radius * R0_of_B0
    @assert target_B0 < max_B0 "TF cannot achieve requested B0 ($target_B0 --> $max_B0)"
    @assert dd.build.oh.max_j < dd.build.oh.critical_j
    @assert dd.build.tf.max_j < dd.build.tf.critical_j
    @assert maximum(dd.solid_mechanics.center_stack.stress.vonmises.oh) < stainless_steel.yield_strength
    @assert maximum(dd.solid_mechanics.center_stack.stress.vonmises.tf) < stainless_steel.yield_strength
    if !operate_at_j_crit
        @assert rel_error(dd.build.oh.flattop_duration, dd.requirements.flattop_duration) <= 0.1 "Relative error on flattop duration is more than 10% ($(dd.build.oh.flattop_duration) --> $(dd.requirements.flattop_duration))"
    end
    @assert rel_error(ϵ, old_ϵ) <= aspect_ratio_tolerance "Plasma aspect ratio changed more than $(aspect_ratio_tolerance) ($old_ϵ --> $ϵ)"

    return actor
end