#= ========== =#
#  HFS sizing  #
#= ========== =#
Base.@kwdef mutable struct FUSEparameters__ActorHFSsizing{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    j_tolerance::Entry{T} = Entry{T}("-", "Tolerance on the OH and TF current limits (overrides ActorFluxSwing.j_tolerance)"; default=0.4)
    stress_tolerance::Entry{T} = Entry{T}("-", "Tolerance on the OH and TF structural stresses limits"; default=0.2)
    error_on_technology::Entry{Bool} = Entry{Bool}("-", "Error if build stresses and current limits are not met"; default=true)
    error_on_performance::Entry{Bool} = Entry{Bool}("-", "Error if requested Bt and flattop duration are not met"; default=true)
    do_plot::Entry{Bool} = Entry{Bool}("-", "Plot"; default=false)
    verbose::Entry{Bool} = Entry{Bool}("-", "Verbose"; default=false)
end

mutable struct ActorHFSsizing{D,P} <: ReactorAbstractActor
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorHFSsizing{P}
    stresses_actor::ActorStresses{D,P}
    fluxswing_actor::ActorFluxSwing{D,P}
end

"""
    ActorHFSsizing(dd::IMAS.dd, act::ParametersAllActors; kw...)

Actor that resizes the High Field Side of the tokamak radial build
* takes into account the OH maximum allowed superconductor current/Field
* takes into account the stresses on the center stack
    
!!! note 
    Manipulates radial build information in `dd.build.layer`
"""
function ActorHFSsizing(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorHFSsizing(dd, act.ActorHFSsizing, act; kw...)
    if actor.par.do_plot
        p = plot(dd.build)
    end
    step(actor)
    finalize(actor)
    if actor.par.do_plot
        display(plot!(p, dd.build; cx=false))
    end
    return actor
end

function ActorHFSsizing(dd::IMAS.dd, par::FUSEparameters__ActorHFSsizing, act::ParametersAllActors; kw...)
    par = act.ActorHFSsizing(kw...)
    fluxswing_actor = ActorFluxSwing(dd, act.ActorFluxSwing)
    stresses_actor = ActorStresses(dd, act.ActorStresses)
    return ActorHFSsizing(dd, par, stresses_actor, fluxswing_actor)
end

function _step(actor::ActorHFSsizing)
    dd = actor.dd
    par = actor.par

    # modify j_tolerance in fluxswing_actor (since actor.fluxswing_actor.par is a copy, this does not affect act.ActorFluxSwing)
    actor.fluxswing_actor.par.j_tolerance = par.j_tolerance

    # Relative error with tolerance used for currents and stresses (not flattop)
    # NOTE: we divide by (abs(target) + 1.0) because critical currents can drop to 0.0!
    # NOTE: we stronlgy penalize going above target, and only gently encourage not going below it (since 
    function target_value(value, target, tolerance)
        tmp = (value .* (1.0 .+ tolerance) .- target) ./ (abs(target) + 1.0)
        if tmp > 0.0
            return exp(tmp) - 1.0
        else
            return -tmp * 1E-3
        end
    end

    function assign_PL_OH_TF(x0)
        # assign optimization arguments
        TFhfs.thickness = TFlfs.thickness = mirror_bound(x0[1], 0.1, 0.9) * CPradius
        OH.thickness = mirror_bound(x0[2], 0.1, 0.9) * (CPradius - TFhfs.thickness - OHTFgap)
        PL.thickness = CPradius - TFhfs.thickness - OH.thickness - OHTFgap
        dd.build.oh.technology.fraction_steel = mirror_bound(x0[3], 0.1, 1.0 - dd.build.oh.technology.fraction_void - 0.1)
        dd.build.tf.technology.fraction_steel = mirror_bound(x0[4], 0.1, 1.0 - dd.build.tf.technology.fraction_void - 0.1)

        # favor relatively small and equal TF and OH
        return norm(((OH.thickness + TFhfs.thickness) / CPradius, abs(OH.thickness - TFhfs.thickness) / CPradius))
    end

    function cost(x0)
        # assign optimization arguments
        c_cst = assign_PL_OH_TF(x0)
        c_cst *= 1E-3

        # evaluate coils currents and stresses
        _step(actor.fluxswing_actor)
        _step(actor.stresses_actor)

        # OH currents and stresses
        if actor.fluxswing_actor.par.operate_oh_at_j_crit && (par.j_tolerance >= 0)
            c_joh = target_value(dd.build.oh.max_j, dd.build.oh.critical_j, par.j_tolerance) # we want max_j to be j_tolerance% below critical_j
        else
            c_joh = 0.0
        end

        if (par.stress_tolerance >= 0)
            c_soh = target_value(maximum(dd.solid_mechanics.center_stack.stress.vonmises.oh), dd.solid_mechanics.center_stack.properties.yield_strength.oh, par.stress_tolerance) # we want stress to be stress_tolerance% below yield_strength
        else
            c_soh = 0.0
        end

        # TF currents and stresses
        if (par.j_tolerance >= 0)
            c_jtf = target_value(dd.build.tf.max_j, dd.build.tf.critical_j, par.j_tolerance) # we want max_j to be j_tolerance% below critical_j
        else
            c_jtf = 0.0
        end

        if (par.stress_tolerance >= 0)
            c_stf = target_value(maximum(dd.solid_mechanics.center_stack.stress.vonmises.tf), dd.solid_mechanics.center_stack.properties.yield_strength.tf, par.stress_tolerance) # we want stress to be stress_tolerance% below yield_strength
        else
            c_stf = 0.0
        end

        # plug stresses
        if !ismissing(dd.solid_mechanics.center_stack.stress.vonmises, :pl) && (par.stress_tolerance >= 0)
            c_spl = target_value(maximum(dd.solid_mechanics.center_stack.stress.vonmises.pl), dd.solid_mechanics.center_stack.properties.yield_strength.pl, par.stress_tolerance)
        else
            c_spl = 0.0
        end

        # flattop
        if actor.fluxswing_actor.par.operate_oh_at_j_crit
            c_flt = abs((dd.build.oh.flattop_duration - dd.requirements.flattop_duration) / dd.requirements.flattop_duration)
        else
            c_flt = 0.0
        end

        if par.verbose
            push!(C_JOH, c_joh)
            push!(C_SOH, c_soh)
            push!(C_JTF, c_jtf)
            push!(C_STF, c_stf)
            push!(C_SPL, c_spl)
            push!(C_FLT, c_flt)
            push!(C_CST, c_cst)
        end

        # total cost
        return norm((norm([c_joh, c_soh]), norm([c_jtf, c_stf]), c_spl, c_flt, c_cst))
    end

    # initialize
    PL = dd.build.layer[1]
    OH = IMAS.get_build_layer(dd.build.layer, type=_oh_)
    BLs = IMAS.get_build_layers(dd.build.layer, type=_blanket_, fs=_hfs_)
    if isempty(BLs)
        BL = PL
    else
        BL = BLs[1]
    end
    old_BL_thickness = BL.thickness
    TFhfs = IMAS.get_build_layer(dd.build.layer, type=_tf_, fs=_hfs_)
    TFlfs = IMAS.get_build_layer(dd.build.layer, type=_tf_, fs=_lfs_)
    plasma = IMAS.get_build_layer(dd.build.layer, type=_plasma_)

    CPradius = TFhfs.end_radius
    OHTFgap = CPradius - TFhfs.thickness - OH.thickness - PL.thickness

    target_B0 = maximum(abs, dd.equilibrium.vacuum_toroidal_field.b0)
    a = (plasma.end_radius - plasma.start_radius) / 2.0
    old_R0 = (plasma.end_radius + plasma.start_radius) / 2.0

    if par.verbose
        C_JOH = Float64[]
        C_SOH = Float64[]
        C_JTF = Float64[]
        C_STF = Float64[]
        C_SPL = Float64[]
        C_FLT = Float64[]
        C_CST = Float64[]
    end

    # optimization
    old_build = deepcopy(dd.build)
    res = Optim.optimize(
        x0 -> cost(x0),
        [TFhfs.thickness / CPradius, OH.thickness / (CPradius - TFhfs.thickness - OHTFgap), dd.build.oh.technology.fraction_steel, dd.build.tf.technology.fraction_steel],
        Optim.NelderMead(),
        Optim.Options(iterations=1000, g_tol=1e-6);
        autodiff=:forward
    )
    assign_PL_OH_TF(res.minimizer)
    step(actor.fluxswing_actor)
    step(actor.stresses_actor)

    R0 = (plasma.start_radius + plasma.end_radius) / 2.0

    function check_technology()
        @assert dd.build.tf.max_j .* (1.0 .+ par.j_tolerance * 0.9) < dd.build.tf.critical_j "TF exceeds critical current: $(dd.build.tf.max_j .* (1.0 .+ par.j_tolerance) / dd.build.tf.critical_j * 100)%"
        @assert dd.build.oh.max_j .* (1.0 .+ par.j_tolerance * 0.9) < dd.build.oh.critical_j "OH exceeds critical current: $(dd.build.oh.max_j .* (1.0 .+ par.j_tolerance) / dd.build.oh.critical_j * 100)%"
        if !ismissing(dd.solid_mechanics.center_stack.stress.vonmises, :pl)
            @assert maximum(dd.solid_mechanics.center_stack.stress.vonmises.pl) .* (1.0 .+ par.stress_tolerance * 0.9) < dd.solid_mechanics.center_stack.properties.yield_strength.pl "PL stresses are too high: $(maximum(dd.solid_mechanics.center_stack.stress.vonmises.pl) .* (1.0 .+ par.stress_tolerance) / dd.solid_mechanics.center_stack.properties.yield_strength.pl * 100)%"
        end
        @assert maximum(dd.solid_mechanics.center_stack.stress.vonmises.oh) .* (1.0 .+ par.stress_tolerance * 0.9) < dd.solid_mechanics.center_stack.properties.yield_strength.oh "OH stresses are too high: $(maximum(dd.solid_mechanics.center_stack.stress.vonmises.oh) .* (1.0 .+ par.stress_tolerance) / dd.solid_mechanics.center_stack.properties.yield_strength.oh * 100)%"
        @assert maximum(dd.solid_mechanics.center_stack.stress.vonmises.tf) .* (1.0 .+ par.stress_tolerance * 0.9) < dd.solid_mechanics.center_stack.properties.yield_strength.tf "TF stresses are too high: $(maximum(dd.solid_mechanics.center_stack.stress.vonmises.tf) .* (1.0 .+ par.stress_tolerance) / dd.solid_mechanics.center_stack.properties.yield_strength.tf * 100)%"
    end

    function check_performance()
        max_B0 = dd.build.tf.max_b_field / TFhfs.end_radius * R0
        @assert target_B0 < max_B0 "TF cannot achieve requested B0 ($target_B0 instead of $max_B0)"
        @assert dd.build.oh.flattop_duration .* (1.0 .+ par.j_tolerance * 0.9) > dd.requirements.flattop_duration "OH cannot achieve requested flattop ($(dd.build.oh.flattop_duration) insted of $(dd.requirements.flattop_duration))"
    end

    function print_details()
        print(res)

        if par.verbose
            p = plot(yscale=:log10, legend=:topright)
            if sum(C_JOH) > 0.0
                plot!(p, C_JOH, label="cost Jcrit OH")
            end
            plot!(p, C_JTF, label="cost Jcrit TF")
            if !ismissing(dd.solid_mechanics.center_stack.stress.vonmises, :pl)
                plot!(p, C_SPL, label="cost stress PL")
            end
            plot!(p, C_SOH, label="cost stress OH")
            plot!(p, C_STF, label="cost stress TF")
            if sum(C_FLT) > 0.0
                plot!(p, C_FLT, label="cost flattop")
            end
            plot!(p, C_CST, label="small TF & OH")
            display(p)
        end

        @show [PL.thickness]
        @show [OH.thickness, dd.build.oh.technology.fraction_steel]
        @show [TFhfs.thickness, dd.build.tf.technology.fraction_steel]
        println()
        @show target_B0
        @show dd.build.tf.max_b_field * TFhfs.end_radius / R0
        println()
        @show dd.build.oh.flattop_duration
        @show dd.requirements.flattop_duration
        println()
        @show dd.build.oh.max_j
        @show dd.build.oh.critical_j
        println()
        @show dd.build.tf.max_j
        @show dd.build.tf.critical_j
        if !ismissing(dd.solid_mechanics.center_stack.stress.vonmises, :pl)
            println()
            @show maximum(dd.solid_mechanics.center_stack.stress.vonmises.pl)
            @show dd.solid_mechanics.center_stack.properties.yield_strength.pl
        end
        println()
        @show maximum(dd.solid_mechanics.center_stack.stress.vonmises.oh)
        @show dd.solid_mechanics.center_stack.properties.yield_strength.oh
        println()
        @show maximum(dd.solid_mechanics.center_stack.stress.vonmises.tf)
        @show dd.solid_mechanics.center_stack.properties.yield_strength.tf
        println()
        @show old_R0 / a
        @show R0 / a
    end

    try
        if par.error_on_technology
            check_technology()
        end
        if par.error_on_performance
            check_performance()
        end
        if par.verbose
            print_details()
        end
    catch e
        print_details()
        dd.build = old_build
        rethrow(e)
    end

    return actor
end
