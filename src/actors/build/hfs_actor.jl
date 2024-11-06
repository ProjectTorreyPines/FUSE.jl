#= ========== =#
#  HFS sizing  #
#= ========== =#
Base.@kwdef mutable struct FUSEparameters__ActorHFSsizing{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    error_on_technology::Entry{Bool} = Entry{Bool}("-", "Error if build stresses and current limits are not met"; default=true)
    error_on_performance::Entry{Bool} = Entry{Bool}("-", "Error if requested Bt and flattop duration are not met"; default=true)
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
    verbose::Entry{Bool} = act_common_parameters(; verbose=false)
end

mutable struct ActorHFSsizing{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorHFSsizing{P}
    stresses_actor::ActorStresses{D,P}
    fluxswing_actor::ActorFluxSwing{D,P}
end

"""
    ActorHFSsizing(dd::IMAS.dd, act::ParametersAllActors; kw...)

Actor that resizes the High Field Side of the tokamak radial build

  - takes into account the OH maximum allowed superconductor current/Field
  - takes into account the stresses on the center stack

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
    cs = dd.solid_mechanics.center_stack
    eqt = dd.equilibrium.time_slice[]

    # Relative error with tolerance used for currents and stresses (not flattop)
    # NOTE: we divide by (abs(target) + 1.0) because critical currents can drop to 0.0!
    # NOTE: we stronlgy penalize going above target, and only gently encourage not going below it (since
    function target_value(value, target, tolerance)
        tmp = (value .* (1.0 .+ tolerance) - target) ./ (abs(target) + 1.0)
        return sign(tmp) * tmp^2
    end

    function assign_PL_OH_TF(x0)
        # assign optimization arguments
        # x0[1]: normalized TF coil thickness
        # x0[2]: normalized OH coil thickness
        # x0[3]: OH steel fraction
        # x0[4]: TF steel fraction
        # x0[5]: TF nose thickness
        TFhfs.thickness = TFlfs.thickness = x0[1] * CPradius
        OH.thickness = x0[2] * (CPradius - TFhfs.thickness - OHTFgap)
        PL.thickness = CPradius - TFhfs.thickness - OH.thickness - OHTFgap
        dd.build.oh.technology.fraction_steel = x0[3]
        dd.build.tf.technology.fraction_steel = x0[4]
        dd.build.tf.nose_hfs_fraction = x0[5]

        # want smallest possible TF and OH
        # keeping them of similar size is a hint for the optimizer to achieve convergence
        # for all things being equal, maximizing steel is good to keep the cost of the magnets down
        return 1E-3 * (
            ((OH.thickness + TFhfs.thickness) / CPradius)^2 + # favor small OH and TF thicknesses
            0.1 * ((OH.thickness - TFhfs.thickness) / CPradius)^2 + # favor OH and TF similiar thickness
            0.1 * (1.0 - dd.build.oh.technology.fraction_steel)^2 + # favor steel over superconductor in OH coil
            0.1 * (1.0 - dd.build.tf.technology.fraction_steel)^2)  # favor steel over superconductor in TF coil
    end

    function cost(x0)
        # assign optimization arguments
        c_cst = assign_PL_OH_TF(x0)

        # evaluate coils currents and stresses
        finalize(step(actor.fluxswing_actor))
        finalize(step(actor.stresses_actor))

        # OH currents and stresses
        if actor.fluxswing_actor.par.operate_oh_at_j_crit && (dd.requirements.coil_j_margin >= 0)
            c_joh = target_value(dd.build.oh.max_j, dd.build.oh.critical_j, dd.requirements.coil_j_margin) # we want max_j to be coil_j_margin% below critical_j
        else
            c_joh = 0.0
        end

        if (dd.requirements.coil_stress_margin >= 0)
            c_soh = target_value(maximum(cs.stress.vonmises.oh), cs.properties.yield_strength.oh, dd.requirements.coil_stress_margin) # we want stress to be coil_stress_margin% below yield_strength
        else
            c_soh = 0.0
        end

        # TF currents and stresses
        if (dd.requirements.coil_j_margin >= 0)
            c_jtf = target_value(dd.build.tf.max_j, dd.build.tf.critical_j, dd.requirements.coil_j_margin) # we want max_j to be coil_j_margin% below critical_j
        else
            c_jtf = 0.0
        end

        if (dd.requirements.coil_stress_margin >= 0)
            c_stf = target_value(maximum(cs.stress.vonmises.tf), cs.properties.yield_strength.tf, dd.requirements.coil_stress_margin) # we want stress to be coil_stress_margin% below yield_strength
        else
            c_stf = 0.0
        end

        # plug stresses
        if !ismissing(cs.stress.vonmises, :pl) && (dd.requirements.coil_stress_margin >= 0)
            c_spl = target_value(maximum(cs.stress.vonmises.pl), cs.properties.yield_strength.pl, dd.requirements.coil_stress_margin)
        else
            c_spl = 0.0
        end

        # flattop
        # Additional 10% of flattop duration for shape control
        if actor.fluxswing_actor.par.operate_oh_at_j_crit
            c_flt = target_value(dd.requirements.flattop_duration, dd.build.oh.flattop_duration, 0.1)
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

        # total cost and constraints
        return c_cst, [c_joh, c_soh, c_flt, c_jtf, c_stf, c_spl], [0.0]
    end

    # initialize
    old_build = deepcopy(dd.build)
    PL = dd.build.layer[1]
    OH = IMAS.get_build_layer(dd.build.layer; type=_oh_)
    TFhfs = IMAS.get_build_layer(dd.build.layer; type=_tf_, fs=_hfs_)
    TFlfs = IMAS.get_build_layer(dd.build.layer; type=_tf_, fs=_lfs_)
    plasma = IMAS.get_build_layer(dd.build.layer; type=_plasma_)
    CPradius = TFhfs.end_radius
    OHTFgap = CPradius - TFhfs.thickness - OH.thickness - PL.thickness
    R0, B0 = eqt.global_quantities.vacuum_toroidal_field.r0, eqt.global_quantities.vacuum_toroidal_field.b0
    target_B0 = abs(B0)

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
    old_logging = actor_logging(dd, false)
    res = nothing
    try
        bounds = (
            [0.1, 0.1, 0.1, 0.1, 0.0],
            [0.9, 0.9, 1.0 - dd.build.oh.technology.fraction_void - 0.1, 1.0 - dd.build.tf.technology.fraction_void - 0.1, 0.9])
        options = Metaheuristics.Options(; seed=1, iterations=50)
        algorithm = Metaheuristics.ECA(; N=50, options)
        IMAS.refreeze!(dd.core_profiles.profiles_1d[], :conductivity_parallel)
        res = Metaheuristics.optimize(cost, bounds, algorithm)
        IMAS.empty!(dd.core_profiles.profiles_1d[], :conductivity_parallel)
        assign_PL_OH_TF(Metaheuristics.minimizer(res))
    finally
        actor_logging(dd, old_logging)
    end
    finalize(step(actor.fluxswing_actor))
    finalize(step(actor.stresses_actor))

    function print_details()
        if res !== nothing
            println(cost(Metaheuristics.minimizer(res)))
            print(res)
        end

        if par.verbose
            p = plot(; yscale=:log10, legend=:topright)
            plot!(p, C_JOH ./ (C_JOH .> 0.0); label="Jcrit OH constraint")
            plot!(p, C_JTF ./ (C_JTF .> 0.0); label="Jcrit TF constraint")
            plot!(p, C_SPL ./ (C_SPL .> 0.0); label="stress PL constraint")
            plot!(p, C_SOH ./ (C_SOH .> 0.0); label="stress OH constraint")
            plot!(p, C_STF ./ (C_STF .> 0.0); label="stress TF constraint")
            plot!(p, C_FLT ./ (C_FLT .> 0.0); label="flattop constraint")
            plot!(p, C_CST ./ (C_CST .> 0.0); label="small TF & OH cost", color=:black)
            display(p)
        end

        @show [PL.thickness]
        @show [OH.thickness, dd.build.oh.technology.fraction_steel]
        @show [TFhfs.thickness, dd.build.tf.technology.fraction_steel]
        @show [dd.build.tf.nose_hfs_fraction]
        println()
        @show target_B0
        @show dd.build.tf.max_b_field * TFhfs.end_radius / R0
        println()
        @show dd.build.oh.flattop_duration
        @show dd.requirements.flattop_duration
        @show dd.build.oh.flattop_duration / dd.requirements.flattop_duration
        println()
        @show dd.build.oh.max_j
        @show dd.build.oh.critical_j
        @show dd.build.oh.critical_j / dd.build.oh.max_j
        println()
        @show dd.build.tf.max_j
        @show dd.build.tf.critical_j
        @show dd.build.tf.critical_j / dd.build.tf.max_j
        if !ismissing(cs.stress.vonmises, :pl)
            println()
            @show maximum(cs.stress.vonmises.pl)
            @show cs.properties.yield_strength.pl
        end
        println()
        @show maximum(cs.stress.vonmises.oh)
        @show cs.properties.yield_strength.oh
        @show cs.properties.yield_strength.oh / maximum(cs.stress.vonmises.oh)
        println()
        @show maximum(cs.stress.vonmises.tf)
        @show cs.properties.yield_strength.tf
        @show cs.properties.yield_strength.tf / maximum(cs.stress.vonmises.tf)
    end

    try
        success = true
        # technology checks
        success = assert_conditions(
            dd.build.tf.max_j .* (1.0 .+ dd.requirements.coil_j_margin * 0.9) < dd.build.tf.critical_j,
            "TF exceeds critical current: $(dd.build.tf.max_j .* (1.0 .+ dd.requirements.coil_j_margin) / dd.build.tf.critical_j * 100)%",
            par.error_on_technology,
            success
        )
        success = assert_conditions(
            dd.build.oh.max_j .* (1.0 .+ dd.requirements.coil_j_margin * 0.9) < dd.build.oh.critical_j,
            "OH exceeds critical current: $(dd.build.oh.max_j .* (1.0 .+ dd.requirements.coil_j_margin) / dd.build.oh.critical_j * 100)%",
            par.error_on_technology,
            success
        )
        if !ismissing(cs.stress.vonmises, :pl)
            success = assert_conditions(
                maximum(cs.stress.vonmises.pl) .* (1.0 .+ dd.requirements.coil_stress_margin * 0.9) < cs.properties.yield_strength.pl,
                "PL stresses are too high: $(maximum(cs.stress.vonmises.pl) .* (1.0 .+ dd.requirements.coil_stress_margin) / cs.properties.yield_strength.pl * 100)%",
                par.error_on_technology,
                success
            )
        end
        success = assert_conditions(
            maximum(cs.stress.vonmises.oh) .* (1.0 .+ dd.requirements.coil_stress_margin * 0.9) < cs.properties.yield_strength.oh,
            "OH stresses are too high: $(maximum(cs.stress.vonmises.oh) .* (1.0 .+ dd.requirements.coil_stress_margin) / cs.properties.yield_strength.oh * 100)%",
            par.error_on_technology,
            success
        )
        success = assert_conditions(
            maximum(cs.stress.vonmises.tf) .* (1.0 .+ dd.requirements.coil_stress_margin * 0.9) < cs.properties.yield_strength.tf,
            "TF stresses are too high: $(maximum(cs.stress.vonmises.tf) .* (1.0 .+ dd.requirements.coil_stress_margin) / cs.properties.yield_strength.tf * 100)%",
            par.error_on_technology,
            success
        )

        # performance checks
        max_B0 = dd.build.tf.max_b_field / TFhfs.end_radius * R0
        success = assert_conditions(target_B0 < max_B0, "TF cannot achieve requested B0 ($target_B0 [T] instead of $max_B0 [T])", par.error_on_performance, success)
        success = assert_conditions(
            dd.build.oh.flattop_duration > dd.requirements.flattop_duration * 0.9,
            "OH cannot achieve requested flattop ($(dd.build.oh.flattop_duration) [s] insted of $(dd.requirements.flattop_duration) [s])",
            par.error_on_performance,
            success
        )

        @assert success "HFS sizing cannot be done within constraints"

        if par.verbose
            print_details()
        end

    catch e
        print_details()
        plot(eqt; cx=true)
        plot!(old_build)
        display(plot!(dd.build; cx=false))
        dd.build = old_build
        rethrow(e)
    end

    return actor
end

function assert_conditions(passing_condition::Bool, message::String, do_raise::Bool, success::Bool)
    if !passing_condition
        @warn message
    end
    return success && !(do_raise && !passing_condition)
end