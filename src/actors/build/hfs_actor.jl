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
    par::OverrideParameters{P,FUSEparameters__ActorHFSsizing{P}}
    act::ParametersAllActors{P}
    stresses_actor::ActorStresses{D,P}
    fluxswing_actor::ActorFluxSwing{D,P}
end

"""
    ActorHFSsizing(dd::IMAS.dd, act::ParametersAllActors; kw...)

Actor that resizes the High Field Side of the tokamak radial build.
It changes the radial build of the center stack (plug, OH, and TF)
accounting for stresses, superconductors critical currents, flux swing, and field requirements.

!!! note

    Manipulates radial build information in `dd.build.layer`
"""
function ActorHFSsizing(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorHFSsizing(dd, act.ActorHFSsizing, act; kw...)
    if actor.par.do_plot
        p = plot(dd.build.layer)
    end
    step(actor)
    finalize(actor)
    if actor.par.do_plot
        display(plot!(p, dd.build; cx=false))
        display(plot(dd.solid_mechanics.center_stack.stress))
    end
    return actor
end

function ActorHFSsizing(dd::IMAS.dd, par::FUSEparameters__ActorHFSsizing, act::ParametersAllActors; kw...)
    logging_actor_init(ActorHFSsizing)
    par = OverrideParameters(par; kw...)
    fluxswing_actor = ActorFluxSwing(dd, act.ActorFluxSwing)
    stresses_actor = ActorStresses(dd, act.ActorStresses)
    return ActorHFSsizing(dd, par, act, stresses_actor, fluxswing_actor)
end

function _step(actor::ActorHFSsizing)
    dd = actor.dd
    par = actor.par
    cs = dd.solid_mechanics.center_stack
    eqt = dd.equilibrium.time_slice[]

    # run optimization on a lower resolution grid
    original_n_points = actor.stresses_actor.par.n_points
    actor.stresses_actor.par.n_points = 5

    # Relative error with tolerance used for currents and stresses (not flattop)
    # NOTE: we divide by (abs(target) + 1.0) because critical currents can drop to 0.0!
    function target_value(value, target, tolerance)
        return (value .* (1.0 .+ tolerance) - target) ./ (abs(target) + 1.0)
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
        if length(x0) == 5
            dd.build.tf.nose_hfs_fraction = x0[5]
        end
    end

    function cost(x0)
        # assign optimization arguments
        assign_PL_OH_TF(x0)

        # evaluate coils currents and stresses
        finalize(step(actor.fluxswing_actor))
        finalize(step(actor.stresses_actor))

        # c_XXX are done such that <0 is OK, and >0 is BAD

        # OH currents
        if dd.requirements.coil_j_margin >= 0
            c_joh = (1.0 + dd.requirements.coil_j_margin) - dd.build.oh.critical_j / dd.build.oh.max_j # we want max_j to be coil_j_margin% below critical_j
            c_joh += 1.0 / (0.1 + dd.build.oh.critical_j) # critical_j can go to 0.0! stay away from current quench conditions
        else
            c_joh = 0.0
        end

        # TF currents
        if dd.requirements.coil_j_margin >= 0
            c_jtf = (1.0 + dd.requirements.coil_j_margin) - dd.build.tf.critical_j / dd.build.tf.max_j # we want max_j to be coil_j_margin% below critical_j
            c_jtf += 1.0 / (0.1 + dd.build.tf.critical_j) # critical_j can go to 0.0! stay away from current quench conditions
        else
            c_jtf = 0.0
        end

        # OH stresses
        if dd.requirements.coil_stress_margin >= 0
            c_soh = (1.0 + dd.requirements.coil_stress_margin) - cs.properties.yield_strength.oh / maximum(cs.stress.vonmises.oh) # we want stress to be coil_stress_margin% below yield_strength
        else
            c_soh = 0.0
        end

        # OH stresses
        if dd.requirements.coil_stress_margin >= 0
            c_stf = (1.0 + dd.requirements.coil_stress_margin) - cs.properties.yield_strength.tf / maximum(cs.stress.vonmises.tf) # we want stress to be coil_stress_margin% below yield_strength
        else
            c_stf = 0.0
        end

        # plug stresses
        if (dd.requirements.coil_stress_margin >= 0) && !ismissing(cs.stress.vonmises, :pl)
            c_spl = (1.0 + dd.requirements.coil_stress_margin) - cs.properties.yield_strength.pl / maximum(cs.stress.vonmises.pl) # we want stress to be coil_stress_margin% below yield_strength
        else
            c_spl = 0.0
        end

        # flattop
        if (dd.requirements.coil_j_margin >= 0) && !ismissing(dd.requirements, :flattop_duration)
            c_flt = (1.0 + dd.requirements.coil_j_margin) - dd.build.oh.flattop_duration / dd.requirements.flattop_duration
            c_flt_abs = 0.0
        else
            c_flt = 0.0
            c_flt_abs = - dd.build.oh.flattop_duration
        end

        # margins
        margins = [c_joh, c_soh, c_jtf, c_stf, c_spl, c_flt]
        c_mgn = norm(max.(margins, 0.0))
        c_Δmn = (maximum(margins) - minimum(margins))^2

        # want smallest possible TF and OH
        # for all things being equal, maximizing steel is good to keep the cost of the magnets down
        if nose
            c_scs = norm((
                OH.thickness / CPradius * (1.0 - dd.build.oh.technology.fraction_steel),
                TFhfs.thickness / CPradius * (1.0 - dd.build.tf.nose_hfs_fraction) * (1.0 - dd.build.tf.technology.fraction_steel)
            ))
        else
            c_scs = norm((
                OH.thickness / CPradius * (1.0 - dd.build.oh.technology.fraction_steel),
                TFhfs.thickness / CPradius * (1.0 - dd.build.tf.technology.fraction_steel)
            ))
        end

        push!(C_JOH, c_joh)
        push!(C_SOH, c_soh)
        push!(C_JTF, c_jtf)
        push!(C_STF, c_stf)
        push!(C_SPL, c_spl)
        push!(C_FLT, c_flt)
        push!(C_SCS, c_scs)
        push!(C_MGN, c_mgn)
        push!(C_ΔMG, c_Δmn)

        # total cost and constraints
        return norm([c_scs, c_mgn, c_Δmn]) + c_flt_abs, [maximum(max.(0.0, margins))], [0.0]
    end

    # initialize
    old_build = deepcopy(dd.build)
    PL = dd.build.layer[1]
    OH = IMAS.get_build_layer(dd.build.layer; type=_oh_)
    TFhfs = IMAS.get_build_layer(dd.build.layer; type=_tf_, fs=_hfs_)
    TFlfs = IMAS.get_build_layer(dd.build.layer; type=_tf_, fs=_lfs_)
    CPradius = TFhfs.end_radius
    OHTFgap = CPradius - TFhfs.thickness - OH.thickness - PL.thickness
    R0, B0 = eqt.global_quantities.vacuum_toroidal_field.r0, eqt.global_quantities.vacuum_toroidal_field.b0
    target_B0 = abs(B0)

    C_JOH = Float64[]
    C_SOH = Float64[]
    C_JTF = Float64[]
    C_STF = Float64[]
    C_SPL = Float64[]
    C_FLT = Float64[]
    C_SCS = Float64[]
    C_MGN = Float64[]
    C_ΔMG = Float64[]

    if Bool(dd.solid_mechanics.center_stack.bucked)
        nose = false
        dd.build.tf.nose_hfs_fraction = 0.0
    else
        nose = true
    end

    # optimization
    old_logging = actor_logging(dd, false)
    res = nothing
    try
        if nose
            bounds = (
                [0.01, 0.01, 0.1, 0.1, 0.01],
                [0.9, 0.9, 1.0 - dd.build.oh.technology.fraction_void - 0.1, 1.0 - dd.build.tf.technology.fraction_void - 0.1, 0.99])
        else
            bounds = (
                [0.01, 0.01, 0.1, 0.1],
                [0.99, 0.99, 1.0 - dd.build.oh.technology.fraction_void - 0.1, 1.0 - dd.build.tf.technology.fraction_void - 0.1])
        end

        options = Metaheuristics.Options(; seed=1, iterations=100)
        algorithm = Metaheuristics.ECA(; N=50, options)
        IMAS.refreeze!(dd.core_profiles.profiles_1d[], :conductivity_parallel)
        res = Metaheuristics.optimize(cost, bounds, algorithm)
        IMAS.empty!(dd.core_profiles.profiles_1d[], :conductivity_parallel)
        assign_PL_OH_TF(Metaheuristics.minimizer(res))
    finally
        actor_logging(dd, old_logging)
    end

    # final value (and stresses_actor on higher fidelity grid)
    finalize(step(actor.fluxswing_actor))
    actor.stresses_actor.par.n_points = original_n_points
    finalize(step(actor.stresses_actor))

    function print_details(; do_plot::Bool=par.do_plot)
        if res !== nothing
            print(res)
        end

        if do_plot
            p = plot(; yscale=:log10, legend=:bottomleft)
            plot!(p, C_MGN ./ (C_MGN .> 0.0); label="satisfy tolerances", alpha=0.9)
            plot!(p, C_SCS ./ (C_SCS .> 0.0); label="minimize superconductor", alpha=0.9)
            plot!(p, C_ΔMG ./ (C_ΔMG .> 0.0); label="clear tolerances by an equal amount", alpha=0.9)
            display(p)

            p = plot(; yscale=:log10, legend=:bottomleft)
            scatter!(p, C_JOH ./ (C_JOH .> 0.0); label="Jcrit OH constraint", alpha=0.25)
            scatter!(p, C_JTF ./ (C_JTF .> 0.0); label="Jcrit TF constraint", alpha=0.25)
            scatter!(p, C_SPL ./ (C_SPL .> 0.0); label="stress PL constraint", alpha=0.25)
            scatter!(p, C_SOH ./ (C_SOH .> 0.0); label="stress OH constraint", alpha=0.25)
            scatter!(p, C_STF ./ (C_STF .> 0.0); label="stress TF constraint", alpha=0.25)
            scatter!(p, C_FLT ./ (C_FLT .> 0.0); label="flattop constraint", alpha=0.25)
            display(p)
        end

        @show Bool(dd.solid_mechanics.center_stack.bucked)
        @show Bool(dd.solid_mechanics.center_stack.noslip)
        @show Bool(dd.solid_mechanics.center_stack.plug)
        println()
        @show [PL.thickness]
        @show [OH.thickness, dd.build.oh.technology.fraction_steel]
        @show [TFhfs.thickness, dd.build.tf.technology.fraction_steel]
        @show [dd.build.tf.nose_hfs_fraction]
        println()
        @show dd.requirements.coil_stress_margin
        @show dd.requirements.coil_j_margin
        println()
        @show target_B0
        @show dd.build.tf.max_b_field * TFhfs.end_radius / R0
        println()
        @show dd.build.oh.flattop_duration
        if !ismissing(dd.requirements, :flattop_duration)
            @show dd.requirements.flattop_duration
            @show dd.build.oh.flattop_duration / dd.requirements.flattop_duration
        end
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
            @show cs.properties.yield_strength.pl / maximum(cs.stress.vonmises.pl)
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
            dd.build.tf.max_j < dd.build.tf.critical_j,
            "TF exceeds critical current: $(dd.build.tf.max_j / dd.build.tf.critical_j * 100)%",
            par.error_on_technology,
            success)
        success = assert_conditions(
            dd.build.oh.max_j < dd.build.oh.critical_j,
            "OH exceeds critical current: $(dd.build.oh.max_j / dd.build.oh.critical_j * 100)%",
            par.error_on_technology,
            success)
        if !ismissing(cs.stress.vonmises, :pl)
            success = assert_conditions(
                maximum(cs.stress.vonmises.pl) < cs.properties.yield_strength.pl,
                "PL stresses are too high: $(maximum(cs.stress.vonmises.pl) / cs.properties.yield_strength.pl * 100)%",
                par.error_on_technology,
                success)
        end
        success = assert_conditions(
            maximum(cs.stress.vonmises.oh) < cs.properties.yield_strength.oh,
            "OH stresses are too high: $(maximum(cs.stress.vonmises.oh) / cs.properties.yield_strength.oh * 100)%",
            par.error_on_technology,
            success)
        success = assert_conditions(
            maximum(cs.stress.vonmises.tf) < cs.properties.yield_strength.tf,
            "TF stresses are too high: $(maximum(cs.stress.vonmises.tf) / cs.properties.yield_strength.tf * 100)%",
            par.error_on_technology,
            success)

        # performance checks
        max_B0 = dd.build.tf.max_b_field / TFhfs.end_radius * R0
        success = assert_conditions(target_B0 < max_B0,
            "TF cannot achieve requested B0 ($target_B0 [T] instead of $max_B0 [T])",
            par.error_on_performance,
            success)
        if ismissing(dd.requirements, :flattop_duration)
            success = assert_conditions(
                dd.build.oh.flattop_duration > 0.0,
                "OH cannot achieve flattop ($(dd.build.oh.flattop_duration) [s])",
                par.error_on_performance,
                success)
        else
            success = assert_conditions(
                dd.build.oh.flattop_duration > dd.requirements.flattop_duration,
                "OH cannot achieve requested flattop ($(dd.build.oh.flattop_duration) [s] insted of $(dd.requirements.flattop_duration) [s])",
                par.error_on_performance,
                success)
        end

        @assert success "HFS sizing cannot be done within constraints"

        if par.verbose
            print_details()
        end

    catch e
        print_details(; do_plot=true)
        plot(eqt; cx=true)
        plot!(old_build)
        display(plot!(dd.build; cx=false))
        display(plot(dd.solid_mechanics.center_stack.stress))
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