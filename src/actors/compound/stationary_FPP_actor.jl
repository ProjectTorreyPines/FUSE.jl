#= ===================== =#
#  ActorStationaryFPP  #
#= ===================== =#
# Define the parameters for the FPP workflow
Base.@kwdef mutable struct FUSEparameters__ActorStationaryFPP{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    convergence_error::Entry{T} = Entry{T}("-", "Convergence error threshold for FPP workflow"; default=5E-2)
    fix_ip::Entry{Bool} = Entry{Bool}("-", "Fix Ip for true) or adjust Ip for false"; default=true)
    betaN_target::Entry{T} = Entry{T}("-", "Target normalized beta (betaN) for FPP workflow"; default=3.6)
    f_Gr_target::Entry{T} = Entry{T}("-", "Target Greenwald fraction for FPP workflow"; default=1.3)
    diff_max_betaN::Entry{T} = Entry{T}("-", "Maximum allowable difference in betaN"; default=0.01)
    diff_max_density::Entry{T} = Entry{T}("-", "Maximum allowable difference in density"; default=0.01)
    coef_source_decrease::Entry{T} = Entry{T}("-", "Coefficient for source decrease adjustment"; default=1.0)
    coef_diff_betaN::Entry{T} = Entry{T}("-", "Coefficient for betaN adjustment"; default=1.0)
    num_iterations::Entry{Int} = Entry{Int}("-", "Number of main workflow iterations"; default=5)
    transport_num::Entry{Int} = Entry{Int}("-", "Maximum number of transport iterations per workflow"; default=5)
    json_num::Entry{Int} = Entry{Int}("-", "Starting iteration number"; default=0)
    
    # Default settings for actor parameters
    tglf_sat_rule::Entry{Symbol} = Entry{Symbol}("-", "Saturation rule for TGLF"; default=:sat0)
    tglf_lump_ions::Entry{Bool} = Entry{Bool}("-", "Lump ions for TGLF"; default=false)
    tglf_electromagnetic::Entry{Bool} = Entry{Bool}("-", "Enable electromagnetic effects for TGLF"; default=false)
    flux_matcher_evolve_pedestal::Entry{Bool} = Entry{Bool}("-", "Evolve pedestal for Flux Matcher"; default=false)
    flux_matcher_optimizer_algorithm::Entry{Symbol} = Entry{Symbol}("-", "Optimizer algorithm for Flux Matcher"; default=:simple)
    flux_matcher_max_iterations::Entry{Int} = Entry{Int}("-", "Maximum iterations for Flux Matcher"; default=100)
    flux_matcher_rho_transport::Entry{AbstractRange} = Entry{AbstractRange}("-", "Rho transport range for Flux Matcher"; default=0.2:0.1:0.8)
    flux_matcher_step_size::Entry{Float64} = Entry{Float64}("-", "Step size for Flux Matcher"; default=0.3)
    flux_matcher_relax::Entry{Float64} = Entry{Float64}("-", "Relaxation factor for Flux Matcher"; default=0.5)
    flux_matcher_evolve_densities::Entry{Symbol} = Entry{Symbol}("-", "Density evolution mode for Flux Matcher"; default=:flux_match)
    tglf_model::Entry{Symbol} = Entry{Symbol}("-", "TGLF model"; default=:TJLF)
    neoclassical_model::Entry{Symbol} = Entry{Symbol}("-", "Neoclassical model"; default=:hirshmansigmar)
    equilibrium_model::Entry{Symbol} = Entry{Symbol}("-", "Equilibrium model"; default=:TEQUILA)
    tequila_relax::Entry{Float64} = Entry{Float64}("-", "Relaxation factor for TEQUILA"; default=0.4)
    tequila_number_of_iterations::Entry{Int} = Entry{Int}("-", "Number of iterations for TEQUILA"; default=3000)
    
    #== display and debugging parameters ==#
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=true)
    verbose::Entry{Bool} = act_common_parameters(; verbose=true)
    save_path::Entry{String} = Entry{String}("-", "Directory path for saving iteration results"; default=joinpath(homedir(), "FPP"))
end

# Define the FPP actor
mutable struct ActorStationaryFPP{D, P} <: CompoundAbstractActor{D, P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorStationaryFPP{P}
    act::ParametersAllActors
    actor_hc::ActorHCD{D, P}
    actor_tr::ActorCoreTransport{D, P}
    actor_eq::ActorEquilibrium{D, P}
    actor_jt::ActorCurrent{D, P}
end

# Constructor for the FPP actor

"""
Constructs an `ActorStationaryFPP` with the specified `dd`, `par`, and `act`, and runs the workflow.

# Arguments
- `dd`: The plasma state (`IMAS.dd` object).
- `par`: Parameters for the FPP workflow (`FUSEparameters__ActorStationaryFPP`).
- `act`: Actor parameter configuration (`ParametersAllActors`).

# Keyword Arguments
Allows overriding of default `par` values via keyword arguments.

# Returns
A constructed `ActorStationaryFPP` instance.
"""

function ActorStationaryFPP(dd::IMAS.dd, par::FUSEparameters__ActorStationaryFPP, act::ParametersAllActors; kw...)
    # Input validation
    if isnothing(dd) || isnothing(act)
        error("`dd` and `act` must be provided and initialized.")
    end
    logging_actor_init(ActorStationaryFPP)
    # Initialize parameters with defaults and apply keyword overrides
    #par = FUSEparameters__ActorStationaryFPP(kw...)
    par = par(kw...)

    # Apply default or overridden values to actor-specific parameters
    act.ActorTGLF.sat_rule = par.tglf_sat_rule
    act.ActorTGLF.lump_ions = par.tglf_lump_ions
    act.ActorTGLF.electromagnetic = par.tglf_electromagnetic
    act.ActorFluxMatcher.evolve_pedestal = par.flux_matcher_evolve_pedestal
    act.ActorFluxMatcher.optimizer_algorithm = par.flux_matcher_optimizer_algorithm
    act.ActorFluxMatcher.max_iterations = par.flux_matcher_max_iterations
    act.ActorFluxMatcher.rho_transport = par.flux_matcher_rho_transport
    act.ActorFluxMatcher.step_size = par.flux_matcher_step_size
    act.ActorFluxMatcher.relax = par.flux_matcher_relax
    act.ActorFluxMatcher.evolve_densities = par.flux_matcher_evolve_densities
    act.ActorTGLF.model = par.tglf_model
    act.ActorNeoclassical.model = par.neoclassical_model
    act.ActorEquilibrium.model = par.equilibrium_model
    act.ActorTEQUILA.relax = par.tequila_relax
    act.ActorTEQUILA.number_of_iterations = par.tequila_number_of_iterations

    actor_hc = ActorHCD(dd, act.ActorHCD, act)
    actor_tr = ActorCoreTransport(dd, act.ActorCoreTransport, act)
    actor_eq = ActorEquilibrium(dd, act.ActorEquilibrium, act)
    actor_jt = ActorCurrent(dd, act.ActorCurrent, act)

    actor = ActorStationaryFPP(dd, par, act, actor_hc, actor_tr, actor_eq, actor_jt)
    step(actor)
    #finalize(actor)
    return actor
end

# Main step function for the FPP workflow
function _step(actor::ActorStationaryFPP)
    dd = actor.dd
    par = actor.par

    # Extract parameters from `par`
    betaN_target = par.betaN_target.value
    f_Gr_target = par.f_Gr_target.value
    diff_max_betaN = par.diff_max_betaN.value
    diff_max_density = par.diff_max_density.value
    coef_source_decrease = par.coef_source_decrease.value
    coef_diff_betaN = par.coef_diff_betaN.value
    num_iterations = par.num_iterations.value
    transport_num = par.transport_num.value
    fix_ip = par.fix_ip.value
    json_num = par.json_num  # Starting iteration number

    # Workflow logic
    for iteration in 1:num_iterations
        println("========== Start of Iteration $iteration ==========")

        num = 0
        diff_betaN = 1.0
        scale_power = 1.0
        scale_pellet = 1.0
        dd_beam = dd.pulse_schedule.ec.beam[1].power_launched.reference[1]
        dd_pellet = dd.pulse_schedule.pellet.launcher[1].frequency.reference[1]
        scale_power_record = Float64[]
        scale_pellet_record = Float64[]

        while abs(diff_betaN) > diff_max_betaN
            if num == transport_num
                println("The number of transport iterations exceeded the limit. Num = $num")
                break
            end

            num += 1
            dd0 = deepcopy(dd)
            FUSE.ActorHCD(dd, actor.act)
            FUSE.ActorFluxMatcher(dd, actor.act; verbose=par.verbose.value, do_plot=par.do_plot.value)

            betaN = dd.core_profiles.global_quantities.beta_tor_norm
            diff_betaN = -(betaN[1] - betaN_target) / betaN_target

            # Fusion power thresholds
            fusion_power_threshold_high = ifelse(dd0.summary.fusion.power.value[1] * 5 >= 5e8, 1.1, 1.3)
            fusion_power_threshold_low = ifelse(dd0.summary.fusion.power.value[1] * 5 >= 5e8, 0.9, 0.7)
            previous_fusion_power = dd0.summary.fusion.power.value[1] * 5
            current_fusion_power = dd.summary.fusion.power.value[1] * 5

            println("Fusion power comparison = ", current_fusion_power / previous_fusion_power)

            # Adjust power and pellet scaling
            if (current_fusion_power / previous_fusion_power) > fusion_power_threshold_high
                println("Fusion power exceeds high threshold.")
                reset_dd_if_needed!(dd, dd0, num, transport_num, diff_betaN, diff_max_betaN)
                adjust_scaling!(dd, dd0, num, dd_beam, dd_pellet, scale_power, scale_pellet, 0.9, 0.9, scale_power_record, scale_pellet_record)
            elseif (current_fusion_power / previous_fusion_power) < fusion_power_threshold_low
                println("Fusion power falls below low threshold.")
                reset_dd_if_needed!(dd, dd0, num, transport_num, diff_betaN, diff_max_betaN)
                adjust_scaling!(dd, dd0, num, dd_beam, dd_pellet, scale_power, scale_pellet, 1.1, 1.1, scale_power_record, scale_pellet_record)
            else
                println("Fusion power within acceptable range.")
                reset_dd_if_needed!(dd, dd0, num, transport_num, diff_betaN, diff_max_betaN)
                adjust_to_target!(dd, dd0, num, diff_betaN, dd_beam, dd_pellet, scale_power, scale_pellet, betaN_target, f_Gr_target, coef_diff_betaN, coef_source_decrease, diff_max_density)
            end
        end

        # Final adjustments and save results
        finalize_scaling!(dd, dd_beam, dd_pellet, scale_power, scale_pellet)
        run_equilibrium_and_current!(dd, actor.act, fix_ip)
        # Save iteration
        json_num += 1
        println("Saving results for Iteration $json_num")
        IMAS.imas2json(dd, joinpath(par.save_path, "Iteration_$(json_num)"); freeze=false)
        println("========== End of Iteration $iteration ==========\n")
        sleep(10)
    end
    # Update `json_num` in `par` to reflect the latest state
    par.json_num = json_num
end


# Helper functions
# 1. Reset the current plasma state to the previous state if transport iterations are within the limit and the betaN deviation exceeds the threshold
"""
    reset_dd_if_needed!(dd, dd0, num, transport_num, diff_betaN, diff_max_betaN)

Reset the current plasma state (`dd`) to the previous state (`dd0`) if the number of transport iterations is within the 
allowable limit and the betaN deviation exceeds the specified threshold.

# Arguments
- `dd`: The current plasma state (mutable `IMAS.dd` object).
- `dd0`: A deep copy of the previous plasma state for fallback.
- `num`: The current number of transport iterations.
- `transport_num`: Maximum allowable number of transport iterations.
- `diff_betaN`: The normalized difference between the current betaN and the target betaN.
- `diff_max_betaN`: Maximum allowable normalized betaN difference.

"""
function reset_dd_if_needed!(dd, dd0, num, transport_num, diff_betaN, diff_max_betaN)
    if num < transport_num && abs(diff_betaN) > diff_max_betaN
        println("Resetting `dd` to previous state due to large betaN difference.")
        dd = deepcopy(dd0)
    end
end
# 2. Adjust power and pellet scaling for fusion power thresholds
"""
    adjust_scaling!(dd, dd0, num, dd_beam, dd_pellet, scale_power, scale_pellet, power_coef, pellet_coef, scale_power_record, scale_pellet_record)

Adjust the power and pellet scaling factors based on the provided coefficients, and log the scaling adjustments.

# Arguments
- `dd`: The current plasma state (mutable `IMAS.dd` object).
- `dd0`: A deep copy of the previous plasma state.
- `num`: Current transport iteration number.
- `dd_beam`: Initial beam power reference.
- `dd_pellet`: Initial pellet frequency reference.
- `scale_power`: Current power scaling factor.
- `scale_pellet`: Current pellet scaling factor.
- `power_coef`: Coefficient for adjusting power scaling.
- `pellet_coef`: Coefficient for adjusting pellet scaling.
- `scale_power_record`: Record of power scaling adjustments.
- `scale_pellet_record`: Record of pellet frequency scaling adjustments.

"""
function adjust_scaling!(dd, dd0, num, dd_beam, dd_pellet, scale_power, scale_pellet, power_coef, pellet_coef, scale_power_record, scale_pellet_record)
    scale_power *= power_coef
    scale_pellet *= pellet_coef
    dd.pulse_schedule.ec.beam[1].power_launched.reference .= dd_beam * scale_power
    dd.pulse_schedule.ec.beam[2].power_launched.reference .= dd_beam * scale_power
    dd.pulse_schedule.pellet.launcher[1].frequency.reference .= dd_pellet * scale_pellet
    push!(scale_power_record, power_coef)
    push!(scale_pellet_record, pellet_coef)
    println("Number of Transport Iteration = ",num)
    println("scale_power = ",scale_power_record)
    println("scale_pellet = ",scale_pellet_record)
    sleep(10)
end

# 3. Adjust power and density to meet targets
"""
    adjust_to_target!(dd, dd0, num, diff_betaN, dd_beam, dd_pellet, scale_power, scale_pellet, betaN_target, f_Gr_target, coef_diff_betaN, coef_source_decrease, diff_max_density)

Adjust power and density scaling to meet the target betaN and Greenwald fraction. Logs the adjustments and intermediate values.

# Arguments
- `dd`: The current plasma state (mutable `IMAS.dd` object).
- `dd0`: A deep copy of the previous plasma state.
- `num`: Current transport iteration number.
- `diff_betaN`: The normalized difference between current and target betaN.
- `dd_beam`: Initial beam power reference.
- `dd_pellet`: Initial pellet frequency reference.
- `scale_power`: Current power scaling factor.
- `scale_pellet`: Current pellet scaling factor.
- `betaN_target`: Target normalized beta.
- `f_Gr_target`: Target Greenwald fraction.
- `coef_diff_betaN`: Coefficient for betaN adjustment.
- `coef_source_decrease`: Coefficient for density source adjustment.
- `diff_max_density`: Maximum allowable normalized density difference.

"""
function adjust_to_target!(dd, dd0, num, diff_betaN, dd_beam, dd_pellet, scale_power, scale_pellet, betaN_target, f_Gr_target, coef_diff_betaN, coef_source_decrease, diff_max_density)
    Paux_scale = 1.0 + coef_diff_betaN * diff_betaN
    scale_power *= Paux_scale
    f_Gr = IMAS.greenwald_fraction(dd)
    density_ratio = f_Gr / f_Gr_target
    ext_particle_coef = 1.0 - (density_ratio - 1.0) * coef_source_decrease
    scale_pellet *= ext_particle_coef
    dd.pulse_schedule.ec.beam[1].power_launched.reference .= dd_beam * scale_power
    dd.pulse_schedule.ec.beam[2].power_launched.reference .= dd_beam * scale_power
    if abs(density_ratio - 1.0) > diff_max_density
        dd.pulse_schedule.pellet.launcher[1].frequency.reference .= dd_pellet * scale_pellet
    end
    println("Number of Transport Iteration = ",num)
    println("betaN = ", dd.core_profiles.global_quantities.beta_tor_norm)
    println("diff_betaN = ",diff_betaN)
    println("Paux_scale = ",Paux_scale)
    println("f_Gr = ", f_Gr)
    println("ext_particle_coef = ",ext_particle_coef)
    push!(scale_power_record,Paux_scale)
    push!(scale_pellet_record,ext_particle_coef)
    println("scale_power = ",scale_power_record)
    println("scale_pellet = ",scale_pellet_record)
    sleep(10)
end
# 4. Finalize power and pellet scaling
"""
    finalize_scaling!(dd, dd_beam, dd_pellet, scale_power, scale_pellet)

Clamp the final power and pellet scaling factors to a predefined range and update the plasma state. Logs the final scaling adjustments.

# Arguments
- `dd`: The current plasma state (mutable `IMAS.dd` object).
- `dd_beam`: Initial beam power reference.
- `dd_pellet`: Initial pellet frequency reference.
- `scale_power`: Final power scaling factor.
- `scale_pellet`: Final pellet scaling factor.

"""
function finalize_scaling!(dd, dd_beam, dd_pellet, scale_power, scale_pellet)
    scale_power = clamp(scale_power, 0.95, 1.05)
    scale_pellet = clamp(scale_pellet, 0.95, 1.05)
    dd.pulse_schedule.ec.beam[1].power_launched.reference .= dd_beam * scale_power
    dd.pulse_schedule.ec.beam[2].power_launched.reference .= dd_beam * scale_power
    dd.pulse_schedule.pellet.launcher[1].frequency.reference .= dd_pellet * scale_pellet
    FUSE.ActorHCD(dd, actor.act)
    println("scale_power = ",scale_power)
    println("Paux_power[1] = ",dd.pulse_schedule.ec.beam[1].power_launched.reference[1])
    println("Paux_power[2] = ",dd.pulse_schedule.ec.beam[2].power_launched.reference[1])
    println("scale_pellet = ",scale_pellet)
    println("pellet_launcher.frequency = ",dd.pulse_schedule.pellet.launcher[1].frequency.reference[1])
end

# 5. Run equilibrium and current evolution
"""
    run_equilibrium_and_current!(dd, act, fix_ip)

Runs the equilibrium and current evolution based on the `fix_ip` parameter.

# Arguments
- `dd`: The current plasma state (mutable `IMAS.dd` object).
- `act`: The `ParametersAllActors` object holding actor-specific parameters.
- `fix_ip`: Boolean indicating whether to run fix Ip option or adjust Ip option.
"""
function run_equilibrium_and_current!(dd, act, fix_ip::Bool)
    if fix_ip
        println("Running Fix Ip mode. This mode is typically used in running the Standard H-mode scenario or during the initial phase of the High-betaP scenario, when the plasma current (Ip) is well-defined or reasonably estimated.")
        for i in 1:3
            FUSE.ActorEquilibrium(dd, act; model=:TEQUILA, ip_from=:pulse_schedule)
            FUSE.ActorCurrent(dd, act; ip_from=:pulse_schedule)
        end
        FUSE.ActorEquilibrium(dd, act; model=:TEQUILA, ip_from=:pulse_schedule, do_plot=true)
        println("betaN = ", dd.core_profiles.global_quantities.beta_tor_norm)
        println("Ip = ", dd.core_profiles.global_quantities.ip)
    else
        println("Running Adjust Ip mode. This mode is typically used in the High-betaP scenario when the Ip doesn't change sgnificantly.")
        oh_com = 0.03
        foh = 1.0
        target_foh = 0.2
        relax = 0.2
        count = 0
        max_iterations = 3

        while abs(foh - target_foh) > oh_com && count < max_iterations
            FUSE.ActorEquilibrium(dd, act; ip_from=:pulse_schedule)
            FUSE.ActorCurrent(dd, act; ip_from=:pulse_schedule)
            curoh = dd.core_profiles.global_quantities.ip[1] - dd.core_profiles.global_quantities.current_non_inductive[1]
            foh = curoh / dd.core_profiles.global_quantities.ip[1]
            println("foh = ", foh)

            if abs(foh - target_foh) > oh_com
                dd.pulse_schedule.flux_control.i_plasma.reference *= (1.0 - relax * (foh - target_foh))
                FUSE.ActorEquilibrium(dd, act; ip_from=:pulse_schedule)
                println("foh - target_foh = ", foh - target_foh)
                println("i_plasma = ", dd.pulse_schedule.flux_control.i_plasma.reference[1])
            end
            count += 1
        end

        FUSE.ActorEquilibrium(dd, act; ip_from=:pulse_schedule)
        println("betaN = ", dd.core_profiles.global_quantities.beta_tor_norm)
        println("Ip = ", dd.core_profiles.global_quantities.ip)
    end
end

