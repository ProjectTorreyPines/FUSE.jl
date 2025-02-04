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
    diff_max_betaN::Entry{T} = Entry{T}("-", "Maximum allowable difference in betaN"; default=0.05)
    diff_max_density::Entry{T} = Entry{T}("-", "Maximum allowable difference in density"; default=0.05)
    coef_source_decrease::Entry{T} = Entry{T}("-", "Coefficient for source decrease adjustment"; default=1.0)
    coef_diff_betaN::Entry{T} = Entry{T}("-", "Coefficient for betaN adjustment"; default=1.0)
    num_iterations::Entry{Int} = Entry{Int}("-", "Number of main workflow iterations"; default=5)
    transport_num::Entry{Int} = Entry{Int}("-", "Maximum number of transport iterations per iteration"; default=5)
    json_num::Entry{Int} = Entry{Int}("-", "Starting iteration number"; default=0)
    
    # Default settings for actor parameters
    tglf_sat_rule::Entry{Symbol} = Entry{Symbol}("-", "Saturation rule for TGLF"; default=:sat0)
    tglf_lump_ions::Entry{Bool} = Entry{Bool}("-", "Lump ions for TGLF"; default=false)
    tglf_electromagnetic::Entry{Bool} = Entry{Bool}("-", "Enable electromagnetic effects for TGLF"; default=false)
    flux_matcher_evolve_pedestal::Entry{Bool} = Entry{Bool}("-", "Evolve pedestal for Flux Matcher"; default=false)
    flux_matcher_optimizer_algorithm ::Entry{Symbol} = Entry{Symbol}("-", "Optimizer algorithm for Flux Matcher"; default=:simple)
    flux_matcher_max_iterations::Entry{Int} = Entry{Int}("-", "Maximum iterations for Flux Matcher"; default=100)
    flux_matcher_rho_transport::Entry{AbstractRange} = Entry{AbstractRange}("-", "Rho transport range for Flux Matcher"; default=0.2:0.1:0.8)
    flux_matcher_step_size::Entry{Float64} = Entry{Float64}("-", "Step size for Flux Matcher"; default=0.2)
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
    do_plot_workflow::Entry{Bool} = Entry{Bool}("-", "Enable or disable workflow-level plotting"; default=true)
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

function ActorStationaryFPP(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorStationaryFPP
    if isnothing(dd) || isnothing(act)
        error("`dd` and `act` must be provided and initialized.")
    end
    #logging_actor_init(ActorStationaryFPP)
    par = par(kw...)
    ensure_save_path(par)

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

    actor_hc = FUSE.ActorHCD(dd, act.ActorHCD, act)
    actor_eq = FUSE.ActorEquilibrium(dd, act.ActorEquilibrium, act; ip_from=:core_profiles)
    actor_tr = FUSE.ActorCoreTransport(dd, act.ActorCoreTransport, act)
    actor_jt = FUSE.ActorCurrent(dd, act.ActorCurrent, act; ip_from=:pulse_schedule)

    actor = ActorStationaryFPP(dd, par, act, actor_hc, actor_tr, actor_eq, actor_jt)

    # Run the workflow with logging
    log_file_path = joinpath(par.save_path, "workflow.log")
    #step(actor)
    log_to_file(log_file_path) do
        step(actor)
    end
    return actor
end


# Main step function for the FPP workflow
function _step(actor::ActorStationaryFPP)
    dd = actor.dd
    par = actor.par

    # Extract parameters from `par`
    betaN_target = par.betaN_target
    f_Gr_target = par.f_Gr_target
    diff_max_betaN = par.diff_max_betaN
    diff_max_density = par.diff_max_density
    coef_source_decrease = par.coef_source_decrease
    coef_diff_betaN = par.coef_diff_betaN
    num_iterations = par.num_iterations
    transport_num = par.transport_num
    fix_ip = par.fix_ip
    json_num = par.json_num  # Starting iteration number
    convergence_threshold = par.convergence_error  # Threshold for convergence
    n_recent = 5  # Number of recent iterations to check for convergence

    # Initialize variables for convergence tracking
    fusion_power_history = Float64[]  # Track fusion power
    betaN_history = Float64[]  # Track betaN values

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
        dd0 = deepcopy(dd)
        while abs(diff_betaN) > diff_max_betaN
            if num == transport_num
                println("The number of transport iterations exceeded the limit. Num = $num")
                break
            end
            num += 1
            FUSE.ActorHCD(dd, actor.act)
            FUSE.ActorFluxMatcher(dd, actor.act; verbose=par.verbose, do_plot=par.do_plot)

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
                println("betaN = ", dd.core_profiles.global_quantities.beta_tor_norm)
                if num < transport_num && abs(diff_betaN) > diff_max_betaN
                    println("Resetting `dd` to previous state due to large betaN difference.")
                    dd = deepcopy(dd0)
                end
                scale_power, scale_pellet = adjust_scaling!(dd, num, dd_beam, dd_pellet, scale_power, scale_pellet, 0.8, 0.8, scale_power_record, scale_pellet_record)
            elseif (current_fusion_power / previous_fusion_power) < fusion_power_threshold_low
                println("Fusion power falls below low threshold.")
                println("betaN = ", dd.core_profiles.global_quantities.beta_tor_norm)
                if num < transport_num && abs(diff_betaN) > diff_max_betaN
                    println("Resetting `dd` to previous state due to large betaN difference.")
                    dd = deepcopy(dd0)
                end    
                scale_power, scale_pellet = adjust_scaling!(dd, num, dd_beam, dd_pellet, scale_power, scale_pellet, 1.2, 1.2, scale_power_record, scale_pellet_record)
            else
                println("Fusion power within acceptable range.")
                println("betaN = ", dd.core_profiles.global_quantities.beta_tor_norm)
                if num < transport_num && abs(diff_betaN) > diff_max_betaN
                    println("Resetting `dd` to previous state due to large betaN difference.")
                    dd = deepcopy(dd0)
                end   
                scale_power, scale_pellet = adjust_to_target!(dd, num, diff_betaN, dd_beam, dd_pellet, scale_power, scale_pellet, betaN_target, f_Gr_target, coef_diff_betaN, coef_source_decrease, diff_max_density, scale_power_record, scale_pellet_record)
            end
        end

        # Final adjustments and save results
        finalize_scaling!(dd, actor.act, dd_beam, dd_pellet, scale_power, scale_pellet)
        run_equilibrium_and_current!(dd, actor.act, fix_ip)

        # Save iteration
        json_num += 1
        file_path = joinpath(par.save_path, "Iteration_$(json_num)")
        println("Saving to path: ", file_path)
        IMAS.imas2json(dd, file_path; freeze=false)
        
        # Append convergence tracking values
        push!(fusion_power_history, dd.summary.fusion.power.value[1] * 5)
        push!(betaN_history, dd.core_profiles.global_quantities.beta_tor_norm[1])
        println("Fusion power history : ", fusion_power_history)
        println("betaN history : ", betaN_history)
        # Convergence check
        if length(fusion_power_history) >= n_recent
            #using Statistics
            avg_power, std_dev_power = calculate_standard_deviation(fusion_power_history, n_recent)
            avg_betaN, std_dev_betaN = calculate_standard_deviation(betaN_history, n_recent)

            println("Average fusion power over last $n_recent iterations: ", avg_power)
            println("Standard deviation of fusion power: ", std_dev_power)
            println("Average betaN over last $n_recent iterations: ", avg_betaN)
            println("Standard deviation of betaN: ", std_dev_betaN)

            if (std_dev_power / avg_power < convergence_threshold) && (std_dev_betaN / avg_betaN < convergence_threshold)
                println("The system has stabilized. Ending workflow.")
                println("The convergences are : ",std_dev_power / avg_power, std_dev_betaN / avg_betaN )
                par.json_num = json_num  # Update parameter to reflect latest iteration
                dd_set = read_json_files_to_dd_set(par.save_path)
                values_dict = read_values_from_dd_set(dd_set, [ :Paux_tot, :βn, :βn_MHD, :ip, :fGW, :Pfusion, :Qfusion, :ip_bs_aux_ohm, :ip_ohm, :ip_bs, :Te0, :Ti0, :ne0])
                plot_results(values_dict; do_plot=true)
                return actor
            else
                println("The system is still changing.")
                println("The convergences are : ",std_dev_power / avg_power, std_dev_betaN / avg_betaN )
            end
        end
        println("========== End of Iteration $iteration ==========\n")
    end
    # Update `json_num` in `par` to reflect the latest state
    par.json_num = json_num

    # Read and plot results
    if par.do_plot_workflow
        #println("Generating workflow-level plots...")
        dd_set = read_json_files_to_dd_set(par.save_path)
        values_dict = read_values_from_dd_set(dd_set, [ :Paux_tot, :βn, :βn_MHD, :ip, :fGW, :Pfusion, :Qfusion, :ip_bs_aux_ohm, :ip_ohm, :ip_bs, :Te0, :Ti0, :ne0])
        plot_results(values_dict; do_plot=true)
    else
        println("Workflow-level plotting is disabled.")
    end
    return actor
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
    return dd
end
# 2. Adjust power and pellet scaling for fusion power thresholds
"""
    adjust_scaling!(dd, num, dd_beam, dd_pellet, scale_power, scale_pellet, power_coef, pellet_coef, scale_power_record, scale_pellet_record)

Adjust the power and pellet scaling factors based on the provided coefficients, and log the scaling adjustments.

# Arguments
- `dd`: The current plasma state (mutable `IMAS.dd` object).
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
function adjust_scaling!(dd, num, dd_beam, dd_pellet, scale_power, scale_pellet, power_coef, pellet_coef, scale_power_record, scale_pellet_record)
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
    return scale_power, scale_pellet
end

# 3. Adjust power and density to meet targets
"""
    adjust_to_target!(dd, num, diff_betaN, dd_beam, dd_pellet, scale_power, scale_pellet, betaN_target, f_Gr_target, coef_diff_betaN, coef_source_decrease, diff_max_density)

Adjust power and density scaling to meet the target betaN and Greenwald fraction. Logs the adjustments and intermediate values.

# Arguments
- `dd`: The current plasma state (mutable `IMAS.dd` object).
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
function adjust_to_target!(dd, num, diff_betaN, dd_beam, dd_pellet, scale_power, scale_pellet, betaN_target, f_Gr_target, coef_diff_betaN, coef_source_decrease, diff_max_density, scale_power_record, scale_pellet_record)
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
    println("diff_betaN = ",diff_betaN)
    println("Paux_scale = ",Paux_scale)
    println("f_Gr = ", f_Gr)
    println("ext_particle_coef = ",ext_particle_coef)
    push!(scale_power_record,Paux_scale)
    push!(scale_pellet_record,ext_particle_coef)
    println("scale_power = ",scale_power_record)
    println("scale_pellet = ",scale_pellet_record)
    return scale_power, scale_pellet
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
function finalize_scaling!(dd, act, dd_beam, dd_pellet, scale_power, scale_pellet)
    scale_power = clamp(scale_power, 0.9, 1.1)
    scale_pellet = clamp(scale_pellet, 0.9, 1.1)
    dd.pulse_schedule.ec.beam[1].power_launched.reference .= dd_beam * scale_power
    dd.pulse_schedule.ec.beam[2].power_launched.reference .= dd_beam * scale_power
    dd.pulse_schedule.pellet.launcher[1].frequency.reference .= dd_pellet * scale_pellet
    FUSE.ActorHCD(dd, act)
    f_Gr = IMAS.greenwald_fraction(dd)
    println("betaN = ", dd.core_profiles.global_quantities.beta_tor_norm)
    println("f_Gr = ", f_Gr)
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

#6.plot_results(values_dict)
"""
Plots key plasma simulation parameters over the course of iterations using the extracted values.
"""
function plot_results(values_dict; do_plot = true)
    if !do_plot
        println("Skipping plots as `do_plot` is set to false.")
        return
    end
    # Extract individual arrays from the dictionary
    Paux_tot_values = values_dict[:Paux_tot]
    βn_values = values_dict[:βn]
    βnMHD_values = values_dict[:βn_MHD]
    ip_values = values_dict[:ip]
    fGW_values = values_dict[:fGW]
    Pfusion_values = values_dict[:Pfusion]
    Qfusion_values = values_dict[:Qfusion]
    ip_bs_aux_ohm_values = values_dict[:ip_bs_aux_ohm]
    ip_ohm_values = values_dict[:ip_ohm]
    ip_bs_values = values_dict[:ip_bs]
    fbs_values = ip_bs_values ./ ip_bs_aux_ohm_values
    foh_values = ip_ohm_values ./ ip_bs_aux_ohm_values
    Te0_values = values_dict[:Te0]
    Ti0_values = values_dict[:Ti0]
    ne0_values = values_dict[:ne0] ./ 1.e20

    # Plotting
    p = plot(layout=(2, 4), size=(1200, 400), left_margin=5Plots.mm, bottom_margin=5Plots.mm)

    # βn and βn_MHD
    plot!(p, 1:length(βn_values), βn_values, subplot=1, title="βn and βn_MHD", ylabel="βn", legend=false, ylims=(0, 6), label="βn")
    plot!(p, 1:length(βnMHD_values), βnMHD_values, subplot=1, color=:red, label="βn_MHD")

    # Pfusion and Paux
    plot!(p, 1:length(Pfusion_values), Pfusion_values, subplot=2, title="Pfusion and Paux", ylabel="Power (MW)", legend=:topleft, ylims=(0, 800), color=:blue, label="Pfusion")
    plot!(p, 1:length(Paux_tot_values), Paux_tot_values, subplot=2, color=:green, label="Paux")

    # Plasma current (Ip)
    plot!(p, 1:length(ip_values), ip_values, subplot=3, title="Plasma Current (Ip)", ylabel="Ip (MA)", legend=false, ylims=(0, 10))

    # Greenwald Fraction (fGW)
    plot!(p, 1:length(fGW_values), fGW_values, subplot=4, title="Greenwald Fraction", ylabel="fGW", legend=false, ylims=(0, 2))

    # Electron and Ion Temperatures (Te0 and Ti0)
    plot!(p, 1:length(Te0_values), Te0_values, subplot=5, title="Core Temperatures (Te0, Ti0)", ylabel="Temperature (keV)", legend=:topleft, color=:blue, label="Te0")
    plot!(p, 1:length(Ti0_values), Ti0_values, subplot=5, color=:red, label="Ti0")

    # Electron density (ne0)
    plot!(p, 1:length(ne0_values), ne0_values, subplot=6, title="Core Electron Density", ylabel="ne0 (10^20 m^-3)", legend=false, ylims=(1.0, 2.5))

    # Bootstrap and Ohmic Fractions (fbs and foh)
    plot!(p, 1:length(fbs_values), fbs_values, subplot=7, title="Bootstrap and Ohmic Fractions", ylabel="Fraction", legend=:topleft, color=:blue, label="fbs")
    plot!(p, 1:length(foh_values), foh_values, subplot=7, color=:red, label="foh")

    # Fusion gain (Qfusion)
    plot!(p, 1:length(Qfusion_values), Qfusion_values, subplot=8, title="Fusion Gain (Qfusion)", ylabel="Q", legend=false, ylims=(0, 20))

    display(p)
end

# 7. calculate_standard_deviation(values, n)
"""
Calculates the average and standard deviation of the last `n` values in the given vector.

# Arguments
- `values`: A vector of numerical values.
- `n`: The number of recent values to include in the calculation.

# Returns
A tuple `(avg, std_dev)` where:
- `avg` is the average of the last `n` values.
- `std_dev` is the standard deviation of the last `n` values.
"""
function calculate_standard_deviation(values::Vector{<:Real}, n::Int)
    # Ensure there are enough iterations
    if length(values) < n
        error("Not enough iterations to calculate standard deviation.")
    end

    # Convert the values to `Float64` explicitly to handle mixed types
    recent_values = Float64.(values[end-n+1:end])

    # Calculate the average and standard deviation
    avg = sum(recent_values) / length(recent_values)
    std_dev = sqrt(sum((x - avg)^2 for x in recent_values) / (length(recent_values) - 1))

    return avg, std_dev
end

# 8. ensure the save path exists 
function ensure_save_path(par::FUSEparameters__ActorStationaryFPP)
    if !isdir(par.save_path)
        mkpath(par.save_path)
        println("Created directory: ", par.save_path)
    end
end

# 9. Function to read values from dd_set
function read_values_from_dd_set(dd_set, keys_to_extract)
    values_dict = Dict{Symbol, Vector{Any}}()
    
    # Initialize the dictionary with empty arrays for each key
    for key in keys_to_extract
        values_dict[key] = []
    end

    # Extract the values for each key
    for dd in dd_set
        ex = IMAS.extract(dd)
        for key in keys_to_extract
            key_index = findfirst(x -> x == key, ex.keys)
            if key_index !== nothing
                key_val = ex.vals[key_index]
                if hasfield(typeof(key_val), :value)
                    push!(values_dict[key], getfield(key_val, :value))
                else
                    push!(values_dict[key], key_val)
                end
            else
                push!(values_dict[key], nothing)  # Add nothing if the key doesn't exist
            end
        end
    end
    return values_dict
end

# 10. Function to read JSON files and convert them to IMAS objects
function read_json_files_to_dd_set(folder_path)
    json_files = filter(x -> occursin("Iteration_", x), readdir(folder_path, join=true))
    
    # Sort files numerically based on the iteration number
    sorted_files = sort(json_files, by = x -> parse(Int, match(r"Iteration_(\d+)", basename(x)).captures[1]))
    
    # Convert JSON files to IMAS objects
    dd_set = [IMAS.json2imas(file) for file in sorted_files]
    return dd_set
end

# 11. Function to redirect 'stdout' to the log file
function log_to_file(func::Function, file_path::String)
    # Ensure the directory for the log file exists
    dir_path = dirname(file_path)
    if !isdir(dir_path)
        mkpath(dir_path)
        println("Created directory: $dir_path")
    end

    # Open the log file and redirect `stdout`
    open(file_path, "w") do log_file
        original_stdout = Base.stdout  # Save the original `stdout`
        try
            redirect_stdout(log_file)  # Redirect `stdout` to the log file
            func()  # Execute the provided function
        finally
            redirect_stdout(original_stdout)  # Restore `stdout`
        end
    end
end
