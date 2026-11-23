#= ========= =#
#  ActorWPED  #
#= ========= =#
@actor_parameters_struct ActorWPED{T} begin
    #== common pedestal parameters==#
    rho_nml::Entry{T} = Entry{T}("-", "Defines rho at which the no man's land region starts")
    rho_ped::Entry{T} = Entry{T}("-", "Defines rho at which the pedestal region starts") # rho_nml < rho_ped
    T_ratio_pedestal::Entry{T} =
        Entry{T}("-", "Ratio of ion to electron temperatures (or rho at which to sample for that ratio, if negative; or rho_nml-(rho_ped-rho_nml) if 0.0)"; default=1.0)
    Te_sep::Entry{T} = Entry{T}("eV", "Separatrix electron temperature"; default=80.0, check=x -> @assert x > 0 "Te_sep must be > 0")
    Ti_sep::Entry{T} = Entry{T}("eV", "Separatrix ion temperature (if not set, Ti_sep = Te_sep * Ti/Te from `T_ratio_pedestal`)"; check=x -> @assert x > 0 "Ti_sep must be > 0")
    Ti_ped_model::Switch{Symbol} = Switch{Symbol}(
        [:ratio, :increment, :heat_flux, :fingerprint],
        "-",
        "Split of the pedestal pressure between Te_ped and Ti_ped: `:ratio` Ti_ped = r Te_ped; `:increment` Ti_ped-Ti_sep = r (Te_ped-Te_sep) (r from `T_ratio_pedestal`); " *
        "`:heat_flux` pedestal power balance with χi/χe = `chi_i_over_chi_e`; `:fingerprint` pedestal power balance with χi = χi_neo + f χe (f from `fingerprint_mode`)";
        default=:ratio)
    fingerprint_mode::Switch{Symbol} = Switch{Symbol}(
        [:MHD, :MTM, :ETG, :ITG_TEM, :custom],
        "-",
        "Dominant pedestal transport mode for `Ti_ped_model=:fingerprint`, setting the turbulent χi/χe (MHD-like=1, MTM=ETG=0.1, ITG/TEM=2, `:custom`=`chi_i_over_chi_e`)";
        default=:MHD)
    chi_i_over_chi_e::Entry{T} = Entry{T}(
        "-",
        "Pedestal χi/χe for `Ti_ped_model=:heat_flux`, and turbulent χi/χe for `Ti_ped_model=:fingerprint` with `fingerprint_mode=:custom`";
        default=1.0,
        check=x -> @assert x >= 0.0 "chi_i_over_chi_e must be >= 0")
    chi_i_neo_model::Switch{Symbol} = Switch{Symbol}([:changhinton, :none], "-", "Neoclassical ion χ added to the turbulent one for `Ti_ped_model=:fingerprint`"; default=:changhinton)
    Ti_ped_self_consistent::Entry{Bool} = Entry{Bool}(
        "-",
        "For `Ti_ped_model=:heat_flux|:fingerprint` with EPED-like models: solve (bracketed bisection) for the Ti_ped that is consistent with the collisional exchange source and χi_neo recomputed on the resulting profiles. Note that the core profiles are shifted with the pedestal and do not otherwise respond";
        default=false)
    ip_from::Switch{Symbol} = switch_get_from(:ip)
    βn_from::Switch{Symbol} = switch_get_from(:βn)
    ne_from::Switch{Symbol} = switch_get_from(:ne_ped)
    zeff_from::Switch{Symbol} = switch_get_from(:zeff_ped)
    #== actor parameters ==#
    ped_to_core_fraction::Entry{T} = Entry{T}("-", "Ratio of edge (@rho=0.9) to core stored energy [0.05 for L-mode, 0.3 for neg-T plasmas, missing keeps original ratio]")
    #== display and debugging parameters ==#
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
end

mutable struct OptimizationWPED
    α_Te::Float64 # exponent for Te that is used to set inverse scale-length at the boundary
    α_Ti::Float64 # same for Ti
    value_bound::Float64 # value of Te at the boundary
end

mutable struct ActorWPED{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.DD{D}
    par::OverrideParameters{P,FUSEparameters__ActorWPED{P}}
    optimization_guesses::OptimizationWPED
    ti_split::Union{Nothing,NamedTuple} # pedestal Te/Ti split (see `Ti_ped_model`)
end

"""
    ActorWPED(dd::IMAS.DD, act::ParametersAllActors; kw...)

Adjusts pedestal temperature profiles to match a target energy ratio between pedestal and core regions.

The actor optimizes temperature profiles at the plasma edge by finding the optimal boundary temperature 
and slope parameters that produce a specified pedestal-to-core stored energy fraction. It uses optimization 
to iteratively adjust the temperature gradient and boundary conditions until the energy balance criterion is satisfied.

Key operations:
- Calculates current core-to-edge energy ratio if target ratio is not specified
- Optimizes boundary temperature value to match target energy ratio
- Adjusts both electron and ion temperature profiles with consistent Ti/Te ratios
- Preserves separatrix temperature boundary conditions
- Updates temperature profiles through H-mode blending functions
"""
function ActorWPED(dd::IMAS.DD, act::ParametersAllActors; kw...)
    actor = ActorWPED(dd, act.ActorWPED; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorWPED(dd::IMAS.DD, par::FUSEparameters__ActorWPED; kw...)
    logging_actor_init(ActorWPED)
    par = OverrideParameters(par; kw...)
    return ActorWPED(dd, par, OptimizationWPED(0.0, 0.0, 3e3), nothing)
end

"""
    _step(actor::ActorWPED)

Executes the optimization procedure to find optimal pedestal temperature profiles.

The step function performs a nested optimization:
1. Outer loop: optimizes the boundary temperature value to match the target energy ratio
2. Inner loop: for each boundary value, optimizes the temperature gradient parameters (α_Te, α_Ti)
3. Uses Brent's method for robust 1D optimization with relative tolerance of 1E-3
4. Updates both electron and ion temperature profiles maintaining Ti/Te ratio consistency
"""
function _step(actor::ActorWPED{D,P}) where {D<:Real,P<:Real}
    dd = actor.dd
    par = actor.par

    rho09 = 0.9
    cp1d = dd.core_profiles.profiles_1d[]

    if par.do_plot
        ppe = plot(cp1d.electrons, :temperature; label="Te before WPED blending")
        ppi = plot(cp1d, :t_i_average; label="Ti before WPED blending", xlabel="rho")
    end

    if ismissing(par, :ped_to_core_fraction)
        core, edge = IMAS.core_edge_energy(cp1d, rho09)
        ped_to_core_fraction = edge / core
    else
        ped_to_core_fraction = par.ped_to_core_fraction
    end

    Ti_over_Te = ti_te_ratio(cp1d, par.T_ratio_pedestal, par.rho_nml, par.rho_ped)

    # Change the last point of the temperatures profiles since that's what used in cost_WPED_ztarget_pedratio()
    cp1d.electrons.temperature[end] = par.Te_sep
    for ion in cp1d.ion
        if !ismissing(ion, :temperature)
            ion.temperature[end] = separatrix_ti(par, Ti_over_Te)
        end
    end

    if par.Ti_ped_self_consistent && par.Ti_ped_model ∈ (:heat_flux, :fingerprint)
        @warn "ActorWPED: Ti_ped_self_consistent=true is only supported by EPED-like pedestal models (EPED, analytic): collisional exchange and χi_neo are taken from the profiles before the pedestal update" maxlog = 1
    end

    # how Ti at the pedestal top follows from Te at the pedestal top
    split = pedestal_ti_split(dd, par, par.rho_ped, Ti_over_Te)

    Te_orig = cp1d.electrons.temperature
    Ti_orig = cp1d.t_i_average
    res_value_bound = Optim.optimize(
        value_bound -> IMAS.force_float64(cost_WPED_ztarget_pedratio(cp1d, Te_orig, Ti_orig, value_bound, ped_to_core_fraction, par.rho_ped, split)),
        1.0,
        IMAS.force_float64(cp1d.electrons.temperature[1]),
        Optim.Brent();
        rel_tol=IMAS.pedestal_search_tolerance.brent[])

    _, α_Te, α_Ti = wped_apply!(cp1d, res_value_bound.minimizer, ped_to_core_fraction, par.rho_ped, split)
    Te_ped = res_value_bound.minimizer
    actor.ti_split = pedestal_ti_split_result(split, merge((; Te_ped), pedestal_ti(split, Te_ped)))

    if D <: ForwardDiff.Dual
        # The searches above ran on value parts. The boundary value and the two exponents solve three
        # conditions together (inverse scale lengths of Te and Ti continuous at the pedestal boundary,
        # pedestal to core energy ratio); their derivatives follow from the implicit function theorem.
        # (ordered so that each condition is the one that determines the unknown at the same position: when a search
        #  ends without meeting its condition, typically an exponent on a plateau of the mismatch, that unknown is held fixed)
        ustar = Float64[α_Te, α_Ti, res_value_bound.minimizer]
        G = u -> wped_residuals!(cp1d, Te_orig, Ti_orig, u[3], u[1], u[2], ped_to_core_fraction, par.rho_ped, split)
        u = IMAS.implicit_solution(G, ustar, D; atol=1E-2)
        wped_residuals!(cp1d, Te_orig, Ti_orig, u[3], u[1], u[2], ped_to_core_fraction, par.rho_ped, split)
    end

    if par.do_plot
        display(plot!(ppe, cp1d.electrons, :temperature; label="Te after WPED blending"))
        display(plot!(ppi, cp1d, :t_i_average; label="Ti after WPED blending"))
    end

    return actor
end

function cost_WPED_ztarget_pedratio(
    cp1d::IMAS.core_profiles__profiles_1d{T},
    Te_orig::Vector{T},
    Ti_orig::Vector{T},
    value_bound::Real,
    ped_to_core_fraction::Real,
    rho_ped::Real,
    split::NamedTuple) where {T<:Real}

    # restore original profiles
    cp1d.electrons.temperature = deepcopy(Te_orig)
    for ion in cp1d.ion
        if !ismissing(ion, :temperature)
            ion.temperature = deepcopy(Ti_orig)
        end
    end

    cost = cost_WPED_ztarget_pedratio!(cp1d, value_bound, ped_to_core_fraction, rho_ped, split)

    return cost
end

function cost_WPED_ztarget_pedratio!(
    cp1d::IMAS.core_profiles__profiles_1d,
    value_bound::Real,
    ped_to_core_fraction::Real,
    rho_ped::Real,
    split::NamedTuple)

    return wped_apply!(cp1d, value_bound, ped_to_core_fraction, rho_ped, split)[1]
end

"""
    wped_apply!(cp1d, value_bound, ped_to_core_fraction, rho_ped, split)

Sets the temperature profiles for a given pedestal boundary value, finding the exponents of the edge
profiles that keep the inverse scale lengths continuous at the pedestal boundary.

Returns the cost on the pedestal to core energy ratio and the two exponents (Te, Ti)
"""
function wped_apply!(
    cp1d::IMAS.core_profiles__profiles_1d,
    value_bound::Real,
    ped_to_core_fraction::Real,
    rho_ped::Real,
    split::NamedTuple)

    res_α_Te = Optim.optimize(α -> IMAS.force_float64(cost_WPED_α_Te!(cp1d, α, value_bound, rho_ped)), -500, 500, Optim.Brent(); rel_tol=IMAS.pedestal_search_tolerance.brent[])
    cost_WPED_α_Te!(cp1d, res_α_Te.minimizer, value_bound, rho_ped)

    Ti_bound = pedestal_ti(split, value_bound).Ti_ped
    res_α_Ti = Optim.optimize(α -> IMAS.force_float64(cost_WPED_α_Ti!(cp1d, α, Ti_bound, rho_ped)), -500, 500, Optim.Brent(); rel_tol=IMAS.pedestal_search_tolerance.brent[])
    cost_WPED_α_Ti!(cp1d, res_α_Ti.minimizer, Ti_bound, rho_ped)

    rho09 = 0.9
    core, edge = IMAS.core_edge_energy(cp1d, rho09)

    cost = (edge / core .- ped_to_core_fraction)^2

    return cost, res_α_Te.minimizer, res_α_Ti.minimizer
end

"""
    wped_residuals!(cp1d, Te_orig, Ti_orig, value_bound, α_Te, α_Ti, ped_to_core_fraction, rho_ped, split)

Sets the temperature profiles from the original ones for given pedestal boundary value and exponents,
and returns the three conditions that the WPED solution satisfies: the signed mismatch of the inverse
scale lengths of Te and of Ti at the pedestal boundary, and the pedestal to core energy ratio minus its target
"""
function wped_residuals!(
    cp1d::IMAS.core_profiles__profiles_1d{T},
    Te_orig::Vector{T},
    Ti_orig::Vector{T},
    value_bound::Real,
    α_Te::Real,
    α_Ti::Real,
    ped_to_core_fraction::Real,
    rho_ped::Real,
    split::NamedTuple) where {T<:Real}

    rho = cp1d.grid.rho_tor_norm
    Te = deepcopy(Te_orig)
    h_Te = IMAS.residual_WPED_α!(rho, Te, α_Te, value_bound, rho_ped)
    cp1d.electrons.temperature = Te
    Ti = deepcopy(Ti_orig)
    h_Ti = IMAS.residual_WPED_α!(rho, Ti, α_Ti, pedestal_ti(split, value_bound).Ti_ped, rho_ped)
    for ion in cp1d.ion
        if !ismissing(ion, :temperature)
            ion.temperature = Ti
        end
    end
    core, edge = IMAS.core_edge_energy(cp1d, 0.9)
    return T[h_Te, h_Ti, edge / core - ped_to_core_fraction]
end

function cost_WPED_α_Ti!(cp1d::IMAS.core_profiles__profiles_1d, α_Ti::Real, value_bound::Real, rho_ped::Real)
    Ti = cp1d.t_i_average
    rho = cp1d.grid.rho_tor_norm
    cost = IMAS.cost_WPED_α!(rho, Ti, α_Ti, value_bound, rho_ped)
    for ion in cp1d.ion
        if !ismissing(ion, :temperature)
            ion.temperature = Ti
        end
    end
    return cost
end

function cost_WPED_α_Te!(cp1d::IMAS.core_profiles__profiles_1d, α_Te::Real, value_bound::Real, rho_ped::Real)
    Te = cp1d.electrons.temperature
    rho = cp1d.grid.rho_tor_norm
    return IMAS.cost_WPED_α!(rho, Te, α_Te, value_bound, rho_ped)
end
