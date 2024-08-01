#= ========= =#
#  ActorWPED  #
#= ========= =#
Base.@kwdef mutable struct FUSEparameters__ActorWPED{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    #== actor parameters ==#
    ped_to_core_fraction::Entry{T} = Entry{T}("-", "Ratio of edge (@rho=0.9) to core stored energy [0.05 for L-mode, 0.3 for neg-T plasmas]")
    rho_ped::Entry{T} = Entry{T}("-", "Defines rho at which the edge region starts")
    #== data flow parameters ==#
    ne_ped_from::Switch{Symbol} = switch_get_from(:ne_ped)
    zeff_ped_from::Switch{Symbol} = switch_get_from(:zeff_ped)
    #== display and debugging parameters ==#
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
end

mutable struct OptimizationWPED
    α_Te::Float64 # exponent for Te that is used to set inverse scale-length at the boundary
    α_Ti::Float64 # same for Ti
    value_bound::Float64 # value of Te at the boundary
end

mutable struct ActorWPED{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorWPED{P}
    optimization_guesses::OptimizationWPED
end

"""
    ActorWPED(dd::IMAS.dd, act::ParametersAllActors; kw...)

Finds the temperature profile at the edge to match the ped_to_core_fraction of stored energy set in par.ped_to_core_fraction
"""
function ActorWPED(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorWPED(dd, act.ActorWPED; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorWPED(dd::IMAS.dd, par::FUSEparameters__ActorWPED; kw...)
    logging_actor_init(ActorWPED)
    par = par(kw...)
    return ActorWPED(dd, par, OptimizationWPED(0.0, 0.0, 3e3))
end

"""
    _step(actor::ActorWPED)

Finds the temperature profile at the edge to match the ped_to_core_fraction of stored energy set in par.ped_to_core_fraction
"""
function _step(actor::ActorWPED{D,P}) where {D<:Real,P<:Real}
    dd = actor.dd
    par = actor.par

    cp1d = dd.core_profiles.profiles_1d[]

    summary_ped = dd.summary.local.pedestal
    @ddtime summary_ped.position.rho_tor_norm = par.rho_ped
    @ddtime summary_ped.n_e.value = IMAS.get_from(dd, Val{:ne_ped}, par.ne_ped_from, par.rho_ped)
    @ddtime summary_ped.zeff.value = IMAS.get_from(dd, Val{:zeff_ped}, par.zeff_ped_from, par.rho_ped)
    IMAS.blend_core_edge(:L_mode, cp1d, summary_ped, NaN, par.rho_ped; what=:densities)

    rho = cp1d.grid.rho_tor_norm
    rho_bound_idx = argmin(abs.(rho .- par.rho_ped))

    if par.do_plot
        q = plot(cp1d.electrons, :temperature; label="Te before WPED blending")
        qq = plot(cp1d, :t_i_average; label="Ti before WPED blending", xlabel="rho")
    end

    Ti_over_Te = cp1d.t_i_average[rho_bound_idx] / cp1d.electrons.temperature[rho_bound_idx]
    res_value_bound = Optim.optimize(
        value_bound -> cost_WPED_ztarget_pedratio(cp1d, value_bound, par.ped_to_core_fraction, par.rho_ped, Ti_over_Te),
        1.0,
        cp1d.electrons.temperature[1],
        Optim.GoldenSection();
        rel_tol=1E-3)

    cost_WPED_ztarget_pedratio!(cp1d, res_value_bound.minimizer, par.ped_to_core_fraction, par.rho_ped, Ti_over_Te)

    @ddtime summary_ped.t_e.value = IMAS.interp1d(rho, cp1d.electrons.temperature).(par.rho_ped)
    @ddtime summary_ped.t_i_average.value = IMAS.interp1d(rho, cp1d.t_i_average).(par.rho_ped)

    if par.do_plot
        display(plot!(q, cp1d.electrons, :temperature; label="Te after WPED blending"))
        display(plot!(qq, cp1d, :t_i_average; label="Ti after WPED blending"))
    end

    return actor
end

function cost_WPED_ztarget_pedratio(
    cp1d::IMAS.core_profiles__profiles_1d,
    value_bound::Real,
    ped_to_core_fraction::Real,
    rho_ped::Real,
    Ti_over_Te::Real
)
    cp1d_copy = deepcopy(cp1d)
    cost = cost_WPED_ztarget_pedratio!(cp1d_copy, value_bound, ped_to_core_fraction, rho_ped, Ti_over_Te)
    return cost
end

function cost_WPED_ztarget_pedratio!(
    cp1d::IMAS.core_profiles__profiles_1d,
    value_bound::Real,
    ped_to_core_fraction::Real,
    rho_ped::Real,
    Ti_over_Te::Real
)
    res_α_Te = Optim.optimize(α -> cost_WPED_α_Te!(cp1d, α, value_bound, rho_ped), -500, 500, Optim.GoldenSection(); rel_tol=1E-3)
    cost_WPED_α_Te!(cp1d, res_α_Te.minimizer, value_bound, rho_ped)

    res_α_Ti = Optim.optimize(α -> cost_WPED_α_Ti!(cp1d, α, value_bound * Ti_over_Te, rho_ped), -500, 500, Optim.GoldenSection(); rel_tol=1E-3)
    cost_WPED_α_Ti!(cp1d, res_α_Ti.minimizer, value_bound * Ti_over_Te, rho_ped)

    # ped_to_core_fraction is defined at ρ = 0.9
    core, edge = core_and_edge_energy(cp1d, 0.9)
    ratio = edge / core

    cost = ((ratio .- ped_to_core_fraction) / ped_to_core_fraction)^2

    return cost
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

function core_and_edge_energy(cp1d::IMAS.core_profiles__profiles_1d, rho_ped::Real)
    p = cp1d.pressure_thermal
    rho_tor_norm = cp1d.grid.rho_tor_norm
    rho_bound_idx = argmin(abs(rho - rho_ped) for rho in rho_tor_norm)
    pedge = p[rho_bound_idx]
    fedge = (k, x) -> (k <= rho_bound_idx) ? pedge : p[k]
    fcore = (k, x) -> (k <= rho_bound_idx) ? (p[k] - pedge) : 0.0
    volume = cp1d.grid.volume
    core_value = 1.5 * trapz(volume, fcore)
    edge_value = 1.5 * trapz(volume, fedge)
    return core_value, edge_value
end