#= ========= =#
#  ActorWPED  #
#= ========= =#
Base.@kwdef mutable struct FUSEparameters__ActorWPED{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    #== common pedestal parameters==#
    rho_nml::Entry{T} = Entry{T}("-", "Defines rho at which the no man's land region starts")
    rho_ped::Entry{T} = Entry{T}("-", "Defines rho at which the pedestal region starts") # rho_nml < rho_ped
    T_ratio_pedestal::Entry{T} =
        Entry{T}("-", "Ratio of ion to electron temperatures (or rho at which to sample for that ratio, if negative; or rho_nml-(rho_ped-rho_nml) if 0.0)"; default=1.0)
    Te_sep::Entry{T} = Entry{T}("-", "Separatrix electron temperature"; default=80.0, check=x -> @assert x > 0 "Te_sep must be > 0")
    ip_from::Switch{Symbol} = switch_get_from(:ip)
    βn_from::Switch{Symbol} = switch_get_from(:βn)
    ne_from::Switch{Symbol} = switch_get_from(:ne_ped)
    zeff_from::Switch{Symbol} = switch_get_from(:zeff_ped)
    #== actor parameters ==#
    density_factor::Entry{T} = Entry{T}("-", "Scale input density by given factor"; default=1.0)
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

    # Throughout FUSE, the "pedestal" density is the density at rho=0.9
    rho09 = 0.9
    @ddtime summary_ped.n_e.value = IMAS.get_from(dd, Val{:ne_ped}, par.ne_from, rho09) * par.density_factor
    @ddtime summary_ped.zeff.value = IMAS.get_from(dd, Val{:zeff_ped}, par.zeff_from, rho09) # zeff is taken as the average value
    @ddtime summary_ped.position.rho_tor_norm = par.rho_ped

    if par.do_plot
        ppe = plot(cp1d.electrons, :temperature; label="Te before WPED blending")
        ppi = plot(cp1d, :t_i_average; label="Ti before WPED blending", xlabel="rho")
    end

    if ismissing(par, :ped_to_core_fraction)
        core, edge = IMAS.core_edge_energy(cp1d, 0.9)
        ped_to_core_fraction = edge / core
    else
        ped_to_core_fraction = par.ped_to_core_fraction
    end

    Ti_over_Te = ti_te_ratio(cp1d, par.T_ratio_pedestal, par.rho_nml, par.rho_ped)

    # Change the last point of the temperatures profiles since that's what used in cost_WPED_ztarget_pedratio()
    cp1d.electrons.temperature[end] = par.Te_sep
    for ion in cp1d.ion
        if !ismissing(ion, :temperature)
            ion.temperature[end] = cp1d.electrons.temperature[end] * Ti_over_Te
        end
    end

    res_value_bound = Optim.optimize(
        value_bound -> cost_WPED_ztarget_pedratio(cp1d, value_bound, ped_to_core_fraction, par.rho_ped, Ti_over_Te),
        1.0,
        cp1d.electrons.temperature[1],
        Optim.GoldenSection();
        rel_tol=1E-3)

    cost_WPED_ztarget_pedratio!(cp1d, res_value_bound.minimizer, ped_to_core_fraction, par.rho_ped, Ti_over_Te)

    @ddtime summary_ped.t_e.value = IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.electrons.temperature).(par.rho_ped)
    @ddtime summary_ped.t_i_average.value = IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.t_i_average).(par.rho_ped)

    if par.do_plot
        display(plot!(ppe, cp1d.electrons, :temperature; label="Te after WPED blending"))
        display(plot!(ppi, cp1d, :t_i_average; label="Ti after WPED blending"))
    end

    return actor
end

function cost_WPED_ztarget_pedratio(
    cp1d::IMAS.core_profiles__profiles_1d,
    value_bound::Real,
    ped_to_core_fraction::Real,
    rho_ped::Real,
    Ti_over_Te::Real)

    cp1d_copy = deepcopy(cp1d)
    cost = cost_WPED_ztarget_pedratio!(cp1d_copy, value_bound, ped_to_core_fraction, rho_ped, Ti_over_Te)
    return cost
end

function cost_WPED_ztarget_pedratio!(
    cp1d::IMAS.core_profiles__profiles_1d,
    value_bound::Real,
    ped_to_core_fraction::Real,
    rho_ped::Real,
    Ti_over_Te::Real)

    res_α_Te = Optim.optimize(α -> cost_WPED_α_Te!(cp1d, α, value_bound, rho_ped), -500, 500, Optim.GoldenSection(); rel_tol=1E-3)
    cost_WPED_α_Te!(cp1d, res_α_Te.minimizer, value_bound, rho_ped)

    res_α_Ti = Optim.optimize(α -> cost_WPED_α_Ti!(cp1d, α, value_bound * Ti_over_Te, rho_ped), -500, 500, Optim.GoldenSection(); rel_tol=1E-3)
    cost_WPED_α_Ti!(cp1d, res_α_Ti.minimizer, value_bound * Ti_over_Te, rho_ped)

    core, edge = IMAS.core_edge_energy(cp1d, 0.9)

    cost = (edge / core .- ped_to_core_fraction)^2

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
