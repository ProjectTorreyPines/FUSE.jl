#= ========= =#
#  ActorWPED  #
#= ========= =#
Base.@kwdef mutable struct FUSEparameters__ActorWPED{T<:Real} <: ParametersActorPlasma{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    #== actor parameters ==#
    ped_to_core_fraction::Entry{T} = Entry{T}("-", "ratio of pedestal stored energy to core stored energy"; default=0.3)
    rho_nml::Entry{T} = Entry{T}("-", "Defines rho at which the no man's land region starts")
    rho_ped::Entry{T} = Entry{T}("-", "Defines rho at which the pedestal region starts") # rho_nml < rho_ped
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
    rho = cp1d.grid.rho_tor_norm
    rho_bound_idx = argmin(abs.(rho .- par.rho_ped))
    rho_edge = range(par.rho_ped, 1.0, 201)

    Ti_over_Te = cp1d.ion[1].temperature[rho_bound_idx] / cp1d.electrons.temperature[rho_bound_idx]

    if par.do_plot
        q = plot(rho, cp1d.electrons.temperature; label="Te before WPED blending", xlabel="rho", ylabel="Te [eV]")
        qq = plot(rho, cp1d.ion[1].temperature; label="Ti before WPED blending", xlabel="rho", ylabel="Ti [eV]")
    end

    res_value_bound = Optim.optimize(
        value_bound -> cost_WPED_ztarget_pedratio(cp1d, value_bound, par.ped_to_core_fraction, par.rho_ped, rho_edge, Ti_over_Te),
        1.0,
        cp1d.electrons.temperature[1],
        Optim.GoldenSection();
        rel_tol=1E-3)

    cost_WPED_ztarget_pedratio!(cp1d, res_value_bound.minimizer, par.ped_to_core_fraction, par.rho_ped, rho_edge, Ti_over_Te)

    summary_ped = dd.summary.local.pedestal
    nominal_rho_ped = max(0.95, par.rho_ped .+ (1.0 - par.rho_ped) / 2.0)
    @ddtime summary_ped.position.rho_tor_norm = nominal_rho_ped
    @ddtime summary_ped.n_e.value = IMAS.get_from(dd, Val{:ne_ped}, ne_from)
    @ddtime summary_ped.z_eff.value = IMAS.get_from(dd, Val{:zeff_ped}, zeff_ped_from)
    @ddtime summary_ped.t_e.value = interp1d(rho, cp1d.electrons.temperature).(nominal_rho_ped)
    @ddtime summary_ped.t_i_average.value = interp1d(rho, cp1d.ion[1].temperature).(nominal_rho_ped)
    blend_core_edge_Hmode(cp1d, summary_ped, par.rho_nml, par.rho_ped; what=:densities)

    if par.do_plot
        display(plot!(q, rho, cp1d.electrons.temperature; label="Te after"))
        display(plot!(qq, rho, cp1d.ion[1].temperature; label="Ti after"))
    end

    return actor
end

function cost_WPED_ztarget_pedratio(
    cp1d::IMAS.core_profiles__profiles_1d,
    value_bound::Real,
    ped_to_core_fraction::Real,
    rho_ped::Real,
    rho_edge::AbstractVector,
    Ti_over_Te::Real)

    cp1d_copy = deepcopy(cp1d)
    cost = cost_WPED_ztarget_pedratio!(cp1d_copy, value_bound, ped_to_core_fraction, rho_ped, rho_edge, Ti_over_Te)
    return cost
end

function cost_WPED_ztarget_pedratio!(
    cp1d::IMAS.core_profiles__profiles_1d,
    value_bound::Real,
    ped_to_core_fraction::Real,
    rho_ped::Real,
    rho_edge::AbstractVector,
    Ti_over_Te::Real)

    rho = cp1d.grid.rho_tor_norm
    Te = cp1d.electrons.temperature
    Ti = cp1d.ion[1].temperature

    res_α_Te = Optim.optimize(α -> cost_WPED_α_Te(cp1d, α, value_bound, rho_ped, rho_edge), -500, 500, Optim.GoldenSection(); rel_tol=1E-3)
    cost_WPED_α!(rho, Te, res_α_Te.minimizer, value_bound, rho_ped, rho_edge)

    res_α_Ti = Optim.optimize(α -> cost_WPED_α_Ti(cp1d, α, value_bound * Ti_over_Te, rho_ped, rho_edge), -500, 500, Optim.GoldenSection(); rel_tol=1E-3)
    cost_WPED_α!(rho, Ti, res_α_Ti.minimizer, value_bound * Ti_over_Te, rho_ped, rho_edge)
    for ion in cp1d.ion[2:end] # be carefull here to select only the thermal ions
        ion.temperature = Ti
    end

    core, edge = core_and_edge_energy(cp1d, rho_ped)
    ratio = edge / core

    cost = ((ratio .- ped_to_core_fraction) / ped_to_core_fraction)^2

    return cost
end

function cost_WPED_α_Ti(cp1d::IMAS.core_profiles__profiles_1d, α_Ti::Real, value_bound::Real, rho_ped::Real, rho_edge::AbstractVector)
    Ti = deepcopy(cp1d.ion[1].temperature)
    rho = cp1d.grid.rho_tor_norm
    cost = cost_WPED_α!(rho, Ti, α_Ti, value_bound, rho_ped, rho_edge)
    for ion in cp1d.ion[2:end] # be carefull here to select only the thermal ions
        ion.temperature = Ti
    end
    return cost
end

function cost_WPED_α_Te(cp1d::IMAS.core_profiles__profiles_1d, α_Te::Real, value_bound::Real, rho_ped::Real, rho_edge::AbstractVector)
    Te = deepcopy(cp1d.electrons.temperature)
    rho = cp1d.grid.rho_tor_norm
    return cost_WPED_α!(rho, Te, α_Te, value_bound, rho_ped, rho_edge)
end

function cost_WPED_α!(rho, temperature, α::Real, value_bound::Real, rho_ped::Real, rho_edge::AbstractVector)
    rho_bound_idx = argmin(abs.(rho .- rho_ped))

    z_whole_profile = IMAS.calc_z(rho, temperature, :backward)

    temperature[1:rho_bound_idx] = temperature[1:rho_bound_idx] .+ (-temperature[rho_bound_idx] + value_bound)

    profile_ped = IMAS.Edge_profile(rho_edge, rho_ped, value_bound, temperature[end], α)
    temperature[rho_bound_idx:end] .= IMAS.interp1d(rho_edge, profile_ped).(rho[rho_bound_idx:end])

    z_target_Te = z_whole_profile[rho_bound_idx]
    z_profile = IMAS.calc_z(rho_edge, profile_ped, :backward)

    return ((z_profile[1] - z_target_Te) / z_target_Te)^2
end

function core_and_edge_energy(cp1d::IMAS.core_profiles__profiles_1d, rho_ped::Real)
    p = cp1d.pressure_thermal
    rho_bound_idx = argmin(abs.(cp1d.grid.rho_tor_norm .- rho_ped))
    p_edge = deepcopy(p)
    p_edge[1:rho_bound_idx] .= p[rho_bound_idx]
    p_core = p .- p_edge
    core_value = 3 / 2 .* IMAS.cumtrapz(cp1d.grid.volume, p_core)
    edge_value = 3 / 2 .* IMAS.cumtrapz(cp1d.grid.volume, p_edge)
    return core_value[end], edge_value[end]
end

"""
    _finalize(actor::ActorWPED)

Writes results to dd.summary.local.pedestal and possibly updates core_profiles
"""
function _finalize(actor::ActorWPED)
    dd = actor.dd
    par = actor.par

    return actor
end
