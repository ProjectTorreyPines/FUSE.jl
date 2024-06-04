import EPEDNN

#= ============= =#
#  ActorPedestal  #
#= ============= =#
Base.@kwdef mutable struct FUSEparameters__ActorPedestal{T<:Real} <: ParametersActorPlasma{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    #== actor parameters ==#
    rho_nml::Entry{T} = Entry{T}("-", "Defines rho at which the no man's land region starts")
    rho_ped::Entry{T} = Entry{T}("-", "Defines rho at which the pedestal region starts") # rho_nml < rho_ped
    T_ratio_pedestal::Entry{T} =
        Entry{T}("-", "Ratio of ion to electron temperatures (or rho at which to sample for that ratio, if negative; or rho_nml-(rho_ped-rho_nml) if 0.0)"; default=1.0)
    ped_factor::Entry{T} = Entry{T}("-", "Pedestal height multiplier"; default=1.0)
    only_powerlaw::Entry{Bool} = Entry{Bool}("-", "EPED-NN uses power-law pedestal fit (without NN correction)"; default=false)
    #== data flow parameters ==#
    ip_from::Switch{Symbol} = switch_get_from(:ip)
    βn_from::Switch{Symbol} = switch_get_from(:βn)
    ne_ped_from::Switch{Symbol} = switch_get_from(:ne_ped)
    zeff_ped_from::Switch{Symbol} = switch_get_from(:zeff_ped)
    update_core_profiles::Entry{Bool} = Entry{Bool}("-", "Update core_profiles"; default=true)
    #== display and debugging parameters ==#
    warn_nn_train_bounds::Entry{Bool} = Entry{Bool}("-", "EPED-NN raises warnings if querying cases that are certainly outside of the training range"; default=false)
end

mutable struct ActorPedestal{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorPedestal{P}
    epedmod::EPEDNN.EPEDmodel
    inputs::Union{Missing,EPEDNN.InputEPED}
    wped::Union{Missing,Real}
    pped::Union{Missing,Real}
end

"""
    ActorPedestal(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evaluates the pedestal boundary condition (height and width)
"""
function ActorPedestal(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorPedestal(dd, act.ActorPedestal; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorPedestal(dd::IMAS.dd, par::FUSEparameters__ActorPedestal; kw...)
    logging_actor_init(ActorPedestal)
    par = par(kw...)
    epedmod = EPEDNN.loadmodelonce("EPED1NNmodel.bson")
    return ActorPedestal(dd, par, epedmod, missing, missing, missing)
end

"""
    _step(actor::ActorPedestal)

Runs pedestal actor to evaluate pedestal width and height
"""
function _step(actor::ActorPedestal)
    dd = actor.dd
    par = actor.par

    eqt = eq.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]

    sol = run_EPED(eqt, cp1d; par.ne_ped_from, par.zeff_ped_from, par.βn_from, par.ip_from, par.only_powerlaw, par.warn_nn_train_bounds)

    if sol.pressure.GH.H < 1.5 * cp1d.pressure_thermal[end] / 1E6
        actor.pped = 1.5 * cp1d.pressure_thermal[end] / 1E6
        actor.wped = max(sol.width.GH.H, 0.01)
        @warn "EPED-NN output pedestal pressure is lower than separatrix pressure, p_ped = p_edge * 1.5 = $(round(actor.pped)) [MPa] assumed."
    else
        actor.pped = sol.pressure.GH.H
        actor.wped = sol.width.GH.H
    end

    return actor
end

"""
    _finalize(actor::ActorPedestal)

Writes results to dd.summary.local.pedestal and possibly updates core_profiles
"""
function _finalize(actor::ActorPedestal)
    dd = actor.dd
    par = actor.par

    cp1d = dd.core_profiles.profiles_1d[]
    impurity = [ion.element[1].z_n for ion in cp1d.ion if Int(floor(ion.element[1].z_n)) != 1][1]
    zi = sum(impurity) / length(impurity)
    nival = actor.inputs.neped * 1e19 * (actor.inputs.zeffped - 1) / (zi^2 - zi)
    nval = actor.inputs.neped * 1e19 - zi * nival
    nsum = actor.inputs.neped * 1e19 + nval + nival
    tped = (actor.pped * 1e6) / nsum / constants.e

    if par.T_ratio_pedestal == 0.0
        T_ratio_pedestal = IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.ion[1].temperature ./ cp1d.electrons.temperature)(par.rho_nml - (par.rho_ped - par.rho_nml))
    elseif par.T_ratio_pedestal <= 0.0
        T_ratio_pedestal = IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.ion[1].temperature ./ cp1d.electrons.temperature)(par.T_ratio_pedestal)
    else
        T_ratio_pedestal = par.T_ratio_pedestal
    end

    dd_ped = dd.summary.local.pedestal
    @ddtime dd_ped.n_e.value = actor.inputs.neped * 1e19
    @ddtime dd_ped.t_e.value = 2.0 * tped / (1.0 + T_ratio_pedestal) * par.ped_factor
    @ddtime dd_ped.t_i_average.value = @ddtime(dd_ped.t_e.value) * T_ratio_pedestal
    @ddtime dd_ped.position.rho_tor_norm = IMAS.interp1d(cp1d.grid.psi_norm, cp1d.grid.rho_tor_norm).(1 - actor.wped * sqrt(par.ped_factor))

    if par.update_core_profiles
        IMAS.blend_core_edge_Hmode(cp1d, dd_ped, par.rho_nml, par.rho_ped)
    end

    return actor
end

"""
    run_EPED(
        eqt::IMAS.equilibrium__time_slice,
        cp1d::IMAS.core_profiles__profiles_1d;
        ne_ped_from::Symbol,
        zeff_ped_from::Symbol,
        βn_from::Symbol,
        ip_from::Symbol,
        only_powerlaw::Bool,
        warn_nn_train_bounds::Bool)

Sets up and runs EPED from data in IMAS and outputs the EPED solution structure
"""
function run_EPED(
    eqt::IMAS.equilibrium__time_slice,
    cp1d::IMAS.core_profiles__profiles_1d;
    ne_ped_from::Symbol,
    zeff_ped_from::Symbol,
    βn_from::Symbol,
    ip_from::Symbol,
    only_powerlaw::Bool,
    warn_nn_train_bounds::Bool)

    m = Int(round(IMAS.A_effective(cp1d) * 2.0)) / 2.0
    if !(m == 2.0 || m == 2.5)
        @warn "EPED-NN is only trained on m_effective = 2.0 & 2.5 , m_effective = $m"
    end

    neped = IMAS.get_from(dd, Val{:ne_ped}, ne_ped_from)
    zeffped = IMAS.get_from(dd, Val{:zeff_ped}, zeff_ped_from)
    βn = IMAS.get_from(dd, Val{:βn}, βn_from)
    ip = IMAS.get_from(dd, Val{:ip}, ip_from)

    Bt = abs(@ddtime(eqt.global_quantities.vacuum_toroidal_field.b0)) * eqt.global_quantities.vacuum_toroidal_field.r0 / eqt.boundary.geometric_axis.r

    inputs = EPEDNN.InputEPED(
        eqt.boundary.minor_radius,
        βn,
        Bt,
        EPEDNN.effective_triangularity(eqt.boundary.triangularity_lower, eqt.boundary.triangularity_upper),
        abs(ip) / 1e6,
        eqt.boundary.elongation,
        m,
        neped / 1e19,
        eqt.boundary.geometric_axis.r,
        zeffped)

    epedmod = EPEDNN.loadmodelonce("EPED1NNmodel.bson")

    return epedmod(inputs; only_powerlaw, warn_nn_train_bounds)
end