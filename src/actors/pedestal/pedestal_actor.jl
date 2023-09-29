import EPEDNN

#= ============= =#
#  ActorPedestal  #
#= ============= =#
Base.@kwdef mutable struct FUSEparameters__ActorPedestal{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    #== actor parameters ==#
    rho_nml::Entry{T} = Entry{T}("-", "Defines rho at which the no man's land region starts")
    rho_ped::Entry{T} = Entry{T}("-", "Defines rho at which the pedestal region starts") # rho_nml < rho_ped in docstring
    T_ratio_pedestal::Entry{T} = Entry{T}("-", "Ratio of ion to electron temperatures"; default=1.0)
    ped_factor::Entry{T} = Entry{T}("-", "Pedestal height multiplier"; default=1.0)
    only_powerlaw::Entry{Bool} = Entry{Bool}("-", "EPED-NN uses power-law pedestal fit (without NN correction)"; default=false)
    #== data flow parameters ==#
    ip_from::Switch{Symbol} = switch_get_from(:ip)
    βn_from::Switch{Symbol} = switch_get_from(:βn)
    update_core_profiles::Entry{Bool} = Entry{Bool}("-", "Update core_profiles"; default=true)
    #== display and debugging parameters ==#
    warn_nn_train_bounds::Entry{Bool} = Entry{Bool}("-", "EPED-NN raises warnings if querying cases that are certainly outside of the training range"; default=false)
end

mutable struct ActorPedestal{D,P} <: PlasmaAbstractActor
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

    eq = dd.equilibrium
    eqt = eq.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]

    m = Int(round(IMAS.A_effective(cp1d) * 2.0)) / 2.0
    if !(m == 2.0 || m == 2.5)
        @warn "EPED-NN is only trained on m_effective = 2.0 & 2.5 , m_effective = $m"
    end

    neped = @ddtime dd.summary.local.pedestal.n_e.value
    zeffped = @ddtime dd.summary.local.pedestal.zeff.value
    Bt = abs(@ddtime(eq.vacuum_toroidal_field.b0)) * eq.vacuum_toroidal_field.r0 / eqt.boundary.geometric_axis.r

    actor.inputs = EPEDNN.InputEPED(
        eqt.boundary.minor_radius,
        IMAS.get_from(dd, Val{:βn}, par.βn_from),
        Bt,
        EPEDNN.effective_triangularity(eqt.boundary.triangularity_lower, eqt.boundary.triangularity_upper),
        abs(IMAS.get_from(dd, Val{:ip}, par.ip_from) / 1e6),
        eqt.boundary.elongation,
        m,
        neped / 1e19,
        eqt.boundary.geometric_axis.r,
        zeffped)

    sol = actor.epedmod(actor.inputs; par.only_powerlaw, par.warn_nn_train_bounds)

    if sol.pressure.GH.H * 1e6 < cp1d.pressure_thermal[end]
        actor.pped = 1.5 * sol.pressure.GH.H
        actor.wped = maximum(sol.width.GH.H, 0.01)
        @warn "EPED-NN output pedestal pressure is lower than separatrix pressure, p_ped=p_edge * 1.5 = $(round(actor.pped*1e6)) [Pa] assumed "
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

    dd_ped = dd.summary.local.pedestal
    @ddtime dd_ped.t_e.value = 2.0 * tped / (1.0 + par.T_ratio_pedestal) * par.ped_factor
    @ddtime dd_ped.t_i_average.value = @ddtime(dd_ped.t_e.value) * par.T_ratio_pedestal
    @ddtime dd_ped.position.rho_tor_norm = IMAS.interp1d(cp1d.grid.psi_norm, cp1d.grid.rho_tor_norm).(1 - actor.wped * sqrt(par.ped_factor))

    if par.update_core_profiles
        IMAS.blend_core_edge_Hmode(cp1d, dd_ped, par.rho_nml, par.rho_ped)
    end

    return actor
end
