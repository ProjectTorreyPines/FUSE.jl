import EPEDNN

#= ============= =#
#  ActorPedestal  #
#= ============= =#
mutable struct ActorPedestal <: PlasmaAbstractActor
    dd::IMAS.dd
    par::ParametersActor
    epedmod::EPEDNN.EPEDmodel
    inputs::Union{Missing,EPEDNN.InputEPED}
    wped::Union{Missing,Real}
    pped::Union{Missing,Real}
end

Base.@kwdef struct FUSEparameters__ActorPedestal{T} <: ParametersActor where {T<:Real}
    update_core_profiles = Entry(Bool, "", "Update core_profiles"; default=true)
    edge_bound = Entry(Real, "", "Defines rho at which edge starts"; default=0.8)
    temp_pedestal_ratio = Entry(Real, "", "Ratio of ion to electron temperatures"; default=1.0)
    ped_factor = Entry(Real, "", "Pedestal height multiplier"; default=1.0)
    warn_nn_train_bounds = Entry(Bool, "", "EPED-NN raises warnings if querying cases that are certainly outside of the training range"; default=false)
    only_powerlaw = Entry(Bool, "", "EPED-NN uses power-law pedestal fit (without NN correction)"; default=false)
end

"""
    ActorPedestal(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evaluates the pedestal boundary condition (height and width)
"""
function ActorPedestal(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorPedestal(kw...)
    actor = ActorPedestal(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function ActorPedestal(dd::IMAS.dd, par::ParametersActor; kw...)
    logging_actor_init(ActorPedestal)
    par = par(kw...)

    epedmod = EPEDNN.loadmodelonce("EPED1NNmodel.bson")
    return ActorPedestal(dd, par, epedmod, missing, missing, missing)
end

"""
    step(actor::ActorPedestal;
        warn_nn_train_bounds::Bool=actor.par.warn_nn_train_bounds,
        only_powerlaw::Bool=false)

Runs pedestal actor to evaluate pedestal width and height
"""
function _step(actor::ActorPedestal;
    warn_nn_train_bounds::Bool=actor.par.warn_nn_train_bounds,
    only_powerlaw::Bool=false)

    dd = actor.dd
    eq = dd.equilibrium
    eqt = eq.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]

    m = [ion.element[1].a for ion in cp1d.ion if Int(floor(ion.element[1].z_n)) == 1]
    m = sum(m) / length(m)
    if m < 2
        m = 2
    elseif m > 2
        m = 2.5
    end

    neped = @ddtime dd.summary.local.pedestal.n_e.value
    zeffped = @ddtime dd.summary.local.pedestal.zeff.value
    Bt = abs(@ddtime(eq.vacuum_toroidal_field.b0)) * eq.vacuum_toroidal_field.r0 / eqt.boundary.geometric_axis.r
    βn = @ddtime(dd.summary.global_quantities.beta_tor_thermal_norm.value)

    actor.inputs = EPEDNN.InputEPED(
        eqt.boundary.minor_radius,
        βn,
        Bt,
        eqt.boundary.triangularity,
        abs(eqt.global_quantities.ip / 1e6),
        eqt.boundary.elongation,
        m,
        neped / 1e19,
        eqt.boundary.geometric_axis.r,
        zeffped)

    sol = actor.epedmod(actor.inputs; only_powerlaw, warn_nn_train_bounds)

    actor.wped = sol.width.GH.H
    if sol.pressure.GH.H * 1e6 < cp1d.pressure_thermal[end]
        actor.pped = 1.5* sol.pressure.GH.H
        @warn "EPED-NN output pedestal pressure is lower than edge pressure, p_ped=p_edge * 1.5 = $(round(actor.pped*1e6)) [Pa] assumed "

    else
        actor.pped = sol.pressure.GH.H
    end
    return actor
end

"""
    _finalize(actor::ActorPedestal;
        temp_pedestal_ratio::Real=actor.par.temp_pedestal_ratio,
        ped_factor::Real=actor.par.ped_factor,
        edge_bound::Real=actor.par.edge_bound,
        update_core_profiles::Bool=actor.par.update_core_profiles)

Writes results to dd.summary.local.pedestal and possibly updates core_profiles
"""
function _finalize(actor::ActorPedestal;
    temp_pedestal_ratio::Real=actor.par.temp_pedestal_ratio,
    ped_factor::Real=actor.par.ped_factor,
    edge_bound::Real=actor.par.edge_bound,
    update_core_profiles::Bool=actor.par.update_core_profiles)

    dd = actor.dd

    cp1d = dd.core_profiles.profiles_1d[]
    impurity = [ion.element[1].z_n for ion in cp1d.ion if Int(floor(ion.element[1].z_n)) != 1][1]
    zi = sum(impurity) / length(impurity)
    nival = actor.inputs.neped * 1e19 * (actor.inputs.zeffped - 1) / (zi^2 - zi)
    nval = actor.inputs.neped * 1e19 - zi * nival
    nsum = actor.inputs.neped * 1e19 + nval + nival
    tped = (actor.pped * 1e6) / nsum / constants.e

    dd_ped = dd.summary.local.pedestal
    @ddtime dd_ped.t_e.value = 2.0 * tped / (1.0 + temp_pedestal_ratio) * ped_factor
    @ddtime dd_ped.t_i_average.value = @ddtime(dd_ped.t_e.value) * temp_pedestal_ratio
    @ddtime dd_ped.position.rho_tor_norm = 1 - actor.wped * sqrt(ped_factor)

    if update_core_profiles
        IMAS.blend_core_edge_Hmode(cp1d, dd_ped, edge_bound)
    end
end
