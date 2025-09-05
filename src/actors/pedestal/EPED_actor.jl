import EPEDNN

#= ========= =#
#  ActorEPED  #
#= ========= =#
Base.@kwdef mutable struct FUSEparameters__ActorEPED{T<:Real} <: ParametersActor{T}
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
    #== actor parameters==#
    ped_factor::Entry{T} = Entry{T}("-", "Pedestal height multiplier (width is scaled by sqrt of this factor)"; default=1.0, check=x -> @assert x > 0 "ped_factor must be > 0")
    only_powerlaw::Entry{Bool} = Entry{Bool}("-", "EPED-NN uses power-law pedestal fit (without NN correction)"; default=true)
    #== display and debugging parameters ==#
    warn_nn_train_bounds::Entry{Bool} = Entry{Bool}("-", "EPED-NN raises warnings if querying cases that are certainly outside of the training range"; default=false)
end

mutable struct ActorEPED{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorEPED{P}}
    epedmod::EPEDNN.EPEDmodel
    inputs::EPEDNN.InputEPED
    wped::Union{Missing,Real} # pedestal width using EPED definition (1/2 width as fraction of psi_norm)
    pped::Union{Missing,Real} # pedestal height using EPED units (MPa)
end

"""
    ActorEPED(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evaluates the pedestal boundary condition (height and width)
"""
function ActorEPED(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorEPED(dd, act.ActorEPED; kw...)
    step(actor)
    finalize(actor)
    return actor
end


function ActorEPED(dd::IMAS.dd, par::FUSEparameters__ActorEPED; kw...)
    logging_actor_init(ActorEPED)
    par = OverrideParameters(par; kw...)
    epedmod = EPEDNN.loadmodelonce("EPED1NNmodel.bson")
    return ActorEPED(dd, par, epedmod, EPEDNN.InputEPED(), missing, missing)
end

"""
    _step(actor::ActorEPED)

Runs pedestal actor to evaluate pedestal width and height
"""
function _step(actor::ActorEPED{D,P}) where {D<:Real,P<:Real}
    dd = actor.dd
    par = actor.par

    cp1d = dd.core_profiles.profiles_1d[]
    sol = run_EPED!(dd, actor.inputs, actor.epedmod; par.ne_from, par.zeff_from, par.βn_from, par.ip_from, par.only_powerlaw, par.warn_nn_train_bounds)

    if sol.pressure.GH.H < 1.1 * cp1d.pressure_thermal[end] / 1e6
        actor.pped = 1.1 * cp1d.pressure_thermal[end] / 1E6
        actor.wped = max(sol.width.GH.H, 0.005)
        @warn "EPED-NN output pedestal pressure is lower than separatrix pressure, p_ped=p_edge * 1.1 = $(round(actor.pped*1e6)) [Pa] assumed "
    else
        actor.pped = sol.pressure.GH.H
        actor.wped = sol.width.GH.H
    end

    return actor
end

"""
    _finalize(actor::ActorEPED)

Writes results to dd.summary.local.pedestal and updates core_profiles
"""
function _finalize(actor::ActorEPED)
    return __finalize(actor)
end

function __finalize(actor::Union{ActorEPED,ActorAnalyticPedestal})
    dd = actor.dd
    par = actor.par

    cp1d = dd.core_profiles.profiles_1d[]
    rho = cp1d.grid.rho_tor_norm

    # NOTE: Standard EPED uses 1/2 width as fraction of psi_norm
    #       Instead FUSE, IMAS and the Hmode_profiles functions use the full width as a function of rho_tor_norm.
    from_ped_to_full_width = 2.0
    position = IMAS.interp1d(cp1d.grid.psi_norm, rho).(1 - actor.wped * from_ped_to_full_width * sqrt(par.ped_factor))
    w_ped = 1.0 - position
    old_t_i_ped = IMAS.interp1d(rho, cp1d.t_i_average).(position)

    impurity = [IMAS.avgZ(ion.element[1].z_n, old_t_i_ped) for ion in cp1d.ion if !IMAS.is_hydrogenic(ion)][1]
    zi = sum(impurity) / length(impurity)
    @assert actor.inputs.zeffped <= zi
    nival = actor.inputs.neped * 1e19 * (actor.inputs.zeffped - 1) / (zi^2 - zi)
    nval = actor.inputs.neped * 1e19 - zi * nival
    nsum = actor.inputs.neped * 1e19 + nval + nival
    tped = (actor.pped * 1e6) / nsum / IMAS.mks.e

    Ti_over_Te = ti_te_ratio(cp1d, par.T_ratio_pedestal, par.rho_nml, par.rho_ped)
    t_e = 2.0 * tped / (1.0 + Ti_over_Te) * par.ped_factor
    t_i_average = t_e * Ti_over_Te

    # Change the last point of the temperatures profiles since
    # The rest of the profile will be taken care by the blend_core_edge_Hmode() function
    cp1d.electrons.temperature[end] = par.Te_sep
    for ion in cp1d.ion
        if !ismissing(ion, :temperature)
            ion.temperature[end] = cp1d.electrons.temperature[end] * Ti_over_Te
        end
    end

    # blend core_profiles core
    cp1d.electrons.temperature = IMAS.blend_core_edge_Hmode(cp1d.electrons.temperature, rho, t_e, w_ped, par.rho_nml, par.rho_ped)
    ti_avg_new = IMAS.blend_core_edge_Hmode(cp1d.t_i_average, rho, t_i_average, w_ped, par.rho_nml, par.rho_ped)
    for ion in cp1d.ion
        if !ismissing(ion, :temperature)
            ion.temperature = ti_avg_new
        end
    end

    return actor
end

function run_EPED(
    dd::IMAS.dd;
    ne_from::Symbol,
    zeff_from::Symbol,
    βn_from::Symbol,
    ip_from::Symbol,
    only_powerlaw::Bool,
    warn_nn_train_bounds::Bool)

    inputs = EPEDNN.InputEPED()
    epedmod = EPEDNN.loadmodelonce("EPED1NNmodel.bson")
    return run_EPED!(dd, inputs, epedmod; ne_from, zeff_from, βn_from, ip_from, only_powerlaw, warn_nn_train_bounds)
end

"""
    run_EPED!(
        dd::IMAS.dd,
        eped_inputs::EPEDNN.InputEPED,
        epedmod::EPEDNN.EPED1NNmodel;
        ne_from::Symbol,
        zeff_from::Symbol,
        βn_from::Symbol,
        ip_from::Symbol,
        only_powerlaw::Bool,
        warn_nn_train_bounds::Bool)

Runs EPED from dd and outputs the EPED solution as the sol struct
"""
function run_EPED!(
    dd::IMAS.dd,
    eped_inputs::EPEDNN.InputEPED,
    epedmod::EPEDNN.EPED1NNmodel;
    ne_from::Symbol,
    zeff_from::Symbol,
    βn_from::Symbol,
    ip_from::Symbol,
    only_powerlaw::Bool,
    warn_nn_train_bounds::Bool)

    cp1d = dd.core_profiles.profiles_1d[]
    eqt = dd.equilibrium.time_slice[]

    m = round(Int, IMAS.A_effective(cp1d) * 2.0, RoundNearest) / 2.0
    if !(m == 2.0 || m == 2.5)
        @warn "EPED-NN is only trained on m_effective = 2.0 & 2.5 , m_effective = $m"
    end

    w_ped_ne = IMAS.pedestal_tanh_width_half_maximum(cp1d.grid.rho_tor_norm, cp1d.electrons.density_thermal)

    # NOTE: Throughout FUSE, the "pedestal" density is the density at rho=0.9
    # the conversion from ne_ped09 to ne_ped with w_ped is based on this
    # peaked density profile: 1.0/IMAS.Hmode_profiles(0.0, 1.0, 100, 1.0, 1.0, w_ped_ne)[90] ∼ 0.828
    # flat density profile: 1.0/IMAS.Hmode_profiles(0.0, 1.0, 100, 0.0, 1.0, 0.05)[90] ∼ 0.866
    rho09 = 0.9
    tanh_width_to_09_factor = 1.0 / IMAS.Hmode_profiles(0.0, 1.0, 100, 0.5, 1.0, w_ped_ne)[90]
    ne09 = IMAS.get_from(dd, Val(:ne_ped), ne_from, rho09)
    neped = ne09 * tanh_width_to_09_factor
    zeffped = IMAS.get_from(dd, Val(:zeff_ped), zeff_from, rho09)
    βn = IMAS.get_from(dd, Val(:βn), βn_from)
    ip = IMAS.get_from(dd, Val(:ip), ip_from)
    Bt = abs(eqt.global_quantities.vacuum_toroidal_field.b0) * eqt.global_quantities.vacuum_toroidal_field.r0 / eqt.boundary.geometric_axis.r

    # NOTE: EPED results can be very sensitive to δu, δl
    #
    # eqt.boundary can have small changes in κ, δu, δl just due to contouring
    # This issue can be mitigated using higher grid resolutions in the equilibrium solver.
    #
    # Here we use the flux surface right inside of the LCFS, and not the LCFS itself.
    # Not only this avoids these sensitivity issues, but it's actually more correct,
    # since the TOQ equilibrium used by EPED is a fixed boundary equilibrium solver,
    # and as such it cuts out psi at 99% or similar.
    if false
        R = eqt.boundary.geometric_axis.r
        a = eqt.boundary.minor_radius
        κ = eqt.boundary.elongation
        δu = eqt.boundary.triangularity_upper
        δl = eqt.boundary.triangularity_lower
    elseif false
        pr, pz = IMAS.boundary(dd.pulse_schedule.position_control)
        R, a, κ, δu, δl, ζou, ζol, ζil, ζiu = IMAS.miller_R_a_κ_δ_ζ(pr, pz)
    else
        R = (eqt.profiles_1d.r_outboard[end-1] + eqt.profiles_1d.r_inboard[end-1]) / 2.0
        a = (eqt.profiles_1d.r_outboard[end-1] - eqt.profiles_1d.r_inboard[end-1]) / 2.0
        κ = eqt.profiles_1d.elongation[end-1]
        δu = eqt.profiles_1d.triangularity_upper[end-1]
        δl = eqt.profiles_1d.triangularity_lower[end-1]
    end

    eped_inputs.a = a
    @assert !isnan(βn)
    eped_inputs.betan = βn
    eped_inputs.bt = Bt
    eped_inputs.delta = EPEDNN.effective_triangularity(δu, δl)
    eped_inputs.ip = abs(ip) / 1e6
    eped_inputs.kappa = κ
    eped_inputs.m = m
    eped_inputs.neped = neped / 1e19
    eped_inputs.r = R
    eped_inputs.zeffped = zeffped

    return epedmod(eped_inputs; only_powerlaw, warn_nn_train_bounds)
end
