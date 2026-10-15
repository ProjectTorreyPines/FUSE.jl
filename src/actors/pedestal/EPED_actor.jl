import EPEDNN

#= ========= =#
#  ActorEPED  #
#= ========= =#
@actor_parameters_struct ActorEPED{T} begin
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
    nn_model::Switch{Symbol} = Switch{Symbol}([:variable_nesep_ratio, :fixed_nesep_ratio], "-",
        "EPED-NN model: :variable_nesep_ratio (ensemble that takes ne_sep/ne_ped and Te_sep as inputs, and predicts the pedestal height with an uncertainty)" *
        " or :fixed_nesep_ratio (single network, trained at a fixed ne_sep/ne_ped ~ 0.25)"; default=:variable_nesep_ratio)
    ped_factor::Entry{T} = Entry{T}("-", "Pedestal height multiplier (width is scaled by sqrt of this factor)"; default=1.0, check=x -> @assert x > 0 "ped_factor must be > 0")
    only_powerlaw::Entry{Bool} = Entry{Bool}("-", "EPED-NN uses power-law pedestal fit (without NN correction). Only applies to nn_model=:fixed_nesep_ratio"; default=true)
    #== display and debugging parameters ==#
    warn_nn_train_bounds::Entry{Bool} = Entry{Bool}("-", "EPED-NN raises warnings if querying cases that are certainly outside of the training range"; default=false)
end

mutable struct ActorEPED{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.DD{D}
    par::OverrideParameters{P,FUSEparameters__ActorEPED{P}}
    epedmod::EPEDNN.EPEDmodel
    inputs::EPEDNN.InputEPED
    wped::Union{Missing,Real} # pedestal width using EPED definition (1/2 width as fraction of psi_norm)
    pped::Union{Missing,Real} # pedestal height using EPED units (MPa)
    σ_frac::Float64 # uncertainty of the pedestal height, as a fraction of the height itself
end

"""
    ActorEPED(dd::IMAS.DD, act::ParametersAllActors; kw...)

Predicts pedestal pressure and width using the EPED neural network model.

The actor utilizes the EPED (Edge Pedestal Equilibrium and Dynamics) model to predict 
pedestal height and width based on global plasma parameters. EPED combines physics-based 
scaling laws with neural network corrections trained on experimental pedestal data.

Model capabilities:
- Physics-based power-law scaling for robust extrapolation
- Optional neural network corrections for improved accuracy
- Calibrated against experimental data from multiple tokamaks
- Handles both conventional and spherical tokamak geometries

Key inputs (extracted from plasma state):
- Machine geometry (R, a, κ, δ, triangularity)
- Global parameters (βn, Ip, Bt, effective mass)  
- Pedestal conditions (ne_ped, Zeff_ped)

Outputs:
- Pedestal pressure height in MPa (pped)
- Pedestal width as fraction of normalized poloidal flux (wped)
- Fractional uncertainty of the pedestal height (σ_frac)
- Automatic fallback to edge pressure + 10% if EPED prediction is too low

Two EPED-NN models are available, see `act.ActorEPED.nn_model`:
- `:variable_nesep_ratio` (default): a deep ensemble that adds the separatrix conditions
  (`ne_sep/ne_ped` and `Te_sep`) to the inputs and predicts the pedestal height with an
  uncertainty; the pedestal width follows the analytic EPED1 width law
- `:fixed_nesep_ratio`: the original single network, trained at a fixed `ne_sep/ne_ped ~ 0.25`,
  which predicts both the pedestal height and its width
"""
function ActorEPED(dd::IMAS.DD, act::ParametersAllActors; kw...)
    actor = ActorEPED(dd, act.ActorEPED; kw...)
    step(actor)
    finalize(actor)
    return actor
end


function ActorEPED(dd::IMAS.DD, par::FUSEparameters__ActorEPED; kw...)
    logging_actor_init(ActorEPED)
    par = OverrideParameters(par; kw...)
    return ActorEPED(dd, par, eped_model(par.nn_model), EPEDNN.InputEPED(), missing, missing, 0.0)
end

"""
    eped_model(nn_model::Symbol)

Load the EPED-NN model selected by `act.ActorEPED.nn_model`
"""
function eped_model(nn_model::Symbol)
    if nn_model == :variable_nesep_ratio
        return EPEDNN.loadmodelonce("EPED1NNensemble.bson")
    elseif nn_model == :fixed_nesep_ratio
        return EPEDNN.loadmodelonce("EPED1NNmodel.bson")
    else
        error("act.ActorEPED.nn_model can only be :variable_nesep_ratio or :fixed_nesep_ratio")
    end
end

"""
    _step(actor::ActorEPED)

Executes the EPED model prediction and validates results against plasma edge conditions.

The step function calls the EPED neural network model with current plasma parameters 
and checks that the predicted pedestal pressure exceeds the separatrix pressure. 
If the prediction is too low, it applies a 10% safety margin above edge pressure 
while maintaining the predicted width.
"""
function _step(actor::ActorEPED{D,P}) where {D<:Real,P<:Real}
    dd = actor.dd
    par = actor.par

    cp1d = dd.core_profiles.profiles_1d[]
    eqt = dd.equilibrium.time_slice[]

    sol = run_EPED!(dd, actor.inputs, actor.epedmod; par.ne_from, par.zeff_from, par.βn_from, par.ip_from, par.Te_sep, par.only_powerlaw, par.warn_nn_train_bounds)
    pped, actor.σ_frac = EPEDNN.pedestal_height(actor.epedmod, actor.inputs, sol)

    if pped < 1.1 * cp1d.pressure_thermal[end] / 1e6
        actor.pped = 1.1 * cp1d.pressure_thermal[end] / 1E6
        @warn "EPED-NN output pedestal pressure is lower than separatrix pressure, p_ped=p_edge * 1.1 = $(round(actor.pped*1e6)) [Pa] assumed "
    else
        actor.pped = pped
    end

    # NOTE: the width is evaluated on `actor.pped`, so that it stays consistent with the pedestal
    # pressure that is actually used, also when the separatrix pressure sets the pedestal height.
    # βpol_ped is only used by the models that predict the pedestal height alone.
    βpol_ped = IMAS.pedestal_poloidal_beta(eqt, actor.pped * 1e6)
    actor.wped = max(EPEDNN.pedestal_width(actor.epedmod, sol, βpol_ped), 0.005)

    return actor
end

"""
    _finalize(actor::ActorEPED)

Applies EPED predictions to plasma temperature and density profiles.

Writes pedestal pressure and width to dd.summary.local.pedestal and updates 
core_profiles by blending the EPED-predicted pedestal conditions with existing 
core profiles using H-mode profile functions with proper particle balance.
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
    w_ped = (1.0 - position) .* par.ped_factor .^0.5
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

"""
    run_EPED(
        dd::IMAS.DD;
        ne_from::Symbol,
        zeff_from::Symbol,
        βn_from::Symbol,
        ip_from::Symbol,
        Te_sep::Real,
        only_powerlaw::Bool,
        warn_nn_train_bounds::Bool,
        nn_model::Symbol=:variable_nesep_ratio)

Runs EPED-NN from dd and returns the pedestal height `pped` [MPa], its width `wped` (1/2 width as a
fraction of psi_norm) and the uncertainty `σ_frac` of the height (as a fraction of the height)
"""
function run_EPED(
    dd::IMAS.DD;
    ne_from::Symbol,
    zeff_from::Symbol,
    βn_from::Symbol,
    ip_from::Symbol,
    Te_sep::Real,
    only_powerlaw::Bool,
    warn_nn_train_bounds::Bool,
    nn_model::Symbol=:variable_nesep_ratio)

    inputs = EPEDNN.InputEPED()
    epedmod = eped_model(nn_model)
    sol = run_EPED!(dd, inputs, epedmod; ne_from, zeff_from, βn_from, ip_from, Te_sep, only_powerlaw, warn_nn_train_bounds)
    pped, σ_frac = EPEDNN.pedestal_height(epedmod, inputs, sol)
    βpol_ped = IMAS.pedestal_poloidal_beta(dd.equilibrium.time_slice[], pped * 1e6)
    return (pped=pped, wped=EPEDNN.pedestal_width(epedmod, sol, βpol_ped), σ_frac=σ_frac)
end

"""
    run_EPED!(
        dd::IMAS.DD,
        eped_inputs::EPEDNN.InputEPED,
        epedmod::EPEDNN.EPEDmodel;
        ne_from::Symbol,
        zeff_from::Symbol,
        βn_from::Symbol,
        ip_from::Symbol,
        Te_sep::Real,
        only_powerlaw::Bool,
        warn_nn_train_bounds::Bool)

Fills `eped_inputs` from dd, runs EPED-NN and outputs the solution of `epedmod`
(a `PedestalSolution` for a `EPED1NNmodel`, height and uncertainty for a `EPED1NNensemble`)
"""
function run_EPED!(
    dd::IMAS.DD,
    eped_inputs::EPEDNN.InputEPED,
    epedmod::EPEDNN.EPEDmodel;
    ne_from::Symbol,
    zeff_from::Symbol,
    βn_from::Symbol,
    ip_from::Symbol,
    Te_sep::Real,
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
    # separatrix conditions: only used by the EPED1NNensemble model
    eped_inputs.nesep_ratio = cp1d.electrons.density_thermal[end] / neped
    eped_inputs.tesep = Te_sep

    return EPEDNN.run_epednn(epedmod, eped_inputs; only_powerlaw, warn_nn_train_bounds)
end
