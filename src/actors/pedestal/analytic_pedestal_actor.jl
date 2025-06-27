import EPEDNN

#= ===================== =#
#  ActorAnalyticPedestal  #
#= ===================== =#
Base.@kwdef mutable struct FUSEparameters__ActorAnalyticPedestal{T<:Real} <: ParametersActor{T}
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
    model::Switch{Symbol} = Switch{Symbol}([:MAST, :NSTX], "-", "Pedestal width model, w_ped~beta_p,ped^gamma, where gamma=1.0 for :NSTX and 0.5 for :MAST"; default=:MAST)
    width_coefficient::Entry{T} = Entry{T}("-", "Pedestal width coefficient, C2, w_ped=C2*beta_p,ped^gamma"; default=0.0)
    height_coefficient::Entry{T} = Entry{T}("-", "Pedestal height coefficient, C1, beta_p,ped=C1*w_ped^(3/4)/ip_norm^(1/3)"; default=0.0)
    ped_factor::Entry{T} = Entry{T}("-", "Pedestal height multiplier (width is scaled by sqrt of this factor)"; default=1.0, check=x -> @assert x > 0 "ped_factor must be > 0")
end

mutable struct ActorAnalyticPedestal{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorAnalyticPedestal{P}}
    inputs::EPEDNN.InputEPED
    wped::Union{Missing,Real}
    pped::Union{Missing,Real}
end

"""
    ActorAnalyticPedestal(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evaluates the pedestal boundary condition (height and width)
"""
function ActorAnalyticPedestal(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorAnalyticPedestal(dd, act.ActorAnalyticPedestal; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorAnalyticPedestal(dd::IMAS.dd, par::FUSEparameters__ActorAnalyticPedestal; kw...)
    logging_actor_init(ActorAnalyticPedestal)
    par = OverrideParameters(par; kw...)
    return ActorAnalyticPedestal(dd, par, EPEDNN.InputEPED(), missing, missing)
end

# Step function
function _step(actor::ActorAnalyticPedestal{D,P}) where {D<:Real,P<:Real}
    dd = actor.dd
    par = actor.par

    cp1d = dd.core_profiles.profiles_1d[]
    eqt = dd.equilibrium.time_slice[]

    w_ped_ne = IMAS.pedestal_tanh_width_half_maximum(cp1d.grid.rho_tor_norm, cp1d.electrons.density_thermal)

    # NOTE: Throughout FUSE, the "pedestal" density is the density at rho=0.9
    # the conversion from ne_ped09 to ne_ped with w_ped is based on this
    # peaked density profile: 1.0 / IMAS.Hmode_profiles(0.0, 1.0, 100, 1.0, 1.0, w_ped_ne)[90] ∼ 0.828
    # flat density profile: 1.0 / IMAS.Hmode_profiles(0.0, 1.0, 100, 0.0, 1.0, 0.05)[90] ∼ 0.866
    rho09 = 0.9
    tanh_width_to_09_factor = 1.0 / IMAS.Hmode_profiles(0.0, 1.0, 100, 0.5, 1.0, w_ped_ne)[90]
    ne09 = IMAS.get_from(dd, Val{:ne_ped}, par.ne_from, rho09)
    neped = ne09 * tanh_width_to_09_factor
    zeffped = IMAS.get_from(dd, Val{:zeff_ped}, par.zeff_from, rho09)
    βn = IMAS.get_from(dd, Val{:βn}, par.βn_from)
    ip = IMAS.get_from(dd, Val{:ip}, par.ip_from)
    Bt = abs(eqt.global_quantities.vacuum_toroidal_field.b0) * eqt.global_quantities.vacuum_toroidal_field.r0 / eqt.boundary.geometric_axis.r

    R = (eqt.profiles_1d.r_outboard[end-1] + eqt.profiles_1d.r_inboard[end-1]) / 2.0
    a = (eqt.profiles_1d.r_outboard[end-1] - eqt.profiles_1d.r_inboard[end-1]) / 2.0
    κ = eqt.profiles_1d.elongation[end-1]

    # poloidal beta
    ip_norm = (ip / 1e6) / (a * Bt)
    lpol = eqt.global_quantities.length_pol
    Bp = IMAS.mks.μ_0 * ip / lpol

    actor.inputs.a = a
    @assert !isnan(βn)
    actor.inputs.bt = Bt
    actor.inputs.ip = abs(ip) / 1e6
    actor.inputs.kappa = κ
    actor.inputs.neped = neped / 1e19
    actor.inputs.r = R
    actor.inputs.zeffped = zeffped

    if par.model == :MAST
        # MAST like width - w_ped~beta_p,ped^0.5
        if par.height_coefficient == 0.0
            par.height_coefficient = 4.0
        end
        if par.width_coefficient == 0.0
            par.width_coefficient = 0.11
        end
        betap_ped = (par.height_coefficient^(8 / 5) * par.width_coefficient^(6 / 5)) / ip_norm^(8 / 15)
        actor.wped = (par.height_coefficient^(4 / 5) * par.width_coefficient^(8 / 5)) / ip_norm^(4 / 15) # full width
        actor.wped = actor.wped / 2 # eped width definition (1/2 width)
        actor.pped = 1e-6 * betap_ped * Bp^2 / 2 / IMAS.mks.μ_0

    elseif par.model == :NSTX
        # NSTX like width - w_ped~beta_p,ped^1.0
        if par.height_coefficient == 0.0
            par.height_coefficient = 2.0
        end
        if par.width_coefficient == 0.0
            par.width_coefficient = 0.4
        end
        betap_ped = (par.height_coefficient^4 * par.width_coefficient^3) / ip_norm^(4 / 3)
        actor.wped = par.height_coefficient^4.0 * par.width_coefficient^4.0 / ip_norm^(4 / 3) # full width
        actor.wped = actor.wped / 2 # eped width definition (1/2 width)
        actor.pped = 1e-6 * betap_ped * Bp^2 / 2 / IMAS.mks.μ_0
    else
        error("Undefined model! Model should be one of :MAST, :NSTX")
    end

    return actor
end

"""
    _finalize(actor::ActorAnalyticPedestal)

Writes results to dd.summary.local.pedestal and updates core_profiles
"""
function _finalize(actor::ActorAnalyticPedestal)
    return __finalize(actor)
end
