import EPEDNN # CHECK. Currently using this to get some parameters but not sure if this is needed or not

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
    height_coefficient::Entry{T} = Entry{T}("-", "Pedestal height coefficient, C1, beta_p,ped=C1*w_ped^(3/4)/IN^(1/3)"; default=0.0)
    ped_factor::Entry{T} = Entry{T}("-", "Pedestal height multiplier (width is scaled by sqrt of this factor)"; default=1.0, check=x -> @assert x > 0 "ped_factor must be > 0")
    #== display and debugging parameters ==#
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
end

mutable struct ActorAnalyticPedestal{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorAnalyticPedestal{P}}
    inputs::EPEDNN.InputEPED # CHECK - have made it work with this but not sure this is the correct approach?
    wped::Union{Missing,Real} # pedestal width  CHECK WHAT DEFINITION IS BEING USED HERE!!!. using EPED definition (1/2 width as fraction of psi_norm)
    pped::Union{Missing,Real} # pedestal height CHECK, It hink this model is in kPA!!!! using EPED units (MPa)
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
function _step(actor::ActorAnalyticPedestal{D,P}) where {D<:Real, P<:Real}
    dd = actor.dd
    par = actor.par
    cp1d = dd.core_profiles.profiles_1d[]
    eqt = dd.equilibrium.time_slice[]

    w_ped_ne = IMAS.pedestal_tanh_width_half_maximum(cp1d.grid.rho_tor_norm, cp1d.electrons.density_thermal)

    # NOTE: Throughout FUSE, the "pedestal" density is the density at rho=0.9
    # the conversion from ne_ped09 to ne_ped with w_ped is based on this
    # peaked density profile: 1.0/IMAS.Hmode_profiles(0.0, 1.0, 100, 1.0, 1.0, w_ped_ne)[90] ∼ 0.828
    # flat density profile: 1.0/IMAS.Hmode_profiles(0.0, 1.0, 100, 0.0, 1.0, 0.05)[90] ∼ 0.866
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

    actor.inputs.a = a
    @assert !isnan(βn)
    actor.inputs.bt = Bt
    actor.inputs.ip = abs(ip) / 1e6
    actor.inputs.kappa = κ
    actor.inputs.neped = neped / 1e19
    actor.inputs.r = R
    actor.inputs.zeffped = zeffped
    eqt = dd.equilibrium.time_slice[]

    IN = (ip/1e6)/(a*Bt) # Normalized plasma current IN=IP[MA]/(a[m]*BT[T])

    # NSTX like width - w_ped~beta_p,ped^0.5
    lpol = dd.equilibrium.time_slice[].global_quantities.length_pol
    Bp = IMAS.mks.μ_0*ip/lpol
    # MAST like width - w_ped~beta_p,ped^0.5
    if par.model == :MAST
        if par.height_coefficient==0.0
            par.height_coefficient=4.0
        end
        if par.width_coefficient==0.0
            par.width_coefficient=0.11
        end

        betap_ped = (par.height_coefficient^(8/5) * par.width_coefficient^(6/5)) / IN^(8/15)
        actor.wped  = (par.height_coefficient^(4/5) * par.width_coefficient^(8/5)) / IN^(4/15)
        #actor.wped = par.height_coefficient^0.8*par.width_coefficient^1.6/IN^0.25 
        #actor.pped = (1/1000.)*32*par.height_coefficient^1.6*par.width_coefficient^1.2*IN^1.5*actor.inputs.bt^2/(1+actor.inputs.kappa^2) 
        #@show(actor.pped)
    # NSTX like width - w_ped~beta_p,ped^1.0
    elseif par.model == :NSTX
        if par.height_coefficient==0.0
            par.height_coefficient=2.0
        end
        if par.width_coefficient==0.0
            par.width_coefficient=0.4
        end
        betap_ped = (par.height_coefficient^4 * par.width_coefficient^3) / IN^(4/3)
        actor.wped = par.height_coefficient^4.0*par.width_coefficient^4.0/IN^(4/3) 
        #actor.pped = (1/1000.)*32*par.height_coefficient^4.0*par.width_coefficient^3.0*IN^(2/3)*actor.inputs.bt^2/(1+actor.inputs.kappa^2)
        #@show(actor.pped)
    else
        error("Undefined model! Model should be one of :MAST, :NSTX")
    end
    actor.pped = 1e-6*betap_ped*Bp^2/2/IMAS.mks.μ_0
    #@show(actor.pped)

    #if par.do_debug
        # Debug output goes here - WIP
    #end

    # Plotting
    if par.do_plot
        # Plotting goes here - WIP
    end

    return actor
end

# Finalize function
"""
    _finalize(actor::ActorAnalyticPedestal)

Writes results to dd.summary.local.pedestal and updates core_profiles
"""
function _finalize(actor::ActorAnalyticPedestal)
    # zeff, neped, pped
    dd = actor.dd
    par = actor.par

    cp1d = dd.core_profiles.profiles_1d[]
    rho = cp1d.grid.rho_tor_norm

    impurity = [ion.element[1].z_n for ion in cp1d.ion if Int(floor(ion.element[1].z_n)) != 1][1]
    zi = sum(impurity) / length(impurity)
    @assert actor.inputs.zeffped <= zi
    nival = actor.inputs.neped * 1e19 * (actor.inputs.zeffped - 1) / (zi^2 - zi)
    nval = actor.inputs.neped * 1e19 - zi * nival
    nsum = actor.inputs.neped * 1e19 + nval + nival
    tped = (actor.pped * 1e6) / nsum / IMAS.mks.e

    Ti_over_Te = ti_te_ratio(cp1d, par.T_ratio_pedestal, par.rho_nml, par.rho_ped)

    # NOTE: EPED uses 1/2 width as fraction of psi_norm instead FUSE, IMAS, and the Hmode_profiles functions use the full width as a function of rho_tor_norm
    from_ped_to_full_width = 1.0
    t_e = 2.0 * tped / (1.0 + Ti_over_Te) * par.ped_factor
    t_i_average = t_e * Ti_over_Te
    position = IMAS.interp1d(cp1d.grid.psi_norm, rho).(1 - actor.wped * from_ped_to_full_width * sqrt(par.ped_factor))
    w_ped = 1.0 - position

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
