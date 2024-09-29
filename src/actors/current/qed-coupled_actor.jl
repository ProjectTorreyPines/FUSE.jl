import QED

#= ======== =#
#  ActorQEDcoupled  #
#= ======== =#
Base.@kwdef mutable struct FUSEparameters__ActorQEDcoupled{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    Δt::Entry{Float64} = Entry{Float64}("s", "Evolve for Δt")
    Nt::Entry{Int} = Entry{Int}("-", "Number of time steps during evolution"; default=100)
    #== display and debugging parameters ==#
    debug::Entry{Bool} = Entry{Bool}("-", "Turn on QED debugging/plotting"; default=false)
end

mutable struct ActorQEDcoupled{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorQEDcoupled{P}
    QO::Union{Nothing,QED.QED_state}
    build::Union{Nothing,QED.QED_build}
end

"""
    ActorQEDcoupled(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evolves the plasma and coil/vessel currents using the QED current diffusion solver.

!!! note

    This actor operates at "dd.global_time", any time advance must be done outside of the actor

        IMAS.new_timeslice!(dd, dd.global_time + Δt)
        dd.global_time += Δt
        ActorQEDcoupled(dd, act)

!!! note

    Stores data in `dd.core_profiles.profiles_1d[].j_ohmic`
"""
function ActorQEDcoupled(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorQEDcoupled(dd, act.ActorQEDcoupled; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorQEDcoupled(dd::IMAS.dd, par::FUSEparameters__ActorQEDcoupled; kw...)
    logging_actor_init(ActorQEDcoupled)
    par = par(kw...)
    return ActorQEDcoupled(dd, par, nothing, nothing)
end

function _step(actor::ActorQEDcoupled)
    dd = actor.dd
    par = actor.par

    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]

    # non_inductive contribution
    B0 = eqt.global_quantities.vacuum_toroidal_field.b0
    JBni = QED.FE(cp1d.grid.rho_tor_norm, cp1d.j_non_inductive .* B0)

    # initialize QED
    if actor.QO === nothing
        actor.QO = qed_init_from_imas(eqt, cp1d; uniform_rho = 33)
    else
        actor.QO.JBni = JBni
    end

    actor.build = qed_build_from_imas(dd, dd.global_time - par.Δt)

    # current diffusion
    time0 = dd.global_time - 0.5 * par.Δt
    actor.QO = QED.evolve(actor.QO, η_imas(dd.core_profiles.profiles_1d[time0]), actor.build, par.Δt, par.Nt; par.debug)

    return actor
end

function _finalize(actor::ActorQEDcoupled)
    dd = actor.dd

    eqt = dd.equilibrium.time_slice[]
    B0 = eqt.global_quantities.vacuum_toroidal_field.b0

    cp1d = dd.core_profiles.profiles_1d[]
    j_total = QED.JB(actor.QO; ρ=cp1d.grid.rho_tor_norm) ./ B0

    if ismissing(cp1d, :j_non_inductive)
        cp1d.j_ohmic = j_total
    else
        cp1d.j_ohmic = j_total .- cp1d.j_non_inductive
    end

    # NEED TO POPULATE COIL CURRENTS BACK IN dd
    build = actor.build
    current_per_turn = build.Ic

    active_coils = dd.pf_active.coil
    passive_loops = dd.pf_passive.loop
    for (k, coil) in enumerate(active_coils)
        IMAS.@ddtime(coil.current.data = current_per_turn[k])
    end
    Nactive = length(active_coils)
    for (k, loop) in enumerate(passive_loops)
        IMAS.@ddtime(loop.current = current_per_turn[k + Nactive])
    end

    # MAYBE POPULATE RESISTANCE SO IT'S NOT RECALCULATED?

    return actor
end

# utils
"""
    qed_init_from_imas(dd::IMAS.dd)

Setup QED build from data in IMAS `dd`
"""

# function loopelement2coil(loop, element)
#     outline = element.geometry.outline
#     @assert length(outline.r) == 4 "For the time being passive structures must be composed of quadrilateral elements"
#     passive_coil = VacuumFields.QuadCoil(outline.r, outline.z)
#     element_turns = ismissing(element, :turns_with_sign) ? 1 : abs(element.turns_with_sign)
#     loop_turns = sum(ismissing(elm, :turns_with_sign) ? 1 : abs(elm.turns_with_sign) for elm in loop.element)
#     Ic = ismissing(loop, :current) ? 0.0 : loop.current / Nturns)
#     VacuumFields.set_current!(passive_coil, Ic )
#     passive_coil.resistance = VacuumFields.resistance(passive_coil, loop.resistivity)
#     return passive_coil
# end

# elements(), turns(), QuadCoil(), and loop2multi() should all wind up in VacuumFields

function qed_build_from_imas(dd::IMAS.dd{D}, time0::D) where {D <: Real}

    coils = VacuumFields.MultiCoils(dd; load_pf_active=true, load_pf_passive=true);

    # Coil-only quantities
    Mcc = [VacuumFields.mutual(c1, c2) for c1 in coils, c2 in coils]
    Ic = [VacuumFields.current(c) / VacuumFields.turns(c) for c in coils] #c urrent per turn
    Rc = [VacuumFields.resistance(c) for c in coils];
    Vc = zero(Ic);

    eqt = dd.equilibrium.time_slice[time0]
    cp1d = dd.core_profiles.profiles_1d[time0]
    Ip = eqt.global_quantities.ip

    # plasma-coil mutuals
    image = VacuumFields.Image(eqt)
    Mpc = [VacuumFields.mutual(image, coil, Ip) for coil in coils]
    @warn "ActorQEDcoupled doesn't contain dMpc_dt yet"
    dMpc_dt = zero(Mpc) # How Mpc changes in time (like shape)... to test later

    # self inductance
    It = IMAS.cumtrapz(cp1d.grid.area, cp1d.j_tor)
    Wp = 0.5 * IMAS.trapz(cp1d.grid.psi, It)
    Li = 2 * Wp / Ip^2 # internal
    ψb = eqt.profiles_1d.psi[end]
    ψc = sum(Mpc[k] * Ic[k] for k in eachindex(coils))
    Le = (ψb - ψc) / Ip # external
    Lp = Li + Le # total self

    # plasma resistance
    # BCL 9/25/24: from Pohm, which may be wrong
    Pohm = dd.core_sources.source[:ohmic].profiles_1d[].electrons.power_inside[end]
    Ini = IMAS.get_time_array(dd.core_profiles.global_quantities, :current_non_inductive, time0, :linear)
    Iohm = Ip - Ini
    Rp = Pohm / (Ip * Iohm)

    # non-inductive voltage
    Vni = Rp * Ini

    # Waveforms
    # These should come from pulse_schedule
    #V_waveforms = fill(W0, length(coils));
    #W = QED.Waveform{D}(t -> dd.pulse_schedule.voltage[1](t))

    Nc = length(coils)
    V_waveforms = Vector{QED.Waveform{D}}(undef, length(coils))
    Vactive = Vwaveforms_from_pulse_schedule(dd, time0)
    @assert length(Vactive) == length(dd.pf_active.coil)
    Nactive = length(Vactive)
    V_waveforms[1:Nactive] .= Vactive
    V_waveforms[Nactive+1:end] .= fill(QED.Waveform{D}(t -> 0.0), Nc - Nactive)
    return QED.QED_build(Ic, Vc, Rc, Mcc, Vni, Rp, Lp, Mpc, dMpc_dt, V_waveforms)
end

function Vwaveforms_from_pulse_schedule(dd::IMAS.dd{D}, t_start::D=IMAS.global_time(dd)) where {D<:Real}
    Nc = length(dd.pf_active.coil)
    time = dd.pulse_schedule.pf_active.time .- t_start
    WFs = Vector{QED.Waveform{D}}(undef, Nc)
    for (k, supply) in enumerate(dd.pulse_schedule.pf_active.supply)
        if k in (1, 2)
            WFs[k] = QED.Waveform(time, supply.voltage.reference)
        elseif k == 3
            Cu = VacuumFields.MultiCoil(dd.pf_active.coil[3])
            Cl = VacuumFields.MultiCoil(dd.pf_active.coil[4])
            Lu = VacuumFields.mutual(Cu, Cu)
            Ll = VacuumFields.mutual(Cl, Cl)
            M  = VacuumFields.mutual(Cu, Cl)
            Lt = (Lu + 2M + Ll)

            facu = (Lu + M) / Lt
            WFs[3] = QED.Waveform(time , facu .* supply.voltage.reference)

            facl = (Ll + M) / Lt
            WFs[4] = QED.Waveform(time , facl .* supply.voltage.reference)
        elseif k < 12
            coil = dd.pf_active.coil[k+1]
            WFs[k+1] = QED.Waveform(time, supply.voltage.reference)
        else
            Cu = VacuumFields.MultiCoil(dd.pf_active.coil[13])
            Cl = VacuumFields.MultiCoil(dd.pf_active.coil[14])
            Lu = VacuumFields.mutual(Cu, Cu)
            Ll = VacuumFields.mutual(Cl, Cl)
            M  = -VacuumFields.mutual(Cu, Cl) # oppositely connected
            Lt = (Lu + 2M + Ll)

            # positive voltage gives negative current in VS3U
            facu = -(Lu + M) / Lt
            WFs[13] = QED.Waveform(time , facu .* supply.voltage.reference)

            # positive voltage gives positive current in VS3U
            facl = (Ll + M) / Lt
            WFs[14] = QED.Waveform(time , facl .* supply.voltage.reference)
        end
    end
    return WFs
end
