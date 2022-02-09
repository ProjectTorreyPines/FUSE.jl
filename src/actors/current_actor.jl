import QED

mutable struct QEDcurrentActor <: AbstractActor
    dd::IMAS.dd
    QI::QED.QED_state
    η#::Base.Callable
    QO::Union{QED.QED_state,Missing}
    initial_time
    tmax
end

#= =============== =#
#  QEDcurrentActor  #
#= =============== =#
function QEDcurrentActor(dd::IMAS.dd)
    QEDcurrentActor(dd, from_imas(dd), η_imas(dd), missing, @ddtime(dd.equilibrium.time), 0.0)
end

function step(actor::QEDcurrentActor, tmax, Nt, Vedge = nothing, Ip = nothing; resume=false)
    if resume
        if actor.QO !== missing
            actor.QI = actor.QO
            actor.QO = missing
        end
        actor.tmax += tmax
    else
        actor.initial_time = @ddtime(dd.equilibrium.time)
        actor.tmax = tmax
    end
    actor.QO = QED.diffuse(actor.QI, actor.η, tmax, Nt; Vedge, Ip)
end

function finalize(actor::QEDcurrentActor)
    dd = actor.dd
    eqt = dd.equilibrium.time_slice[]

    newtime = actor.initial_time + actor.tmax
    resize!(dd.equilibrium.time_slice, newtime)
    dd.equilibrium.time_slice[newtime] = new = deepcopy(eqt)

    dΡ_dρ = new.profiles_1d.rho_tor[end]
    ρ = new.profiles_1d.rho_tor / dΡ_dρ

    new.profiles_1d.q = 1.0 ./ actor.QO.ι.(ρ)
    new.profiles_1d.j_tor = actor.QO.JtoR.(ρ) ./ new.profiles_1d.gm9
    new.time = newtime

    new
end


# utils
function from_imas(dd::IMAS.dd)
    eqt = dd.equilibrium.time_slice[]

    dΡ_dρ = eqt.profiles_1d.rho_tor[end]
    ρ = eqt.profiles_1d.rho_tor / dΡ_dρ

    B₀ = @ddtime(dd.equilibrium.vacuum_toroidal_field.b0)

    fsa_R⁻² = QED.FE(ρ, eqt.profiles_1d.gm1)
    F = QED.FE(ρ, eqt.profiles_1d.f)

    # Require dV_dρ=0 on-axis
    tmp = dΡ_dρ .* eqt.profiles_1d.dvolume_drho_tor
    tmp[1] = 0.0
    dV_dρ = QED.FE(ρ, tmp)

    ι = QED.FE(ρ, 1.0 ./ eqt.profiles_1d.q)
    JtoR = QED.FE(ρ, eqt.profiles_1d.j_tor .* eqt.profiles_1d.gm9)

    JBni = nothing
    if !ismissing(dd.core_profiles.profiles_1d[], :j_non_inductive)
        JBni = QED.FE(dd.core_profiles.profiles_1d[].grid.rho_tor_norm, dd.core_profiles.profiles_1d[].j_non_inductive .* B₀)
    end

    return QED.QED_state(ρ, dΡ_dρ, B₀, fsa_R⁻², F, dV_dρ, ι, JtoR, JBni = JBni)
end

function η_imas(dd::IMAS.dd)
    rho = dd.core_profiles.profiles_1d[1].grid.rho_tor_norm
    η = 1.0 ./ dd.core_profiles.profiles_1d[1].conductivity_parallel
    return QED.FE(rho, η)
end