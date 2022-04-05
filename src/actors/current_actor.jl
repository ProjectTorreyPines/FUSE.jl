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
    push!(dd.equilibrium.vacuum_toroidal_field.b0,dd.equilibrium.vacuum_toroidal_field.b0[end])

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
    rho_tor = eqt.profiles_1d.rho_tor
    B0 = @ddtime(dd.equilibrium.vacuum_toroidal_field.b0)
    gm1 = eqt.profiles_1d.gm1
    f = eqt.profiles_1d.f
    dvolume_drho_tor = eqt.profiles_1d.dvolume_drho_tor
    q = eqt.profiles_1d.q
    j_tor = eqt.profiles_1d.j_tor
    gm9 = eqt.profiles_1d.gm9

    if !isempty(dd.core_profiles.profiles_1d) && !ismissing(dd.core_profiles.profiles_1d[],:j_non_inductive)
        prof1d = dd.core_profiles.profiles_1d[]
        ρ_j_non_inductive = (prof1d.grid.rho_tor_norm, prof1d.j_non_inductive)
    else
        ρ_j_non_inductive = nothing
    end

    return QED.initialize(rho_tor, B0, gm1, f, dvolume_drho_tor, q, j_tor, gm9; ρ_j_non_inductive)
end

function η_imas(dd::IMAS.dd)
    rho = dd.core_profiles.profiles_1d[].grid.rho_tor_norm
    η = 1.0 ./ dd.core_profiles.profiles_1d[].conductivity_parallel
    return QED.FE(rho, η)
end