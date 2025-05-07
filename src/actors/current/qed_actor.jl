import QED

#= ======== =#
#  ActorQED  #
#= ======== =#
Base.@kwdef mutable struct FUSEparameters__ActorQED{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    Δt::Entry{Float64} = Entry{Float64}("s", "Evolve for Δt (Inf for steady state)"; default=Inf, check=x -> @assert x >= 0 "Δt must be >= 0.0")
    Nt::Entry{Int} = Entry{Int}("-", "Number of time steps during evolution"; default=100, check=x -> @assert x > 0 "Nt must be > 0")
    solve_for::Switch{Symbol} = Switch{Symbol}([:ip, :vloop], "-", "Solve for specified Ip or Vloop"; default=:ip)
    qmin_desired::Entry{Float64} = Entry{Float64}("-", "Keep the minimum magnitude of the q-profile above this value"; default=1.0, check=x -> @assert x >= 0 "qmin_desired >= 0")
    #== data flow parameters ==#
    ip_from::Switch{Symbol} = switch_get_from(:ip)
    vloop_from::Switch{Symbol} = switch_get_from(:vloop)
end

mutable struct ActorQED{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorQED{P}}
    act::ParametersAllActors{P}
    ip_controller::ActorControllerIp{D,P}
    QO::Union{Nothing,QED.QED_state}
end

"""
    ActorQED(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evolves the plasma current using the QED current diffusion solver.

The fundamental quantitiy being solved is `j_total` in `dd.core_profiles.profiles_1d[]`

!!! note

    This actor operates at "dd.global_time", any time advance must be done outside of the actor

        IMAS.new_timeslice!(dd, dd.global_time + Δt)
        dd.global_time += Δt
        ActorQED(dd, act)
"""
function ActorQED(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorQED(dd, act.ActorQED, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorQED(dd::IMAS.dd, par::FUSEparameters__ActorQED, act::ParametersAllActors; kw...)
    logging_actor_init(ActorQED)
    par = OverrideParameters(par; kw...)
    ip_controller = ActorControllerIp(dd, act.ActorControllerIp)
    return ActorQED(dd, par, act, ip_controller, nothing)
end

function _step(actor::ActorQED)
    dd = actor.dd
    par = actor.par

    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]

    B0 = eqt.global_quantities.vacuum_toroidal_field.b0
    # no ohmic, no sawteeth, no time dependent
    j_non_inductive = IMAS.total_sources(dd.core_sources, cp1d; time0=dd.global_time, exclude_indexes=[7, 409, 701], fields=[:j_parallel]).j_parallel

    # initialize QED
    # we must reinitialize to update the equilibrium metrics
    actor.QO = qed_init_from_imas(actor; uniform_rho=501)

    # QED calculates the total current based on q, which goes to infinity at the separatrix
    # this leads to some small but not negligible difference in the total current calculated
    # internally by QED and the one from `dd`.
    ratio = QED.Ip(actor.QO) / IMAS.Ip(cp1d, eqt)

    if par.Nt == 0
        # only initialize, nothing to do

    elseif par.Δt > 0.0 && par.Δt < Inf
        # current diffusion
        t0 = dd.global_time - par.Δt
        t1 = dd.global_time

        if par.solve_for == :vloop && par.vloop_from == :controllers__ip
            # staircase approach to call Ip control at each step of the current ramp: one QED diffuse call for each time step
            # NOTE: `QED.diffuse` is inefficient, would be beneficial to have a `QED.diffuse!` function without allocations
            δt = par.Δt / par.Nt
            No = par.Nt
            Ni = 1
        else
            δt = t1 - t0
            No = 1
            Ni = par.Nt
        end

        for time0 in range(t0, t1, No + 1)[1:end-1]
            if par.solve_for == :ip
                Ip = IMAS.get_from(dd, Val{:ip}, par.ip_from; time0)
                Vedge = nothing
            else
                # run Ip controller if vloop_from == :controllers__ip
                if par.vloop_from == :controllers__ip
                    finalize(step(actor.ip_controller; time0))
                end
                Ip = nothing
                Vedge = IMAS.get_from(dd, Val{:vloop}, par.vloop_from; time0)
            end

            # check where q<1 based on the the q-profile at the previous
            # time-step and change the resisitivity to keep q>1
            qval = 1.0 ./ abs.(actor.QO.ι.(cp1d.grid.rho_tor_norm))
            i_qdes = findlast(qval .< par.qmin_desired)
            if i_qdes === nothing
                i_qdes = 0
            end

            η_jardin, flattened_j_non_inductive = η_JBni_sawteeth(cp1d, j_non_inductive, i_qdes)
            actor.QO.JBni = QED.FE(cp1d.grid.rho_tor_norm, flattened_j_non_inductive .* B0)

            actor.QO = QED.diffuse(actor.QO, η_jardin, δt, Ni; Vedge, Ip, debug=false)

            cp1d.j_total = QED.JB(actor.QO; ρ=cp1d.grid.rho_tor_norm) ./ B0
            cp1d.j_non_inductive = flattened_j_non_inductive
            eqt.profiles_1d.q =  1.0 ./ actor.QO.ι.(eqt.profiles_1d.rho_tor_norm)
#            @ddtime(dd.core_profiles.global_quantities.ip = QED.Ip(actor.QO))
        end

    elseif par.Δt == Inf
        # steady state solution
        if par.solve_for == :ip
            Ip = IMAS.get_from(dd, Val{:ip}, par.ip_from)
            Vedge = nothing
        else
            Ip = nothing
            Vedge = IMAS.get_from(dd, Val{:vloop}, par.vloop_from)
        end

        # we need to run steady state twice, the first time to find the q-profile when the
        # current fully relaxes, and the second time we change the resisitivity to keep q>1
        i_qdes = 0
        for _ in (1, 2)
            η_jardin, flattened_j_non_inductive = η_JBni_sawteeth(cp1d, j_non_inductive, i_qdes)
            actor.QO.JBni = QED.FE(cp1d.grid.rho_tor_norm, flattened_j_non_inductive .* B0)

            actor.QO = QED.steady_state(actor.QO, η_jardin; Vedge, Ip)

            # check where q<1
            qval = 1.0 ./ abs.(actor.QO.ι.(cp1d.grid.rho_tor_norm))
            i_qdes = findlast(qval .< par.qmin_desired)
            if i_qdes === nothing
                break
            end

            cp1d.j_total = QED.JB(actor.QO; ρ=cp1d.grid.rho_tor_norm) ./ B0
            cp1d.j_non_inductive = flattened_j_non_inductive
#            @ddtime(dd.core_profiles.global_quantities.ip = QED.Ip(actor.QO))
        end
    else
        error("act.ActorQED.Δt = $(par.Δt) is not valid")
    end

    return actor
end

# utils
"""
    qed_init_from_imas(actor::ActorQED{D,P}; uniform_rho::Int, j_tor_from::Symbol=:core_profiles, ip_from::Union{Symbol,Real}=j_tor_from) where {D<:Real,P<:Real}

Setup QED from data in IMAS `dd`
"""
function qed_init_from_imas(actor::ActorQED{D,P}; uniform_rho::Int, j_tor_from::Symbol=:core_profiles, ip_from::Union{Symbol,Real}=j_tor_from) where {D<:Real,P<:Real}
    dd = actor.dd
    par = actor.par

    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]
    B0 = eqt.global_quantities.vacuum_toroidal_field.b0

    rho_tor = eqt.profiles_1d.rho_tor
    gm1 = eqt.profiles_1d.gm1
    f = eqt.profiles_1d.f
    dvolume_drho_tor = eqt.profiles_1d.dvolume_drho_tor
    q = eqt.profiles_1d.q
    gm9 = eqt.profiles_1d.gm9

    # DO NOT use the equilibrium j_tor, since it's quality depends on the quality/resolution of the equilibrium solver
    # better to use the j_tor from core_profiles, which is the same quantity that is input in the equilibrium solver
    if j_tor_from === :equilibrium
        j_tor = eqt.profiles_1d.j_tor
    elseif j_tor_from === :core_profiles
        j_tor = IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.j_tor, :cubic).(IMAS.norm01(rho_tor))
    else
        error("j_tor_from must be :equilibrium or :core_profiles")
    end

    if ip_from === :equilibrium
        Ip0 = IMAS.get_from(dd, Val{:ip}, :equilibrium)
    elseif ip_from === :core_profiles
        Ip0 = IMAS.get_from(dd, Val{:ip}, :core_profiles)
    elseif typeof(ip_from) <: Real
        Ip0 = ip_from
    else
        error("ip_from must be :equilibrium, :core_profiles, or a real number")
    end

    if ismissing(cp1d, :j_non_inductive)
        ρ_j_non_inductive = nothing
    else
        i_qdes = findlast(abs.(eqt.profiles_1d.q) .< par.qmin_desired)
        if i_qdes === nothing
            i_qdes = 0
        end
        _, j_non_inductive = η_JBni_sawteeth(cp1d, cp1d.j_non_inductive, i_qdes)
        ρ_j_non_inductive = (cp1d.grid.rho_tor_norm, j_non_inductive)
    end

    ρ_grid = collect(range(0.0, 1.0, uniform_rho))

    return QED.initialize(rho_tor, B0, gm1, f, dvolume_drho_tor, q, j_tor, gm9; ρ_j_non_inductive, ρ_grid, Ip0)
end

"""
    η_JBni_sawteeth(cp1d::IMAS.core_profiles__profiles_1d{T}, j_non_inductive::Vector{T}, i_qdes::Int; use_log::Bool=true) where {T<:Real}

returns

  - resistivity profile using Jardin's model for stationary sawteeth changes the plasma resistivity to raise q>1

  - non-inductive profile with flattening of the current inside of the inversion radius
"""
function η_JBni_sawteeth(cp1d::IMAS.core_profiles__profiles_1d{T}, j_non_inductive::Vector{T}, i_qdes::Int; use_log::Bool=true) where {T<:Real}

    η = 1.0 ./ cp1d.conductivity_parallel

    if i_qdes != 0
        # flattened current resistivity as per Jardin's model
        η[1:i_qdes] .= η[i_qdes]

        # flatten non-inductive current contribution
        rho0 = cp1d.grid.rho_tor_norm[i_qdes]
        width = min(rho0 / 4, 0.05)
        j_non_inductive = IMAS.flatten_profile!(copy(j_non_inductive), cp1d.grid.rho_tor_norm, cp1d.grid.area, rho0, width)
    end

    return QED.η_FE(cp1d.grid.rho_tor_norm, η; use_log), j_non_inductive
end
