"""
    workflow_pedestal_transport(
        ini::ParametersAllInits,
        act::ParametersAllActors,
        w_ped::Real,
        Te_ped::Real;
        nrhos::Int=5)

Initializes a plasma with a given pedestal width `w_ped` and pedestal electron temperature
`Te_ped`, converges equilibrium and current, and evaluates the turbulent fluxes on `nrhos`
radial locations placed across the pedestal (from `1-0.75*w_ped` to `1-0.25*w_ped`).

Returns `(dd, Qe)` where `Qe` is the electron energy flux profile on the transport grid.
"""
function workflow_pedestal_transport(
    ini::ParametersAllInits,
    act::ParametersAllActors,
    w_ped::Real,
    Te_ped::Real;
    nrhos::Int=5)

    ini = deepcopy(ini)
    act = deepcopy(act)

    ini.core_profiles.w_ped = w_ped
    ini.core_profiles.Te_ped = Te_ped

    rhos = 1.0 .- range(0.75 * w_ped, 0.25 * w_ped; length=nrhos)
    act.ActorFluxCalculator.rho_transport = rhos
    act.ActorTGLF.rho_transport = rhos

    dd = init(ini, act)

    act.ActorCurrent.ip_from = :core_profiles
    act.ActorCurrent.vloop_from = :core_profiles

    ActorEquilibrium(dd, act; ip_from=:core_profiles)
    latest_equilibrium_grids!(dd)
    ActorCurrent(dd, act)

    ActorEquilibrium(dd, act; ip_from=:core_profiles)
    latest_equilibrium_grids!(dd)

    ActorFluxCalculator(dd, act)

    Qe = dd.core_transport.model[1].profiles_1d[1].electrons.energy.flux

    return dd, Qe
end

"""
    workflow_pedestal_transport_scan(
        ini::ParametersAllInits,
        act::ParametersAllActors;
        Tepeds::AbstractVector{<:Real},
        w_peds::AbstractVector{<:Real},
        nrhos::Int=5,
        ntasks::Int=Base.Threads.nthreads())

Scans [`workflow_pedestal_transport`](@ref) over the `w_peds` × `Tepeds` grid.

Returns a named tuple `(dds, Qes, w_peds, Tepeds, failures)` where:

  - `dds[iw, it]` is the `dd` of each successful run (`#undef` where the run failed)
  - `Qes[:, iw, it]` is the electron energy flux profile (zeros where the run failed)
  - `failures` lists the `(iw, it, exception)` of the runs that errored out

Runs are dispatched concurrently with `asyncmap`, which overlaps the runs while they are
blocked on external codes (e.g. TGLF); set `ntasks=1` to run serially.
"""
function workflow_pedestal_transport_scan(
    ini::ParametersAllInits,
    act::ParametersAllActors;
    Tepeds::AbstractVector{<:Real},
    w_peds::AbstractVector{<:Real},
    nrhos::Int=5,
    ntasks::Int=Base.Threads.nthreads())

    dds = Matrix{IMAS.dd}(undef, length(w_peds), length(Tepeds))
    Qes = zeros(nrhos, length(w_peds), length(Tepeds))
    failures = Tuple{Int,Int,Any}[]

    params = [(iw, it, w_ped, Te_ped) for (iw, w_ped) in enumerate(w_peds), (it, Te_ped) in enumerate(Tepeds)]

    results = asyncmap(params; ntasks) do (iw, it, w_ped, Te_ped)
        try
            @info "Running iw=$iw/$(length(w_peds)) it=$it/$(length(Tepeds))"
            dd, Qe = workflow_pedestal_transport(ini, act, w_ped, Te_ped; nrhos)
            return (iw, it, dd, Qe, nothing)
        catch e
            @warn "Failed for iw=$iw it=$it: $e"
            return (iw, it, nothing, nothing, e)
        end
    end

    for (iw, it, dd, Qe, err) in results
        if err === nothing
            dds[iw, it] = dd
            Qes[:, iw, it] = Qe
        else
            push!(failures, (iw, it, err))
        end
    end

    return (dds=dds, Qes=Qes, w_peds=collect(w_peds), Tepeds=collect(Tepeds), failures=failures)
end
