#= ========================= =#
#  ActorTearingModeTransport  #
#= ========================= =#
@actor_parameters_struct ActorTearingModeTransport{T} begin
    m::Entry{Vector{Int}} = Entry{Vector{Int}}("-", "Poloidal mode numbers m of the tearing modes"; default=Int[])
    n::Entry{Vector{Int}} = Entry{Vector{Int}}("-", "Toroidal mode numbers n of the tearing modes"; default=Int[])
    rho_res::Entry{Vector{Float64}} = Entry{Vector{Float64}}("-", "rho_tor_norm of the resonant surfaces (empty or NaN = find abs(q) = m/n)"; default=Float64[])
    island_width::Entry{Vector{Float64}} = Entry{Vector{Float64}}("-", "Full island widths in rho_tor_norm (0.0 = mode off)"; default=Float64[])
    w_c::Entry{T} = Entry{T}("-", "Critical island width in rho_tor_norm: smooths the island edges and turns off islands smaller than w_c"; default=0.02, check=x -> @assert x > 0.0 "w_c must be > 0.0")
    chi_e::Entry{T} = Entry{T}("m^2/s", "Electron heat diffusivity inside the islands"; default=20.0, check=x -> @assert x >= 0.0 "chi_e must be >= 0.0")
    chi_i::Entry{T} = Entry{T}("m^2/s", "Ion heat diffusivity inside the islands"; default=20.0, check=x -> @assert x >= 0.0 "chi_i must be >= 0.0")
    D::Entry{T} = Entry{T}("m^2/s", "Particle diffusivity inside the islands"; default=20.0, check=x -> @assert x >= 0.0 "D must be >= 0.0")
    chi_tor::Entry{T} = Entry{T}("m^2/s", "Momentum diffusivity inside the islands"; default=20.0, check=x -> @assert x >= 0.0 "chi_tor must be >= 0.0")
    rho_transport::Entry{AbstractVector{T}} = Entry{AbstractVector{T}}("-", "rho_tor_norm values to compute the island fluxes on"; default=0.25:0.1:0.85)
end

mutable struct ActorTearingModeTransport{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.DD{D}
    par::OverrideParameters{P,FUSEparameters__ActorTearingModeTransport{P}}
    rho_res::Vector{Float64}
    window::Vector{D}
    function ActorTearingModeTransport(dd::IMAS.DD{D}, par::FUSEparameters__ActorTearingModeTransport{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorTearingModeTransport)
        par = OverrideParameters(par; kw...)
        return new{D,P}(dd, par, Float64[], D[])
    end
end

"""
    ActorTearingModeTransport(dd::IMAS.DD, act::ParametersAllActors; kw...)

Models the effect of tearing modes (magnetic islands) on transport with ad hoc diffusivities inside the islands.

Inside an island, parallel transport along the helical field lines short-circuits radial transport.
This is mimicked by adding the diffusivities `chi_e`, `chi_i`, `D` and `chi_tor` over the island width
(see `IMAS.island_window` and `IMAS.tearing_mode_transport!`). When these are much larger than the turbulent
diffusivities, the flux matcher drives the gradients inside the island below the turbulence threshold and
the profiles flatten. The bootstrap hole and the confinement degradation then follow self-consistently.

The fluxes are stored in `dd.core_transport.model[:ntm]` and are evaluated by `ActorFluxCalculator`
together with the turbulent and neoclassical fluxes. Without tearing modes, `dd.core_transport.model[:ntm]` is removed.

!!! note
    `ActorFluxMatcher.rho_transport` must have points inside the islands (eg. at `rho_res` and `rho_res ± w/4`),
    otherwise the islands are not seen by the flux matcher.

# Parameters
- `m`, `n`: Mode numbers, used to locate the resonant surfaces
- `rho_res`: Resonant surface locations. If empty or `NaN`, the outermost surface where abs(q) = m/n is used
- `island_width`: Full island widths in rho_tor_norm (0.0 turns off the mode)
- `w_c`: Critical island width in rho_tor_norm
- `chi_e`, `chi_i`, `D`, `chi_tor`: Diffusivities inside the islands [m²/s]

For example, a 2/1 island of width 0.08 at the q=2 surface: `m=[2], n=[1], island_width=[0.08]`
"""
function ActorTearingModeTransport(dd::IMAS.DD, act::ParametersAllActors; kw...)
    actor = ActorTearingModeTransport(dd, act.ActorTearingModeTransport; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorTearingModeTransport{D,P}) where {D<:Real,P<:Real}
    dd = actor.dd
    par = actor.par
    cp1d = dd.core_profiles.profiles_1d[]
    rho = cp1d.grid.rho_tor_norm

    n_modes = length(par.m)
    @assert length(par.n) == n_modes "n must have the same length as m"
    @assert length(par.island_width) == n_modes "island_width must have the same length as m"
    @assert isempty(par.rho_res) || length(par.rho_res) == n_modes "rho_res must be empty or have the same length as m"

    actor.rho_res = fill(NaN, n_modes)
    actor.window = zeros(D, length(rho))
    for k in 1:n_modes
        rho_res = isempty(par.rho_res) ? NaN : par.rho_res[k]
        if isnan(rho_res)
            rho_q = IMAS.rho_q_surface(rho, cp1d.q, par.m[k] / par.n[k])
            rho_q === nothing && continue
            rho_res = rho_q
        end
        actor.rho_res[k] = rho_res
        actor.window .+= IMAS.island_window(rho, rho_res, par.island_width[k], par.w_c)
    end
    # overlapping islands do not add up
    actor.window .= min.(actor.window, 1.0)

    return actor
end

function _finalize(actor::ActorTearingModeTransport)
    dd = actor.dd
    par = actor.par

    if isempty(par.m)
        index = findfirst(model -> model.identifier.index == IMAS.name_2_index(dd.core_transport.model)[:ntm], dd.core_transport.model)
        index === nothing || deleteat!(dd.core_transport.model, index)
        return actor
    end

    cp1d = dd.core_profiles.profiles_1d[]
    eqt = dd.equilibrium.time_slice[]

    model = resize!(dd.core_transport.model, :ntm; wipe=false)
    model.identifier.name = join(["$(m)/$(n)" for (m, n) in zip(par.m, par.n)], ", ")
    m1d = resize!(model.profiles_1d)
    IMAS.tearing_mode_transport!(m1d, cp1d, par.rho_transport, actor.window; par.chi_e, par.chi_i, par.D, par.chi_tor, R0=eqt.boundary.geometric_axis.r)

    return actor
end
