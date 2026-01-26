#= ===================== =#
#  ActorCoreRadHeatFlux  #
#= ===================== =#
@actor_parameters_struct ActorCoreRadHeatFlux{T} begin
    N::Entry{Int} = Entry{Int}("-", "Number of launched photons"; default=100000)
    r::Entry{Vector{T}} = Entry{Vector{T}}("m", "Vector of r at outermidplane"; default=T[])
    q::Entry{Vector{T}} = Entry{Vector{T}}("W m^-2", "Vector of parallel power density at outer midplane"; default=T[])
    levels::Entry{Union{Int,Vector}} =
        Entry{Union{Int,Vector}}("-", "If Int it defines number of levels in SOL, if vector it corresponds to the psi levels to build SOL"; default=20)
    merge_wall::Entry{Bool} = Entry{Bool}("-", "Merge dd.wall in mesh for the heat flux "; default=true)
    step::Entry{T} = Entry{T}("m", " Step for discretization of the default wall mesh (dd.wall)"; default=0.1)
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
end

mutable struct ActorCoreRadHeatFlux{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorCoreRadHeatFlux{P}}
    wall_heat_flux::Union{Nothing,IMAS.WallHeatFlux}
end

function ActorCoreRadHeatFlux(dd::IMAS.dd{D}, par::FUSEparameters__ActorCoreRadHeatFlux{P}; kw...) where {D<:Real,P<:Real}
    logging_actor_init(ActorCoreRadHeatFlux)
    par = OverrideParameters(par; kw...)
    return ActorCoreRadHeatFlux(dd, par, nothing)
end

"""
    ActorCoreRadHeatFlux(dd::IMAS.dd, act::ParametersAllActors; kw...)

Computes wall heat flux from core plasma radiation using Monte Carlo photon transport.

This actor calculates the heat flux deposited on plasma-facing components by electromagnetic
radiation (bremsstrahlung, line radiation, synchrotron) emitted from the core plasma.
Core radiation represents a significant fraction of total plasma power and directly impacts
wall loading and cooling requirements.

**Physical Modeling:**
- Volume emission from all core radiation sources (impurity line radiation, bremsstrahlung)
- Line-of-sight photon transport from emission location to wall intersection
- Accounts for plasma opacity and absorption effects
- Integration over entire plasma volume and all wall surfaces

**Computational Method:**
1. **Source Definition:** Uses radiation source profiles from `dd.core_sources` including
   contributions from impurity radiation, bremsstrahlung, and other radiative losses

2. **Monte Carlo Transport:** Launches statistical photons from emission locations with:
   - Proper spatial weighting based on local radiation source strength
   - Random directional sampling to capture all emission angles
   - Line-of-sight transport to first wall intersection

3. **Wall Deposition:** Accumulates photon energy deposition on wall mesh elements,
   accounting for geometric factors and surface orientations

**Key Outputs:**
- `q_core_rad`: Heat flux on wall surface from core radiation transport
- Statistical convergence based on number of photons traced (parameter N)
- Wall mesh coordinates for spatial distribution analysis

This calculation provides essential data for first wall and blanket thermal design,
complementing particle heat flux calculations for comprehensive wall loading analysis.
"""
function ActorCoreRadHeatFlux(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorCoreRadHeatFlux(dd, act.ActorCoreRadHeatFlux; kw...)
    step(actor)
    finalize(actor)
    return actor
end

"""    _step(actor::ActorCoreRadHeatFlux)

Executes Monte Carlo photon transport calculation for core radiation heat flux.

The calculation procedure:

1. **Wall Mesh Setup:** Creates the same computational wall mesh used for particle heat flux,
   ensuring consistency between different wall loading calculations.

2. **Radiation Source Processing:** Extracts radiation source profiles from `dd.core_sources`:
   - Computes total volumetric radiation source density (W/m³)
   - Calculates integrated core radiation power for normalization
   - Uses negative source values (losses for core_sources convention)

3. **Monte Carlo Photon Transport:** 
   - Launches N statistical photons from plasma volume using IMAS.core_radiation_heat_flux
   - Each photon carries appropriate statistical weight based on local source strength
   - Traces line-of-sight paths from emission point to wall intersection
   - Accumulates energy deposition on wall mesh elements

4. **Result Assembly:** Creates WallHeatFlux object containing:
   - Wall mesh coordinates (R, Z) and curvilinear distance
   - Core radiation heat flux distribution `q_core_rad`

5. **Optional Visualization:** If enabled, produces plots showing:
   - 1D heat flux profile along the wall perimeter
   - 2D visualization with photon trajectories and SOL boundary
   - Statistical photon launch locations in the plasma

Statistical convergence improves with larger N (number of photons traced).
"""
function _step(actor::ActorCoreRadHeatFlux{D,P}) where {D<:Real,P<:Real}
    dd = actor.dd
    par = actor.par

    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]

    # Compute wall mesh - same for both core radiation and particles (see ActorParticleHeatFlux)
    # (Rwall, Zwall)        wall mesh (with intersections of SOL) -  m
    # s                     curvilinear abscissa computed from (Rwall, Zwall), clockwise starting at OMP - m
    Rwall, Zwall, s, _ = IMAS.mesher_heat_flux(dd; par.r, par.q, par.merge_wall, par.levels, par.step)

    # Parameters for heat flux due to core radiarion
    total_rad_source1d = IMAS.total_radiation_sources(dd.core_sources, cp1d; time0=dd.global_time)
    source_1d = -total_rad_source1d.electrons.energy # minus sign because loss for dd.core_sources
    Prad_core = -total_rad_source1d.electrons.power_inside[end]

    # Compute the heat flux due to the core radiation - q_core_rad - W/m2
    q_core_rad = IMAS.core_radiation_heat_flux(eqt, cp1d.grid.psi_norm, source_1d, par.N, Rwall, Zwall, Prad_core)

    HF = IMAS.WallHeatFlux{D}(; r=Rwall, z=Zwall, q_core_rad, s)
    actor.wall_heat_flux = HF

    #plot
    if par.do_plot
        ll = @layout [a{0.6w,0.9h} b{0.4w}]
        p = plot(; layout=ll, size=(1500, 500))
        plot!(p, HF; which_plot=:oneD, q=:core_radiation, subplot=1)
        sol = IMAS.sol(dd; levels=1)
        photons, _ = IMAS.define_particles(eqt, cp1d.grid.psi_norm, source_1d, par.N)
        plot!(p, photons, actor.dd.equilibrium.time_slice[]; colorbar_entry=false, subplot=2)
        plot!(p, sol; subplot=2, line_z=nothing, color=:black)
        plot!(p, HF; q=:core_radiation, plot_type=:path, subplot=2)
        display(p)
    end

    return actor
end

"""    _finalize(actor::ActorCoreRadHeatFlux)

Finalizes core radiation heat flux calculation for integration with other analyses.

This function is currently a placeholder for future development that will:
1. Store core radiation heat flux results in appropriate IMAS data structures
2. Combine with particle heat flux and other wall loading contributions
3. Provide integrated wall loading data for thermal and structural analysis
4. Enable coupling with neutronics and materials activation calculations

The core radiation heat flux results are currently accessible through the actor's
`wall_heat_flux` field and can be visualized using the built-in plotting capabilities.
"""
function _finalize(actor::ActorCoreRadHeatFlux)
    # Finalize: work in progress - populate dd
    return actor
end