#= ===================== =#
#  ActorParticleHeatFlux  #
#= ===================== =#

@actor_parameters_struct ActorParticleHeatFlux{T} begin
    r::Entry{Vector{T}} = Entry{Vector{T}}("m", "Vector of r at outermidplane"; default=T[])
    q::Entry{Vector{T}} = Entry{Vector{T}}("W m^-2", "Vector of parallel power density at outer midplane"; default=T[])
    levels::Entry{Union{Int,Vector{T}}} =
        Entry{Union{Int,Vector{T}}}("-", "If Int it defines number of levels in SOL, if vector it corresponds to the psi levels to build SOL"; default=20)
    merge_wall::Entry{Bool} = Entry{Bool}("-", "Merge dd.wall in mesh for the heat flux "; default=true)
    step::Entry{T} = Entry{T}("m", " Step for discretization of the default wall mesh (dd.wall)"; default=0.1)
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
end

mutable struct ActorParticleHeatFlux{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorParticleHeatFlux{P}}
    wall_heat_flux::Union{Nothing,IMAS.WallHeatFlux}
end

function ActorParticleHeatFlux(dd::IMAS.dd{D}, par::FUSEparameters__ActorParticleHeatFlux{P}; kw...) where {D<:Real,P<:Real}
    logging_actor_init(ActorParticleHeatFlux)
    par = OverrideParameters(par; kw...)
    return ActorParticleHeatFlux(dd, par, nothing)
end

"""
    ActorParticleHeatFlux(dd::IMAS.dd, act::ParametersAllActors; kw...)

Computes wall heat flux from charged particle transport across the scrape-off layer (SOL).

This actor calculates the heat flux deposited on plasma-facing components by charged particles
(ions and electrons) flowing along magnetic field lines from the plasma edge to the wall.
The calculation accounts for:

**Physical Modeling:**
- Parallel heat conduction along magnetic field lines in the SOL
- Power density decay profiles from outer midplane to divertor/limiters
- Geometric projection effects from field-line-following to wall-normal flux
- Integration over the entire first wall and divertor surfaces

**Computational Approach:**
1. Creates a computational mesh of the first wall incorporating SOL flux surfaces
2. Traces open field lines from the separatrix to wall intersection points  
3. Applies prescribed power decay profiles along field lines
4. Computes both parallel (field-line-following) and perpendicular (wall-normal) heat fluxes

**Key Outputs:**
- `q_part`: Heat flux perpendicular to wall surface from particle transport
- `q_parallel`: Heat flux parallel to magnetic field lines  
- Wall mesh coordinates and curvilinear distance along the wall

The results provide essential input for material design, cooling system requirements,
and component lifetime assessments in tokamak reactor studies.
"""
function ActorParticleHeatFlux(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorParticleHeatFlux(dd, act.ActorParticleHeatFlux; kw...)
    step(actor)
    finalize(actor)
    return actor
end

"""    _step(actor::ActorParticleHeatFlux)

Calculates charged particle heat flux on the wall using SOL field line tracing.

The calculation procedure:

1. **Wall Mesh Generation:** Creates a computational mesh of the first wall that incorporates
   intersections with SOL flux surfaces, using either prescribed power decay profiles or
   default parameters if not specified.

2. **SOL Field Line Tracing:** Traces open magnetic field lines from the separatrix through
   the SOL to their intersection points with the wall mesh.

3. **Heat Flux Calculation:** Computes particle heat flux using the IMAS.particle_heat_flux
   function which:
   - Applies power decay profiles along field lines
   - Accounts for geometric projection from parallel to perpendicular flux components
   - Integrates contributions over all field lines intersecting each wall element

4. **Result Storage:** Creates a WallHeatFlux object containing wall coordinates, heat flux
   distributions, and curvilinear distance along the wall.

5. **Optional Visualization:** If plotting is enabled, generates 1D and 2D visualizations
   showing heat flux distribution along the wall and in the SOL region.

The results provide both perpendicular (wall-normal) and parallel heat flux components.
"""
function _step(actor::ActorParticleHeatFlux{D,P}) where {D<:Real, P<:Real}
    dd = actor.dd
    par = actor.par

    # Compute wall mesh - same for both core radiation and particles (see ActorCoreRadHeatFlux)
    # (Rwall, Zwall)        wall mesh (with intersections of SOL) -  m
    # s                     curvilinear abscissa computed from (Rwall, Zwall), clockwise starting at OMP - m
    # SOL                   list of OpenFieldLines used to compute (Rwall, Zwall) 
    # (r,q)                 Hypothesis of power density decay at omp for definition of SOL
    Rwall, Zwall, s, SOL, r, q = IMAS.mesher_heat_flux(dd; par.r, par.q, par.merge_wall, par.levels, par.step)

    # Compute the heat flux due to the influx of charged particles
    # q_part                 Heat flux due to particles perpendicular to the wall       - W/m2
    # q_parallel             Heat flux due to particles parallel to the magnetic field  - W/m2
    q_part, q_parallel = IMAS.particle_heat_flux(SOL, Rwall, Zwall, r, q)

    HF = IMAS.WallHeatFlux{D}(; r=Rwall, z=Zwall, q_part, q_parallel, s)
    actor.wall_heat_flux = HF

    #plot
    if par.do_plot
        ll = @layout [a{0.6w,0.9h} b{0.4w}]
        p = plot(; layout=ll, size=(1500, 500))
        plot!(p, HF; which_plot=:oneD, q=:particle, subplot=1)
        sol = IMAS.sol(dd; levels=100)
        plot!(p, sol; subplot=2, colorbar_entry=false, color=:grays)
        plot!(p, HF; q=:particle, plot_type=:path, subplot=2)
        display(p)
    end

    return actor
end

# Populate dd
"""    _finalize(actor::ActorParticleHeatFlux)

Finalizes particle heat flux calculation and prepares for integration with other physics.

This function is currently a placeholder for future development that will:
1. Store calculated heat flux results in the appropriate IMAS data structure locations
2. Integrate with other wall loading calculations (e.g., core radiation, neutron heating)
3. Provide data for downstream materials and cooling system analyses

The heat flux results are currently stored in the actor's `wall_heat_flux` field and can be
accessed directly or through the plotting functionality.
"""
function _finalize(actor::ActorParticleHeatFlux)
    # Finalize: work in progress - populate dd
    return actor
end