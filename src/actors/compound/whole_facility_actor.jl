#= ================== =#
#  ActorWholeFacility  #
#= ================== =#
@actor_parameters_struct ActorWholeFacility{T} begin
    update_plasma::Entry{Bool} = Entry{Bool}("-", "Run plasma related actors"; default=true)
    update_build::Entry{Bool} = Entry{Bool}("-", "Optimize tokamak build"; default=true)
end

mutable struct ActorWholeFacility{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorWholeFacility{P}}
    act::ParametersAllActors
    StationaryPlasma::Union{Nothing,ActorStationaryPlasma{D,P}}
    StabilityLimits::Union{Nothing,ActorPlasmaLimits{D,P}}
    FluxSwing::Union{Nothing,ActorFluxSwing{D,P}}
    Stresses::Union{Nothing,ActorStresses{D,P}}
    HFSsizing::Union{Nothing,ActorHFSsizing{D,P}}
    LFSsizing::Union{Nothing,ActorLFSsizing{D,P}}
    CXbuild::Union{Nothing,ActorCXbuild{D,P}}
    PFdesign::Union{Nothing,ActorPFdesign{D,P}}
    PFactive::Union{Nothing,ActorPFactive{D,P}}
    PassiveStructures::Union{Nothing,ActorPassiveStructures{D,P}}
    VerticalStability::Union{Nothing,ActorVerticalStability{D,P}}
    Neutronics::Union{Nothing,ActorNeutronics{D,P}}
    Blanket::Union{Nothing,ActorBlanket{D,P}}
    Divertors::Union{Nothing,ActorDivertors{D,P}}
    BalanceOfPlant::Union{Nothing,ActorBalanceOfPlant{D,P}}
    Costing::Union{Nothing,ActorCosting{D,P}}
end

"""
    ActorWholeFacility(dd::IMAS.dd, act::ParametersAllActors; kw...)

Integrates all physics, engineering, and economic models for complete tokamak facility design.

This compound actor orchestrates the full workflow for designing and analyzing a complete
tokamak power plant. It coordinates plasma physics modeling with engineering design and
economic optimization to provide a comprehensive facility assessment.

Execution workflow:
1. **Plasma Design**: Self-consistent plasma solutions with optimized PF coil positions
2. **Plasma Physics**: Stationary plasma equilibrium with all transport and heating
3. **Engineering Design**: Radial build optimization for all major tokamak systems
4. **Magnetic Systems**: PF coil repositioning and passive structure analysis  
5. **Nuclear Analysis**: Neutronics, blanket optimization, and tritium breeding
6. **Plasma-Facing Components**: Divertor heat loads and thermal management
7. **Power Systems**: Balance of plant and thermal-to-electric conversion
8. **Economics**: Complete facility costing and economic analysis

Key actors executed:
- **ActorStationaryPlasma**: Converged plasma equilibrium and transport solution
- **ActorPlasmaLimits**: Operational stability limits and MHD analysis
- **ActorHFSsizing**: High-field-side radial build optimization (TF, OH, structural)
- **ActorLFSsizing**: Low-field-side TF leg design for ripple and maintenance
- **ActorCXbuild**: 2D cross-sectional geometry generation
- **ActorFluxSwing**: Flux swing analysis for plasma scenarios
- **ActorStresses**: Structural and electromagnetic stress analysis
- **ActorPFdesign**: PF coil system design and optimization
- **ActorPFactive**: PF current optimization for equilibrium control
- **ActorPassiveStructures**: Passive conductor modeling
- **ActorNeutronics**: Neutron flux and nuclear heating analysis
- **ActorBlanket**: Tritium breeding blanket optimization
- **ActorDivertors**: Divertor heat flux and cooling analysis
- **ActorBalanceOfPlant**: Power conversion and auxiliary systems
- **ActorCosting**: Economic analysis and cost optimization

The actor handles complex interdependencies between systems, such as re-optimizing
the PF coil layout after radial build changes and updating the first wall geometry
based on the final equilibrium solution.

!!! note

    Stores comprehensive facility design data across all `dd` structures
"""
function ActorWholeFacility(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorWholeFacility(dd, act.ActorWholeFacility, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorWholeFacility(dd::IMAS.dd, par::FUSEparameters__ActorWholeFacility, act::ParametersAllActors; kw...)
    logging_actor_init(ActorWholeFacility)
    par = OverrideParameters(par; kw...)

    return ActorWholeFacility(dd, par, act,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing)
end

function _step(actor::ActorWholeFacility)
    dd = actor.dd
    par = actor.par
    act = deepcopy(actor.act)

    if !isempty(dd.build.layer) && par.update_build
        # ActorPFdesign optimizing optimizes the coils location, to make sure that ActorStationaryPlasma will be able to keep the equilibrium in place.
        # Internally, ActorPFdesign calls ActorPFactive, which finds the currents that satisfy the (boundary, x-points, strike-points) constraints for a given coils configuration
        actor.PFdesign = ActorPFdesign(dd, act)
    end

    if par.update_plasma
        # ActorStationaryPlasma iterates between ActorCurrent, ActorHCD, ActorCoreTransport, and ActorEquilibrium to find a self-consistent stationary plasma solution
        actor.StationaryPlasma = ActorStationaryPlasma(dd, act)
    end

    if isempty(dd.build.layer)
        @warn "ActorWholeFacility: skipping engineering/costing actors since build is missing"

    else
        if par.update_build
            # ActorHFSsizing changes the radial build of the center stack (plug, OH, and TF)
            # accounting for stresses, superconductors critical currents, flux swing, and field requirements
            actor.HFSsizing = ActorHFSsizing(dd, act)

            # ActorLFSsizing changes the radial build of the outer TF leg to satisfy requirement of ripple and maintenance ports.
            actor.LFSsizing = ActorLFSsizing(dd, act)

            # ActorCXbuild takes the radial build information and generates 2D cross-section outlines of the build layers.
            # At this point we let these outlines ingnore the position of the pf coils, since we'll regenerate them.
            actor.CXbuild = ActorCXbuild(dd, act; layers_aware_of_pf_coils = false)

            # We now re-position the PF coils based on the new 2D build.
            # If we plan on rebuilding the first wall with ActorCXbuild, let the strike-points do whatever they naturally want to do.
            strike_points_weight_bkp = act.ActorPFactive.strike_points_weight
            if act.ActorCXbuild.rebuild_wall
                act.ActorPFactive.strike_points_weight = 0.0
            end
            actor.PFdesign = ActorPFdesign(dd, act)
            act.ActorPFactive.strike_points_weight = strike_points_weight_bkp

            # We now re-solve the equilibrium with the new coils positions
            ActorEquilibrium(dd, act; ip_from=:pulse_schedule)

            # and we update the first wall based on the new equilibrium
            if act.ActorCXbuild.rebuild_wall
                actor.CXbuild = ActorCXbuild(dd, act)
            end

        else
            # Here we just make sure to populate the dd with all the information that creating a new build would give.
            actor.FluxSwing = ActorFluxSwing(dd, act)
            actor.Stresses = ActorStresses(dd, act)
            actor.PFactive = ActorPFactive(dd, act)
        end

        # ActorNeutronics evaluates the neutron flux on the first wall
        actor.Neutronics = ActorNeutronics(dd, act)

        # ActorBlanket optimizes the radial build thickness of the first wall, blanket, shield and Li6 enrichment to achieve a target TBR
        actor.Blanket = ActorBlanket(dd, act)
        # We must re-generate the CX build since we updated the radial build
        actor.CXbuild = ActorCXbuild(dd, act)

        # ActorPassiveStructures populates dd.pf_passive based on the vacuum vessel layer(s)
        actor.PassiveStructures = ActorPassiveStructures(dd, act)

        # ActorDivertors evaluates divertor fluxes using reduced SOL models
        actor.Divertors = ActorDivertors(dd, act)

        # ActorPlasmaLimits evaluates opearational plasma limits, including vertical stability and low-n ideal MHD stability
        actor.StabilityLimits = ActorPlasmaLimits(dd, act)

        # ActorBalanceOfPlant evaluates power needs and optimal thermal-to-electric power conversion efficiency
        actor.BalanceOfPlant = ActorBalanceOfPlant(dd, act)

        # ActorCosting costs the plant
        actor.Costing = ActorCosting(dd, act)
    end

    return actor
end
