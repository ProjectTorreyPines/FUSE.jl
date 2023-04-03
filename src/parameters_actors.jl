mutable struct ParametersActors{T} <: ParametersAllActors where {T<:Real}
    _parent::WeakRef
    _name::Symbol
    ActorCXbuild::FUSEparameters__ActorCXbuild{T}
    ActorFluxSwing::FUSEparameters__ActorFluxSwing{T}
    ActorLFSsizing::FUSEparameters__ActorLFSsizing{T}
    ActorHFSsizing::FUSEparameters__ActorHFSsizing{T}
    ActorStresses::FUSEparameters__ActorStresses{T}
    ActorEquilibrium::FUSEparameters__ActorEquilibrium{T}
    ActorSolovev::FUSEparameters__ActorSolovev{T}
    ActorCHEASE::FUSEparameters__ActorCHEASE{T}
    ActorPFcoilsOpt::FUSEparameters__ActorPFcoilsOpt{T}
    ActorPassiveStructures::FUSEparameters__ActorPassiveStructures{T}
    ActorBlanket::FUSEparameters__ActorBlanket{T}
    ActorBalanceOfPlant::FUSEparameters__ActorBalanceOfPlant{T}
    ActorPowerNeeds::FUSEparameters__ActorPowerNeeds{T}
    ActorThermalCycle::FUSEparameters__ActorThermalCycle{T}
    ActorHeatTransfer::FUSEparameters__ActorHeatTransfer{T}
    ActorQEDcurrent::FUSEparameters__ActorQEDcurrent{T}
    ActorSteadyStateCurrent::FUSEparameters__ActorSteadyStateCurrent{T}
    ActorDivertors::FUSEparameters__ActorDivertors{T}
    ActorNBIsimple::FUSEparameters__ActorNBIsimple{T}
    ActorECsimple::FUSEparameters__ActorECsimple{T}
    ActorICsimple::FUSEparameters__ActorICsimple{T}
    ActorLHsimple::FUSEparameters__ActorLHsimple{T}
    ActorTauenn::FUSEparameters__ActorTauenn{T}
    ActorCosting::FUSEparameters__ActorCosting{T}
    ActorNeutronics::FUSEparameters__ActorNeutronics{T}
    ActorNeoclassical::FUSEparameters__ActorNeoclassical{T}
    ActorPedestal::FUSEparameters__ActorPedestal{T}
    ActorTGLF::FUSEparameters__ActorTGLF{T}
    ActorCoreTransport::FUSEparameters__ActorCoreTransport{T}
    ActorTransportSolver::FUSEparameters__ActorTransportSolver{T}
    ActorEquilibriumTransport::FUSEparameters__ActorEquilibriumTransport{T}
    ActorWholeFacility::FUSEparameters__ActorWholeFacility{T}
    ActorStabilityLimits::FUSEparameters__ActorStabilityLimits{T}
    ActorBetaLimit::FUSEparameters__ActorBetaLimit{T}
    ActorCurrentLimit::FUSEparameters__ActorCurrentLimit{T}
    ActorDensityLimit::FUSEparameters__ActorDensityLimit{T}
end

function ParametersActors{T}() where {T<:Real}
    act = ParametersActors{T}(
        WeakRef(nothing),
        :act,
        FUSEparameters__ActorCXbuild{T}(),
        FUSEparameters__ActorFluxSwing{T}(),
        FUSEparameters__ActorLFSsizing{T}(),
        FUSEparameters__ActorHFSsizing{T}(),
        FUSEparameters__ActorStresses{T}(),
        FUSEparameters__ActorEquilibrium{T}(),
        FUSEparameters__ActorSolovev{T}(),
        FUSEparameters__ActorCHEASE{T}(),
        FUSEparameters__ActorPFcoilsOpt{T}(),
        FUSEparameters__ActorPassiveStructures{T}(),
        FUSEparameters__ActorBlanket{T}(),
        FUSEparameters__ActorBalanceOfPlant{T}(),
        FUSEparameters__ActorPowerNeeds{T}(),
        FUSEparameters__ActorThermalCycle{T}(),
        FUSEparameters__ActorHeatTransfer{T}(),
        FUSEparameters__ActorQEDcurrent{T}(),
        FUSEparameters__ActorSteadyStateCurrent{T}(),
        FUSEparameters__ActorDivertors{T}(),
        FUSEparameters__ActorNBIsimple{T}(),
        FUSEparameters__ActorECsimple{T}(),
        FUSEparameters__ActorICsimple{T}(),
        FUSEparameters__ActorLHsimple{T}(),
        FUSEparameters__ActorTauenn{T}(),
        FUSEparameters__ActorCosting{T}(),
        FUSEparameters__ActorNeutronics{T}(),
        FUSEparameters__ActorNeoclassical{T}(),
        FUSEparameters__ActorPedestal{T}(),
        FUSEparameters__ActorTGLF{T}(),
        FUSEparameters__ActorCoreTransport{T}(),
        FUSEparameters__ActorTransportSolver{T}(),
        FUSEparameters__ActorEquilibriumTransport{T}(),
        FUSEparameters__ActorWholeFacility{T}(),
        FUSEparameters__ActorStabilityLimits{T}(),
        FUSEparameters__ActorBetaLimit{T}(),
        FUSEparameters__ActorCurrentLimit{T}(),
        FUSEparameters__ActorDensityLimit{T}()
    )
    setup_parameters!(act)
    return act
end

function ParametersActors()
    return ParametersActors{Float64}()
end

"""
    act2json(act::ParametersAllActors, filename::AbstractString; kw...)

Save the FUSE parameters to a JSON file with give `filename`
`kw` arguments are passed to the JSON.print function
"""
function act2json(act::ParametersAllActors, filename::AbstractString; kw...)
    return SimulationParameters.par2json(act, filename; kw...)
end

function json2act(filename::AbstractString)
    return SimulationParameters.json2par(filename, ParametersActors())
end
