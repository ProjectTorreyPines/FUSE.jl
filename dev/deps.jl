deps = "ADAS BalanceOfPlantSurrogate BoundaryPlasmaModels CHEASE CoordinateConventions EPEDNN FiniteElementHermite Fortran90Namelists FuseUtils FusionMaterials FuseExchangeProtocol IMAS IMASdd MXHEquilibrium MeshTools MillerExtendedHarmonic NEO NNeutronics QED RABBIT SimulationParameters TEQUILA TGLFNN TJLF VacuumFields XSteam ThermalSystemModels"

txt = String["""
# Dependencies

The FUSE project is built upon multiple Julia packages, many of which reside in the [https://github.com/ProjectTorreyPines](https://github.com/ProjectTorreyPines) organization on GitHub.

![FUSE dependencies](./assets/deps.svg)

"""]
for dep in sort!(split(deps))
    if dep in ("Fortran90Namelists", )
        continue
    end
    push!(txt, "* [$dep](https://projecttorreypines.github.io/$dep.jl/dev/) [[repo](https://github.com/ProjectTorreyPines/$dep.jl)]")
end

open("$(@__DIR__)/deps.md", "w") do io
    return write(io, join(txt, "\n"))
end
