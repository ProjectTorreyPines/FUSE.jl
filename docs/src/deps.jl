# NOTE: to update the deps.svg file, run `make deps_dag` in the FUSE folder.

deps = "ADAS BalanceOfPlantSurrogate BoundaryPlasmaModels CHEASE CoordinateConventions EGGO EPEDNN FiniteElementHermite FRESCO FusionMaterials FuseExchangeProtocol GACODE HelpPlots IMAS IMASdd IMASutils MXHEquilibrium MillerExtendedHarmonic NeoclassicalTransport NNeutronics QED RABBIT SimulationParameters TEQUILA TJLF TORBEAM TroyonBetaNN VacuumFields TurbulentTransport ThermalSystemModels"

txt = String["""
# Ecosystem

The FUSE project is built upon multiple Julia packages, many of which reside in the [https://github.com/ProjectTorreyPines](https://github.com/ProjectTorreyPines) organization on GitHub.

![FUSE dependencies](./assets/deps.svg)

"""]
for dep in sort!(split(deps))
    push!(txt, "* [$dep](https://projecttorreypines.github.io/$dep.jl/) [[repo](https://github.com/ProjectTorreyPines/$dep.jl)]")
end

open("$(@__DIR__)/deps.md", "w") do io
    return write(io, join(txt, "\n"))
end
