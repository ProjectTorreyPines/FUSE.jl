deps = "ADAS BalanceOfPlantSurrogate BoundaryPlasmaModels CHEASE CoordinateConventions EPEDNN FiniteElementHermite IMASutils FusionMaterials FuseExchangeProtocol IMAS IMASdd MXHEquilibrium MeshTools MillerExtendedHarmonic NEO NNeutronics QED RABBIT SimulationParameters TEQUILA TGLFNN TJLF VacuumFields XSteam ThermalSystemModels"

txt = String["""
# Ecosystem

The FUSE project is built upon multiple Julia packages, many of which reside in the [https://github.com/ProjectTorreyPines](https://github.com/ProjectTorreyPines) organization on GitHub.

![FUSE dependencies](./assets/deps.svg)

"""]
for dep in sort!(split(deps))
    if dep in ("XSteam", "BoundaryPlasmaModels")
        continue
        push!(txt, "* $dep [[repo](https://github.com/ProjectTorreyPines/$dep.jl)]")
    elseif dep in ("ThermalSystemModels",)
        push!(txt, "* [$dep](https://projecttorreypines.github.io/$dep.jl/) [[repo](https://github.com/ProjectTorreyPines/$dep.jl)]")
    else
        push!(txt, "* [$dep](https://projecttorreypines.github.io/$dep.jl/dev/) [[repo](https://github.com/ProjectTorreyPines/$dep.jl)]")
    end
end

open("$(@__DIR__)/deps.md", "w") do io
    return write(io, join(txt, "\n"))
end
