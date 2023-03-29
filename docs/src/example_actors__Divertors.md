# Divertors


```@julia
using Revise
using FUSE
using Plots;
FUSE.logging(Logging.Info);
```

### Initialize the ITER case case
[ITER case documentation](https://fuse.help/cases.html#ITER)


```@julia
dd, ini, act = FUSE.init(:ITER, init_from=:ods, do_plot=true);
```

### Run Actors that will be needed for the Divertors


```@julia
FUSE.ActorEquilibriumTransport(dd, act);
FUSE.ActorCXbuild(dd, act)
FUSE.ActorNeutronics(dd, act; do_plot=true);
```

### Running the divertors actor
[ActorDivertors documentation](https://fuse.help/actors.html#Divertors)


```@julia
FUSE.ActorDivertors(dd, act)
dd.divertors
```

### Some divertor metrics that should be implemented


```@julia
# Divertor gasc

IMAS.widthSOL_eich(dd)

# Outputs: 
#        OUT["divertor metrics"] dict:
#        "widthSOL"
#        "PR",
#        "PBR",
#        "PBpR",
#        "heatFluxParallel",
#        "heatFluxPoloidal",
#        "divDeliveredHeatFlux",
#        "qdivPeak",
#        "divRadFraction",

eq = dd.equilibrium
eqt = eq.time_slice[]
eq1d = eqt.profiles_1d
cp1d = dd.core_profiles.profiles_1d[]

major_radius = eqt.boundary.geometric_axis.r
minor_radius = eqt.boundary.minor_radius
aspect_ratio = major_radius / minor_radius
power_SOL = IMAS.total_power_source(IMAS.total_sources(dd))
power_SOL = 219.9e6
Bpol_average = eqt.global_quantities.ip * (4.0 * pi * 1e-7) / eqt.global_quantities.length_pol

widthSOL = 1.35e-3 * (power_SOL / 1e6)^(-0.02) * major_radius^0.04 * Bpol_average^(-0.92) * aspect_ratio^(-0.42) # Eich scaling (NF 53 093031)

divRadFraction = 0.004482781088966957  # is small just took the GASC calculated value

divOBQFraction = 0.8 # Fraction of total power_SOL directed to the outer strike point(s) 0.8 is gasc assumption

if length(dd.divertors.divertor) == 2
    divVertQFraction = 0.5  # Fraction of total power_SOL directed to the upper divertor (assumed to be the divertor with the larger heat flux due to gradB drifts, slight unbalance in dRsep, etc.)
else
    divVertQFraction = 1.0
end

divPoloidalFluxExpansion = 1.0 # "Poloidal flux expansion factor between outboard midplane and divertor targets", gasc_standard = 1.0
divPoloidalFieldLineAngle = 10 # "Angle between poloidal field and divertor plate at footprint" gasc standard = 10 deg
divTotFieldLineAngle = 3.0# "Total tilt angle between edge field and divertor plate at footprint" default = 3.0 deg
diverter_wetted_area = 2 * pi * (major_radius - 0.5 * minor_radius) * widthSOL * divPoloidalFluxExpansion / sin(divPoloidalFieldLineAngle / 180.0 * pi)
divRadArea = 2 * pi * (major_radius - 0.5 * minor_radius) * minor_radius

qdivPeak = power_SOL * divVertQFraction * divOBQFraction * ((1 - divRadFraction) / diverter_wetted_area + divRadFraction / divRadArea)



# calc heat fluxes at outboard midplane
heatFluxPoloidal = power_SOL / (2.0 * pi * (major_radius + minor_radius) * length(dd.divertors.divertor) * widthSOL)
heatFluxParallel = heatFluxPoloidal * eqt.global_quantities.magnetic_axis.b_field_tor / Bpol_average


# calculate unmitigated heat flux at divertor target (no radiation, no flux expansion)
divDeliveredHeatFlux = heatFluxParallel * sin(divTotFieldLineAngle * pi / 180.0)

@show widthSOL, power_SOL / 1e6, qdivPeak / 1e6, divDeliveredHeatFlux / 1e6, "MW/m^2"
@show diverter_wetted_area * divDeliveredHeatFlux / 1e6
"""
     "power_SOL": 219.92776950858612,
     "widthSOL": 0.0009726739483561888,
     "PR": 45.58624844119182,
     "PBR": 214.85581559316367,
     "PBpR": 34.99428003654665,
     "heatFluxParallel": 17809.933014915052,
     "heatFluxPoloidal": 2900.7629215688094,
     "divDeliveredHeatFlux": 932.0998749583548,
     "qdivPeak": 601.7560074909352,
     "divRadFraction": 0.004482781088966957

    "divRadFraction": 0.004482781088966957
     "magneticPoloidalField": 0.7676499214821493,

""";

# divDeliveredHeatFlux, heatFluxParallel, 

```
