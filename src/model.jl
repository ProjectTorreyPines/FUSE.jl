using Pkg
Pkg.activate("./FUSE")
using FUSE
using FUSE.IMAS

ϵ = 0.32
κ = 1.7
δ = 0.33
R0 = 6.2
B0 = 5.3
ip = 15.E6
beta_n = 12.0

resolution = 33

dd=IMAS.dd()
FUSE.init(dd.equilibrium, 0.0; B0, R0, ϵ, δ, κ, beta_n, ip, x_point=false)

eqactor = FUSE.SolovevEquilibriumActor(dd.equilibrium, 0.0, alpha=0, qstar=1.5)
FUSE.step(eqactor)
dd.equilibrium = FUSE.finalize(eqactor, resolution)



# #= ============== =#
# # initialization #
# #= ============== =#

# PP = FUSE.fuse_parameters[:PLASMA_PARAMETERS]
# PM = FUSE.fuse_parameters[:PHYSICS_MODELS]

# core1D = FUSE.core_profiles__profiles_1d()

# n = 11

# core1D.grid.rho_tor_norm = range(0.0, 1.0, length=n)
# core1D.electrons.density = (rho_tor_norm;_...) -> PP[:ne0] .* (1.0 .- rho_tor_norm.^2).^PP[:Sn]
# core1D.electrons.temperature = (rho_tor_norm;_...) -> PP[:Te0] .* (1.0 .- rho_tor_norm.^2).^PP[:St]

# core1D.j_total = (x;_...) -> (1.0 .- x.^2).^PP[:Sj]

# equil = FUSE.equilibrium__time_slice()
# equil.profiles_1d.psi = range(0.0, 1.0, length=n)
# equil.profiles_1d.elongation = (psi;_...) -> psi .* 0.0 .+ PP[:elongation]
# equil.profiles_1d.geometric_axis.r = (psi;_...) -> psi .* 0.0 .+ PP[:R0]

# aspect_ratio(R, a) = R./a

# function aspect_ratio(fds::FDS)
#     R = fds["equilibirum.time_slice[:].profiles_1d.geometric_axis.r"]
#     a = fds["equilibirum.time_slice[:].profiles_1d.a"]
#     aspect_ratio(R,a)
# end

# bootstrapCoefficient = FUSE.collisionless_bootstrap(PM[:bootstrapModel], PP[:elongation], PP[:St], PP[:Sn], PP[:Sj], PP[:Zeff])
# println(bootstrapCoefficient)

# # maxStableElongation(aspectRatio) = 2.43 + 65.0 * exp(-aspectRatio / 0.376)

# # #elongation_fraction = elongation / maxStableElongation