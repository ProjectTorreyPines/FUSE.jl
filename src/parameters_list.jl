function AllParameters()

    params = Parameters()

    equilibrium = params.equilibrium = Parameters()
    equilibrium.B0 = Entry(Real, missing, IMAS.equilibrium__vacuum_toroidal_field, :b0)
    equilibrium.R0 = Entry(Real, missing, IMAS.equilibrium__vacuum_toroidal_field, :r0)
    equilibrium.ϵ = Entry(Real, missing, "", "Plasma aspect ratio")
    equilibrium.δ = Entry(Real, missing, IMAS.equilibrium__time_slice___boundary, :triangularity)
    equilibrium.κ = Entry(Real, missing, IMAS.equilibrium__time_slice___boundary, :elongation)
    equilibrium.βn = Entry(Real, missing, IMAS.equilibrium__time_slice___global_quantities, :beta_normal)
    equilibrium.ip = Entry(Real, missing, IMAS.equilibrium__time_slice___global_quantities, :ip)
    equilibrium.x_point = Entry(Bool, missing, IMAS.equilibrium__time_slice___boundary, :x_point)

    build = params.build = Parameters()
    build.n_oh_coils = Entry(Int, missing, "", "number of OH coils")
    build.n_pf_coils_inside = Entry(Int, missing, "", "number of PF coils inside of the TF")
    build.n_pf_coils_outside = Entry(Int, missing, "", "number of PF coils outside of the TF")

    return params
end
