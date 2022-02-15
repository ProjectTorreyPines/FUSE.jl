top_level_parameters = [:general, :equilibrium, :coil, :build, :gasc, :ods]

function Parameters()
    params = Parameters(Dict{Symbol,Union{Parameter,Parameters}}())
    for item in top_level_parameters
        setproperty!(params, item, Parameters(item))
    end
    return params
end

function Parameters(what::Symbol)
    if what in top_level_parameters
        params = Parameters(Dict{Symbol,Union{Parameter,Parameters}}())

        if what == :general
            params.init_from = Entry(Symbol, missing, "", "Initialize run from") # [:ods, :scalars, :gasc]

        elseif what == :equilibrium
            params.B0 = Entry(Real, missing, IMAS.equilibrium__vacuum_toroidal_field, :b0)
            params.B0 = Entry(Real, missing, IMAS.equilibrium__vacuum_toroidal_field, :b0)
            params.R0 = Entry(Real, missing, IMAS.equilibrium__vacuum_toroidal_field, :r0)
            params.Z0 = Entry(Real, missing, IMAS.equilibrium__vacuum_toroidal_field, :r0)
            params.ϵ = Entry(Real, missing, "", "Plasma aspect ratio")
            params.δ = Entry(Real, missing, IMAS.equilibrium__time_slice___boundary, :triangularity)
            params.κ = Entry(Real, missing, IMAS.equilibrium__time_slice___boundary, :elongation)
            params.βn = Entry(Real, missing, IMAS.equilibrium__time_slice___global_quantities, :beta_normal)
            params.ip = Entry(Real, missing, IMAS.equilibrium__time_slice___global_quantities, :ip)
            params.x_point = Entry(Bool, missing, IMAS.equilibrium__time_slice___boundary, :x_point)
            params.symmetric = Entry(Bool, missing, "", "Is plasma up-down symmetric")
            params.ngrid = Entry(Int, 129, "", "Resolution of the equilibrium grid")

        elseif what == :coil
            params.green_model = Entry(Symbol, :simple, "", "Model to be used for the Greens function table of the PF coils") # [:simple, :....]

        elseif what == :build
            params.n_oh_coils = Entry(Int, missing, "", "Number of OH coils")
            params.n_pf_coils_inside = Entry(Int, missing, "", "Number of PF coils inside of the TF")
            params.n_pf_coils_outside = Entry(Int, missing, "", "Number of PF coils outside of the TF")

            params.is_nuclear_facility = Entry(Bool, missing, "", "Is this a nuclear facility")

        elseif what == :gasc
            params.filename = Entry(String, missing, "", "Output GASC .json file from which data will be loaded")
            params.case = Entry(Int, missing, "", "Number of the GASC run to load")
            params.no_small_gaps = Entry(Bool, true, "", "Remove small gaps from the GASC radial build")

        elseif what == :ods
            params.filename = Entry(String, missing, "", "ODS.json file from which equilibrium is loaded")
        end

    else
        params = Parameters()

        if what == :ITER
            params.ods.filename = joinpath(dirname(abspath(@__FILE__)), "..", "sample", "ITER_eq_ods.json")

            params.equilibrium.R0 = 6.2
            params.equilibrium.ϵ = 0.32
            params.equilibrium.κ = 1.85
            params.equilibrium.δ = 0.485
            params.equilibrium.B0 = 5.3
            params.equilibrium.Z0 = 0.4
            params.equilibrium.ip = 15.E6
            params.equilibrium.βn = 2.0
            params.equilibrium.x_point = true
            params.equilibrium.symmetric = false

            params.general.init_from = :ods

            params.build.is_nuclear_facility = true
            params.build.n_oh_coils = 6
            params.build.n_pf_coils_inside = 0
            params.build.n_pf_coils_outside = 6

        elseif what == :CAT
            params.ods.filename = joinpath(dirname(abspath(@__FILE__)), "..", "sample", "CAT_eq_ods.json")

            params.general.init_from = :ods

            params.build.is_nuclear_facility = false
            params.build.n_oh_coils = 6
            params.build.n_pf_coils_inside = 0
            params.build.n_pf_coils_outside = 6

        elseif what == :D3D
            params.ods.filename = joinpath(dirname(abspath(@__FILE__)), "..", "sample", "D3D_eq_ods.json")

            params.general.init_from = :ods

            params.build.is_nuclear_facility = false
            params.build.n_oh_coils = 20
            params.build.n_pf_coils_inside = 18
            params.build.n_pf_coils_outside = 0

        elseif what == :FPP
            params.gasc.filename = joinpath(dirname(abspath(@__FILE__)), "..", "sample", "FPP_fBS_PBpR_scan.json")
            params.gasc.case = 59
            params.gasc.no_small_gaps = true

            params.general.init_from = :gasc

            params.build.is_nuclear_facility = true
            params.build.n_oh_coils = 6
            params.build.n_pf_coils_inside = 0
            params.build.n_pf_coils_outside = 6
        end
    end

    #    error("Invalid Parameters group: $what")

    return params
end
