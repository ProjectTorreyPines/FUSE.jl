import TGLFNN: run_tglf, run_tglfnn, InputTGLF, flux_solution

#= ============= =#
#  ActorTGLF      #
#= ============= =#
mutable struct ActorTGLF <: PlasmaAbstractActor
    dd::IMAS.dd
    par::ParametersActor
    input_tglfs::AbstractVector{<:InputTGLF}
    flux_solutions::AbstractVector{<:flux_solution}
end

function ParametersActor(::Type{Val{:ActorTGLF}})
    par = ParametersActor(nothing)
    par.nn = Entry(Bool, "", "Use TGLF-NN"; default=true)
    par.sat_rule = Switch([:sat0, :sat0quench, :sat1, :sat1geo, :sat2], "", "Saturation rule"; default=:sat0)
    par.electromagnetic = Entry(Bool, "", "Electromagnetic or electrostatic"; default=false)
    par.rho_transport = Entry(AbstractVector{<:Real}, "", "rho_tor_norm values to compute tglf fluxes on"; default=0.2:0.1:0.8)
    par.warn_nn_train_bounds = Entry(Bool, "", "Raise warnings if querying cases that are certainly outside of the training range"; default=false)
    return par
end

"""
    ActorTGLF(dd::IMAS.dd, act::ParametersAllActors; kw...)

The ActorTGLF evaluates the TGLF predicted turbulence at a set of rho_tor_norm grid points
"""
function ActorTGLF(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorTGLF(kw...)
    actor = ActorTGLF(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function ActorTGLF(dd::IMAS.dd, par::ParametersActor; kw...)
    par = par(kw...)
    input_tglfs = Vector{InputTGLF}(undef, length(par.rho_transport))
    return ActorTGLF(dd, par, input_tglfs, flux_solution[])
end

"""
    step(actor::ActorTGLF)

Runs TGLF actor to evaluate the turbulence flux on a Vector of gridpoints
"""
function _step(actor::ActorTGLF)
    par = actor.par
    dd = actor.dd

    eq1d = dd.equilibrium.time_slice[].profiles_1d
    cp1d = dd.core_profiles.profiles_1d[]
    ix_eq = [argmin(abs.(eq1d.rho_tor_norm .- rho)) for rho in par.rho_transport]
    ix_cp = [argmin(abs.(cp1d.grid.rho_tor_norm .- rho)) for rho in par.rho_transport]
    for (k, (gridpoint_eq, gridpoint_cp)) in enumerate(zip(ix_eq, ix_cp))
        actor.input_tglfs[k] = inputtglf(dd, gridpoint_eq, gridpoint_cp, par.sat_rule, par.electromagnetic)
    end

    tglf_model = string(par.sat_rule) * "_" * (par.electromagnetic ? "em" : "es")

    anomalous_index = IMAS.name_2_index(dd.core_transport.model)[:anomalous]
    model = resize!(dd.core_transport.model, "identifier.index" => anomalous_index)
    model.identifier.name = (par.nn ? "TGLF-NN" : "TGLF") * " " * tglf_model
    m1d = resize!(model.profiles_1d)
    m1d.grid_flux.rho_tor_norm = par.rho_transport

    if par.nn
        actor.flux_solutions = map(input_tglf -> run_tglfnn(input_tglf, par.warn_nn_train_bounds; model_filename=tglf_model), actor.input_tglfs)
    else
        actor.flux_solutions = asyncmap(input_tglf -> run_tglf(input_tglf), actor.input_tglfs)
    end

    return actor
end

"""
    finalize(actor::ActorTGLF)

Writes results to dd.core_transport
"""
function finalize(actor::ActorTGLF)
    dd = actor.dd
    par = actor.par
    cp1d = dd.core_profiles.profiles_1d[]
    eqt = dd.equilibrium.time_slice[]

    model = findfirst(:anomalous, actor.dd.core_transport.model)
    m1d = model.profiles_1d[]
    m1d.electrons.energy.flux = zeros(length(par.rho_transport))
    m1d.total_ion_energy.flux = zeros(length(par.rho_transport))
    m1d.electrons.particles.flux = zeros(length(par.rho_transport))
    m1d.momentum_tor.flux = zeros(length(par.rho_transport))
    for (tglf_idx, rho) in enumerate(par.rho_transport)
        rho_transp_idx = findfirst(i -> i == rho, m1d.grid_flux.rho_tor_norm)
        rho_cp_idx = argmin(abs.(cp1d.grid.rho_tor_norm .- rho))
        m1d.electrons.energy.flux[rho_transp_idx] = actor.flux_solutions[tglf_idx].ENERGY_FLUX_e * IMAS.gyrobohm_energy_flux(cp1d, eqt)[rho_cp_idx] # W / m^2
        m1d.total_ion_energy.flux[rho_transp_idx] = actor.flux_solutions[tglf_idx].ENERGY_FLUX_i * IMAS.gyrobohm_energy_flux(cp1d, eqt)[rho_cp_idx] # W / m^2
        m1d.electrons.particles.flux[rho_transp_idx] = actor.flux_solutions[tglf_idx].PARTICLE_FLUX_e * IMAS.gyrobohm_particle_flux(cp1d, eqt)[rho_cp_idx] # 1 / m^2 / s
        m1d.momentum_tor.flux[rho_transp_idx] = actor.flux_solutions[tglf_idx].STRESS_TOR_i * IMAS.gyrobohm_momentum_flux(cp1d, eqt)[rho_cp_idx] #
    end
    return actor
end

"""
    inputtglf(dd, gridpoint_eq, gridpoint_cp)::InputTGLF

Evaluate TGLF input parameters at given radii
"""
function inputtglf(dd::IMAS.dd, gridpoint_eq::Integer, gridpoint_cp::Integer, sat::Symbol=:sat0, electromagnetic::Bool=false)
    e = IMAS.gacode_units.e
    k = IMAS.gacode_units.k
    me = IMAS.gacode_units.me
    m_to_cm = IMAS.gacode_units.m_to_cm
    T_to_Gauss = IMAS.gacode_units.T_to_Gauss

    eq = dd.equilibrium
    eqt = eq.time_slice[]
    eq1d = eqt.profiles_1d
    cp1d = dd.core_profiles.profiles_1d[]
    ions = cp1d.ion
    if length(ions) > 2
        ions = lump_ions_as_bulk_and_impurity!(deepcopy(ions))
    end
    mi = ions[1].element[1].a * 1.6726e-24

    Rmaj = IMAS.interp1d(eq1d.rho_tor_norm, eq1d.gm8 * m_to_cm).(cp1d.grid.rho_tor_norm)
    Rmaj = IMAS.interp1d(eq1d.rho_tor_norm, m_to_cm * 0.5 * (eq1d.r_outboard .+ eq1d.r_inboard)).(cp1d.grid.rho_tor_norm)

    rmin = IMAS.r_min_core_profiles(cp1d, eqt)

    q_profile = IMAS.interp1d(eq1d.rho_tor_norm, eq1d.q).(cp1d.grid.rho_tor_norm)
    kappa = IMAS.interp1d(eq1d.rho_tor_norm, eq1d.elongation).(cp1d.grid.rho_tor_norm)
    delta = IMAS.interp1d(eq1d.rho_tor_norm, 0.5 * (eq1d.triangularity_lower + eq1d.triangularity_upper)).(cp1d.grid.rho_tor_norm)
    zeta = IMAS.interp1d(eq1d.rho_tor_norm, 0.25 * (eq1d.squareness_lower_inner .+ eq1d.squareness_lower_outer .+ eq1d.squareness_upper_inner .+ eq1d.squareness_upper_outer)).(cp1d.grid.rho_tor_norm)

    a = rmin[end]
    q = q_profile[gridpoint_cp]

    Te = cp1d.electrons.temperature
    dlntedr = -IMAS.calc_z(rmin, Te)
    Te = Te[gridpoint_cp]
    dlntedr = dlntedr[gridpoint_cp]

    ne = cp1d.electrons.density_thermal .* 1e-6
    dlnnedr = -IMAS.calc_z(rmin, ne)
    ne = ne[gridpoint_cp]
    dlnnedr = dlnnedr[gridpoint_cp]

    Bt = @ddtime eq.vacuum_toroidal_field.b0
    bunit = IMAS.bunit(eqt)[gridpoint_eq] * T_to_Gauss

    input_tglf = InputTGLF()

    signb = sign(Bt)
    signq = sign(q)
    input_tglf.SIGN_BT = signb
    input_tglf.SIGN_IT = signb * signq

    input_tglf.NS = length(ions) + 1
    input_tglf.MASS_1 = me / mi
    input_tglf.TAUS_1 = 1.0
    input_tglf.AS_1 = 1.0
    input_tglf.ZS_1 = -1

    input_tglf.RLNS_1 = a .* dlnnedr
    input_tglf.RLTS_1 = a .* dlntedr

    c_s = IMAS.c_s(cp1d)[gridpoint_cp]
    w0 = -1 * cp1d.rotation_frequency_tor_sonic
    if any(i -> i == 0, w0)
        w0p = zeros(length(w0))
    else
        w0p = IMAS.gradient(rmin, w0)
    end
    gamma_p = -Rmaj[gridpoint_cp] * w0p[gridpoint_cp]
    gamma_e = -rmin[gridpoint_cp] / q * w0p[gridpoint_cp]
    mach = Rmaj[gridpoint_cp] * w0[gridpoint_cp] / c_s
    input_tglf.VPAR_1 = -input_tglf.SIGN_IT * mach
    input_tglf.VPAR_SHEAR_1 = -1 * input_tglf.SIGN_IT * (a / c_s) * gamma_p
    input_tglf.VEXB_SHEAR = 1 * gamma_e * (a / c_s)

    for iion in 1:length(ions)
        species = iion + 1
        setfield!(input_tglf, Symbol("MASS_$species"), ions[iion].element[1].a / ions[1].element[1].a)
        setfield!(input_tglf, Symbol("ZS_$species"), Int(floor(ions[iion].element[1].z_n / ions[1].element[1].z_n)))

        Ti = ions[iion].temperature
        dlntidr = -IMAS.calc_z(rmin, Ti)
        Ti = Ti[gridpoint_cp]
        dlntidr = dlntidr[gridpoint_cp]

        ni = ions[iion].density_thermal .* 1e-6
        dlnnidr = -IMAS.gradient(rmin, ni) ./ ni
        ni = ni[gridpoint_cp]
        dlnnidr = dlnnidr[gridpoint_cp]

        setfield!(input_tglf, Symbol("TAUS_$species"), Ti / Te)
        setfield!(input_tglf, Symbol("AS_$species"), ni / ne)
        setfield!(input_tglf, Symbol("VPAR_$species"), input_tglf.VPAR_1)
        setfield!(input_tglf, Symbol("VPAR_SHEAR_$species"), input_tglf.VPAR_SHEAR_1)
        setfield!(input_tglf, Symbol("RLNS_$species"), a * dlnnidr)
        setfield!(input_tglf, Symbol("RLTS_$species"), a * dlntidr)
    end

    input_tglf.BETAE = 8.0 * pi * ne * k * Te / bunit^2
    loglam = 24.0 - log(sqrt(ne) / (Te))
    input_tglf.XNUE = a / c_s * sqrt(ions[1].element[1].a) * e^4 * pi * ne * loglam / (sqrt(me) * (k * Te)^1.5)
    input_tglf.ZEFF = cp1d.zeff[gridpoint_cp]
    rho_s = IMAS.rho_s(cp1d, eqt)[gridpoint_cp]
    input_tglf.DEBYE = 7.43e2 * sqrt(Te / (ne)) / rho_s
    input_tglf.RMIN_LOC = rmin[gridpoint_cp] / a
    input_tglf.RMAJ_LOC = Rmaj[gridpoint_cp] / a
    input_tglf.ZMAJ_LOC = 0
    input_tglf.DRMINDX_LOC = 1.0

    drmaj = IMAS.gradient(rmin, Rmaj)

    input_tglf.DRMAJDX_LOC = drmaj[gridpoint_cp]
    input_tglf.DZMAJDX_LOC = 0.0

    input_tglf.Q_LOC = abs(q)

    input_tglf.KAPPA_LOC = kappa[gridpoint_cp]

    skappa = rmin .* IMAS.gradient(rmin, kappa) ./ kappa
    sdelta = rmin .* IMAS.gradient(rmin, delta)
    szeta = rmin .* IMAS.gradient(rmin, zeta)

    input_tglf.S_KAPPA_LOC = skappa[gridpoint_cp]
    input_tglf.DELTA_LOC = delta[gridpoint_cp]
    input_tglf.S_DELTA_LOC = sdelta[gridpoint_cp]
    input_tglf.ZETA_LOC = zeta[gridpoint_cp]
    input_tglf.S_ZETA_LOC = szeta[gridpoint_cp]

    press = cp1d.pressure_thermal
    Pa_to_dyn = 10.0

    dpdr = IMAS.gradient(rmin, press .* Pa_to_dyn)[gridpoint_cp]
    input_tglf.P_PRIME_LOC = abs(q) / (rmin[gridpoint_cp] / a)^2 * rmin[gridpoint_cp] / bunit^2 * dpdr

    dqdr = IMAS.gradient(rmin, q_profile)[gridpoint_cp]
    s = rmin[gridpoint_cp] / q * dqdr
    input_tglf.Q_PRIME_LOC = q^2 * a^2 / rmin[gridpoint_cp]^2 * s

    # saturation rules
    input_tglf.ALPHA_ZF = 1.0 # 1=default, -1=low ky cutoff kypeak search
    input_tglf.USE_MHD_RULE = true # default is true
    input_tglf.NKY = 12 # 17 for NN, 12 for validation
    if sat == :sat2
        input_tglf.UNITS = "CGYRO"
        input_tglf.SAT_RULE = 2
        input_tglf.KYGRID_MODEL = 4
        input_tglf.NMODES = 5
        input_tglf.NBASIS_MAX = 6
        input_tglf.USE_AVE_ION_GRID = true
        input_tglf.ALPHA_QUENCH = 0 # 0=spectral shift, 1=quench
    else
        input_tglf.UNITS = "GYRO"
        if sat == :sat1
            input_tglf.SAT_RULE = 1
            input_tglf.ALPHA_QUENCH = 0
        elseif sat == :sat1geo
            input_tglf.SAT_RULE = 1
            input_tglf.ALPHA_QUENCH = 0
            input_tglf.UNITS = "CGYRO"
        elseif sat == :sat0
            input_tglf.SAT_RULE = 0
            input_tglf.ALPHA_QUENCH = 0
        elseif sat == :sat0quench
            input_tglf.SAT_RULE = 0
            input_tglf.ALPHA_QUENCH = 1
        end
        input_tglf.KYGRID_MODEL = 1
        input_tglf.NMODES = 2
        input_tglf.NBASIS_MAX = 4
        input_tglf.USE_AVE_ION_GRID = false # default is false
    end

    # electrostatic/electromagnetic
    if electromagnetic
        input_tglf.USE_BPER = true
        input_tglf.USE_BPAR = true
    else
        input_tglf.USE_BPER = false
        input_tglf.USE_BPAR = false
    end

    return input_tglf
end

"""
    lump_ions_as_bulk_and_impurity!(ions::IMAS.IDSvector{<:IMAS.core_profiles__profiles_1d___ion})

Changes core_profiles.ion to 2 species, bulk specie (H, D, T) and combined impurity specie by weigthing masses and densities 
"""
function lump_ions_as_bulk_and_impurity!(ions::IMAS.IDSvector{<:IMAS.core_profiles__profiles_1d___ion})
    if length(ions) < 2
        error("TAUENN requires two ion species to run")
    elseif any(!ismissing(ion, :density_fast) for ion in ions)
        error("lump_ions_as_bulk_and_impurity! is not setup for handling fast ions")
    elseif length(ions) == 2
        return ions
    end
    zs = [ion.element[1].z_n for ion in ions]
    as = [ion.element[1].a for ion in ions]

    bulk_index = findall(zs .== 1)
    impu_index = findall(zs .!= 1)

    ratios = zeros(length(ions[1].density_thermal), length(ions))
    for index in [bulk_index, impu_index]
        ntot = zeros(length(ions[1].density_thermal))
        for ix in index
            ratios[:, ix] = ions[ix].density_thermal * zs[ix]
            ntot .+= ions[ix].density_thermal * zs[ix]
        end
        for ix in index
            ratios[:, ix] ./= ntot
        end
    end

    # bulk ions
    push!(ions, IMAS.core_profiles__profiles_1d___ion())
    bulk = ions[end]
    resize!(bulk.element, 1)
    bulk.label = "bulk"

    # impurity ions
    push!(ions, IMAS.core_profiles__profiles_1d___ion())
    impu = ions[end]
    resize!(impu.element, 1)
    impu.label = "impurity"

    # weight different ion quantities based on 
    for (index, ion) in [(bulk_index, bulk), (impu_index, impu)]
        ion.element[1].z_n = 0.0
        ion.element[1].a = 0.0
        for ix in index # z_average is tricky since it's a single constant for the whole profile
            ion.element[1].z_n += sum(zs[ix] .* ratios[:, ix]) / length(ratios[:, ix])
            ion.element[1].a += sum(as[ix] .* ratios[:, ix]) / length(ratios[:, ix])
        end
        for item in [:density_thermal, :temperature, :rotation_frequency_tor]
            if typeof(getfield(ions[index[1]], item)) <: Union{Vector{T},T} where {T<:AbstractFloat}
                setfield!(ion, item, getfield(ions[index[1]], item) .* 0.0)
                for ix in index
                    setfield!(ion, item, getfield(ion, item) .+ getfield(ions[ix], item) .* ratios[:, ix])
                end
            end
        end
    end

    for k in reverse(1:length(ions)-2)
        deleteat!(ions, k)
    end

    return ions
end