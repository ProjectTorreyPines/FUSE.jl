import TGLFNN:run_tglf,run_tglfnn
#= ============= =#
#  ActorTGLF      #
#= ============= =#
mutable struct ActorTGLF <: PlasmaAbstractActor
    dd::IMAS.dd
    par::ParametersActor
    #tglfmod::TGLFNN.TGLFmodel
    InputTGLF::InputTGLF
    flux_solution::flux_solution
end

function ParametersActor(::Type{Val{:ActorTGLF}})
    par = ParametersActor(nothing)
    par.tglf_model = Switch([:tglf_sat0,:tglfnn], "", "TGLF model to run"; default=:tglfnn)
    par.rho_tglf = Entry(Real, "", "rho value to compute tglf fluxes"; default=0.6)
    par.warn_nn_train_bounds = Entry(Bool, "", "Raise warnings if querying cases that are certainly outside of the training range"; default=false)
    return par
end

"""
    ActorTGLF(dd::IMAS.dd, act::ParametersAllActors; kw...)

The ActorTGLF evaluates the TGLF predicted turbulence at single rho 
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
    return ActorTGLF(dd, par,InputTGLF,flux_solution)
end
"""
    step(actor::ActorTGLF;
        warn_nn_train_bounds::Bool=actor.par.warn_nn_train_bounds,
        only_powerlaw::Bool=false)

Runs pedestal actor to evaluate pedestal width and height
"""
function step(actor::ActorTGLF)
    
    eq = dd.equilibrium
    eqt = eq.time_slice[]
    eq1d = eqt.profiles_1d
    cp1d = dd.core_profiles.profiles_1d[]

    rho_cp = cp1d.grid.rho_tor_norm
    rho_eq = eq1d.rho_tor_norm
    gridpoint_eq = argmin(abs.(rho_eq .- rho_tglf))
    gridpoint_cp = argmin(abs.(rho_cp .- rho_tglf))

    input_tglf = inputtglf(dd, gridpoint_eq, gridpoint_cp)

    if par.tglf_model == :tglfnn
        sol = run_tglfnn(input_tglf, warn_nn_train_bounds)
    elseif par.tglf_model == :tglf_sat0
        sol = run_tglf(input_tglf)
    end

    actor.flux_solution = sol

    return actor
end

"""
    finalize(actor::ActorTGLF;
        temp_pedestal_ratio::Real=actor.par.temp_pedestal_ratio,
        eped_factor::Real=actor.par.eped_factor)

Writes results to dd.summary.local.pedestal
"""
function finalize(actor::ActorTGLF)

    dd = actor.dd
    PARTICLE_FLUX_e
    STRESS_TOR_i
    ENERGY_FLUX_e
    ENERGY_FLUX_i

    @ddtime dd.core_transport.model[0].profiles_1d[].electrons.energy.flux = actor.flux_solution.ENERGY_FLUX_e
    @ddtime dd.core_transport.model[0].profiles_1d[].electrons.particle.flux = actor.flux_solution.PARTICLE_FLUX_e
    @ddtime dd.core_transport.model[0].profiles_1d[].momentum_tor.flux =  actor.flux_solution.STRESS_TOR_i
    @ddtime core_transport.model[0].profiles_1d[].total_ion_energy  = actor.flux_solution.ENERGY_FLUX_i
    return actor
end

"""
    inputtglf(dd, gridpoint_eq, gridpoint_cp)::InputTGLF

Evaluate TGLF input parameters at given radii
"""
function inputtglf(dd::IMAS.dd, gridpoint_eq::Integer, gridpoint_cp::Integer)::InputTGLF
    e = 4.8032e-10
    k = 1.6022e-12
    me = 9.1094e-28
    c = 2.9979e10
    m_to_cm = 1e2
    T_to_Gauss = 1e4

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
    rmin_eq = m_to_cm * 0.5 * (eq1d.r_outboard - eq1d.r_inboard)
    rmin = IMAS.interp1d(eq1d.rho_tor_norm, rmin_eq).(cp1d.grid.rho_tor_norm)

    q = eq1d.q[gridpoint_eq]
    kappa = IMAS.interp1d(eq1d.rho_tor_norm, eq1d.elongation).(cp1d.grid.rho_tor_norm)
    delta = IMAS.interp1d(eq1d.rho_tor_norm, 0.5 * (eq1d.triangularity_lower + eq1d.triangularity_lower)).(cp1d.grid.rho_tor_norm)

    a = rmin[end]

    Te = cp1d.electrons.temperature
    dlntedr = -IMAS.gradient(rmin, Te) ./ Te
    Te = Te[gridpoint_cp]
    dlntedr = dlntedr[gridpoint_cp]

    ne = cp1d.electrons.density_thermal .* 1e-6
    dlnnedr = -IMAS.gradient(rmin, ne) ./ ne
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

    c_s = sqrt(k * Te / mi)
    w0 = -1*cp1d.rotation_frequency_tor_sonic
    w0p = IMAS.gradient(w0, rmin)
    gamma_p = -Rmaj[gridpoint_cp]*w0p[gridpoint_cp]
    gamma_e = -rmin[gridpoint_cp]/q*w0p[gridpoint_cp]
    mach = Rmaj[gridpoint_cp] * w0[gridpoint_cp] / c_s

    input_tglf.VPAR_1 = -input_tglf.SIGN_IT * mach
    input_tglf.VPAR_SHEAR_1 = -1 * input_tglf.SIGN_IT * (a / c_s) * gamma_p
    input_tglf.VEXB_SHEAR = -1*gamma_e * (a / c_s)

    for iion in 1:length(ions)
        species = iion + 1
        setfield!(input_tglf, Symbol("MASS_$species"), ions[iion].element[1].a / ions[1].element[1].a)
        setfield!(input_tglf, Symbol("ZS_$species"), Int(floor(ions[iion].element[1].z_n / ions[1].element[1].z_n)))

        Ti = ions[iion].temperature
        dlntidr = -IMAS.gradient(rmin, Ti) ./ Ti
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
    rho_s = c_s / (e * abs(bunit) / (mi * c))
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

    input_tglf.S_KAPPA_LOC = skappa[gridpoint_cp]
    input_tglf.DELTA_LOC = delta[gridpoint_cp]
    input_tglf.S_DELTA_LOC = sdelta[gridpoint_cp]
    input_tglf.ZETA_LOC = 0.0
    input_tglf.S_ZETA_LOC = 0.0

    press = cp1d.pressure_thermal
    Pa_to_dyn = 10.0

    dpdr = IMAS.gradient(rmin, press .* Pa_to_dyn)[gridpoint_cp]
    input_tglf.P_PRIME_LOC = abs(q) / (rmin[gridpoint_cp] / a)^2 * rmin[gridpoint_cp] / bunit^2 * dpdr

    dqdr = IMAS.gradient(rmin_eq, eq1d.q)[gridpoint_eq]
    s = rmin[gridpoint_cp] / q * dqdr
    input_tglf.Q_PRIME_LOC = q^2 * a^2 / rmin[gridpoint_cp]^2 * s

    input_tglf._Qgb = ne * k * Te * c_s * (rho_s / a)^2 * 1.0e-7

    return input_tglf
end
