import TGLFNN: run_tglf, run_tglfnn, InputTGLF, flux_solution
#= ============= =#
#  ActorTGLF      #
#= ============= =#
mutable struct ActorTGLF <: PlasmaAbstractActor
    dd::IMAS.dd
    par::ParametersActor
    input_tglf::Union{InputTGLF,Missing}
    flux_solution::Union{flux_solution,Missing}
end

function ParametersActor(::Type{Val{:ActorTGLF}})
    par = ParametersActor(nothing)
    par.tglf_model = Switch([:tglf_sat0, :tglfnn], "", "TGLF model to run"; default=:tglfnn)
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
    eq1d = dd.equilibrium.time_slice[].profiles_1d
    cp1d = dd.core_profiles.profiles_1d[]

    rho_cp = cp1d.grid.rho_tor_norm
    rho_eq = eq1d.rho_tor_norm
    gridpoint_eq = argmin(abs.(rho_eq .- par.rho_tglf))
    gridpoint_cp = argmin(abs.(rho_cp .- par.rho_tglf))

    input_tglf = inputtglf(dd, gridpoint_eq, gridpoint_cp)

    return ActorTGLF(dd, par, input_tglf, missing)
end

"""
    step(actor::ActorTGLF)

Runs TGLF actor to evaluate the turbulence flux on a gridpoint
"""
function step(actor::ActorTGLF)
    par = actor.par

    if par.tglf_model == :tglfnn
        sol = run_tglfnn(input_tglf, par.warn_nn_train_bounds)
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
    cp1d = dd.core_profiles.profiles_1d[]
    eqt = dd.equilibrium.time_slice[]

    model = dd.core_transport.model[[idx for idx in keys(actor.dd.core_transport.model) if actor.dd.core_transport.model[idx].identifier.name == string(actor.par.tglf_model)][1]]
    rho_transp_idx = findfirst(i -> i == actor.par.rho_tglf, model.profiles_1d[].grid_flux.rho_tor_norm)
    rho_cp_idx = findfirst(i -> i == actor.par.rho_tglf, cp1d.grid.rho_tor_norm)
    if isnothing(rho_transp_idx)
        error("The transport grid doesn't have the tglf flux matching gridpoint")
    end
    model.profiles_1d[].electrons.energy.flux[rho_transp_idx] = actor.flux_solution.ENERGY_FLUX_e * IMAS.gyrobohm_energy_flux(cp1d, eqt)[rho_cp_idx] # W / m^2
    model.profiles_1d[].total_ion_energy.flux[rho_transp_idx] = actor.flux_solution.ENERGY_FLUX_i * IMAS.gyrobohm_energy_flux(cp1d, eqt)[rho_cp_idx] # W / m^2
    model.profiles_1d[].electrons.particles.flux[rho_transp_idx] = actor.flux_solution.PARTICLE_FLUX_e * IMAS.gyrobohm_particle_flux(cp1d, eqt)[rho_cp_idx] # 1 / m^2 / s
    model.profiles_1d[].momentum_tor.flux[rho_transp_idx] = actor.flux_solution.STRESS_TOR_i * IMAS.gyrobohm_momentum_flux(cp1d, eqt)[rho_cp_idx] #
    return actor
end

"""
    inputtglf(dd, gridpoint_eq, gridpoint_cp)::InputTGLF

Evaluate TGLF input parameters at given radii
"""
function inputtglf(dd::IMAS.dd, gridpoint_eq::Integer, gridpoint_cp::Integer)
    e = IMAS.gacode_units.e
    k = IMAS.gacode_units.k
    me = IMAS.gacode_units.me
    c = IMAS.gacode_units.c
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
    rmin = IMAS.r_min_core_profiles(cp1d, eqt)

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

    c_s = IMAS.c_s(cp1d)[gridpoint_cp]
    w0 = -1 * cp1d.rotation_frequency_tor_sonic
    if any(i -> i == 0, w0)
        w0p = zeros(length(w0))
    else
        w0p = IMAS.gradient(w0, rmin)
    end
    gamma_p = -Rmaj[gridpoint_cp] * w0p[gridpoint_cp]
    gamma_e = -rmin[gridpoint_cp] / q * w0p[gridpoint_cp]
    mach = Rmaj[gridpoint_cp] * w0[gridpoint_cp] / c_s
    input_tglf.VPAR_1 = -input_tglf.SIGN_IT * mach
    input_tglf.VPAR_SHEAR_1 = -1 * input_tglf.SIGN_IT * (a / c_s) * gamma_p
    input_tglf.VEXB_SHEAR = -1 * gamma_e * (a / c_s)

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

    input_tglf.S_KAPPA_LOC = skappa[gridpoint_cp]
    input_tglf.DELTA_LOC = delta[gridpoint_cp]
    input_tglf.S_DELTA_LOC = sdelta[gridpoint_cp]
    input_tglf.ZETA_LOC = 0.0
    input_tglf.S_ZETA_LOC = 0.0

    press = cp1d.pressure_thermal
    Pa_to_dyn = 10.0

    dpdr = IMAS.gradient(rmin, press .* Pa_to_dyn)[gridpoint_cp]
    input_tglf.P_PRIME_LOC = abs(q) / (rmin[gridpoint_cp] / a)^2 * rmin[gridpoint_cp] / bunit^2 * dpdr

    dqdr = IMAS.gradient(m_to_cm .* 0.5 .* (eq1d.r_outboard .- eq1d.r_inboard), eq1d.q)[gridpoint_eq]
    s = rmin[gridpoint_cp] / q * dqdr
    input_tglf.Q_PRIME_LOC = q^2 * a^2 / rmin[gridpoint_cp]^2 * s

    return input_tglf
end

"""
    tglf_multi(dd::IMAS.dd, model::Symbol, rho_gridpoints::Vector{<:Real})

Sets up the transport grid and runs ActorTGLF on all the transport gird points serially
"""
function tglf_multi(dd::IMAS.dd, model::Symbol, rho_gridpoints::Vector{<:Real})
    model = resize!(dd.core_transport.model, "identifier.name" => string(model))
    IMAS.setup_transport_grid!(model, rho_gridpoints)
    return map(rho_tglf -> FUSE.ActorTGLF(dd, act, rho_tglf=rho_tglf), rho_gridpoints)

end
"""
    lump_ions_as_bulk_and_impurity!(ions::IMAS.IDSvector{IMAS.core_profiles__profiles_1d___ion})

Changes core_profiles.ion to 2 species, bulk specie (H, D, T) and combined impurity specie by weigthing masses and densities 
"""
function lump_ions_as_bulk_and_impurity!(ions::IMAS.IDSvector{IMAS.core_profiles__profiles_1d___ion})
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