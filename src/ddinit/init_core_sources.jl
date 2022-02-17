import NumericalIntegration: cumul_integrate

function init_core_sources(dd::IMAS.dd, par::Parameters)
    init_from = par.general.init_from

    if init_from == :gasc
        gasc = GASC(par.gasc.filename, par.gasc.case)
        init_core_sources(dd, gasc)

    elseif init_from == :ods
        dd1 = IMAS.json2imas(par.ods.filename)
        if length(keys(dd1.core_sources)) > 0
            dd.core_sources = dd1.core_sources
        else
            init_from = :scalars
        end
    end

    if init_from == :scalars
        init_core_sources(dd.core_sources, dd.equilibrium)
    end

    return dd
end

function init_core_sources(dd::IMAS.dd, gasc::GASC)
    gasc = gasc.solution

    heating_power = gasc["OUTPUTS"]["current drive"]["powerAux"]
    cd_power = heating_power / gasc["current drive"]["auxCDPowerFactor"]
    plug = heating_power / gasc["INPUTS"]["power efficiency"]["efficiencyAux"]

    cd_powers = Dict()
    for system in ["NNB", "NB", "LH", "FW", "EC", "HI"]
        cd_powers[system] = cd_power * gasc["INPUTS"]["current drive"]["$(system)CDFraction"]
    end

    init_core_sources(dd.core_sources, dd.equilibrium)
    return dd
end


function init_core_sources(cs::IMAS.core_sources, eq::IMAS.equilibrium; Paux_e::Real = 0.0, Paux_i::Real = 0.0, ngrid::Int = 51)
    empty!(cs)

    resize!(cs.source, 1)
    cs.source[1].identifier.name = "arb"
    cs.source[1].identifier.index = 901
    cs.source[1].identifier.description = "Arbitrary source from FUSE transport initialization"

    resize!(cs.source[1].profiles_1d)
    cs1d = cs.source[1].profiles_1d[]
    cs1d.grid.rho_tor_norm = rho = LinRange(0.0, 1.0, ngrid)

    rho_eq = eq.time_slice[].profiles_1d.rho_tor_norm
    cs1d.grid.volume = IMAS.interp(rho_eq, eq.time_slice[].profiles_1d.volume)[rho]

    auxHeatingProfile = exp.(-4.0 * rho)
    pow_prof = cumul_integrate(cs1d.grid.volume, auxHeatingProfile)
    pow_prof = pow_prof ./ pow_prof[end]

    cs1d.electrons.power_inside = pow_prof .* Paux_e
    cs1d.total_ion_power_inside = pow_prof .* Paux_i

    return cs
end


