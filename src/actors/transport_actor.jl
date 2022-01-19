using NumericalIntegration
using AD_TAUENN

function Hmode_profiles(edge::Real, ped::Real, core::Real, ngrid::Int, expin::Real, expout::Real, widthp::Real)
    w_E1 = 0.5 * widthp  # width as defined in eped
    xphalf = 1.0 - w_E1

    xped = xphalf - w_E1
 
    pconst = 1.0 - tanh((1.0 - xphalf) / w_E1)
    a_t = 2.0 * (ped - edge) / (1.0 + tanh(1.0) - pconst)

    coretanh = 0.5 * a_t * (1.0 - tanh(-xphalf / w_E1) - pconst) + edge

    xpsi = LinRange(1e-50, 1, ngrid)

    val = [0.5 * a_t * (1. - tanh((xpsi[i] - xphalf) / w_E1) - pconst) + edge * 1. for i in 1:ngrid]

    xtoped = xpsi / xped
    grid = LinRange(0, 1,ngrid)
    for (i,ival) in enumerate(grid)
        if xtoped[i]<0 
            @inbounds val[i] = val[i] + (core - coretanh)
        elseif xtoped[i] ^ expin < 1.0
            @inbounds val[i] = val[i] + (core - coretanh) * (1.0 - xtoped[i] ^ expin) ^ expout
        end
    end

    return val
end

function init(cs::IMAS.core_sources, eq::IMAS.equilibrium; Paux_e::Real, Paux_i::Real)
    resize!(cs.source, 1)
    resize!(cs.source[1].profiles_1d,1)
    for time_index in 1:length(eq.time_slice)
        cs1d = cs.source[1].profiles_1d[time_index]
        cs.source[1].identifier.name = "arb"
        cs.source[1].identifier.index = 901
        cs.source[1].identifier.description = "Arbitrary source from transport initialization"
        cs1d.grid.rho_tor_norm = rho = eq.time_slice[time_index].profiles_1d.rho_tor_norm
        cs1d.grid.volume = vol = eq.time_slice[time_index].profiles_1d.volume

        auxHeatingProfile =  exp.(-4.0 * rho)
        pow_prof = cumul_integrate(vol, auxHeatingProfile)
        pow_prof = pow_prof ./ pow_prof[end]

        cs1d.electrons.power_inside = pow_prof .* Paux_e
        cs1d.total_ion_power_inside = pow_prof .* Paux_i
    end
    return cs
end

function init(cp::IMAS.core_profiles, eqt::IMAS.equilibrium__time_slice; neped::Real, ne_peaking::Real, Te_ped::Real, Te_peaking::Real, w_ped::Real, zeff::Real, Paux::Real,T_ratio=1.0, n_points=101)
    resize!(cp.profiles_1d,1)
    cpt = cp.profiles_1d[1]
    cpt.grid.rho_tor_norm = rho =  LinRange(0, 1, n_points)
    cpt.zeff = ones(n_points) .* zeff
    cpt.rotation_frequency_tor_sonic = 5e3 * abs(Paux/1e6 * 1.0 + 0.5) .* (1.0 .- rho)

    # Set ions
    resize!(cpt.ion, 2)
    cpt.ion[1].label = "D"
    resize!(cpt.ion[1].element, 1)
    cpt.ion[1].element[1].z_n = 1
    cpt.ion[1].element[1].a = 2
    resize!(cpt.ion[2].element, 1)
    cpt.ion[2].label = "C"
    cpt.ion[2].element[1].z_n = 6
    cpt.ion[2].element[1].a = 12

    # Set densities
    ne_core = ne_peaking * neped
    #ne = TAUENN_AD.Hmode_profiles(0.5 * neped, neped, ne_core, length(rho), 0.9, 0.9, w_ped)
    ne = Hmode_profiles(0.5 * neped, neped, ne_core, length(rho), 0.9, 0.9, w_ped)

    cpt.electrons.density = ne
    zimp1 = 6.0
    niFraction = zeros(2)
    niFraction[2] = (zeff - 1.0) / (zimp1 * (zimp1 - 1.0))
    niFraction[1] = 1.0 - zimp1 * niFraction[2]
    for i in 1:length(cpt.ion)
        cpt.ion[i].density = ni = ne .* niFraction[i]
    end

    # Set temperatures
    betaN = eqt.global_quantities.beta_normal
    Bt = eqt.global_quantities.magnetic_axis.b_field_tor
    Ip = eqt.global_quantities.ip
    a = eqt.boundary.minor_radius

    Te_core = 10. * betaN * abs(Bt * (Ip/1e6)) / a / (ne_core/1e20) / (2.0 * 1.6e1 * 4.0 * pi * 1.0e-4)
#    Te = TAUENN_AD.Hmode_profiles(80., teped, tcore, length(rho), Te_peaking, Te_peaking, w_ped)
    Te = Hmode_profiles(80., Te_ped, Te_core, length(rho), Te_peaking, Te_peaking, w_ped)
    cpt.electrons.temperature = Te
    cpt.ion[1].temperature = Te ./ T_ratio
    cpt.ion[2].temperature = Te ./ T_ratio

    # to be done as an expression
    # prof1d.pressure_thermal = 1.6e-19 .* ne .* Te
    # for i in 1:length(prof1d.ion)
    #     ni = prof1d.ion[i].density
    #     prof1d.ion[i].temperature = Ti = tval .* inputs.Tratio
    #     pion =  1.6e-19 .* ni .* Ti
    #     prof1d.pressure_thermal = prof1d.pressure_thermal .+ pion
    # end
    return cp
end

#= ================ =#
#     TAUENN actor   #
#= ================ =#

mutable struct TaueNNactor <: AbstractActor
    cp::IMAS.core_profiles
    eqt::IMAS.equilibrium__time_slice
    cs::IMAS.core_sources
    rho_fluxmatch::Real
    eped_factor::Real
    temp_shape::Real
end

function TaueNNactor(cp::IMAS.core_profiles, eq::IMAS.equilibrium, cs::IMAS.core_sources; rho_fluxmatch=0.4, eped_factor=1.0, temp_shape=1.8)
    time_index = argmax([is_missing(eqt.global_quantities,:ip) ? 0.0 : abs(eqt.global_quantities.ip) for eqt in eq.time_slice])
    return TaueNNactor(cp, eq.time_slice[time_index], cs)
end

function TaueNNactor(dd::IMAS.dd; rho_fluxmatch=0.4, eped_factor=1.0, Tshape=1.8)
    time_index = argmax([is_missing(eqt.global_quantities,:ip) ? 0.0 : abs(eqt.global_quantities.ip) for eqt in dd.equilibrium.time_slice])
    return TaueNNactor(dd.core_profiles, dd.equilibrium.time_slice[time_index], dd.core_sources)
end


# step
function step(tauennactor::TaueNNactor)
    # run tauenn
    tauenn(tauennactor.cp, tauenator.eqt, tauennactor.cs, rho_fluxmatch, eped_factor, temp_shape)
    print("step of tauennactor")
end