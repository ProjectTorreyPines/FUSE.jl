function LH_analysis(dd::IMAS.dd, transition_start::Float64, density_transition_end::Float64, temperature_transition_end::Float64=density_transition_end; do_plot::Bool=true)
    rho = dd.core_profiles.profiles_1d[1].grid.rho_tor_norm
    index09 = argmin_abs(rho, 0.9)
    time = dd.core_profiles.time

    tau_n = density_transition_end - transition_start
    tau_t = temperature_transition_end - transition_start

    # density timescale
    ne = [IMAS.ne_line(dd.equilibrium.time_slice[cp1d.time], cp1d) for cp1d in dd.core_profiles.profiles_1d]
    ne0 = [cp1d.electrons.density[1] for cp1d in dd.core_profiles.profiles_1d]
    ne09 = [cp1d.electrons.density[index09] for cp1d in dd.core_profiles.profiles_1d]

    # zeff timescale
    zeff = [sum(cp1d.zeff) / length(cp1d.grid.rho_tor_norm) for cp1d in dd.core_profiles.profiles_1d]
    zeff0 = [cp1d.zeff[1] for cp1d in dd.core_profiles.profiles_1d]
    zeff09 = [cp1d.zeff[index09] for cp1d in dd.core_profiles.profiles_1d]

    # temperature timescale
    te = [IMAS.trapz(rho, cp1d.electrons.temperature) for cp1d in dd.core_profiles.profiles_1d]
    te0 = [cp1d.electrons.temperature[1] for cp1d in dd.core_profiles.profiles_1d]
    te09 = [cp1d.electrons.temperature[index09] for cp1d in dd.core_profiles.profiles_1d]

    # smooth density based on L-H transition timescale 
    fs = 1 / (time[2] - time[1])
    band = 1 / tau_n
    ne_smoothed = IMAS.lowpassfilter(ne, fs, band)
    zeff_smoothed = IMAS.lowpassfilter(zeff09, fs, band)

    index_L = (time .< transition_start - tau_n / 2)
    index_no_LH = (time .< transition_start - tau_n / 2) .|| (time .> density_transition_end + tau_n / 2)
    index_H = (time .> density_transition_end + tau_n / 2)

    ne_Lvalue = IMAS.interp1d(time, ne_smoothed).(transition_start - tau_n / 4)
    ne_Hvalue = IMAS.interp1d(time, ne_smoothed).(density_transition_end + tau_n / 4)

    zeff_Lvalue = IMAS.interp1d(time, zeff_smoothed).(transition_start - tau_n / 4)
    zeff_Hvalue = IMAS.interp1d(time, zeff_smoothed).(density_transition_end + tau_n / 4)

    # smoothed H-mode density trace
    ne_H = deepcopy(ne_smoothed)
    ne_H[index_L] .*= ne_Hvalue / ne_Lvalue
    ne_H = IMAS.interp1d(time[index_no_LH], ne_H[index_no_LH], :pchip).(time)

    # smoothed H-mode zeff trace
    zeff_H = deepcopy(zeff_smoothed)
    zeff_H[index_L] .*= zeff_Hvalue / zeff_Lvalue
    zeff_H = IMAS.interp1d(time[index_no_LH], zeff_H[index_no_LH], :pchip).(time)

    # smoothed L-mode density trace
    ne_L = deepcopy(ne_smoothed)
    ne_L[index_H] .*= ne_Lvalue / ne_Hvalue
    ne_L = IMAS.interp1d(time[index_no_LH], ne_L[index_no_LH], :pchip).(time)

    # smoothed L-mode density trace
    zeff_L = deepcopy(zeff_smoothed)
    zeff_L[index_H] .*= zeff_Lvalue / zeff_Hvalue
    zeff_L = IMAS.interp1d(time[index_no_LH], zeff_L[index_no_LH], :pchip).(time)

    # WPED
    ped_to_core(W) = W[2] / W[1]
    W_ped_to_core = [ped_to_core(IMAS.core_edge_energy(cp1d, 0.9)) for cp1d in dd.core_profiles.profiles_1d]
    W_ped_to_core_fraction = sum(W_ped_to_core[index_L]) / sum(index_L)

    if do_plot
        p = plot(; layout=(4, 1), size=(1000, 1000), link=:x)

        # Density
        subplot = 1
        plot!(; title="Electron density", subplot, legend_position=:bottomright, ylim=(0, Inf))
        plot!(time, ne0; label="On axis density", lw=2, subplot)
        plot!(time, ne; label="Line average density", lw=2, subplot)
        plot!(time, ne09; label="Density @rho=0.9", lw=2, subplot)
        plot!(time, ne_L; label="Line average density L-mode", lw=2, subplot)
        plot!(time, ne_H; label="Line average density H-mode", lw=2, subplot)
        vline!([transition_start, density_transition_end]; label="Density transition time", subplot)

        # Zeff
        subplot = 2
        plot!(; title="Zeff", subplot, link=:x, legend_position=:bottomright, ylim=(0, Inf))
        plot!(time, zeff0; label="On axis Zeff", lw=2, subplot)
        plot!(time, zeff; label="Average Zeff", lw=2, subplot)
        plot!(time, zeff09; label="Zeff @rho=0.9", lw=2, subplot)
        plot!(time, zeff_L; label="Zeff L-mode", lw=2, subplot)
        plot!(time, zeff_H; label="Zeff H-mode", lw=2, subplot)
        vline!([transition_start, density_transition_end]; label="Density transition time", subplot)

        # Temperature
        subplot = 3
        plot!(; title="Electron temperature", subplot, link=:x, legend_position=:bottomright, ylim=(0, Inf))
        plot!(time, te0; label="On axis temperature", lw=2, subplot)
        plot!(time, te; label="Average temperature", lw=2, subplot)
        plot!(time, te09; label="Temperature @rho=0.9", lw=2, subplot)
        vline!([transition_start, temperature_transition_end]; label="Temperature transition time", subplot)

        # L-H power threshold based on scaling law
        subplot = 4
        plot!(; title="L-H power thresold (no radiation)", subplot, xlabel="Time [s]", link=:x, legend_position=:bottomright, ylim=(0, Inf))
        injected_power = IMAS.total_power(dd.pulse_schedule, time; tau_smooth=time[2] - time[1])
        scaling_power = [IMAS.scaling_L_to_H_power(dd; time0) for time0 in time]
        plot!(time, injected_power ./ scaling_power; label="Power / Power threshold", lw=2, subplot)
        vline!([transition_start, temperature_transition_end, density_transition_end]; label="Transition times", subplot)
        #        scatter!([dd.global_time], [IMAS.L_H_threshold(dd)], subplot, label="Power threshold (with radiation)")
        #plot!(time, 1.0./W_ped_to_core; subplot, label="Core to pedestal stored energy")

        display(p)
    end

    return (time=time, tau_n=tau_n, tau_t=tau_t,
        mode_transitions=Dict{Float64,Symbol}(0.0 => :L_mode, transition_start => :H_mode),
        ne_L=ne_L, ne_H=ne_H, ne_L_over_H=ne_Lvalue / ne_Hvalue,
        zeff_L=zeff_L, zeff_H=zeff_H, zeff_L_over_H=zeff_Lvalue / zeff_Hvalue,
        W_ped_to_core_fraction=W_ped_to_core_fraction)
end