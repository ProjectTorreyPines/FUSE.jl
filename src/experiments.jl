function ne_line_without_LH_transition(dd::IMAS.dd, transition_start::Float64, density_transition_end::Float64, temperature_transition_end::Float64; do_plot::Bool=true)
    rho = dd.core_profiles.profiles_1d[1].grid.rho_tor_norm
    index09 = argmin(abs.(rho .- 0.9))
    time = dd.core_profiles.time

    # density timescale
    ne = IMAS.interp1d(dd.pulse_schedule.density_control.time, dd.pulse_schedule.density_control.n_e_line.reference).(time)
    ne0 = [cp1d.electrons.density[1] for cp1d in dd.core_profiles.profiles_1d]
    ne09 = [cp1d.electrons.density[index09] for cp1d in dd.core_profiles.profiles_1d]
    tau_n = density_transition_end - transition_start

    # temperature timescale
    te = [IMAS.trapz(rho, cp1d.electrons.temperature) for cp1d in dd.core_profiles.profiles_1d]
    te0 = [cp1d.electrons.temperature[1] for cp1d in dd.core_profiles.profiles_1d]
    te09 = [cp1d.electrons.temperature[index09] for cp1d in dd.core_profiles.profiles_1d]
    tau_t = temperature_transition_end - transition_start

    # smooth density based on L-H transition timescale 
    fs = 1 / (time[2] - time[1])
    band = 1 / tau_n
    ne_smoothed = IMAS.lowpassfilter(ne, fs, band)

    Lvalue = IMAS.interp1d(time, ne_smoothed).(transition_start)
    Hvalue = IMAS.interp1d(time, ne_smoothed).(density_transition_end)

    index_L = (time .< transition_start - tau_n / 2)
    index_no_LH = (time .< transition_start - tau_n / 2) .|| (time .> density_transition_end + tau_n / 2)
    index_H = (time .> density_transition_end + tau_n / 2)

    # smoothed H-mode density trace
    ne_H = deepcopy(ne_smoothed)
    ne_H[index_L] .*= Hvalue / Lvalue
    ne_H = IMAS.interp1d(time[index_no_LH], ne_H[index_no_LH], :pchip).(time)

    # smoothed L-mode density trace
    ne_L = deepcopy(ne_smoothed)
    ne_L[index_H] .*= Lvalue / Hvalue
    ne_L = IMAS.interp1d(time[index_no_LH], ne_L[index_no_LH], :pchip).(time)

    # WPED
    ped_to_core(W) = W[2]/W[1]
    W_ped_to_core = [ped_to_core(IMAS.core_edge_energy(cp1d, 0.9)) for cp1d in dd.core_profiles.profiles_1d]
    W_ped_to_core_fraction = sum(W_ped_to_core[index_L]) / sum(index_L)

    if do_plot
        p = plot(; layout=(3, 1), size=(1000, 1000), link=:x)

        # Density
        plot!(; title="Electron density", subplot=1, legend_position=:bottomright, ylim=(0,Inf))
        plot!(dd.core_profiles.time, ne0; label="On axis density", lw=2, subplot=1)
        plot!(time, ne; label="Line average density", lw=2, subplot=1)
        plot!(dd.core_profiles.time, ne09; label="Density @rho=0.9", lw=2, subplot=1)
        plot!(time, ne_L; label="Line average density L-mode", lw=2, subplot=1)
        plot!(time, ne_H; label="Line average density H-mode", lw=2, subplot=1)
        vline!([transition_start, density_transition_end]; label="Density transition time", subplot=1)

        # Temperature
        plot!(; title="Electron temperature", subplot=2, link=:x, legend_position=:bottomright, ylim=(0,Inf))
        plot!(dd.core_profiles.time, te0; label="On axis temperature", lw=2, subplot=2)
        plot!(dd.core_profiles.time, te; label="Average temperature", lw=2, subplot=2)
        plot!(dd.core_profiles.time, te09; label="Temperature @rho=0.9", lw=2, subplot=2)
        vline!([transition_start, temperature_transition_end]; label="Temperature transition time", subplot=2)

        # L-H power threshold based on scaling law
        plot!(; title="L-H power thresold (no radiation)", subplot=3, xlabel="Time [s]", link=:x, legend_position=:bottomright, ylim=(0,Inf))
        injected_power = IMAS.total_power(dd.pulse_schedule, dd.core_profiles.time; tau_smooth=dd.core_profiles.time[2] - dd.core_profiles.time[1])
        scaling_power = [IMAS.scaling_L_to_H_power(dd; time0) for time0 in dd.core_profiles.time]
        plot!(dd.core_profiles.time, injected_power ./ scaling_power; label="Power / Power threshold", lw=2, subplot=3)
        vline!([transition_start, temperature_transition_end, density_transition_end]; label="Transition times", subplot=3)
        #        scatter!([dd.global_time], [IMAS.L_H_threshold(dd)], subplot=3, label="Power threshold (with radiation)")
        plot!(dd.core_profiles.time, 1.0./W_ped_to_core; subplot=3, label="Core to pedestal stored energy")

        display(p)
    end

    return (ne_L=ne_L, ne_H=ne_H, time=time, ne_L_over_H=Lvalue / Hvalue, tau_n=tau_n, tau_t=tau_t, W_ped_to_core_fraction = W_ped_to_core_fraction, mode_transitions=Dict{Float64,Symbol}(0.0 => :L_mode, transition_start => :H_mode))
end