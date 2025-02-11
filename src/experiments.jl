function ne_line_without_LH_transition(dd::IMAS.dd, transition_start::Float64, density_transition_end::Float64, temperature_transition_end::Float64; do_plot::Bool=true)
    rho = dd.core_profiles.profiles_1d[1].grid.rho_tor_norm
    time = dd.core_profiles.time

    # density timescale
    ne = IMAS.interp1d(dd.pulse_schedule.density_control.time, dd.pulse_schedule.density_control.n_e_line.reference).(time)
    tau_n = density_transition_end - transition_start

    # temperature timescale
    te = [IMAS.trapz(rho, cp1d.electrons.temperature) for cp1d in dd.core_profiles.profiles_1d]
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

    if do_plot
        p = plot(; layout=(3, 1), size=(1000, 1000), link=:x)

        # Density
        plot!(; title="Line average density", subplot=1)
        plot!(time, ne; label="Density data", lw=2, subplot=1)
        plot!(time, ne_L; label="Density L-mode", lw=2, subplot=1)
        plot!(time, ne_H; label="Density H-mode", lw=2, subplot=1)
        vline!([transition_start, density_transition_end]; label="Density transition time", subplot=1)
        vline!([dd.global_time]; label="dd.global_time", ls=:dash, subplot=1)
        plot!(dd.pulse_schedule.density_control.n_e_line, :reference; subplot=1, xlabel="")

        # Temperature
        plot!(; title="Average electron temperature", subplot=2, link=:x)
        plot!(dd.core_profiles.time, te; label="Temperature data", lw=2, subplot=2)
        vline!([transition_start, temperature_transition_end]; label="Temperature transition time", subplot=2)
        vline!([dd.global_time]; label="dd.global_time", ls=:dash, subplot=2)

        # L-H power threshold based on scaling law
        plot!(; title="L-H power thresold", subplot=3, xlabel="Time [s]", link=:x)
        injected_power = IMAS.total_power(dd.pulse_schedule, dd.core_profiles.time; tau_smooth=dd.core_profiles.time[2] - dd.core_profiles.time[1])
        scaling_power = [IMAS.scaling_L_to_H_power(dd; time0) for time0 in dd.core_profiles.time]
        plot!(dd.core_profiles.time, injected_power ./ scaling_power; label="Power threshold", lw=2, subplot=3)
        vline!([transition_start]; label="Transition start time", subplot=3)
        vline!([dd.global_time]; label="dd.global_time", ls=:dash, subplot=3)
        display(p)
    end

    return (L=ne_L, H=ne_H, time=time, L_over_H=Lvalue / Hvalue, tau_n=tau_n, tau_t=tau_t)
end