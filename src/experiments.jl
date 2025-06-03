"""
    LH_analysis(dd::IMAS.dd, transition_start::Float64=0.0, tau_n::Float64=0.4, tau_t::Float64=tau_n; do_plot::Bool=true)

This function takes the experimental data and returns time traces for the H-mode electron density and Zeff time traces
that would go in pulse_schedule. It also returns scalar factors `ne_L_over_H` and `zeff_L_over_H` that can be used to
obtain the L-mode traces by multiplying the H-mode traces.

Named tuple with:

  - `time::Vector{Float64}`: Time vector from core profiles.
  - `tau_n::Float64`: Duration [s] of the density transition.
  - `tau_t::Float64`: Duration [s] of the temperature transition.
  - `mode_transitions::Dict{Float64,Symbol}`: Dictionary of key transition times with labels `:L_mode` and `:H_mode`.
  - `ne_L::Vector{Float64}`: Smoothed L-mode line-averaged electron density.
  - `ne_H::Vector{Float64}`: Smoothed H-mode line-averaged electron density.
  - `ne_L_over_H::Float64`: Ratio of L-mode to H-mode density values.
  - `zeff_L::Vector{Float64}`: Smoothed L-mode Zeff.
  - `zeff_H::Vector{Float64}`: Smoothed H-mode Zeff.
  - `zeff_L_over_H::Float64`: Ratio of L-mode to H-mode Zeff values.
  - `W_ped_to_core_fraction::Float64`: Average ratio of pedestal to core stored energy during the L-mode phase.
"""
function LH_analysis(dd::IMAS.dd; scale_LH::Real=1.0, transition_start::Real=0.0, tau_n::Real=0.4, tau_t::Real=tau_n * 0.75, do_plot::Bool=true)
    rho = dd.core_profiles.profiles_1d[1].grid.rho_tor_norm
    index09 = argmin_abs(rho, 0.9)
    time = dd.core_profiles.time

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

    # powers for L-H threshold
    injected_power = IMAS.total_power(dd.pulse_schedule, time; tau_smooth=time[2] - time[1])
    scaling_power = [IMAS.scaling_L_to_H_power(dd; time0) for time0 in time]

    # LH detection
    is_H_mask = [false for time0 in time]

    # function to take the ratio of pedestal to core stored energy
    ped_to_core(W) = W[2] / W[1]

    # auto-detect mode transitions
    mode_transitions = Dict{Float64,Symbol}()
    mode_transitions[0.0] = :L_mode
    n = 1 # maximum number of H-mode transitions allowed
    for k in eachindex(time)
        if k == 1 || (n > 0 && is_H_mask[k-1] == true && injected_power[k] < scaling_power[k] * scale_LH * 1.1)
            is_H_mask[k] = false
        elseif n > 0 && is_H_mask[k-1] == false && injected_power[k] > scaling_power[k] * scale_LH * 0.9
            is_H_mask[k] = true
            n -= 1
        else
            is_H_mask[k] = is_H_mask[k-1]
        end
    end
    for (h_start, h_end) in get_true_ranges(is_H_mask)
        mode_transitions[time[h_start]] = :H_mode
        if h_end < length(time)
            mode_transitions[time[h_end]] = :L_mode
        end
    end

    if transition_start < 0.0 || length(mode_transitions) == 1
        ne_L = ne_H = ne
        zeff_L = zeff_H = zeff
        W_ped_to_core_fraction = sum([ped_to_core(IMAS.core_edge_energy(dd.core_profiles.profiles_1d[time0], 0.9)) for time0 in time]) / length(time)
        ne_L_over_H = 1.0
        zeff_L_over_H = 1.0
        transition_start = -1.0
    else
        if transition_start == 0.0
            transition_start = sort!(collect(keys(mode_transitions)))[2]
        elseif transition_start > 0.0
            mode_transitions[transition_start] = :H_mode
            is_H_mask = [false for time0 in time]
            is_H_mask[time.>=transition_start] .= true
        end
        smoothed_H_mask =
            (IMAS.smooth_beam_power(time, Float64.(is_H_mask), tau_n / 2) .+ reverse(IMAS.smooth_beam_power(-reverse(time), Float64.(reverse(is_H_mask)), tau_n / 2))) / 2.0
        index_H = smoothed_H_mask .>= 0.95
        index_L = smoothed_H_mask .<= 0.05
        index_H[1] = false
        index_L[1] = true

        if sum(index_L) == 0 || sum(index_H) == 0
            display(plot(time, injected_power ./ scaling_power; title="L-H power thresold (no radiation)", xlabel="Time [s]", xlim=(0.0, Inf)))
        end

        density_transition_end = transition_start + tau_n
        temperature_transition_end = transition_start + tau_t

        # smooth density based on L-H transition timescale 
        fs = 1 / (time[2] - time[1])
        band = 1 / tau_n
        ne_smoothed = IMAS.lowpassfilter(ne, fs, band)
        zeff_smoothed = IMAS.lowpassfilter(zeff09, fs, band)

        # ne_L/ne_H as a ne_L_over_H factor
        ne_L = IMAS.interp1d(time[index_L], ne_smoothed[index_L], :constant).(time)
        ne_L_over_H = sum(ne_L[index_H] / ne_smoothed[index_H]) / sum(index_H)
        index = sortperm([time[index_L]; time[index_H]])
        ne_H = IMAS.interp1d([time[index_L]; time[index_H]][index], [ne_smoothed[index_L] ./ ne_L_over_H; ne_smoothed[index_H]][index], :pchip).(time)
        ne_L = ne_H .* ne_L_over_H

        # zeff_L/zeff_H as a ne_L_over_H factor
        zeff_L = IMAS.interp1d(time[index_L], zeff_smoothed[index_L], :constant).(time)
        zeff_L_over_H = sum(zeff_L[index_H] / zeff_smoothed[index_H]) / sum(index_H)
        index = sortperm([time[index_L]; time[index_H]])
        zeff_H = IMAS.interp1d([time[index_L]; time[index_H]][index], [zeff_smoothed[index_L] ./ zeff_L_over_H; zeff_smoothed[index_H]][index], :pchip).(time)
        zeff_L = zeff_H .* zeff_L_over_H

        # WPED
        W_ped_to_core = [ped_to_core(IMAS.core_edge_energy(cp1d, 0.9)) for cp1d in dd.core_profiles.profiles_1d]
        W_ped_to_core_fraction = sum(W_ped_to_core[index_L]) / sum(index_L)
    end

    if do_plot
        p = plot(; layout=(4, 1), size=(1000, 1000), link=:x)

        background_color_legend = PlotUtils.Colors.RGBA(1.0, 1.0, 1.0, 0.6)

        # Density
        subplot = 1
        plot!(; title="Electron density", subplot, legend_position=:bottomright, xlim=(0.0, Inf), ylim=(0, Inf), background_color_legend)
        plot!(time, ne0; label="On axis density", lw=2, subplot)
        plot!(time, ne; label="Line average density", lw=2, subplot)
        plot!(time, ne09; label="Density @rho=0.9", lw=2, subplot)
        plot!(time, ne_L; label="Line average density L-mode", lw=2, subplot)
        plot!(time, ne_H; label="Line average density H-mode", lw=2, subplot)
        if transition_start >= 0.0
            vline!([transition_start, density_transition_end]; label="Density transition time", subplot)
        end

        # Zeff
        subplot = 2
        plot!(; title="Zeff", subplot, link=:x, legend_position=:bottomright, xlim=(0.0, Inf), ylim=(0, Inf), background_color_legend)
        plot!(time, zeff0; label="On axis Zeff", lw=2, subplot)
        plot!(time, zeff; label="Average Zeff", lw=2, subplot)
        plot!(time, zeff09; label="Zeff @rho=0.9", lw=2, subplot)
        plot!(time, zeff_L; label="Zeff L-mode", lw=2, subplot)
        plot!(time, zeff_H; label="Zeff H-mode", lw=2, subplot)
        if transition_start >= 0.0
            vline!([transition_start, density_transition_end]; label="Density transition time", subplot)
        end

        # Temperature
        subplot = 3
        plot!(; title="Electron temperature", subplot, link=:x, legend_position=:bottomright, xlim=(0.0, Inf), ylim=(0, Inf), background_color_legend)
        plot!(time, te0; label="On axis temperature", lw=2, subplot)
        plot!(time, te; label="Average temperature", lw=2, subplot)
        plot!(time, te09; label="Temperature @rho=0.9", lw=2, subplot)
        if transition_start >= 0.0
            vline!([transition_start, temperature_transition_end]; label="Temperature transition time", subplot)
        end

        # L-H power threshold based on scaling law
        subplot = 4
        plot!(; title="L-H power thresold (no radiation)", subplot, xlabel="Time [s]", link=:x, legend_position=:bottomright, xlim=(0.0, Inf), ylim=(0, Inf), background_color_legend)
        plot!(time, injected_power ./ scaling_power; label="Power / Power threshold", lw=2, subplot)
        if transition_start >= 0.0
            vline!([transition_start, temperature_transition_end, density_transition_end]; label="Transition times", subplot)
        end
        hline!([scale_LH * 0.9, scale_LH * 1.1]; subplot, label="")
        true_ranges = get_true_ranges(is_H_mask)
        for (start_idx, end_idx) in true_ranges
            vspan!([time[start_idx-1], time[end_idx-1]]; color=:lightgray, alpha=0.5, label=false, subplot)
        end

        display(p)
    end

    return (time=time, tau_n=tau_n, tau_t=tau_t,
        mode_transitions=mode_transitions,
        ne_L=ne_L, ne_H=ne_H, ne_L_over_H=ne_L_over_H,
        zeff_L=zeff_L, zeff_H=zeff_H, zeff_L_over_H=zeff_L_over_H,
        W_ped_to_core_fraction=W_ped_to_core_fraction)
end

function get_true_ranges(mask)
    ranges = []
    i = 1
    while i <= length(mask)
        if mask[i]
            start_idx = i
            while i <= length(mask) && mask[i]
                i += 1
            end
            end_idx = i
            push!(ranges, (start_idx, end_idx))
        else
            i += 1
        end
    end
    return ranges
end