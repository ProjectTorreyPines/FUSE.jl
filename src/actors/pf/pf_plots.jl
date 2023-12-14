#= ======== =#
#  plotting  #
#= ======== =#
"""
    plot_ActorPF_cx(
        actor::Union{ActorPFactive{D,P},ActorPFcoilsOpt{D,P}};
        time_index=nothing,
        equilibrium=true,
        build=true,
        coils_flux=false,
        rail=false,
        plot_r_buffer=1.6) where {D<:Real,P<:Real}

Plot recipe for ActorPFcoilsOpt and ActorPFactive
"""
@recipe function plot_ActorPF_cx(
    actor::Union{ActorPFactive{D,P},ActorPFcoilsOpt{D,P}};
    time_index=nothing,
    equilibrium=true,
    build=true,
    coils_flux=false,
    rail=false,
    plot_r_buffer=1.6) where {D<:Real,P<:Real}

    @assert typeof(time_index) <: Union{Nothing,Integer}
    @assert typeof(equilibrium) <: Bool
    @assert typeof(build) <: Bool
    @assert typeof(coils_flux) <: Bool
    @assert typeof(rail) <: Bool
    @assert typeof(plot_r_buffer) <: Real

    dd = actor.dd
    par = actor.par

    if time_index === nothing
        time_index = findfirst(x -> x.time == dd.global_time, actor.eq_out.time_slice)
    end
    time0 = actor.eq_out.time_slice[time_index].time

    # if there is no equilibrium then treat this as a field_null plot
    eqt2d = findfirst(:rectangular, actor.eq_out.time_slice[time_index].profiles_2d)
    field_null = false
    if eqt2d === nothing || ismissing(eqt2d, :psi)
        coils_flux = equilibrium
        field_null = true
    end

    # when plotting coils_flux the build is not visible anyways
    if coils_flux
        build = false
    end

    # setup plotting area
    xlim = [0.0, maximum(dd.build.layer[end].outline.r)]
    ylim = [minimum(dd.build.layer[end].outline.z), maximum(dd.build.layer[end].outline.z)]
    xlim --> xlim * plot_r_buffer
    ylim --> ylim
    aspect_ratio --> :equal

    # plot build
    if build
        @series begin
            exclude_layers --> [:oh]
            dd.build
        end
    end

    # plot coils_flux
    if coils_flux
        ngrid = 129
        R = range(xlim[1], xlim[2], ngrid)
        Z = range(ylim[1], ylim[2], Int(ceil(ngrid * (ylim[2] - ylim[1]) / (xlim[2] - xlim[1]))))

        coils = GS_IMAS_pf_active__coil{D,D}[]
        for coil in dd.pf_active.coil
            if IMAS.is_ohmic_coil(coil)
                coil_tech = dd.build.oh.technology
            else
                coil_tech = dd.build.pf_active.technology
            end
            coil = GS_IMAS_pf_active__coil(coil, coil_tech, par.green_model)
            coil.time_index = time_index
            push!(coils, coil)
        end

        # ψ coil currents
        ψbound = actor.eq_out.time_slice[time_index].global_quantities.psi_boundary
        ψ = [sum(VacuumFields.ψ(coil, r, z; Bp_fac=2π) for coil in coils) for r in R, z in Z]

        ψmin = minimum(x -> isnan(x) ? Inf : x, ψ)
        ψmax = maximum(x -> isnan(x) ? -Inf : x, ψ)
        ψabsmax = maximum(x -> isnan(x) ? -Inf : x, abs.(ψ))

        if field_null
            clims = (-ψabsmax / 10 + ψbound, ψabsmax / 10 + ψbound)
        else
            clims = (ψmin, ψmax)
        end

        @series begin
            seriestype --> :contourf
            c --> :diverging
            colorbar_entry --> false
            levels --> range(clims[1], clims[2], 21)
            linewidth --> 0.0
            R, Z, transpose(ψ)
        end

        if field_null
            @series begin
                seriestype --> :contour
                colorbar_entry --> false
                levels --> [ψbound]
                linecolor --> :gray
                R, Z, transpose(ψ)
            end
        end

        @series begin
            wireframe --> true
            exclude_layers --> [:oh]
            dd.build
        end
    end

    # plot equilibrium
    if equilibrium
        if field_null
            pc = dd.pulse_schedule.position_control
            @series begin
                cx := true
                label --> "Field null region"
                seriescolor --> :red
                IMAS.boundary(pc, 1)
            end
        else
            @series begin
                cx := true
                label --> "Final"
                seriescolor --> :red
                actor.eq_out.time_slice[time_index]
            end
            @series begin
                cx := true
                label --> "Target"
                seriescolor --> :blue
                lcfs --> true
                linestyle --> :dash
                actor.dd.equilibrium.time_slice[time_index]
            end
        end
    end

    # plot pf_active coils
    @series begin
        time0 --> time0
        dd.pf_active
    end

    # plot optimization rails
    if rail
        @series begin
            label --> (build ? "Coil opt. rail" : "")
            dd.build.pf_active.rail
        end
    end

end

"""
    plot_ActorPFcoilsOpt_trace(
        trace::PFcoilsOptTrace,
        what::Symbol=:cost;
        start_at=1)

Plot recipe for ActorPFcoilsOpt optimization trace

Attributes:

  - what::Symbol=:cost or :currents or individual fields of the PFcoilsOptTrace structure
  - start_at=::Int=1 index of the first element of the trace to start plotting
"""
@recipe function plot_ActorPFcoilsOpt_trace(
    trace::PFcoilsOptTrace,
    what::Symbol=:cost;
    start_at=1)

    @assert typeof(start_at) <: Integer

    start_at = minimum([start_at, length(trace.cost_total)])
    x = start_at:length(trace.cost_total)
    legend --> :bottomleft
    if what == :cost
        if sum(trace.cost_lcfs[start_at:end]) > 0.0
            data = trace.cost_lcfs[start_at:end]
            index = data .> 0.0
            @series begin
                label --> "ψ"
                yscale --> :log10
                x[index], data[index]
            end
        end
        if sum(trace.cost_currents[start_at:end]) > 0.0
            data = trace.cost_currents[start_at:end]
            index = data .> 0.0
            @series begin
                label --> "currents"
                yscale --> :log10
                x[index], data[index]
            end
        end
        if sum(trace.cost_oh[start_at:end]) > 0.0
            data = trace.cost_oh[start_at:end]
            index = data .> 0.0
            @series begin
                label --> "oh"
                yscale --> :log10
                x[index], data[index]
            end
        end
        if sum(trace.cost_1to1[start_at:end]) > 0.0
            data = trace.cost_1to1[start_at:end]
            index = data .> 0.0
            @series begin
                label --> "1to1"
                yscale --> :log10
                x[index], data[index]
            end
        end
        if sum(trace.cost_spacing[start_at:end]) > 0.0
            data = trace.cost_spacing[start_at:end]
            index = data .> 0.0
            @series begin
                label --> "spacing"
                yscale --> :log10
                x[index], data[index]
            end
        end
        @series begin
            label --> "total"
            yscale --> :log10
            linestyle --> :dash
            color --> :black
            # ylim --> [minimum(trace.cost_total[start_at:end]) / 10,maximum(trace.cost_total[start_at:end])]
            x, trace.cost_total[start_at:end]
        end

    elseif what == :params
        nparams = length(getfield(trace, what)[1]) - 1

        for k in 1:nparams
            @series begin
                label --> "#$k"
                x, [getfield(trace, what)[i][k] for i in eachindex(trace.cost_total)][start_at:end]
            end
        end

    else
        @series begin
            if occursin("cost_", String(what))
                yscale --> :log10
            end
            label --> String(what)
            x, getfield(trace, what)[start_at:end]
        end
    end
end
