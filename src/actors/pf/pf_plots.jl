#= ======== =#
#  plotting  #
#= ======== =#
@recipe function plot_ActorPFdesign_cx(actor::ActorPFdesign)
    @series begin
        actor.actor_pf
    end
end

"""
    plot_ActorPF_cx(
        actor::ActorPFactive{D,P};
        equilibrium=true,
        build=true,
        coils_flux=false,
        rails=false,
        plot_r_buffer=1.6) where {D<:Real,P<:Real}

Plot recipe for ActorPFdesign and ActorPFactive
"""
@recipe function plot_ActorPFactive_cx(
    actor::ActorPFactive{D,P};
    equilibrium=true,
    build=true,
    coils_flux=false,
    rails=true,
    control_points=true,
    plot_r_buffer=1.6) where {D<:Real,P<:Real}

    id = IMAS.recipe_id_for_help_plot(actor)
    HelpPlots.assert_type_and_record_argument(id, Bool, "Show equilibrium"; equilibrium)
    HelpPlots.assert_type_and_record_argument(id, Bool, "Show build"; build)
    HelpPlots.assert_type_and_record_argument(id, Bool, "Show coils flux"; coils_flux)
    HelpPlots.assert_type_and_record_argument(id, Bool, "Show rails"; rails)
    HelpPlots.assert_type_and_record_argument(id, Bool, "Show control points"; control_points)
    HelpPlots.assert_type_and_record_argument(id, Float64, "How much to buffer R axis of the plot to fit the legend"; plot_r_buffer)

    dd = actor.dd
    par = actor.par

    # if there is no equilibrium then treat this as a field_null plot
    eqt2d = findfirst(:rectangular, actor.eqt_out.profiles_2d)
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
            alpha --> 0.25
            label := false
            dd.build.layer
        end
        @series begin
            exclude_layers --> [:oh]
            wireframe := true
            dd.build.layer
        end
    end

    # plot coils_flux
    if coils_flux
        ngrid = 129
        R = range(xlim[1], xlim[2], ngrid)
        Z = range(ylim[1], ylim[2], Int(ceil(ngrid * (ylim[2] - ylim[1]) / (xlim[2] - xlim[1]))))

        coils = VacuumFields.GS_IMAS_pf_active__coil{D,D}[]
        for coil in dd.pf_active.coil
            if IMAS.is_ohmic_coil(coil)
                coil_tech = dd.build.oh.technology
            else
                coil_tech = dd.build.pf_active.technology
            end
            coil = VacuumFields.GS_IMAS_pf_active__coil(coil, coil_tech, par.green_model)
            push!(coils, coil)
        end

        # ψ coil currents
        ψbound = actor.eqt_out.global_quantities.psi_boundary
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

    # plot optimization rails
    if rails
        @series begin
            label --> (build ? "Coil opt. rail" : "")
            alpha --> 0.5
            dd.build.pf_active.rail
        end
    end

    # plot pf_active coils
    @series begin
        time0 --> actor.eqt_out.time
        dd.pf_active
    end

    # plot equilibrium
    if equilibrium
        if field_null
            pc = dd.pulse_schedule.position_control
            @series begin
                cx := true
                label --> "Field null region"
                color --> :red
                r, z = IMAS.boundary(pc, 1)
                r, z
            end
        else
            @series begin
                cx := true
                label --> "Final (λ_reg=$(round(log10(actor.λ_regularize);digits=1)))"
                color --> :red
                actor.eqt_out
            end
            @series begin
                cx := true
                label --> "Original"
                color --> :gray
                lcfs --> true
                lw := 1
                actor.dd.equilibrium.time_slice[]
            end
        end
    end

    # plot control points
    if control_points
        if !isempty(actor.iso_control_points)
            @series begin
                color := :blue
                linestyle := :dash
                linewidth := 1.5
                label := "Iso-flux constraints"
                [cpt.R1 for cpt in actor.iso_control_points], [cpt.Z1 for cpt in actor.iso_control_points]
            end
        end
        if !isempty(actor.flux_control_points)
            @series begin
                color := :blue
                seriestype := scatter
                markerstrokewidth := 0
                label := "Flux constraints"
                [cpt.R for cpt in actor.flux_control_points], [cpt.Z for cpt in actor.flux_control_points]
            end
        end
        if !isempty(actor.saddle_control_points)
            @series begin
                color := :blue
                seriestype := scatter
                markerstrokewidth := 0
                marker := :star
                label := "Saddle constraints"
                [cpt.R for cpt in actor.saddle_control_points], [cpt.Z for cpt in actor.saddle_control_points]
            end
        end
    end
end
