import TurbulentTransport: TurbulenceMode, AbstractModeIdentification, TJLFModeIdentification, NNModeIdentification, identify_modes, run_modeid_nn, MODE_COLORS, MODE_LABELS
using LaTeXStrings

#= =========== =#
#  ActorModeID  #
#= =========== =#
@actor_parameters_struct ActorModeID{T} begin
    model::Switch{Symbol} = Switch{Symbol}([:TJLF, :ModeIDNN], "-", "Mode identification method: :TJLF (full quasilinear analysis) or :ModeIDNN (neural network, < 1 ms)"; default=:TJLF)
    modeid_model::Entry{String} = Entry{String}("-", "ModeID NN model filename (BSON), only used when model=:ModeIDNN"; default="modeid_qlgyro_sat3_azf-1")
    rho_transport::Entry{AbstractVector{T}} = Entry{AbstractVector{T}}("-", "ρ transport grid for mode identification"; default=0.25:0.1:0.85)
    sat_rule::Switch{Symbol} = Switch{Symbol}([:sat0, :sat0quench, :sat1, :sat1geo, :sat2, :sat3], "-", "Saturation rule"; default=:sat3)
    electromagnetic::Entry{Bool} = Entry{Bool}("-", "Electromagnetic or electrostatic"; default=true)
    lump_ions::Entry{Bool} = Entry{Bool}("-", "Lumps the fuel species (D,T) as well as the impurities together"; default=true)
    MXH_modes::Entry{Int} = Entry{Int}("-", "Number of MXH harmonics"; default=1)
    warn_nn_train_bounds::Entry{Bool} = Entry{Bool}("-", "Warn if NN inputs are outside training bounds (only for ModeIDNN)"; default=false)
    em_threshold::Entry{T} = Entry{T}("-", "EM/ES QL weight ratio threshold for MTM/KBM vs ITG/TEM/ETG classification (only for TJLF)"; default=0.5)
    ion_electron_threshold::Entry{T} = Entry{T}("-", "Ion/electron ES QL weight ratio threshold for TEM vs ETG classification (only for TJLF)"; default=0.5)
    ky_etg::Entry{T} = Entry{T}("-", "ky threshold above which ES electron-direction modes are classified as ETG (only for TJLF)"; default=2.0)
end

mutable struct ActorModeID{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorModeID{P}}
    labeled_dds::Vector{Pair{String,IMAS.dd{D}}}
    results::Vector{Pair{String,Vector{<:AbstractModeIdentification{D}}}}
end

"""
    ActorModeID(dd::IMAS.dd, act::ParametersAllActors; kw...)
    ActorModeID(labeled_dds::AbstractVector{<:Pair{String,<:IMAS.dd}}, act::ParametersAllActors; kw...)

Identifies the dominant turbulence mode (ITG, TEM, KBM, ETG, MTM) driving heat flux
at each radial location.

Supports two methods:
- `model=:TJLF` (default): Full TJLF quasilinear analysis — runs TJLF at each radial point,
  classifies modes from QL weight ratios and frequencies. Accurate but slow (~seconds per ρ).
- `model=:ModeIDNN`: Neural network classification — directly predicts the dominant mode from
  TGLF input parameters in < 1 ms for all radial points. Uses ensemble of residual networks
  trained on QLGYRO SAT3 mode identification data.

Accepts either a single `dd` or a vector of `"label" => dd` pairs for side-by-side comparison
of mode identification across different equilibria (e.g. initial vs flux-matched).

Results are visualized with `plot(actor)`, which shows color-coded markers
(ITG=green, TEM=orange, KBM=violet, ETG=blue, MTM=red) at each radial location,
stacked vertically when comparing multiple equilibria.

# Examples

```julia
# Single dd with TJLF (default, slow)
actor = ActorModeID(dd, act)

# Single dd with ModeID NN (fast)
act.ActorModeID.model = :ModeIDNN
actor = ActorModeID(dd, act)

# Compare initial vs flux-matched
actor = ActorModeID(["initial" => dd, "flux-matched" => dd_tjlf], act)
```
"""
function ActorModeID(dd::IMAS.dd{D}, act::ParametersAllActors; kw...) where {D<:Real}
    actor = ActorModeID(dd, act.ActorModeID; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorModeID(labeled_dds::AbstractVector{<:Pair{String,<:IMAS.dd{D}}}, act::ParametersAllActors; kw...) where {D<:Real}
    dd_first = first(labeled_dds).second
    logging_actor_init(ActorModeID)
    par = OverrideParameters(act.ActorModeID; kw...)
    entries = Pair{String,IMAS.dd{D}}[l => d for (l, d) in labeled_dds]
    actor = ActorModeID(dd_first, par, entries, Pair{String,Vector{<:AbstractModeIdentification{D}}}[])
    step(actor)
    finalize(actor)
    return actor
end

function ActorModeID(dd::IMAS.dd{D}, par::FUSEparameters__ActorModeID; kw...) where {D<:Real}
    logging_actor_init(ActorModeID)
    par = OverrideParameters(par; kw...)
    entries = Pair{String,IMAS.dd{D}}["" => dd]
    return ActorModeID(dd, par, entries, Pair{String,Vector{<:AbstractModeIdentification{D}}}[])
end

"""
    _step(actor::ActorModeID)

Identifies the dominant turbulence mode at each `rho_transport` location.

When `model=:TJLF`: constructs TJLF inputs, runs full quasilinear analysis, classifies from QL weights.
When `model=:ModeIDNN`: constructs TGLF inputs, runs neural network forward pass (< 1 ms).
"""
function _step(actor::ActorModeID{D,P}) where {D<:Real,P<:Real}
    par = actor.par
    empty!(actor.results)

    if par.model == :TJLF
        for (label, dd) in actor.labeled_dds
            input_tglfs = InputTGLF(dd, par.rho_transport, par.sat_rule, par.electromagnetic, par.lump_ions; MXH_modes=par.MXH_modes)
            input_tjlfs = [InputTJLF{D}(input_tglfs[k]) for k in eachindex(par.rho_transport)]

            mode_ids = identify_modes(input_tjlfs;
                em_threshold=par.em_threshold,
                ion_electron_threshold=par.ion_electron_threshold,
                ky_etg=par.ky_etg)

            push!(actor.results, label => mode_ids)
        end
    elseif par.model == :ModeIDNN
        for (label, dd) in actor.labeled_dds
            mode_ids = run_modeid_nn(dd, par.rho_transport;
                model_filename=par.modeid_model,
                warn_nn_train_bounds=par.warn_nn_train_bounds,
                MXH_modes=par.MXH_modes)

            push!(actor.results, label => mode_ids)
        end
    end

    return actor
end

const _MODE_MARKER_SHAPES = [:circle, :diamond, :square, :utriangle, :dtriangle, :hexagon, :star5]

function _mode_color(mode::TurbulenceMode, alpha::Real=0.5)
    c = Plots.Colors.parse(Plots.Colors.Colorant, MODE_COLORS[mode])
    return Plots.Colors.RGBA(Plots.Colors.red(c), Plots.Colors.green(c), Plots.Colors.blue(c), alpha)
end

@recipe function plot_ActorModeID(actor::ActorModeID)
    par = actor.par
    rho = collect(par.rho_transport)
    n_results = length(actor.results)

    size --> (900, 150 + 100 * n_results)
    xlabel --> L"\rho"
    title --> "Dominant Turbulence Mode"
    xlims --> (0.0, 1.0)
    legend --> :outerright
    grid --> false
    framestyle --> :box

    ylims --> (0.3, n_results + 0.7)
    yticks --> (1:n_results, [p.first == "" ? "dd" : p.first for p in actor.results])

    legend_shown = Set{TurbulenceMode}()

    for (row, (label, mode_ids)) in enumerate(actor.results)
        for mode in instances(TurbulenceMode)
            mask = [mid.dominant_mode == mode for mid in mode_ids]
            any(mask) || continue
            show_label = mode ∉ legend_shown
            if show_label
                push!(legend_shown, mode)
            end
            @series begin
                seriestype := :scatter
                label --> (show_label ? MODE_LABELS[mode] : "")
                markercolor --> _mode_color(mode)
                markersize --> 8
                markershape --> :circle
                markerstrokewidth --> 1
                markerstrokecolor --> :white
                rho[mask], fill(row, sum(mask))
            end
        end
    end
end

"""
    plot(actor::ActorModeID, Val(:profiles))

Four-panel plot of electron temperature, electron density, ion temperature, and toroidal
rotation with mode-colored scatter markers superimposed at each `rho_transport` location.

When multiple `dd` objects are provided, profile lines are differentiated by line style and
mode markers by marker shape. Color always encodes the dominant turbulence mode.
"""
function Plots.plot(actor::ActorModeID, ::Val{:profiles}; kwargs...)
    par = actor.par
    rho_transport = collect(par.rho_transport)

    line_styles = [:solid, :dash, :dot, :dashdot, :dashdotdot]
    ylabels = [L"T_e~[\mathrm{eV}]", L"n_e~[\mathrm{m}^{-3}]", L"T_i~[\mathrm{eV}]", L"\omega_\mathrm{tor}~[\mathrm{rad/s}]"]

    p = plot(layout=(2, 2), size=(1000, 700), margin=5 * Plots.Measures.mm; kwargs...)

    for (dd_idx, (label, dd)) in enumerate(actor.labeled_dds)
        cp1d = dd.core_profiles.profiles_1d[]
        rho = cp1d.grid.rho_tor_norm
        lbl = label == "" ? "dd" : label
        ls = line_styles[mod1(dd_idx, length(line_styles))]

        Te = ismissing(cp1d.electrons, :temperature) ? fill(NaN, length(rho)) : cp1d.electrons.temperature
        ne = ismissing(cp1d.electrons, :density) ? fill(NaN, length(rho)) : cp1d.electrons.density
        Ti = (!isempty(cp1d.ion) && !ismissing(cp1d.ion[1], :temperature)) ? cp1d.ion[1].temperature : fill(NaN, length(rho))
        vtor = (!isempty(cp1d.ion) && !ismissing(cp1d.ion[1], :rotation_frequency_tor)) ? cp1d.ion[1].rotation_frequency_tor : fill(NaN, length(rho))

        panels = [Te, ne, Ti, vtor]

        for (sp, profile) in enumerate(panels)
            plot!(p, rho, profile;
                subplot=sp,
                ylabel=ylabels[sp],
                xlabel=L"\rho",
                label=lbl,
                linestyle=ls,
                xlims=(0.0, 1.0),
                ylims=(0.0, Inf))
        end
    end

    legend_shown = Set{TurbulenceMode}()

    for (res_idx, ((label, mode_ids), (_, dd))) in enumerate(zip(actor.results, actor.labeled_dds))
        cp1d = dd.core_profiles.profiles_1d[]
        rho = cp1d.grid.rho_tor_norm

        Te = ismissing(cp1d.electrons, :temperature) ? fill(NaN, length(rho)) : cp1d.electrons.temperature
        ne = ismissing(cp1d.electrons, :density) ? fill(NaN, length(rho)) : cp1d.electrons.density
        Ti = (!isempty(cp1d.ion) && !ismissing(cp1d.ion[1], :temperature)) ? cp1d.ion[1].temperature : fill(NaN, length(rho))
        vtor = (!isempty(cp1d.ion) && !ismissing(cp1d.ion[1], :rotation_frequency_tor)) ? cp1d.ion[1].rotation_frequency_tor : fill(NaN, length(rho))

        profiles = [Te, ne, Ti, vtor]

        for mode in instances(TurbulenceMode)
            mask = [mid.dominant_mode == mode for mid in mode_ids]
            any(mask) || continue
            show_label = mode ∉ legend_shown
            if show_label
                push!(legend_shown, mode)
            end

            rho_vals = rho_transport[mask]
            nearest_idxs = [argmin_abs(rho, r) for r in rho_vals]

            for (sp, profile) in enumerate(profiles)
                y_vals = profile[nearest_idxs]
                scatter!(p, rho_vals, y_vals;
                    subplot=sp,
                    label=(sp == 1 && show_label ? MODE_LABELS[mode] : ""),
                    markercolor=_mode_color(mode),
                    markersize=6,
                    markershape=:circle,
                    markerstrokewidth=1,
                    markerstrokecolor=:white)
            end
        end
    end

    return p
end

function Base.show(io::IO, ::MIME"text/plain", actor::ActorModeID)
    println(io, "ActorModeID ($(actor.par.model)):")
    if isempty(actor.results)
        println(io, "  (not yet run)")
        return
    end
    rho = collect(actor.par.rho_transport)
    frac_label = actor.par.model == :ModeIDNN ? "probability" : "of |flux|"
    for (label, mode_ids) in actor.results
        header = label == "" ? "dd" : label
        println(io, "  [$header]")
        for (k, mid) in enumerate(mode_ids)
            frac = round(mid.dominant_mode_fraction * 100; digits=1)
            println(io, "    ρ=$(round(rho[k]; digits=3)): $(MODE_LABELS[mid.dominant_mode]) ($(frac)% $frac_label)")
        end
    end
end
