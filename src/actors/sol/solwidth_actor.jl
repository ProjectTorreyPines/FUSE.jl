#= ============= =#
#  ActorSOLWidth  #
#= ============= =#
@actor_parameters_struct ActorSOLWidth{T} begin
    model::Switch{Symbol} = Switch{Symbol}([:HD, :HD2T, :SSF, :SSF2T, :SSF2TFLUX], "-", "solwidth model (HD2T, SSF2T, SSF2TFLUX solve for Ti_sep; HD and SSF assume Ti=2Te)"; default=:SSF2T)
    solwidth_path::Entry{String} = Entry{String}("-", "Directory holding the `solwidth` package (e.g. CTCModeling/modules/sol-width-model), added to Python's sys.path. Leave empty if solwidth is pip-installed"; default="")
    Psol::Entry{T} = Entry{T}("W", "Power across the separatrix (defaults to IMAS.power_sol when not set)")
    ne_sep::Entry{T} = Entry{T}("m^-3", "Separatrix electron density given to the model (defaults to core_profiles ne at rho=1 when not set)")
    topology::Switch{Symbol} = Switch{Symbol}([:auto, :LSN, :USN, :DN], "-", "Magnetic topology (:auto derives it from the x-points lying on the LCFS)"; default=:auto)
    error_on_failure::Entry{Bool} = Entry{Bool}("-", "Throw an error (instead of a warning) when the solwidth solver does not converge"; default=false)
    verbose::Entry{Bool} = act_common_parameters(; verbose=false)
end

mutable struct ActorSOLWidth{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorSOLWidth{P}}
    inputs::Dict{String,Any}
    result::Dict{String,Any}
end

"""
    ActorSOLWidth(dd::IMAS.dd, act::ParametersAllActors; kw...)

Couples to the `solwidth` reduced SOL-width models (CTCModeling/modules/sol-width-model:
HD, HD2T, SSF, SSF2T, SSF2TFLUX) through PythonCall.

Global plasma parameters (Rmaj, Rmin, κ, δ, Ip, Btor, Psol and its electron fraction, Zeff, Aeff,
topology) and the separatrix density are taken from `dd` and solved by `solwidth.solve_one`.
The model returns the upstream (separatrix) and target temperatures and the SOL widths (λq, λn, ...).

Note that `ne_sep` is an *input* of the solwidth models: they predict `Te_sep`, `Ti_sep` at that density.

Requires `using PythonCall` (which loads the FUSE `PythonCallExt` extension) with a Python that has numpy and scipy, e.g.:

    ENV["JULIA_CONDAPKG_BACKEND"] = "Null"
    ENV["JULIA_PYTHONCALL_EXE"] = "/path/to/python"
    using PythonCall

!!! note

    Stores `Te_sep`, `Ti_sep` and `ne_sep` in `dd.edge_profiles.profiles_1d[].electrons/t_i_average` at
    `rho_pol_norm=1.0`. The full solwidth result is kept in `actor.result`.
    To feed it back to the pedestal set `act.ActorPedestal.Te_sep = actor.result["Te_sep"]`.
"""
function ActorSOLWidth(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorSOLWidth(dd, act.ActorSOLWidth; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorSOLWidth(dd::IMAS.dd{D}, par::FUSEparameters__ActorSOLWidth{P}; kw...) where {D<:Real,P<:Real}
    logging_actor_init(ActorSOLWidth)
    par = OverrideParameters(par; kw...)
    return ActorSOLWidth(dd, par, Dict{String,Any}(), Dict{String,Any}())
end

"""
    _solwidth_solve(model::AbstractString, pars::AbstractDict, solwidth_path::AbstractString)

Solve one solwidth parameter set and return the ModelResult as a Dict. Implemented in the PythonCallExt extension.
"""
function _solwidth_solve end

"""
    _solwidth_topology(eqt::IMAS.equilibrium__time_slice; tol::Float64=0.01)

Topology string for solwidth from the x-points that lie on the LCFS (within `tol` [m])
"""
function _solwidth_topology(eqt::IMAS.equilibrium__time_slice; tol::Float64=0.01)
    pr = eqt.boundary.outline.r
    pz = eqt.boundary.outline.z
    on_lcfs = [xp.z for xp in eqt.boundary.x_point if minimum(hypot.(pr .- xp.r, pz .- xp.z)) < tol]
    lower = any(z -> z < 0.0, on_lcfs)
    upper = any(z -> z > 0.0, on_lcfs)
    if lower && upper
        return "DN"
    elseif upper
        return "USN"
    else
        return "LSN"
    end
end

function _step(actor::ActorSOLWidth)
    dd = actor.dd
    par = actor.par

    if Base.get_extension(@__MODULE__, :PythonCallExt) === nothing
        error("ActorSOLWidth requires PythonCall: run `using PythonCall` (with a Python that has numpy and scipy) before using this actor")
    end

    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]

    Psol_par = getproperty(par, :Psol, missing)
    Psol = ismissing(Psol_par) ? IMAS.power_sol(dd.core_sources, cp1d) : Float64(Psol_par)

    total_source = IMAS.total_sources(dd.core_sources, cp1d; time0=dd.global_time, fields=[:power_inside, :total_ion_power_inside])
    power_electrons = total_source.electrons.power_inside[end]
    power_ions = total_source.total_ion_power_inside[end]
    Pfrac = clamp(power_electrons / (power_electrons + power_ions), 0.0, 1.0)

    ne_sep_par = getproperty(par, :ne_sep, missing)
    ne_sep = ismissing(ne_sep_par) ? cp1d.electrons.density_thermal[end] : Float64(ne_sep_par)

    r_out = eqt.profiles_1d.r_outboard[end]
    r_in = eqt.profiles_1d.r_inboard[end]

    topology = par.topology == :auto ? _solwidth_topology(eqt) : string(par.topology)

    actor.inputs = Dict{String,Any}(
        "name" => "FUSE t=$(dd.global_time)",
        "Rmaj" => (r_out + r_in) / 2.0,
        "Rmin" => (r_out - r_in) / 2.0,
        "kappa" => eqt.boundary.elongation,
        "delta" => eqt.boundary.triangularity,
        "Ip" => abs(eqt.global_quantities.ip),
        "Btor" => abs(IMAS.B0_geo(eqt)),
        "Psol" => Psol,
        "Pfrac" => Pfrac,
        "nsep" => ne_sep,
        "Zeff" => cp1d.zeff[end],
        "Aeff" => IMAS.A_effective(cp1d),
        "magnetic_topology" => topology
    )

    actor.result = _solwidth_solve(string(par.model), actor.inputs, par.solwidth_path)

    if !actor.result["success"]
        msg = "ActorSOLWidth: solwidth $(par.model) did not converge [$(strip(actor.result["message"]))], residual norm = $(actor.result["residual_norm"])"
        par.error_on_failure ? error(msg) : @warn msg
    end

    if par.verbose
        @info "ActorSOLWidth inputs: $(actor.inputs)"
        @info @sprintf("ActorSOLWidth %s: Te_sep=%.1f eV  Ti_sep=%.1f eV  Te_targ=%.1f eV  λq=%.2f mm  residual=%.2e",
            par.model, actor.result["Te_sep"], actor.result["Ti_sep"], actor.result["Te_targ"], actor.result["lambda_q"] * 1e3, actor.result["residual_norm"])
    end

    return actor
end

function _finalize(actor::ActorSOLWidth)
    dd = actor.dd
    result = actor.result
    ne_sep = actor.inputs["nsep"]

    ep1d = resize!(dd.edge_profiles.profiles_1d)
    if ismissing(ep1d.grid, :rho_pol_norm) || length(ep1d.grid.rho_pol_norm) == 1
        ep1d.grid.rho_pol_norm = [1.0]
        ep1d.electrons.density = [ne_sep]
        ep1d.electrons.temperature = [result["Te_sep"]]
        ep1d.t_i_average = [result["Ti_sep"]]
    else
        @assert ep1d.grid.rho_pol_norm[1] == 1.0 "ep1d.grid.rho_pol_norm = $(ep1d.grid.rho_pol_norm)"
        ep1d.electrons.density[1] = ne_sep
        ep1d.electrons.temperature[1] = result["Te_sep"]
        ep1d.t_i_average[1] = result["Ti_sep"]
    end

    return actor
end
