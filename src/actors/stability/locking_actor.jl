import Roots
using Plots
using LaTeXStrings
import ModeLocking
using ModeLocking: ODEparams, NNparams, LockingNNModel, LockingResults
#import FUSE: coordinates


#================== =#
#  ActorLockingProbability  #
#================== =#
Base.@kwdef mutable struct FUSEparameters__ActorLocking{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    m_pol::Entry{Int} = Entry{Int}("_", "poloidal mode number of the mode"; default=2)
    n_tor::Entry{Int} = Entry{Int}("_", "toroidal mode number of the mode"; default=1)
    grid_size::Entry{Int} = Entry{Int}("-", "grid resolution for control space"; default=100)
    t_final::Entry{Float64} = Entry{Float64}("-", "Final integration time in units of tearing time (~ms)"; default=100.)
    time_steps::Entry{Int} = Entry{Int}("-", "number of time steps for the ODE integration"; default=200)
    overwrite_params::Entry{Bool} = Entry{Bool}("-", "Whether to overwrite ODE parameters to reproduce PoP2024 results"; default=false)
    control_type::Switch{Symbol} = Switch{Symbol}([:EF, :LinStab, :NLsaturation], # EF: error field
        "-",                                                            # LinStab: vary stability_index,
        "Use a user specified Control case to run the locking models"; default=:EF) # NLsaturation: vary NL saturation
    task::Switch{Symbol} = Switch{Symbol}(
        [:solve_system, :single_case, :calc_prob, :Monte_Carlo, :eval_prob, :transfer_learning],
        "-",
        "Solve the system on the full grid (:solve_system), run a single case (:single_case), build the probability model (:calc_prob — engine chosen by prob_method), evaluate it at operating points (:eval_prob), fine-tune it (:transfer_learning), or Monte Carlo (NOT implemented)";
        default = :solve_system
    )
    application::Switch{String} = Switch{String}(
        ["RP-RW", "RP-IW", "RP-RP-RW", "RP-RW-IW"],
        "-", 
        "Type of application: 'RP-RW' for resistive plasma with a single rational surface interacting with a resistive wall; 
                              'RP-IW' for resistive plasma with a single rational surface interacting with an ideal wall;
                              'RP-RP-RW' for resistive plasma with two rational surfaces interacting with a resistive wall;
                              'RP-RW-IW' for resistive plasma with a single rational surface interacting with both walls"; 
        default="RP-RW"
    )
    EF_phase::Entry{Float64} = Entry{Float64}("-", "Phase of the applied error field (degrees)"; default=0.)
    NL_saturation::Entry{Bool} = Entry{Bool}("-", "Nonlinear saturation parameter for the mode"; default=false)
    RPRW_stability_index::Entry{Float64} = Entry{Float64}(
        "-", 
        "Stability index of the system (set to Neg. value for now)"; default=-0.5)
    b0::Entry{Float64} = Entry{Float64}(
        "Tesla", 
        "Scale for magnetic perturbations, usually ~10Gauss"; default=1.e-3)
    t0::Entry{Float64} = Entry{Float64}(
        "seconds", 
        "Characteristic time scale for normalization , usually TM/RW growth rate"; default=1.e-3)
    r0::Entry{Float64} = Entry{Float64}(
        "meter", 
        "Length scale for the integration , usually minor radius"; default=1.)
    source_torque::Entry{Float64} = Entry{Float64}("Newton.meter", "NBI torque"; default=5.)
    error_field::Entry{Float64} = Entry{Float64}(
        "Gauss",
        "Fixed n=1 error field. Used as the ODE right-hand side's ε whenever the error " *
        "field is NOT the swept control (control_type=:LinStab/:NLsaturation); for " *
        "control_type=:EF the swept Control2 supplies ε instead and this is ignored. " *
        "This is the only place the error field is set — ode_params.error_field is " *
        "derived from it and holds the dimensionless form ModeLocking reads";
        default=10.0)
    plot_orientation::Switch{Symbol} = Switch{Symbol}([:portrait, :landscape], "-",
        "Tile-plot layout: :portrait = 3×2 (paper), :landscape = 2×3 (slides)"; default=:portrait)
    save_plots::Entry{Bool} = Entry{Bool}(
        "-",
        "Write the locking figures to disk at the end of any task that produces results " *
        "(:solve_system, :calc_prob, :eval_prob, :transfer_learning). Which figures appear " *
        "depends on what the task produced: the probability map needs a trained model, and " *
        "the operating-point markers need :eval_prob";
        default=false)
    plots_dir::Entry{String} = Entry{String}(
        "-",
        "Directory for save_plots output; empty = current working directory";
        default="")
    op_times::Entry{Vector{Float64}} = Entry{Vector{Float64}}(
        "s",
        "Times at which to evaluate locking probability; empty = use current dd.global_time";
        default=Float64[])
    op_C1::Entry{Float64} = Entry{Float64}(
        "kHz",
        "Operating-point rotation frequency for :single_case; NaN = use source_torque default";
        default=NaN)
    op_C2::Entry{Union{Float64,Vector{Float64}}} = Entry{Union{Float64,Vector{Float64}}}(
        "-",
        "Operating-point Control2 — scalar or one entry per op_times (a scalar is broadcast). " *
        "Units follow control_type: Gauss (error field) for :EF; Gauss (n=1 br amplitude at the " *
        "rational surface, inverted to Δt) for :LinStab; native α for :NLsaturation. " *
        "NaN = use ode_params default";
        default=NaN)
    P_locked_threshold::Entry{Float64} = Entry{Float64}(
        "-",
        "Locking probability above which the limit is considered breached";
        default=0.5)
    prob_method::Switch{Symbol} = Switch{Symbol}([:nn, :conv, :kde], "-",
        "Engine that turns the classified ODE grid into P(locked): neural net, windowed convolution, or KDE";
        default=:nn)
    conv_window_C1::Entry{Int} = Entry{Int}(
        "-",
        "Window size along the C1 (Ω₀) axis for prob_method=:conv/:kde (must be odd); ignored by :nn";
        default=5)
    conv_window_C2::Entry{Int} = Entry{Int}(
        "-",
        "Window size along the C2 (EF/Δ'/α) axis for prob_method=:conv/:kde (must be odd); ignored by :nn";
        default=5)
end


mutable struct ActorLocking{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorLocking{P}}
    ode_params::Union{Nothing, ODEparams}
    results::Union{Nothing, LockingResults}
    nn_params::NNparams
    # Operating points resolved by :eval_prob — `nothing` until it runs.
    # Grouped rather than kept as three parallel vectors so they cannot disagree
    # in length, and so the names do not collide with par.op_C1/op_C2, which are
    # user *inputs* in different units (par.op_C1 is a scalar in kHz).
    #   times : evaluated times [s]
    #   C1    : dimensionless rotation (f·t₀) at each time, computed from dd
    #   C2    : Control2 in user-facing units (Gauss for :EF and :LinStab)
    eval::Union{Nothing,@NamedTuple{times::Vector{Float64}, C1::Vector{Float64}, C2::Vector{Float64}}}

    function ActorLocking(
        dd::IMAS.dd{D},
        par::FUSEparameters__ActorLocking{P};
        ode_params = nothing,
        nn_params  = NNparams(),
        kw...
    ) where {D<:Real,P<:Real}

        logging_actor_init(ActorLocking)

        # Apply standard FUSE parameter overrides
        par = OverrideParameters(par; kw...)

        # Handle ODE params. par.error_field (Gauss) is the single user-facing
        # error field; ode_params.error_field is derived from it and holds the
        # dimensionless form ModeLocking's resolve_control reads.
        ode = if ode_params === nothing
            ODEparams()
        elseif ode_params isa NamedTuple
            ODEparams(; ode_params...)
        elseif ode_params isa ODEparams
            ode_params
        else
            error("ode_params must be nothing, ODEparams, or NamedTuple")
        end
        # Safe on every path, including a reused ODEparams: the value is assigned
        # from par.error_field rather than scaled in place, so it cannot compound
        _normalize_error_field!(ode, par)

        return new{D,P}(dd, par, ode, nothing, nn_params, nothing)
    end
end


"""
    ActorLocking(dd::IMAS.dd, act::ParametersAllActors; kw...)

Provides a common interface to run different equilibrium actors
"""
function ActorLocking(dd::IMAS.dd, act::ParametersAllActors; kw...)
    #par = act.ActorLocking(kw...)
    actor = ActorLocking(dd, act.ActorLocking; kw...)
    step(actor)
    finalize(actor)
    return actor
end



function _step(actor::ActorLocking)
    dd = actor.dd
    par = actor.par
    task = par.task
    application = par.application

    # Populate the physical parameters needed to solve the ODEs. ode_params is
    # always set by the constructor (with error_field already normalized), so
    # building a bare ODEparams here would silently reintroduce a Gauss-valued
    # error_field into the model.
    actor.ode_params === nothing && error("actor.ode_params is unset — construct the actor through ActorLocking(dd, par; ode_params=...)")
    actor.ode_params = set_up_ode_params!(dd, par, actor.ode_params)

    # Main driver routine
    # Time evolve the ODEs, calculate locking probability, or load/evaluate model
    if task == :single_case
        # Default the operating rotation to what the plasma is actually doing at
        # dd.global_time, same as :eval_prob; an explicit op_C1 (kHz) overrides
        C1 = if isnan(par.op_C1)
            c1 = _rotation_at_rat_surface(dd, par, actor.ode_params.rat_surface)
            @info "op_C1 not set — using rotation from dd at ρ=$(round(actor.ode_params.rat_surface; sigdigits=4)): C1=$(round(c1; sigdigits=4))"
            c1
        else
            par.op_C1 * 1e3 * par.t0   # kHz → dimensionless f·t₀
        end
        c2_user = _op_C2_scalar(par)   # first entry when op_C2 was given as a vector
        C2 = isnan(c2_user) ? nothing : _op_C2_to_control2(actor, c2_user)
        @info "Solving one case: C1=$(C1)  C2=$(something(c2_user, "default"))"
        solve_one_case(par, actor.ode_params, application; C1, C2)
    
    elseif task == :solve_system
        @info "Solving the ODEs for $(application) system"
        compute_br_max(actor)   # log br_max implied by current Drw + max(EF)
        
	# Solve the ODE system on the whole control grid, normalize, and
        # classify (prob=nothing; filled in-place by train_locking_nn below)
        actor.results = _solve_grid_and_classify(actor, application)

        # Checkpoint ODE results to disk so calc_prob can run in a future session
        save_ode_results(actor)

    elseif task == :calc_prob
        # Build the probability model — ODEs are NOT re-solved.
        # Use in-memory results if available; otherwise load from disk.
        if actor.results === nothing
            @info "No in-memory ODE results — loading from disk"
            load_ode_results!(actor)
        end

        # Both windowed engines centre their window on a grid point, so an even
        # size has no centre — reject it here rather than silently off-setting.
        if par.prob_method in (:conv, :kde)
            iseven(par.conv_window_C1) && error("conv_window_C1 = $(par.conv_window_C1) must be odd for prob_method=$(repr(par.prob_method))")
            iseven(par.conv_window_C2) && error("conv_window_C2 = $(par.conv_window_C2) must be odd for prob_method=$(repr(par.prob_method))")
        end

        # Build the probability model with the engine par.prob_method selects.
        # Both windows are available to both windowed engines.
        if par.prob_method == :nn
            @info "Training NN classifier to calculate probability of locking as a function of Controls"
            train_locking_nn(actor)

        elseif par.prob_method == :conv
            @info "Computing convolution probability (window = $(par.conv_window_C1)×$(par.conv_window_C2))"
            conv_locking_probability(actor)
        elseif par.prob_method == :kde
            @info "Computing KDE probability (window = $(par.conv_window_C1)×$(par.conv_window_C2))"
            kde_locking_probability(actor)
        else
            error("Unknown prob_method: $(repr(par.prob_method)) — supported engines are :nn, :conv, :kde")
        end

        save_prob_model(actor)

    elseif task == :eval_prob
        # Load the saved probability model named by par.prob_method.
        if actor.results === nothing
            @info "No in-memory ODE results — loading from disk"
            load_ode_results!(actor)
        end
        # Prefer the model already in memory (e.g. just trained by :calc_prob);
        # only hit disk when nothing has been loaded yet this session.
        # No silent fallback across engines: loading a different model than the
        # one asked for produced plausible-looking numbers from the wrong source.
        if actor.results.prob === nothing
            actor.results.prob = load_prob_model(actor; method=par.prob_method)
        end
        prob_model = actor.results.prob

        # Determine times to evaluate — default to current dd.global_time
        times = isempty(par.op_times) ? [dd.global_time] : collect(par.op_times)

        # op_C2 is per-time: a scalar is broadcast across all op_times, a vector
        # must match op_times one-for-one. Values stay in user units here; the
        # mapping into the model's Control2 axis happens per point below (and is
        # re-derived from the same helper in _finalize).
        C2_user = if par.op_C2 isa AbstractVector
            length(par.op_C2) == length(times) ||
                error("op_C2 has $(length(par.op_C2)) entries but op_times has $(length(times)); " *
                      "pass a scalar to broadcast, or one op_C2 per op_time")
            collect(Float64, par.op_C2)
        else
            fill(Float64(par.op_C2), length(times))
        end

        rt = actor.ode_params.rat_surface

        # Accumulate locally, then publish as one unit so the three vectors
        # cannot be left inconsistent if the loop throws part-way through
        eval_C1 = Float64[]
        eval_C2 = Float64[]

        t_orig = dd.global_time
        for (t, c2_user) in zip(times, C2_user)
            dd.global_time  = t
            C1_op       = _rotation_at_rat_surface(dd, par, rt)
            C2_ctrl     = _op_C2_to_control2(actor, c2_user)
            P_locked    = prob_model(C1_op, C2_ctrl)
            push!(eval_C1, C1_op)
            push!(eval_C2, c2_user)
            @info "  t=$(round(t*1e3; digits=1)) ms  C1=$(round(C1_op; sigdigits=4))  C2=$(round(c2_user; sigdigits=4)) [user] → $(round(C2_ctrl; sigdigits=4)) [Control2]  P=$(round(P_locked; sigdigits=4))"
        end
        dd.global_time = t_orig
        actor.eval = (times=times, C1=eval_C1, C2=eval_C2)
        @info "eval_prob: $(length(times)) point(s)  model=$(typeof(prob_model))"

    elseif task == :transfer_learning
        # Solve the (typically focused/sparse, via ode_params.Control1_min/max
        # and Control2_min/max + a smaller grid_size) grid for the new dd,
        # normalize, and classify — same pipeline as :solve_system.
        actor.results = _solve_grid_and_classify(actor, application)

        # Fine-tune the saved base NN model on the new equilibrium's data,
        # freezing all but the last layer.
        @info "Loading base NN model for transfer learning"
        base_model = load_locking_nn(; control_type=par.control_type)

        X_new, y_new = ModeLocking.prepare_nn_data(actor.results.locking_labels,
                                        actor.ode_params.Control1, actor.ode_params.Control2)

        actor.results.prob = ModeLocking.transfer_learn_locking_nn(base_model, X_new, y_new; nn_params=actor.nn_params)
        save_locking_nn(actor; filename="nn_model_$(par.control_type)_TL.bson")

    elseif task == :Monte_Carlo
        error("Monte-Carlo not implemented yet")

    else
        error("Unknown task: $(task)")

    end

    # Figures are opt-in. Gated on results rather than on task: :single_case
    # displays its own traces and leaves results unset, and every task that does
    # produce results is worth plotting. save_locking_plots itself decides which
    # figures apply — the probability map only when a model exists, the
    # operating-point markers only when :eval_prob populated actor.eval.
    if par.save_plots && actor.results !== nothing
        save_locking_plots(actor; dir=isempty(par.plots_dir) ? pwd() : par.plots_dir)
    end

    return actor
end

"""
    _locking_mode_entry(dd, par) → toroidal_mode

Create (or fetch) the `mhd_linear` toroidal_mode entry for this m/n at the
current `dd.global_time` and fill its identity fields.  `stability_metric` is
deliberately not touched here — see `_finalize`.
"""
function _locking_mode_entry(dd::IMAS.dd, par)
    mhd_ts = resize!(dd.mhd_linear.time_slice; wipe=false)
    mode = resize!(mhd_ts.toroidal_mode,
        "perturbation_type.name" => "Island locking m=$(par.m_pol)/n=$(par.n_tor)",
        "n_tor" => par.n_tor)
    mode.perturbation_type.description = "Tearing-mode locking hazard (ActorLocking)"
    mode.m_pol_dominant = Float64(par.m_pol)
    return mode
end


"""
    _lower_fold_C2(bb, c1, c2, C1q; nl) → Union{Nothing, Tuple{Float64,Bool}}

Locate the lower (small-`Control2`) edge of the bistable band at rotation `C1q`,
returning `(C2_fold, bracketed)` in dimensionless `Control2` units.

`bb` is `results.bifurcation_bounds`, laid out with **rows = Control1 (rotation)
and columns = Control2** — a consequence of the column-major flattening in
`set_control_parameters!` (`Control1` varies fastest) combined with the square `reshape`
in `ModeLocking.calculate_bifurcation_bounds`.  `bb < 0` marks bistability in
both branches: the linear branch stores the cubic discriminant (D < 0 → three
real roots), the NL branch stores -1.0 where ≥ 2 positive real roots exist.

`bracketed=true` flags the NL branch, where `bb ∈ {-1,+1}` carries no within-cell
information and the fold is only resolved to one grid column.

Returns `nothing` — meaning *unresolved*, never *safe* — when the boundary was
not computed, when `C1q` falls outside the scanned rotation range, or when the
scan floor is already bistable (the fold lies at or below `Control2_min`, so any
value reported would be an artifact of the scan bounds).
"""
function _lower_fold_C2(bb::Union{AbstractMatrix,Nothing},
                        c1::AbstractVector, c2::AbstractVector, C1q::Real; nl::Bool)
    bb === nothing && return nothing
    @assert size(bb) == (length(c1), length(c2)) "bifurcation_bounds layout changed: expected (Control1, Control2) = $((length(c1), length(c2))), got $(size(bb))"

    (C1q < first(c1) || C1q > last(c1)) && return nothing
    row = @view bb[argmin(abs.(c1 .- C1q)), :]   # nearest scanned rotation

    j = findfirst(<(0), row)                     # first bistable column
    (j === nothing || j == 1) && return nothing

    nl && return (Float64(c2[j]), true)          # ±1 only: no sub-grid refinement

    v1, v2 = row[j-1], row[j]                    # interpolate the D = 0 crossing
    return (Float64(c2[j-1] + v1 / (v1 - v2) * (c2[j] - c2[j-1])), false)
end


function _finalize(actor::ActorLocking)
    dd  = actor.dd
    par = actor.par
    results    = actor.results
    ode_params = actor.ode_params

    results === nothing && return actor

    # Fraction of the scanned grid that classified as locked. This moves with
    # Control1/Control2_min/max, so it characterises the scan box rather than the
    # plasma — logged, never written to dd.
    frac_locked = count(==(2), results.locking_labels) / length(results.locking_labels)
    @info "Locked fraction of scanned grid: $(round(frac_locked; sigdigits=3))"

    # An operating point needs a probability model and time-aligned coordinates
    # Lengths are equal by construction now that eval is built as one unit
    ev = actor.eval
    have_op = results.prob !== nothing && ev !== nothing && !isempty(ev.C1)

    if !have_op
        # Record the mode's identity, but leave stability_metric unfilled rather
        # than substituting a grid statistic under the same name.
        _locking_mode_entry(dd, par)
        return actor
    end

    # ev.C2 holds user units (Gauss for :EF and :LinStab); the model axis is
    # dimensionless Control2. Reuse the mapping _step applied so the two can't diverge.
    c2_dim(c2) = _op_C2_to_control2(actor, c2)

    c1_axis = unique(ode_params.Control1)
    c2_axis = unique(ode_params.Control2)

    # Ascending order: resize! at global_time can only append forward in time
    order = sortperm(ev.times)

    # Bistability is written all-or-nothing across op times: set_time_array fills
    # gaps by repeating the previous value (or NaN-padding a late-created array),
    # so a ragged `fraction` would read back as a real limit at times where the
    # fold was never resolved.
    folds = [_lower_fold_C2(results.bifurcation_bounds, c1_axis, c2_axis, ev.C1[k];
                            nl=par.NL_saturation) for k in order]
    write_bistability = all(!isnothing, folds)
    if !write_bistability && any(!isnothing, folds)
        @info "Bistability limit skipped: fold unresolved at $(count(isnothing, folds))/$(length(folds)) operating times"
    end

    # The `C2_op / C2_fold > 1` convention assumes a positive Control2 axis. For
    # :LinStab the axis is Δt^RW < 0, so the ratio
    # of two negatives loses the "further past the fold" ordering — skip rather
    # than publish a limit whose direction is undefined.
    if write_bistability
        c2_ops = [c2_dim(ev.C2[k]) for k in order]
        if any(<=(0), c2_ops) || any(f -> f[1] <= 0, folds)
            @info "Bistability limit skipped: Control2 axis is not strictly positive (control_type=$(par.control_type)), so the fraction convention does not apply"
            write_bistability = false
        end
    end

    prob_model = resize!(dd.limits.model,
        "identifier.name" => "TM locking m=$(par.m_pol)/n=$(par.n_tor) probability")
    prob_model.identifier.description =
        "Locking probability at operating point < $(par.P_locked_threshold) (ActorLocking)"

    if write_bistability
        bist_model = resize!(dd.limits.model,
            "identifier.name" => "Island locking m=$(par.m_pol)/n=$(par.n_tor) bistability")
        bist_model.identifier.description =
            "Operating point below bistability onset (ActorLocking)" *
            (any(f -> f[2], folds) ? " [grid-bracketed: NL saturation]" : "")
    end

    t_orig = dd.global_time
    for (i, k) in enumerate(order)
        dd.global_time = ev.times[k]
        c2 = c2_dim(ev.C2[k])

        # stability_metric carries exactly one quantity: P(locked) at the op point
        P_op = results.prob(ev.C1[k], c2)
        _locking_mode_entry(dd, par).stability_metric = P_op
        @ddtime(prob_model.fraction = P_op / par.P_locked_threshold)   # > 1 → locked

        if write_bistability
            @ddtime(bist_model.fraction = c2 / folds[i][1])            # > 1 → bistable
        end
    end
    dd.global_time = t_orig

    return actor
end


# ─────────────────────────────────────────────────────────────────────────────
#  Setup: build ODE parameters from dd (called once per actor run)
# ─────────────────────────────────────────────────────────────────────────────

function set_up_ode_params!(dd::IMAS.dd, par, ode_params::ODEparams)
    """
    Initialize the ODE parameters for the locking actor.
    
    Args:
        dd: IMAS data structure
        par: Parameters for the simulation
    Returns:
        ode_params: Initialized ODE parameters
    """
    
    # find the normalized radius of the q=2 surface
    q_surf = par.m_pol / par.n_tor
    q_prof = dd.equilibrium.time_slice[].profiles_1d.q
    rho = dd.equilibrium.time_slice[].profiles_1d.rho_tor_norm
    ode_params.rat_surface = find_rat_surface(q_prof, rho, q_surf)

    # calculate the stability indices and mutual inductances
    ode_params = calculate_stability_index!(dd, par, ode_params)

    # Set physical parameters in dimensionless form
    ode_params = set_phys_params!(dd, par, ode_params)
    
    # Overwrite params to reproduce PoP2024
    if par.overwrite_params
        @info "Overwriting ODE parameters to reproduce PoP2024 results"
        ode_params.mu = 0.1
        ode_params.Inertia = 1
        ode_params.rat_surface = 0.67
        par.EF_phase = -90.0
    end

    
    # Prepare control parameters based on the control type
    ode_params = set_control_parameters!(dd, par, ode_params)

    return ode_params
end




function find_rat_surface(q_prof::Vector{Float64}, rho::Vector{Float64}, rat_surface::Float64)
    q_interp = IMAS.interp1d(rho, q_prof)
    f = x -> abs(q_interp(x)) - rat_surface
    roots = Roots.find_zeros(f, rho[begin], rho[end])
    isempty(roots) && error("No q=$rat_surface surface found in rho ∈ [$(rho[begin]), $(rho[end])]")
    rho_rat = maximum(roots)
    @info "Found q=$rat_surface surface at: $rho_rat ($(length(roots)) crossing(s))"
    return rho_rat
end



function calculate_stability_index!(dd::IMAS.dd, par, ode_params::ODEparams)
    rt = ode_params.rat_surface  
    rw = ode_params.res_wall
    rc = ode_params.control_surf
    m0 = par.m_pol
    
    rat21 = (rw / rt)^m0
    rat12 = rat21^(-1)
    rat32 = (rc / rw)^m0
    rat23 = rat32^(-1)

    ode_params.l12 = (2 * m0 / rw) / (rat21 - rat12)
    ode_params.l21 = (2 * m0 / rt) / (rat21 - rat12)
    ode_params.l32 = (2 * m0 / rw) / (rat32 - rat23)

    DeltaWl = (m0 / rw) * (rat21 + rat12) / (rat21 - rat12)
    DeltaWr = -(m0 / rw) * (rat32 + rat23) / (rat32 - rat23)

    ode_params.DeltaW = DeltaWr - DeltaWl
    ode_params.stability_index = par.RPRW_stability_index + ode_params.l21 * ode_params.l12 / ode_params.DeltaW

    # For :LinStab the swept Control2 supplies Δt directly (ModeLocking.resolve_control
    # takes Deltat = C2), so stability_index — and therefore RPRW_stability_index — has
    # no effect on the solved system. It survives only in the br_max diagnostic log and
    # in the _Drw output filename tag.
    if par.control_type == :LinStab
        @warn @sprintf(
            "control_type=:LinStab — par.RPRW_stability_index (%.4g) is unused: Δ_RW is the swept quantity, set with ode_params.Control2_min/max (in Δ_RW)",
            par.RPRW_stability_index)
    end

    # NOTE: the :LinStab sweep bounds are NOT clamped here. Control2_min/max are
    # user input and stay exactly as supplied; set_control_parameters! derives the
    # effective range into locals when it builds the grid. Read the range that was
    # actually used off `extrema(ode_params.Control2)`.

    return ode_params
end


function set_phys_params!(dd::IMAS.dd, par, ode_params::ODEparams)
    """
    Set the physical parameters in dimensionless form
    
    Args:
        dd: IMAS data structure
        par: Parameters for the simulation
        ode_params: ODE parameters to be set
    Returns:
        ode_params: Updated ODE parameters with physical constants set
    """

    # Define some constants
    #      Also NEED Zeff
    cp1d = dd.core_profiles.profiles_1d[]
    eqp1d = dd.equilibrium.time_slice[].profiles_1d
    mu0_val = IMAS.mks.μ_0
    
    # Set the scales for non-dimensionalization
    psi0 = par.b0 * par.r0
    U0 = psi0^2 * par.r0 / mu0_val
    
    mass_ion = cp1d.ion[1].element[1].a * IMAS.mks.m_p  # kg, mass of the main ion species

    rt = ode_params.rat_surface  # dimensionless, scaled by r0
    rho_cp = cp1d.grid.rho_tor_norm
    rho_eq = eqp1d.rho_tor_norm

    # Mass density at the rational surface: prefer ion density_thermal,
    # fall back to electron density / Z_i (quasi-neutrality) if ion density
    # is missing or zero.
    ion_dens = cp1d.ion[1].density_thermal
    if !ismissing(ion_dens) && !isempty(ion_dens) && any(!iszero, ion_dens)
        mass_dens_atq2 = mass_ion * IMAS.interp1d(rho_cp, ion_dens)(rt)
    else
        Z_i = cp1d.ion[1].element[1].z_n
        ne = cp1d.electrons.density_thermal
        ni_from_ne = ne ./ Z_i
        mass_dens_atq2 = mass_ion * IMAS.interp1d(rho_cp, ni_from_ne)(rt)
        @warn "ion density_thermal unavailable — using n_e/Z_i (quasi-neutrality) for mass density"
    end
 
    # flux-surface quantities at the rational surface
    R0           = IMAS.interp1d(rho_eq, eqp1d.gm8)(rt)     # <R>  [m]
    R2_avg       = IMAS.interp1d(rho_eq, eqp1d.gm10)(rt)    # <R²> [m²]
    area_rt = IMAS.interp1d(rho_eq, eqp1d.surface)(rt)  # S = 2π∫R dl [m²]
    gm2_rt       = IMAS.interp1d(rho_eq, eqp1d.gm2)(rt)     # <|∇ρ_tor|²/R²> [m⁻²]

    # rotation: core (ρ=0.1) and at the rational surface
    rot_profile = cp1d.ion[1].rotation_frequency_tor
    rot_core  = IMAS.interp1d(rho_cp, rot_profile)(0.1)
    rot_at_rs = IMAS.interp1d(rho_cp, rot_profile)(rt)
    
    # NBI torque at the rational surface:
    # Prefer the FUSE NBI actor output (dd.core_sources) if it has been populated upstream;
    # fall back to the user-specified par.source_torque scalar if not.
    torque_at_rat_surf = try
        src1d         = IMAS.total_sources(dd.core_sources, cp1d;
                            time0=dd.global_time, fields=[:torque_tor_inside])
        rho_src       = src1d.grid.rho_tor_norm
        torque_inside = src1d.torque_tor_inside
        if isempty(rho_src) || all(iszero, torque_inside)
            error("empty")
        end
        torque_val = IMAS.interp1d(rho_src, torque_inside)(rt)
        @info "NBI torque at rational surface from dd.core_sources: $(round(torque_val; sigdigits=4)) N·m"
        torque_val
    catch
        @warn "dd.core_sources torque not available — falling back to par.source_torque = $(par.source_torque) N·m"
        par.source_torque
    end

    # Calculate the drag coefficient in SI — must be positive (viscous drag, not drive)
    muSI = torque_at_rat_surf / rot_at_rs
    if muSI < 0
        @warn "Negative viscous drag μ_SI = $(round(muSI; sigdigits=4)) (torque and rotation have opposite signs); using |μ|"
        muSI = abs(muSI)
    end

    # Moment of inertia of toroidal shell: I = ρ · δ · S · <R²>
    # S = flux surface area (m²), δ = physical layer width (m)
    delta_phys = ode_params.layer_width * par.r0
    inertia = mass_dens_atq2 * delta_phys * area_rt * R2_avg

    # set nonlinear saturation
    if par.NL_saturation == false
        ode_params.saturation_param = 0.
    end

    # Set the dimensionless quantities
    ode_params.mu          = muSI / (U0 * par.t0)
    ode_params.Inertia     = inertia / (U0 * par.t0^2)
    ode_params.tor_fac      = area_rt * gm2_rt

    tau_mu = inertia / muSI
    rot_core_kHz = rot_core / (2π * 1e3)
    rot_rs_kHz   = rot_at_rs / (2π * 1e3)

    @info """Physical parameters set:
      torque_at_rat_surf = $(round(torque_at_rat_surf; sigdigits=4)) N·m
      rot_core(ρ≈0.1)   = $(round(rot_core_kHz; sigdigits=4)) kHz
      rot(q=2, ρ=$(round(rt; digits=3))) = $(round(rot_rs_kHz; sigdigits=4)) kHz
      μ_SI               = $(round(muSI; sigdigits=4)) N·m·s
      inertia            = $(round(inertia; sigdigits=4)) kg·m²
      tau_μ (viscous time) = $(round(tau_mu; sigdigits=4)) s
      mass_dens(q=2)     = $(round(mass_dens_atq2; sigdigits=4)) kg/m³
      ⟨R⟩(q=2)   = gm8  = $(round(R0; sigdigits=4)) m
      ⟨R²⟩(q=2)  = gm10 = $(round(R2_avg; sigdigits=4)) m²
      S(q=2)             = $(round(area_rt; sigdigits=4)) m²
      gm2(q=2)           = $(round(gm2_rt; sigdigits=4)) m⁻²
      S·gm2              = $(round(area_rt * gm2_rt; sigdigits=4))
      δ_phys             = $(round(delta_phys; sigdigits=4)) m
      ψ₀ = b0·r0         = $(round(psi0; sigdigits=4)) T·m
      U0 = ψ₀²·r0/μ₀    = $(round(U0; sigdigits=4)) N·m
      μ (dimensionless)  = $(round(ode_params.mu; sigdigits=4))
      I (dimensionless)  = $(round(ode_params.Inertia; sigdigits=4))
      Δ'                 = $(round(ode_params.stability_index; sigdigits=4))
      ΔW                 = $(round(ode_params.DeltaW; sigdigits=4))
      l12                = $(round(ode_params.l12; sigdigits=4))
      l21                = $(round(ode_params.l21; sigdigits=4))
      l32                = $(round(ode_params.l32; sigdigits=4))
      error_field        = $(par.error_field) Gauss  →  $(ode_params.error_field) (dimensionless EF flux function)
      EF_phase           = $(par.EF_phase)°"""

    return ode_params
end


function set_control_parameters!(dd::IMAS.dd, par, ode_params::ODEparams) 
    """
    Prepare the control parameters based on the control type.
    
    Args:
        control_type: Symbol indicating the type of control
        ode_params: ODE parameters to be modified
        par: Parameters for the simulation
    Returns:
        ode_params: Updated ODE parameters with control settings
    """

    control_type = par.control_type
    N = par.grid_size
    M = par.grid_size
    b0 = par.b0
    m_pol = par.m_pol

    l21 = ode_params.l21
    l12 = ode_params.l12
    DeltaW = ode_params.DeltaW
    rt = ode_params.rat_surface
    rc    = ode_params.control_surf
    c2min = ode_params.Control2_min
    c2max = ode_params.Control2_max

    # Rotation at the rational surface in kHz (from dd)
    rho_rot = dd.core_sources.source[1].profiles_1d[1].grid.rho_tor_norm
    rot_prof = dd.core_profiles.profiles_1d[1].rotation_frequency_tor_sonic
    rot_at_rs_rads = IMAS.interp1d(rho_rot, rot_prof)(rt)   # rad/s
    rot_at_rs_kHz  = rot_at_rs_rads / (2π * 1e3)            # kHz

    # Control1_min/max are in kHz and stay as the user supplied them; the grid is
    # built from a local raised, if needed, to cover the actual rotation
    c1min = ode_params.Control1_min
    c1max = max(ode_params.Control1_max, rot_at_rs_kHz)
    if c1max != ode_params.Control1_max
        @info @sprintf("Control1_max raised from %.4g to %.4g kHz to cover the rotation at the rational surface",
                       ode_params.Control1_max, c1max)
    end
    @info("Control1 (f₀) range: [$(c1min), $(round(c1max; sigdigits=4))] kHz")

    # Build grid in kHz, convert to dimensionless frequency units: f_dim = f_kHz * 1e3 * t0
    # (the 2π relating f to ω lives explicitly in the RHS)
    C1_kHz_vals = range(c1min, c1max, length=N) |> collect
    C1_dim_vals = C1_kHz_vals .* (1e3 * par.t0)
    ode_params.Control1 = vec(repeat(C1_dim_vals, 1, M))


    # Initialize the other control parameter based on the control type
    if control_type == :EF
        # For :EF, Control2_min/max are interpreted directly as the EF amplitude
        # in Gauss. Convert to the flux-equivalent perturbation psi_eps (still
        # referred to as "Eps" in the code), which carries units of b0*r0 (same
        # as psi0):
        #   EF_Tesla = EF_Gauss * 1e-4
        #   psi_eps  = -i * r_c * EF_Tesla / m_pol   (magnitude = r_c*EF_Tesla/m_pol)
        # The -i indicates psi_eps is -90deg out of phase with the true EF — not
        # yet propagated as an actual phase offset in the RHS (TODO, left as-is).
        EF_Gauss_vals = range(c2min, c2max, length=M) |> collect
        EF_Tesla_vals = EF_Gauss_vals .* 1.e-4 / b0   # Gauss -> Tesla -> dimensionless
        Control2_vals = rc .* EF_Tesla_vals ./ m_pol  # psi_eps, units of b0*r0
        @info("Maximum error field is $(c2max) Gauss")
    elseif control_type == :LinStab
        # Control2_min/max are Δ_RW here, not Δt. Δ_RW is the physically meaningful
        # quantity (Δt = Δ_RW + Δt_crit with Δt_crit purely geometric), and the
        # stability requirement reads directly as Δ_RW < 0 instead of being encoded
        # as an offset from Δt_crit. The grid itself holds Δt, which is what
        # ModeLocking.resolve_control expects from Control2.
        Δt_crit = l21 * l12 / DeltaW
        c2max < 0.0 || error("control_type=:LinStab — Control2_max is Δ_RW and must be < 0 for the RP-RW system to be weakly stable; got $(c2max)")
        c2min < c2max || error("control_type=:LinStab — Control2_min ($(c2min)) must be below Control2_max ($(c2max)) in Δ_RW")
        Control2_vals = collect(range(c2min, c2max, length=M)) .+ Δt_crit
        @info @sprintf("RP-RW stability sweep: Δ_RW ∈ [%.4g, %.4g]  →  Δt ∈ [%.4g, %.4g]  (Δt_crit=%.4g)",
                       c2min, c2max, first(Control2_vals), last(Control2_vals), Δt_crit)
        @info("Fixed EF in ODEs: $(par.error_field) Gauss → $(ode_params.error_field) (dimensionless)")

    else
        Control2_vals = range(c2min, c2max, length=M) |> collect
        @info("Fixed EF in ODEs: $(par.error_field) Gauss → $(ode_params.error_field) (dimensionless)")
    end

    # NOTE: error_field is NOT converted here. It is normalized once, in the
    # ActorLocking constructor (see _normalize_error_field!), because this
    # function runs on every `step` and an in-place Gauss→dimensionless rescale
    # here compounded silently whenever an ODEparams was reused across actors.

    Control2 = vec(repeat(Control2_vals', N, 1))
    ode_params.Control2 = Control2

    return ode_params
end


# ─────────────────────────────────────────────────────────────────────────────
#  Thin wrappers around ModeLocking — keep backward-compatible call sites
# ─────────────────────────────────────────────────────────────────────────────

"""
    solve_one_case(par, ode_params::ODEparams, application::String)

Solve and plot a single trajectory (task = :single_case). Delegates the ODE
solve to `ModeLocking.simulate_one_case` and the plotting to
`ModeLocking.plot_time_traces`.
"""
function solve_one_case(par, ode_params::ODEparams, application::String;
                        C1::Union{Float64,Nothing}=nothing, C2::Union{Float64,Nothing}=nothing)
    sol, norm_t = ModeLocking.simulate_one_case(ode_params, application, par.n_tor, deg2rad(par.EF_phase), par.control_type,
                                                  par.source_torque, par.t_final, par.time_steps; C1, C2)

    # Time-dependent figures: TM (ψ_tN), RWM (ψ_wN, RP-RW only), and Ω_tN vs time
    fig = ModeLocking.plot_time_traces(norm_t, sol.t; t0=par.t0)

    # Same traces, but in physical units (Gauss for magnetic amplitudes, kHz for rotation)
    fig_phys = ModeLocking.plot_time_traces(sol, [par.b0, par.t0, par.r0])

    # Stack both figures in a single combined plot so the second doesn't
    # overtake/replace the first on display
    fig_combined = plot(fig_phys, fig; layout=(2, 1), size=(800, 800))
    display(fig_combined)

    return fig
end


"""
    _solve_grid_and_classify(actor, application) → LockingResults

Solve the ODE system over the current `actor.ode_params.Control1/Control2`
grid, normalize the results, and classify them via k-means into
locked/unlocked. When NL saturation is off, also computes the analytic
bifurcation boundary.

Returns a `LockingResults` with `prob = nothing` — callers fill that in
(e.g. via `train_locking_nn` for `:solve_system`, or
`transfer_learn_locking_nn` for `:transfer_learning`).
"""
function _solve_grid_and_classify(actor::ActorLocking, application::String)
    par = actor.par
    return ModeLocking.solve_and_classify(actor.ode_params, application, par.n_tor, deg2rad(par.EF_phase), par.control_type,
                                           par.t_final, par.time_steps, par.NL_saturation, par.grid_size)
end


"""
    train_locking_nn(actor) → LockingNNModel

Train a binary NN classifier (C1, C2) → P(locked) on the k-means labels in
`actor.results`, using `actor.nn_params`.  Updates `actor.results.prob`
in-place and returns the model.

The stored model is callable: `actor.results.prob(C1, C2)` ∈ [0, 1].
To search for better hyperparameters first, call `tune_locking_nn(actor)`.
"""
function train_locking_nn(actor::ActorLocking)
    actor.results === nothing && error("No results — run the actor with task=:solve_system first")
    return ModeLocking.train_locking_nn(actor.results, actor.ode_params, actor.nn_params)
end


"""
    tune_locking_nn(actor; n_trials=20, n_folds=3, rng=GLOBAL_RNG) → NNparams

Random hyperparameter search using k-fold CV on the current actor results.
After finding the best configuration, retrains the full model with those
hyperparameters and updates `actor.results.prob` in-place.

Returns the best `NNparams` (useful for Task 2 transfer learning).
"""
function tune_locking_nn(actor::ActorLocking; kwargs...)
    actor.results === nothing && error("No results — run the actor with task=:solve_system first")
    return ModeLocking.tune_locking_nn(actor.results, actor.ode_params; kwargs...)
end


"""
    conv_locking_probability(actor; window_C1=par.conv_window_C1, window_C2=par.conv_window_C2)

Compute locking probability via 2D box-filter convolution and store in
`actor.results.prob`.  The result is callable as `actor.results.prob(C1, C2)`.

Call after `:solve_system` or `load_ode_results!(actor)`.
"""
function conv_locking_probability(actor::ActorLocking;
                                   window_C1::Int=actor.par.conv_window_C1,
                                   window_C2::Int=actor.par.conv_window_C2)
    actor.results === nothing && error("No results — run the actor with task=:solve_system first")
    actor.results.prob = ModeLocking.conv_locking_probability(
        actor.results, actor.ode_params, actor.par.grid_size;
        window_C1, window_C2)
    return actor.results.prob
end


"""
    kde_locking_probability(actor; window_C1=par.conv_window_C1, window_C2=par.conv_window_C2)

Compute locking probability via Gaussian KDE and store in
`actor.results.prob`.  The result is callable as `actor.results.prob(C1, C2)`.

Uses the same window parameters as `conv_locking_probability`;
the Gaussian σ along each axis is `window / 4`.
"""
function kde_locking_probability(actor::ActorLocking;
                                  window_C1::Int=actor.par.conv_window_C1,
                                  window_C2::Int=actor.par.conv_window_C2)
    actor.results === nothing && error("No results — run the actor with task=:solve_system first")
    actor.results.prob = ModeLocking.kde_locking_probability(
        actor.results, actor.ode_params, actor.par.grid_size;
        window_C1, window_C2)
    return actor.results.prob
end


"""
    save_ode_results(actor; filename="ode_results.bson", dir=ModeLocking.LOCKING_RESULTS_DIR) → path

Save the ODE grid results to disk so that `task=:calc_prob` can be run in a
future session without re-solving the ODEs.
"""
function save_ode_results(actor::ActorLocking; kwargs...)
    actor.results === nothing && error("No ODE results — run task=:solve_system first")
    return ModeLocking.save_ode_results(actor.results, actor.ode_params;
        control_type=actor.par.control_type, kwargs...)
end


"""
    load_ode_results!(actor; filename="ode_results.bson", dir=ModeLocking.LOCKING_RESULTS_DIR)

Load previously saved ODE results from disk into `actor.results` and
`actor.ode_params.Control1/Control2`.  Called by `task=:calc_prob` when
`actor.results` is nothing (fresh session).
"""
function load_ode_results!(actor::ActorLocking; kwargs...)
    results, Control1, Control2 = ModeLocking.load_ode_results(;
        control_type=actor.par.control_type, kwargs...)
    actor.ode_params === nothing && (actor.ode_params = ODEparams())
    actor.ode_params.Control1 = Control1
    actor.ode_params.Control2 = Control2
    actor.results = results
    return actor.results
end


"""
    save_prob_model(actor; kwargs...) → path

Save the current probability model (NN, convolution, or KDE) to disk.
Dispatches based on the type of `actor.results.prob`.
  - `LockingNNModel` → `nn_model_<control_type>.bson`
  - `ConvProbModel`  → `conv<w1>x<w2>_model_<control_type>.bson`
  - `KDEProbModel`   → `kde<w1>x<w2>_model_<control_type>.bson`
"""
function save_prob_model(actor::ActorLocking; kwargs...)
    actor.results === nothing      && error("No results — run task=:solve_system first")
    actor.results.prob === nothing && error("No trained model — run train_locking_nn, conv_, or kde_locking_probability first")
    ct  = actor.par.control_type
    nl  = actor.par.NL_saturation
    prob = actor.results.prob
    if prob isa LockingNNModel
        return ModeLocking.save_locking_nn(prob; control_type=ct, nl_sat=nl, kwargs...)
    elseif prob isa ModeLocking.ConvProbModel
        return ModeLocking.save_conv_prob(prob; control_type=ct, nl_sat=nl, window_C1=actor.par.conv_window_C1, window_C2=actor.par.conv_window_C2, kwargs...)
    elseif prob isa ModeLocking.KDEProbModel
        return ModeLocking.save_kde_prob(prob; control_type=ct, nl_sat=nl, window_C1=actor.par.conv_window_C1, window_C2=actor.par.conv_window_C2, kwargs...)
    else
        error("Unknown prob model type: $(typeof(prob))")
    end
end
save_locking_nn(actor::ActorLocking; kwargs...) = save_prob_model(actor; kwargs...)

"""
    load_prob_model(actor; method=:nn, kwargs...)

Load a saved probability model from disk.
  - `method=:nn`   → loads `nn_model_<control_type>.bson`
  - `method=:conv` → loads `conv<w1>x<w2>_model_<control_type>.bson`
  - `method=:kde`  → loads `kde<w1>x<w2>_model_<control_type>.bson`
"""
function load_prob_model(actor::ActorLocking; method::Symbol=:nn, kwargs...)
    ct  = actor.par.control_type
    nl  = actor.par.NL_saturation
    if method == :nn
        return ModeLocking.load_locking_nn(; control_type=ct, nl_sat=nl, kwargs...)
    elseif method == :conv
        return ModeLocking.load_conv_prob(; control_type=ct, nl_sat=nl, window_C1=actor.par.conv_window_C1, window_C2=actor.par.conv_window_C2, kwargs...)
    elseif method == :kde
        return ModeLocking.load_kde_prob(; control_type=ct, nl_sat=nl, window_C1=actor.par.conv_window_C1, window_C2=actor.par.conv_window_C2, kwargs...)
    else
        error("Unknown method: $method — use :nn, :conv, or :kde")
    end
end
load_locking_nn(; kwargs...) = ModeLocking.load_locking_nn(; kwargs...)


# ─────────────────────────────────────────────────────────────────────────────
#  Plotting wrappers
# ─────────────────────────────────────────────────────────────────────────────

"""
    plot_sols(actor) → (fig1, fig2, fig3)

Calls plot_scatter, plot_phase_diagrams, and (when a trained NN model is
available) plot_probability.  Returns all three handles; fig3 is nothing
when no model has been trained yet.
"""
plot_sols(actor::ActorLocking) = ModeLocking.plot_sols(actor.results, actor.ode_params, actor.par.grid_size, actor.par.control_type;
    b0=actor.par.b0, t0=actor.par.t0, m_pol=Float64(actor.par.m_pol), orientation=actor.par.plot_orientation,
    shot_label=_shot_label(actor.dd; overwrite=actor.par.overwrite_params))

function _shot_label(dd::IMAS.dd; overwrite::Bool=false, with_time::Bool=true)
    lbl = try
        pulse = dd.dataset_description.data_entry.pulse
        time  = dd.global_time
        t_ms  = round(time * 1e3; digits=1)
        if pulse > 0
            with_time ? "$(pulse), t = $(t_ms) ms" : "$(pulse)"
        else
            "default"
        end
    catch
        "default"
    end
    overwrite ? lbl * " (OW)" : lbl
end

"plot_scatter(actor) → Figure 1 — multi-panel scatter of state variables"
plot_scatter(actor::ActorLocking) = ModeLocking.plot_scatter(actor.results;
    orientation=actor.par.plot_orientation, shot_label=_shot_label(actor.dd; overwrite=actor.par.overwrite_params))

"plot_phase_diagrams(actor) → Figure 2 — pcolor of Ω_n and ψ_tn over control space"
plot_phase_diagrams(actor::ActorLocking) = ModeLocking.plot_phase_diagrams(actor.results, actor.ode_params, actor.par.grid_size, actor.par.control_type;
    b0=actor.par.b0, t0=actor.par.t0, m_pol=Float64(actor.par.m_pol), orientation=actor.par.plot_orientation,
    shot_label=_shot_label(actor.dd; overwrite=actor.par.overwrite_params))

"plot_probability(actor) → Figure 3 — locking probability with optional operating-point markers"
function plot_probability(actor::ActorLocking;
                          op_C1::Vector{Float64}=actor.eval === nothing ? Float64[] : actor.eval.C1,
                          op_C2::Vector{Float64}=actor.eval === nothing ? Float64[] : actor.eval.C2,
                          op_label::Vector{String}=actor.eval === nothing ? String[] :
                              ["t = $(round(t*1e3; digits=1)) ms" for t in actor.eval.times],
                          show_shot_label::Bool=true)
    drw      = round(actor.par.control_type == :LinStab ?
                     actor.ode_params.Control2_max : actor.par.RPRW_stability_index; sigdigits=3)
    base_lbl = _shot_label(actor.dd; overwrite=actor.par.overwrite_params, with_time=false)
    shot_lbl = show_shot_label ?
                   latexstring("\\mathrm{$(base_lbl)},\\;\\Delta^t_{rw} = $(drw)") : ""
    # ModeLocking places the operating-point markers on the DISPLAY axis. That
    # axis is Gauss for :EF (op_C2 is already Gauss) and 1/α for :NLsaturation
    # (ModeLocking inverts α itself), but Δt for :LinStab — where op_C2 is the br
    # amplitude in Gauss. Map it so the markers land on the axis being drawn.
    op_C2_disp = actor.par.control_type == :LinStab ?
                 [_op_C2_to_control2(actor, c2; verbose=false) for c2 in op_C2] : op_C2

    ModeLocking.plot_probability(actor.results, actor.ode_params, actor.par.control_type;
        b0=actor.par.b0, t0=actor.par.t0, m_pol=Float64(actor.par.m_pol),
        shot_label=shot_lbl, op_label=op_label,
        op_C1=op_C1, op_C2=op_C2_disp)
end

"""
    save_locking_plots(actor; dir, format=:png)

Save all available locking plots to `dir` with descriptive filenames encoding
the run metadata: control type, grid size, application, shot, time.

Example filenames:
  scatter_EF_100x100_RPRW_175060.2500_Drw0.5.png       (shot 175060, t=2500ms)
  NN_probability_all_op_points_EF_100x100_RPRW_175060.2500_Drw0.5.png  (all op points)
  NN_probability_op_EF_100x100_RPRW_175060.450_Drw0.5.png              (single snapshot, t=450ms)
"""
function save_locking_plots(actor::ActorLocking; dir::String, format::Symbol=:png)
    actor.results === nothing && error("No results to plot — run the actor first")
    isdir(dir) || mkpath(dir)

    par = actor.par
    dd  = actor.dd

    ow = par.overwrite_params ? "_OW" : ""
    shot_tag = try
        pulse = dd.dataset_description.data_entry.pulse
        time  = dd.global_time
        t_ms  = round(Int, time * 1e3)
        (pulse > 0) ? "_$(pulse).$(t_ms)$(ow)" : "_default$(ow)"
    catch
        "_default$(ow)"
    end

    drw     = round(abs(par.control_type == :LinStab ?
                        actor.ode_params.Control2_max : par.RPRW_stability_index); sigdigits=3)
    drw_tag = "_Drw$(drw)"
    tag = "$(par.control_type)_$(par.grid_size)x$(par.grid_size)_$(replace(par.application, "-" => ""))$(shot_tag)$(drw_tag)"
    ext = string(format)

    Plots.savefig(plot_scatter(actor),       joinpath(dir, "scatter_$(tag).$(ext)"))
    Plots.savefig(plot_phase_diagrams(actor), joinpath(dir, "phase_$(tag).$(ext)"))
    if actor.results.prob !== nothing
        prob = actor.results.prob
        prob_prefix = if prob isa LockingNNModel
            "NN_probability"
        elseif prob isa ModeLocking.ConvProbModel
            "Conv$(par.conv_window_C1)x$(par.conv_window_C2)_probability"
        elseif prob isa ModeLocking.KDEProbModel
            "KDE$(par.conv_window_C1)x$(par.conv_window_C2)_probability"
        else
            "probability"
        end
        # Clean plot (no operating-point overlay)
        Plots.savefig(plot_probability(actor; op_C1=Float64[], op_C2=Float64[]),
                      joinpath(dir, "$(prob_prefix)_$(tag).$(ext)"))
        # Duplicate with every operating point overlaid at once
        if actor.eval !== nothing && !isempty(actor.eval.C1)
            Plots.savefig(plot_probability(actor; show_shot_label=false),
                          joinpath(dir, "$(prob_prefix)_all_op_points_$(tag).$(ext)"))
            # One figure per snapshot with a single yellow star;
            # filename encodes the snapshot time in place of dd.global_time
            for i in eachindex(actor.eval.C1)
                t_snap = round(actor.eval.times[i] * 1e3; digits=1)
                snap_shot_tag = try
                    pulse = dd.dataset_description.data_entry.pulse
                    ow_str = par.overwrite_params ? "_OW" : ""
                    (pulse > 0) ? "_$(pulse).$(t_snap)$(ow_str)" : "_default$(ow_str)"
                catch
                    "_default$(par.overwrite_params ? "_OW" : "")"
                end
                snap_tag = "$(par.control_type)_$(par.grid_size)x$(par.grid_size)_$(replace(par.application, "-" => ""))$(snap_shot_tag)$(drw_tag)"
                lbl_i = ["t = $(t_snap) ms"]
                Plots.savefig(
                    plot_probability(actor;
                        op_C1=[actor.eval.C1[i]], op_C2=[actor.eval.C2[i]],
                        op_label=lbl_i, show_shot_label=false),
                    joinpath(dir, "$(prob_prefix)_op_$(snap_tag).$(ext)"))
            end
        end
    end
    @info "Saved locking plots to $dir"
    return nothing
end


# ─────────────────────────────────────────────────────────────────────────────
#  br diagnostic: max normal b-field at r_t ↔ RPRW_stability_index
# ─────────────────────────────────────────────────────────────────────────────

"""
    _eps_ref(actor) -> ε

Reference error-field flux perturbation (dimensionless psi_eps) for the br ↔ Δ_RW
relation.  Which quantity carries ε depends on `control_type`: for `:EF` the
swept `Control2` *is* the error field, so the peak of the scan is used; for
`:LinStab` and `:NLsaturation` the swept axis is Δt / α and ε is the fixed
`ode_params.error_field` (made dimensionless once by `_normalize_error_field!`).
"""
function _eps_ref(actor::ActorLocking)
    op  = actor.ode_params
    par = actor.par
    if par.control_type == :EF
        return isempty(op.Control2) ?
               op.Control2_max * 1e-4 / par.b0 * op.control_surf / Float64(par.m_pol) :
               maximum(op.Control2)
    else
        # Already dimensionless — _normalize_error_field! converted it once, in the
        # constructor. Re-normalizing here would apply (1e-4/b0)·r_c/m a second time.
        return op.error_field
    end
end

"""
    br_drw(actor; direction, drw, br_Gauss, eps, verbose) -> Float64

Locked-state relation between the peak normal b-field at the rational surface and
the effective RP-RW stability index Δ_RW (`RPRW_stability_index`):

    psit_max = l₂₁ · l₃₂ · ε / (Δw · Δ_RW)
    br_max   = (m / rt) · psit_max

Derived from the RP-RW steady state (`dydt[1]`/`dydt[4]` with Ω=0, α=0) and
independently reproduced by the Q=0 limit of `_nl_bifurcation_indicator`.  Note
Δ_RW = Δt − l₂₁·l₁₂/Δw, so the Δt dependence is carried entirely by Δ_RW.

`direction=:forward`  — given `drw` (Δ_RW), return br in Gauss.
`direction=:backward` — given `br_Gauss`, return Δ_RW.

`eps` defaults to `_eps_ref(actor)`.

Δ_RW comes out negative on its own: Δw < 0 while l₂₁, l₃₂, ε and ψt are positive.
No sign is imposed, so a positive result means the inputs are unphysical and is
reported rather than silently flipped.
"""
function br_drw(actor::ActorLocking; direction::Symbol,
                drw::Float64 = actor.par.control_type == :LinStab ?
                               actor.ode_params.Control2_max :
                               actor.par.RPRW_stability_index,
                br_Gauss::Float64=NaN,
                eps::Float64=NaN,
                verbose::Bool=true)
    op = actor.ode_params
    m  = Float64(actor.par.m_pol)
    b0 = actor.par.b0
    ε  = isnan(eps) ? _eps_ref(actor) : eps

    if direction === :forward
        psit_max = op.l21 * op.l32 * ε / (op.DeltaW * drw)
        br_Gauss = (m / op.rat_surface) * psit_max * b0 * 1e4   # b_r = (m/r)·ψ
        verbose && @info @sprintf(
            "br_drw[→]: Δw=%.4g  l₂₁=%.4g  l₃₂=%.4g  ε=%.4g  Δ_RW=%.4g  →  ψt_max=%.4g  br=%.4g G",
            op.DeltaW, op.l21, op.l32, ε, drw, psit_max, br_Gauss)
        return br_Gauss

    elseif direction === :backward
        isnan(br_Gauss) && error("br_drw(direction=:backward) requires br_Gauss")
        psit_max = (br_Gauss / (b0 * 1e4)) * op.rat_surface / m  # ψ = r·b_r/m
        out = op.l21 * op.l32 * ε / (op.DeltaW * psit_max)
        out >= 0 && error(@sprintf(
            "br_drw[←]: Δ_RW = %.4g ≥ 0 for br = %.4g G — the RP-RW system is not weakly stable at this amplitude (check ε=%.4g and Δw=%.4g)",
            out, br_Gauss, ε, op.DeltaW))
        verbose && @info @sprintf(
            "br_drw[←]: br=%.4g G  →  ψt_max=%.4g   (Δw=%.4g  l₂₁=%.4g  l₃₂=%.4g  ε=%.4g)  →  Δ_RW=%.4g",
            br_Gauss, psit_max, op.DeltaW, op.l21, op.l32, ε, out)
        return out

    else
        error("br_drw: direction must be :forward or :backward, got $(repr(direction))")
    end
end

"""
Forward br relation — returns `(br_norm, br_Gauss)`.  Thin wrapper over `br_drw`.

`drw` defaults to the Δ_RW the run actually uses. For `:LinStab` that is the top
of the sweep, `ode_params.Control2_max` (Δ_RW there), which gives the largest br
in the scan since br ∝ 1/|Δ_RW|. For `:EF`/`:NLsaturation` Δ_RW is fixed and comes
from `par.RPRW_stability_index`, which `:LinStab` ignores entirely.
"""
function compute_br_max(actor::ActorLocking;
                        drw::Float64 = actor.par.control_type == :LinStab ?
                                       actor.ode_params.Control2_max :
                                       actor.par.RPRW_stability_index,
                        eps_max::Float64=NaN)
    br_Gauss = br_drw(actor; direction=:forward, drw, eps=eps_max)
    return br_Gauss / (actor.par.b0 * 1e4), br_Gauss
end

"Backward br relation — returns Δ_RW.  Thin wrapper over `br_drw`."
compute_drw_from_br(actor::ActorLocking, br_Gauss::Float64; eps_max::Float64=NaN) =
    br_drw(actor; direction=:backward, br_Gauss, eps=eps_max)


"""
    _normalize_error_field!(ode_params, par) -> ODEparams

Derive `ode_params.error_field` from `par.error_field` (Gauss) into the
dimensionless psi_eps normalization the model uses throughout:

    errF = EF_Gauss · 1e-4 / b0 · r_c / m_pol

the same scaling `:EF` applies when building the Control2 grid, so a fixed error
field and a swept one land in identical units.  `ModeLocking.resolve_control`
reads `ode_params.error_field` expecting this normalization, which is why the
field cannot simply hold Gauss.

The value is *assigned* from `par`, never accumulated onto whatever the field
already held, so calling this repeatedly is a no-op — the compounding that made
error_field drift by ~1000× is impossible by construction.

`control_surf` is only ever read at run time, never reassigned, so its value is
already final here.
"""
function _normalize_error_field!(ode_params::ODEparams, par)
    norm = (1e-4 / par.b0) * ode_params.control_surf / Float64(par.m_pol)
    want = par.error_field * norm
    # ODEparams' own default is in Gauss; anything else means the caller set it
    # via ode_params, which par.error_field now overrides
    had = ode_params.error_field
    if had != ODEparams().error_field && !isapprox(had, want; rtol=1e-8)
        @warn @sprintf(
            "ode_params.error_field=%.4g is ignored — par.error_field=%.4g Gauss is authoritative (→ %.4g dimensionless)",
            had, par.error_field, want)
    end
    ode_params.error_field = want
    return ode_params
end

"""
    _rotation_at_rat_surface(dd, par, rt) -> Float64

Dimensionless rotation (f·t₀) at the rational surface `rt`, interpolated from
`core_profiles` at the current `dd.global_time`.  This is the Control1 value the
plasma actually sits at, as opposed to a user-supplied `par.op_C1`.
"""
function _rotation_at_rat_surface(dd::IMAS.dd, par, rt::Real)
    rho_cp      = dd.core_profiles.profiles_1d[].grid.rho_tor_norm
    rot_profile = dd.core_profiles.profiles_1d[].ion[1].rotation_frequency_tor
    rot_rs_rads = IMAS.interp1d(rho_cp, rot_profile)(rt)
    return (rot_rs_rads / (2π)) * par.t0
end

"""
    _op_C2_scalar(par) -> Float64

Single `op_C2` value for the scalar-only tasks (`:single_case`), taking the first
entry when `op_C2` was supplied as a vector.
"""
function _op_C2_scalar(par)
    par.op_C2 isa AbstractVector || return Float64(par.op_C2)
    return isempty(par.op_C2) ? NaN : Float64(first(par.op_C2))
end

"""
    _op_C2_to_control2(actor, c2_user) -> Float64

Map an operating-point `op_C2` from user units onto the dimensionless `Control2`
axis the probability model was trained on:

  - `:EF`           — Gauss (error field) → psi_eps = EF·1e-4/b0 · r_c/m
  - `:LinStab`      — Gauss (n=1 br amplitude at the rational surface) → Δt, by
                      inverting the locked-state br relation for Δ_RW and then
                      applying Δt = Δ_RW + l₂₁·l₁₂/Δw (the inverse of the
                      `stability_index` assignment in `calculate_stability_index!`)
  - `:NLsaturation` — α, already native

Warns when a `:LinStab` inversion lands outside the scanned Control2 range, since
the probability model is then extrapolating.
"""
function _op_C2_to_control2(actor::ActorLocking, c2_user::Real; verbose::Bool=true)
    par = actor.par
    op  = actor.ode_params

    if par.control_type == :EF
        return Float64(c2_user) * 1e-4 / par.b0 * op.control_surf / Float64(par.m_pol)

    elseif par.control_type == :LinStab
        # br (Gauss) → ψt → Δ_RW → Δt.  Δt_crit is the marginal value where
        # Δ_RW = 0 and the RP-RW system stops being weakly stable; Δt must stay
        # strictly below it. br_drw already errors on Δ_RW ≥ 0, so this is a
        # belt-and-braces check that also names the limit in the message.
        drw     = br_drw(actor; direction=:backward, br_Gauss=Float64(c2_user), verbose)
        Δt_crit = op.l21 * op.l12 / op.DeltaW
        Δt      = drw + Δt_crit
        Δt >= Δt_crit && error(@sprintf(
            "op_C2: br=%.4g G → Δt=%.4g ≥ Δt_crit=%.4g — RP-RW system would not be weakly stable",
            Float64(c2_user), Δt, Δt_crit))
        c2min, c2max = isempty(op.Control2) ? (op.Control2_min, op.Control2_max) : extrema(op.Control2)
        if !(c2min <= Δt <= c2max)
            @warn @sprintf(
                "op_C2: br=%.4g G → Δt=%.4g is outside the scanned Control2 range [%.4g, %.4g] (Δt_crit=%.4g) — probability model is extrapolating",
                Float64(c2_user), Δt, c2min, c2max, Δt_crit)
        end
        return Δt

    else
        return Float64(c2_user)
    end
end
