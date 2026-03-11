#= ============= =#
#  ActorPFdesign  #
#= ============= =#
@actor_parameters_struct ActorPFdesign{T} begin
    symmetric::Entry{Bool} = Entry{Bool}("-", "Force PF coils location to be up-down symmetric"; default=true)
    update_equilibrium::Entry{Bool} = Entry{Bool}("-", "Overwrite target equilibrium with the one that the coils can actually make"; default=false)
    model::Switch{Symbol} = Switch{Symbol}([:none, :uniform, :optimal], "-", "Coil placement strategy"; default=:optimal)
    reset_rails::Entry{Bool} = Entry{Bool}("-", "Reset PF coils rails"; default=true)
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
    verbose::Entry{Bool} = act_common_parameters(; verbose=false)
end

mutable struct ActorPFdesign{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorPFdesign{P}}
    act::ParametersAllActors{P}
    actor_pf::ActorPFactive{D,P}
end

"""
    ActorPFdesign(dd::IMAS.dd, act::ParametersAllActors; kw...)

Optimizes poloidal field (PF) coil placement and sizing to achieve target equilibrium configurations.
Uses optimization algorithms to find coil positions that minimize current requirements while satisfying
magnetic flux and force balance constraints for plasma control.

# Design approaches (controlled by `model`)
- `:none`: Uses existing coil configuration without modification
- `:uniform`: Places coils uniformly on predefined rails
- `:optimal`: Optimizes both coil placement and currents using least-squares minimization

# Optimization process (for `:optimal` mode)
1. Initializes PF coil rails based on radial build geometry
2. Uses nested optimization: outer loop optimizes coil positions, inner loop finds currents
3. Minimizes flux/saddle constraint violations while penalizing large currents
4. Includes coil spacing penalties to prevent overlapping
5. Sizes coils based on current requirements and engineering margins

# Key constraints
- Flux control points (boundary and internal flux surfaces)
- Saddle control points (X-points and magnetic axis)
- Current density limits and coil spacing requirements
- Up-down symmetry (if `symmetric=true`)

# Key outputs
- Optimized PF coil positions and sizes (`dd.pf_active.coil[].element[].geometry`)
- Coil currents that achieve target equilibrium
- Updated equilibrium solution (if `update_equilibrium=true`)

!!! note

    Manupulates data in `dd.pf_active`
"""
function ActorPFdesign(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorPFdesign(dd, act.ActorPFdesign, act; kw...)
    finalize(step(actor))
    return actor
end

function ActorPFdesign(dd::IMAS.dd, par::FUSEparameters__ActorPFdesign, act::ParametersAllActors; kw...)
    logging_actor_init(ActorPFdesign)
    par = OverrideParameters(par; kw...)
    actor_pf = ActorPFactive(dd, act.ActorPFactive; par.update_equilibrium)
    return ActorPFdesign(dd, par, act, actor_pf)
end

"""
    _step(actor::ActorPFdesign)

Executes PF coil design optimization based on the selected model:
- For `:optimal`: Uses nested optimization with coil placement (outer) and current finding (inner)
- Includes cost penalties for coil spacing and current magnitudes
- Applies engineering constraints (current limits, spacing) and symmetry requirements
- Sizes coils based on current requirements after optimization completes
- Always runs ActorPFactive at the end to update equilibrium with final coil configuration
"""
function _step(actor::ActorPFdesign{T}) where {T<:Real}
    dd = actor.dd
    par = actor.par
    eqt = dd.equilibrium.time_slice[]

    if par.model == :none
        # leave untouched

    elseif par.model in [:uniform, :optimal]
        # reset pf coil rails
        if par.reset_rails
            init_pf_active!(dd.pf_active, dd.build, eqt)
        end

        # optimize coil placement
        if par.model == :optimal
            actor.actor_pf.λ_regularize = -1.0
            _step(actor.actor_pf)

            function placement_cost(packed::Vector{Float64}; prog::Any)
                # update dd.pf_active from packed vector
                optim_coils = actor.actor_pf.setup_cache.optim_coils
                actor.actor_pf.λ_regularize = unpack_rail!(packed, optim_coils, par.symmetric, dd.build)

                # find currents
                _step(actor.actor_pf)

                # make coils that are close to one another more expensive
                cost_spacing = 0.0
                if length(optim_coils) > 0
                    for (k1, c1) in enumerate(optim_coils)
                        for (k2, c2) in enumerate(optim_coils)
                            if k1 < k2
                                d = sqrt((c1.r - c2.r)^2 + (c1.z - c2.z)^2)
                                s = sqrt((c1.width + c2.width)^2 + (c1.height + c2.height)^2)
                                if !(IMAS.is_ohmic_coil(VacuumFields.imas(c1)) && IMAS.is_ohmic_coil(VacuumFields.imas(c2)))
                                    cost_spacing = max(cost_spacing, (s - d) / d)
                                end
                            end
                        end
                    end
                end

                coils = (coil for coil in vcat(actor.actor_pf.setup_cache.fixed_coils, actor.actor_pf.setup_cache.pinned_coils, actor.actor_pf.setup_cache.optim_coils))
                cost_currents = norm((coil.current_per_turn * coil.turns for coil in coils)) / eqt.global_quantities.ip

                cost = norm((actor.actor_pf.cost, 0.1 * cost_spacing))^2 * (1 .+ cost_currents)

                if prog !== nothing
                    ProgressMeter.next!(prog; showvalues=[("constraints", actor.actor_pf.cost), ("spacing", cost_spacing), ("currents", cost_currents)])
                end

                return cost
            end

            old_logging = actor_logging(dd, false)
            par.verbose && ProgressMeter.ijulia_behavior(:clear)
            prog = ProgressMeter.ProgressUnknown(;dt=0.1, desc="Calls:", enabled=par.verbose)
            try
                packed, bounds = pack_rail(dd.build, actor.actor_pf.λ_regularize, par.symmetric)
                res = Optim.optimize(x -> placement_cost(x; prog), packed, Optim.NelderMead())#, Optim.Options(; g_tol=1E-6))
                packed = res.minimizer
                actor.actor_pf.λ_regularize = unpack_rail!(packed, actor.actor_pf.setup_cache.optim_coils, par.symmetric, dd.build)
                if par.verbose
                    println(res)
                end
                # size the PF coils based on the currents they are carrying
                size_pf_active(actor.actor_pf.setup_cache.optim_coils, eqt; min_size=1.0, tolerance=dd.requirements.coil_j_margin, par.symmetric)
            finally
                actor_logging(dd, old_logging)
            end
        end
    end

    _step(actor.actor_pf) # must run actor_pf to update the equilibrium accordingly

    return actor
end

"""
    _finalize(actor::ActorPFdesign)

Update actor.eqt2d_out 2D equilibrium PSI based on coils currents
"""
function _finalize(actor::ActorPFdesign{D,P}) where {D<:Real,P<:Real}
    par = actor.par
    finalize(actor.actor_pf)
    if par.do_plot
        display(plot(actor.actor_pf; rails=true))
    end
    return actor
end
