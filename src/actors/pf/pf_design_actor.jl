#= ============= =#
#  ActorPFdesign  #
#= ============= =#
Base.@kwdef mutable struct FUSEparameters__ActorPFdesign{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    symmetric::Entry{Bool} = Entry{Bool}("-", "Force PF coils location to be up-down symmetric"; default=true)
    do_plot::Entry{Bool} = Entry{Bool}("-", "Plot"; default=false)
    verbose::Entry{Bool} = Entry{Bool}("-", "Verbose"; default=false)
end

mutable struct ActorPFdesign{D,P} <: ReactorAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorPFdesign{P}
    actor_pf::ActorPFactive{D,P}
end

"""
    ActorPFdesign(dd::IMAS.dd, act::ParametersAllActors; kw...)

Optimize PF coil locations to achieve desired equilibrium

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
    par = par(kw...)
    actor_pf = ActorPFactive(dd, act.ActorPFactive)
    return ActorPFdesign(dd, par, actor_pf)
end

"""
    _step(actor::ActorPFdesign)

Find currents that satisfy boundary and flux/saddle constraints in a least-square sense
"""
function _step(actor::ActorPFdesign{T}) where {T<:Real}
    dd = actor.dd
    par = actor.par

    # reset pf coil rails
    n_coils = [rail.coils_number for rail in dd.build.pf_active.rail]
    init_pf_active!(dd.pf_active, dd.build, dd.equilibrium.time_slice[], n_coils)

    actor.actor_pf.λ_regularize = -1.0
    step(actor.actor_pf)

    function placement_cost(packed::Vector{Float64}; prog::Any)
        # update dd.pf_active from packed vector
        optim_coils = actor.actor_pf.setup_cache.optim_coils
        actor.actor_pf.λ_regularize = unpack_rail!(packed, optim_coils, par.symmetric, dd.build)

        # find currents
        _step(actor.actor_pf)
        # size_pf_active(actor.actor_pf.setup_cache.optim_coils)

        # make coils that are close to one another more expensive
        cost_spacing = 0.0
        if length(optim_coils) > 0
            for (k1, c1) in enumerate(optim_coils)
                for (k2, c2) in enumerate(optim_coils)
                    if k1 < k2
                        d = sqrt((c1.r - c2.r)^2 + (c1.z - c2.z)^2)
                        s = sqrt((c1.width + c2.width)^2 + (c1.height + c2.height)^2)
                        if !(IMAS.is_ohmic_coil(imas(c1)) && IMAS.is_ohmic_coil(imas(c2)))
                            cost_spacing = max(cost_spacing, (s - d) / d)
                        end
                    end
                end
            end
        end

        cost = norm([actor.actor_pf.cost, cost_spacing])^2

        if prog !== nothing
            ProgressMeter.next!(prog; showvalues=[("constraints", actor.actor_pf.cost), ("spacing", cost_spacing)])
        end

        return cost
    end

    old_logging = actor_logging(dd, false)
    prog = ProgressMeter.ProgressUnknown(; desc="Calls:", enabled=par.verbose)
    try
        packed, bounds = pack_rail(dd.build, actor.actor_pf.λ_regularize, par.symmetric)
        res = Optim.optimize(x -> placement_cost(x; prog), packed, Optim.NelderMead())#, Optim.Options(; g_tol=1E-6))
        packed = res.minimizer
        actor.actor_pf.λ_regularize = unpack_rail!(packed, actor.actor_pf.setup_cache.optim_coils, par.symmetric, dd.build)
        if par.verbose
            println(res)
        end
    finally
        actor_logging(dd, old_logging)
    end

    size_pf_active(actor.actor_pf.setup_cache.optim_coils)

    return actor
end

"""
    _finalize(actor::ActorPFdesign)

Update actor.eq_out 2D equilibrium PSI based on coils currents
"""
function _finalize(actor::ActorPFdesign{D,P}) where {D<:Real,P<:Real}
    par = actor.par
    finalize(actor.actor_pf)
    if par.do_plot
        display(plot(actor.actor_pf; rails=true))
    end
    return actor
end
