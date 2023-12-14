#= ============= =#
#  ActorPFdesign  #
#= ============= =#
Base.@kwdef mutable struct FUSEparameters__ActorPFdesign{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    dr_sep::Entry{T} = Entry{T}("m", "Distance between primary and secondary separatrix evaluated at the midplane")
end

mutable struct ActorPFdesign{D,P} <: ReactorAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorPFdesign{P}
    actor_pf::ActorPFactive{D,P}
end

"""
    ActorPFdesign(dd::IMAS.dd, act::ParametersAllActors; kw...)

Optimize pf_active to achieve desired equilibrium

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
    eqt=dd.equilibrium.time_slice[]

    # reset pf coil rails
    n_coils = [rail.coils_number for rail in dd.build.pf_active.rail]
    init_pf_active!(dd.pf_active, dd.build, dd.equilibrium.time_slice[], n_coils)

    actor.actor_pf.λ_regularize = -1.0
    step(actor.actor_pf)

    # Get coils (as GS3_IMAS_pf_active__coil) organized by their function and initialize them
    fixed_coils, pinned_coils, optim_coils = fixed_pinned_optim_coils(actor.actor_pf, :rail)

    symmetric = false
    packed, bounds = pack_rail(dd.build, actor.actor_pf.λ_regularize, symmetric)

    function placement_cost(packed::Vector{Float64})
        # update dd.pf_active from packed vector
        actor.actor_pf.λ_regularize = unpack_rail!(packed, optim_coils, symmetric, dd.build)

        # find currents
        _step(actor.actor_pf)

        # make nearby coils more expensive
        const_spacing = 0.0
        if length(optim_coils) > 0
            for (k1, c1) in enumerate(optim_coils)
                for (k2, c2) in enumerate(optim_coils)
                    if k1 < k2
                        d = sqrt((c1.r - c2.r)^2 + (c1.z - c2.z)^2)
                        s = sqrt((c1.width + c2.width)^2 + (c1.height + c2.height)^2)
                        if IMAS.is_ohmic_coil(imas(c1)) && IMAS.is_ohmic_coil(imas(c2))
                        else
                            const_spacing = max(const_spacing, s - d)
                        end
                    end
                end
            end
        end

        @show actor.actor_pf.cost
        return norm([actor.actor_pf.cost, const_spacing])
        #return [actor.actor_pf.cost], [const_spacing], [0.0]
    end

    old_logging = actor_logging(dd, false)

    # options = Metaheuristics.Options(; seed=1, iterations=50)
    # algorithm = Metaheuristics.DE(; N=20, options)
    # res = Metaheuristics.optimize(x -> placement_cost(x), bounds, algorithm)
    # packed = Metaheuristics.minimizer(res)

    res = Optim.optimize(placement_cost, packed, Optim.NelderMead(), Optim.Options(g_tol=1E-6))
    packed = res.minimizer

    #if verbose
    println(res)
    #end
    actor.actor_pf.λ_regularize = unpack_rail!(packed, optim_coils, symmetric, dd.build)
    actor_logging(dd, old_logging)

    return actor
end

"""
    _finalize(actor::ActorPFdesign)

Update actor.eq_out 2D equilibrium PSI based on coils currents
"""
function _finalize(actor::ActorPFdesign{D,P}) where {D<:Real,P<:Real}
    finalize(actor.actor_pf)
    return actor
end
