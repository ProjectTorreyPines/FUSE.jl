#= ============== =#
#  OH TF stresses  #
#= ============== =#
Base.@kwdef mutable struct FUSEparameters__ActorStresses{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    do_plot::Entry{Bool} = Entry(Bool, "-", "plot"; default=false)
    n_points::Entry{Int} = Entry(Int, "-", "Number of grid points"; default=5)
end

mutable struct ActorStresses <: ReactorAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorStresses
    function ActorStresses(dd::IMAS.dd, par::FUSEparameters__ActorStresses; kw...)
        logging_actor_init(ActorStresses)
        par = par(kw...)
        return new(dd, par)
    end
end

"""
    ActorStresses(dd::IMAS.dd, act::ParametersAllActors; kw...)

Estimates mechanical stresses on the center stack

!!! note 
    Stores data in `dd.solid_mechanics`
"""
function ActorStresses(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorStresses
    actor = ActorStresses(dd, par; kw...)
    step(actor; par.n_points)
    finalize(actor)
    if par.do_plot
        display(plot(actor.dd.solid_mechanics.center_stack.stress))
    end
    return actor
end

function _step(actor::ActorStresses; n_points::Integer=5)
    eq = actor.dd.equilibrium
    bd = actor.dd.build
    sm = actor.dd.solid_mechanics

    R_tf_in = IMAS.get_build(bd, type=_tf_, fs=_hfs_).start_radius
    R_tf_out = IMAS.get_build(bd, type=_tf_, fs=_hfs_).end_radius
    R0 = (R_tf_in + R_tf_out) / 2.0
    B0 = maximum(eq.vacuum_toroidal_field.b0)
    Bz_oh = bd.oh.max_b_field
    R_oh_in = IMAS.get_build(bd, type=_oh_).start_radius
    R_oh_out = IMAS.get_build(bd, type=_oh_).end_radius
    f_struct_tf = bd.tf.technology.fraction_stainless
    f_struct_oh = bd.oh.technology.fraction_stainless

    bucked = sm.center_stack.bucked == 1
    noslip = sm.center_stack.noslip == 1
    plug = sm.center_stack.plug == 1
    empty!(sm.center_stack)

    for oh_on in [true, false]
        solve_1D_solid_mechanics!(
            sm.center_stack,
            R0,
            B0,
            R_tf_in,
            R_tf_out,
            oh_on ? Bz_oh : 0.0,
            R_oh_in,
            R_oh_out;
            bucked=bucked,
            noslip=noslip,
            plug=plug,
            f_struct_tf=f_struct_tf,
            f_struct_oh=f_struct_oh,
            f_struct_pl=1.0,
            n_points=n_points,
            verbose=false
        )
    end

    return actor
end

@recipe function plot_ActorStresses(actor::ActorStresses)
    @series begin
        actor.dd.solid_mechanics.center_stack.stress
    end
end