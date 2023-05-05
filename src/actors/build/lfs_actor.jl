#= ========== =#
#  LFS sizing  #
#= ========== =#
Base.@kwdef mutable struct FUSEparameters__ActorLFSsizing{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    do_plot::Entry{Bool} = Entry(Bool, "-", "Plot"; default=false)
    verbose::Entry{Bool} = Entry(Bool, "-", "Verbose"; default=false)
end

mutable struct ActorLFSsizing <: ReactorAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorLFSsizing
    function ActorLFSsizing(dd::IMAS.dd, par::FUSEparameters__ActorLFSsizing; kw...)
        logging_actor_init(ActorLFSsizing)
        par = par(kw...)
        return new(dd, par)
    end
end

"""
    ActorLFSsizing(dd::IMAS.dd, act::ParametersAllActors; kw...)

Actor that resizes the Low Field Side of the tokamak radial build
* Places TF outer leg at radius required to meet the dd.build.tf.ripple requirement
* Other low-field side layers are scaled proportionally

!!! note 
    Manipulates radial build information in `dd.build.layer`
"""
function ActorLFSsizing(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorLFSsizing(dd, act.ActorLFSsizing; kw...)
    if actor.par.do_plot
        plot(dd.build)
    end
    step(actor)
    finalize(actor)
    if actor.par.do_plot
        display(plot!(dd.build; cx=false))
    end
    return actor
end

function _step(actor::ActorLFSsizing)
    dd = actor.dd
    par = actor.par

    new_TF_radius = IMAS.R_tf_ripple(IMAS.get_build(dd.build, type=_plasma_).end_radius, dd.build.tf.ripple, dd.build.tf.coils_n)

    itf = IMAS.get_build(dd.build, type=_tf_, fs=_lfs_, return_index=true) - 1
    iplasma = IMAS.get_build(dd.build, type=_plasma_, return_index=true) + 1

    # resize layers proportionally
    # start from the vacuum gaps before resizing the material layers
    for vac in [true, false]
        old_TF_radius = IMAS.get_build(dd.build, type=_tf_, fs=_lfs_).start_radius
        delta = new_TF_radius - old_TF_radius
        if par.verbose
            println("TF outer leg radius changed by $delta [m]")
        end
        thicknesses = [dd.build.layer[k].thickness for k in iplasma:itf if !vac || lowercase(dd.build.layer[k].material) == "vacuum"]
        for k in iplasma:itf
            if !vac || lowercase(dd.build.layer[k].material) == "vacuum"
                dd.build.layer[k].thickness *= (1 + delta / sum(thicknesses))
                hfs_thickness = IMAS.get_build(dd.build, identifier=dd.build.layer[k].identifier, fs=_hfs_).thickness
                if dd.build.layer[k].thickness < hfs_thickness
                    dd.build.layer[k].thickness = hfs_thickness
                end
            end
        end
    end

    return actor
end