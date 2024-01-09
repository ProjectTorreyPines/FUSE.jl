#= ========== =#
#  LFS sizing  #
#= ========== =#
Base.@kwdef mutable struct FUSEparameters__ActorLFSsizing{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    do_plot::Entry{Bool} = Entry{Bool}("-", "Plot"; default=false)
    verbose::Entry{Bool} = Entry{Bool}("-", "Verbose"; default=false)
    maintenance::Switch{Symbol} = Switch{Symbol}([:vertical, :horizontal, :none], "-", "Scheme for installation/removal of in-vessel components"; default=:none)
    tor_modularity::Entry{Int} = Entry{Int}("-", "Number of toroidal modules of blanket normalized to number of TF coils (must be >= 1)"; default=2)
    pol_modularity::Entry{Int} = Entry{Int}("-", "Number of poloidal modules of each toroidal blanket sector (must be >= 1)"; default=1)
end

mutable struct ActorLFSsizing{D,P} <: ReactorAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorLFSsizing{P}
    function ActorLFSsizing(dd::IMAS.dd{D}, par::FUSEparameters__ActorLFSsizing{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorLFSsizing)
        par = par(kw...)
        return new{D,P}(dd, par)
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

    TF_lfs_layer = IMAS.get_build_layer(dd.build.layer, type=_tf_, fs=_lfs_)
    vessel_lfs_layer = IMAS.get_build_layer(dd.build.layer, type=_vessel_, fs=_lfs_)

    # calculate TF leg radius that gives required TF field ripple
    ripple_TF_radius = IMAS.R_tf_ripple(IMAS.get_build_layer(dd.build.layer, type=_plasma_).end_radius, dd.build.tf.ripple, dd.build.tf.coils_n)

    # calculate vacuum port geometry for maintenance
    if par.maintenance != :none
        rVP_hfs_ib, rVP_hfs_ob, rVP_lfs_ib, rVP_lfs_ob = IMAS.vertical_maintenance(dd.build; par.tor_modularity, par.pol_modularity)

        if par.tor_modularity == 1
            # if tor_modularity is 1, then lfs vessel end_radius must coincide with rVP_lfs_ob
            maintenance_TF_radius = rVP_lfs_ob + (TF_lfs_layer.start_radius - vessel_lfs_layer.end_radius)
        else            
            # if tor_modularity is 2, then TF coil end_radius can be flush with rVP_lfs_ob
            maintenance_TF_radius = rVP_lfs_ob - (TF_lfs_layer.thickness)
        end
    else
        maintenance_TF_radius = 0.0
    end

    # new TF radius is the larger of the two constraints
    new_TF_radius = maximum([ripple_TF_radius, maintenance_TF_radius])

    # resize layers proportionally
    # start from the vacuum gaps before resizing the material layers
    old_TF_radius = TF_lfs_layer.start_radius
    delta = new_TF_radius - old_TF_radius
    if par.verbose
        if par.maintenance != :none
            println("$(par.maintenance) maintenance (tor_modularity=$(par.tor_modularity), pol_modularity=$(par.pol_modularity)) requires lfs port wall at $rVP_lfs_ob [m]")
            println("ripple_TF_radius = $ripple_TF_radius [m], maintenance_TF_radius = $maintenance_TF_radius [m]")
        end
        println("old_TF_radius = $old_TF_radius [m], new_TF_radius = $new_TF_radius [m]")
        println("TF outer leg radius changed by $delta [m]")
    end

    itf = IMAS.get_build_index(dd.build.layer, type=_tf_, fs=_lfs_) - 1
    iplasma = IMAS.get_build_index(dd.build.layer, type=_plasma_) + 1
    for vac in (true) #, false) ## Removed the 'false' iteration due to double-counting, not sure why this was implimented (DBW)
        thicknesses = [dd.build.layer[k].thickness for k in iplasma:itf if !vac || lowercase(dd.build.layer[k].material) == "vacuum"]
        for k in iplasma:itf
            if !vac || lowercase(dd.build.layer[k].material) == "vacuum"
                dd.build.layer[k].thickness *= (1 + delta / sum(thicknesses))
                hfs_thickness = IMAS.get_build_layer(dd.build.layer, identifier=dd.build.layer[k].identifier, fs=_hfs_).thickness
                if dd.build.layer[k].thickness < hfs_thickness
                    dd.build.layer[k].thickness = hfs_thickness
                end
            end
        end
    end
    
    return actor
end