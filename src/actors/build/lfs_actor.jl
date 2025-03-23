#= ========== =#
#  LFS sizing  #
#= ========== =#
Base.@kwdef mutable struct FUSEparameters__ActorLFSsizing{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    maintenance::Switch{Symbol} = Switch{Symbol}([:vertical, :horizontal, :none], "-", "Scheme for installation/removal of in-vessel components"; default=:none)
    tor_modularity::Entry{Int} =
        Entry{Int}("-", "Number of toroidal modules of blanket normalized to number of TF coils"; default=2, check=x -> @assert x > 0 "must be: tor_modularity > 0")
    pol_modularity::Entry{Int} =
        Entry{Int}("-", "Number of poloidal modules of each toroidal blanket sector"; default=1, check=x -> @assert x in [1, 2] "must be: pol_modularity in [1,2]")
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
    verbose::Entry{Bool} = act_common_parameters(; verbose=false)
end

mutable struct ActorLFSsizing{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorLFSsizing{P}}
    function ActorLFSsizing(dd::IMAS.dd{D}, par::FUSEparameters__ActorLFSsizing{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorLFSsizing)
        par = OverrideParameters(par; kw...)
        return new{D,P}(dd, par)
    end
end

"""
    ActorLFSsizing(dd::IMAS.dd, act::ParametersAllActors; kw...)

Actor that resizes the Low Field Side of the tokamak radial build
It changes the location of the outer TF leg by takig into account requirement of ripple and maintenance ports.

!!! note

    Manipulates radial build information in `dd.build.layer`
"""
function ActorLFSsizing(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorLFSsizing(dd, act.ActorLFSsizing; kw...)
    if actor.par.do_plot
        plot(dd.build.layer)
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

    TF_lfs_layer = IMAS.get_build_layer(dd.build.layer; type=_tf_, fs=_lfs_)

    # calculate TF leg radius that gives required TF field ripple
    ripple_TF_radius = IMAS.R_tf_ripple(IMAS.get_build_layer(dd.build.layer; type=_plasma_).end_radius, dd.build.tf.ripple, dd.build.tf.coils_n)

    # calculate vacuum port geometry for maintenance
    if par.maintenance != :none
        rVP_hfs_ib, rVP_hfs_ob, rVP_lfs_ib, rVP_lfs_ob = IMAS.vertical_maintenance(dd.build; par.tor_modularity, par.pol_modularity)

        if par.tor_modularity == 1
            # if tor_modularity is 1, then lfs vessel end_radius must coincide with rVP_lfs_ob
            vessel_lfs_layer = IMAS.get_build_layers(dd.build.layer; type=_vessel_, fs=_lfs_)[1]
            @assert TF_lfs_layer.start_radius > vessel_lfs_layer.end_radius
            maintenance_TF_radius = rVP_lfs_ob + (TF_lfs_layer.start_radius - vessel_lfs_layer.end_radius)
        else
            # if tor_modularity is 2, then TF coil end_radius can be flush with rVP_lfs_ob
            maintenance_TF_radius = rVP_lfs_ob - TF_lfs_layer.thickness
        end
    else
        maintenance_TF_radius = 0.0
    end

    # new TF radius is the larger of the two constraints
    new_TF_radius = maximum([ripple_TF_radius, maintenance_TF_radius])
    old_TF_radius = TF_lfs_layer.start_radius
    delta = new_TF_radius - old_TF_radius

    if par.verbose
        if par.maintenance != :none
            println("$(par.maintenance) maintenance (tor_modularity=$(par.tor_modularity), pol_modularity=$(par.pol_modularity)) requires lfs port wall at $rVP_lfs_ob [m]")
            println("ripple_TF_radius = $ripple_TF_radius [m], maintenance_TF_radius = $maintenance_TF_radius [m]")
        end
        println("old_TF_radius = $old_TF_radius [m], new_TF_radius = $new_TF_radius [m]")
    end

    ivessel = IMAS.get_build_indexes(dd.build.layer; type=_vessel_, fs=_lfs_)[1] - 1
    iplasma = IMAS.get_build_index(dd.build.layer; type=_plasma_) + 1

    # resize first vacuum gap between VV and plasma
    for k in ivessel:-1:iplasma
        if lowercase(dd.build.layer[k].material) == "vacuum"
            old_thickness = dd.build.layer[k].thickness
            dd.build.layer[k].thickness += delta
            hfs_thickness = IMAS.get_build_layer(dd.build.layer; dd.build.layer[k].identifier, fs=_hfs_).thickness
            if dd.build.layer[k].thickness < hfs_thickness
                dd.build.layer[k].thickness = hfs_thickness
            end
            if par.verbose
                println("TF outer leg radius changed by $(dd.build.layer[k].thickness - old_thickness) [m]")
            end
            break
        end
    end

    return actor
end
