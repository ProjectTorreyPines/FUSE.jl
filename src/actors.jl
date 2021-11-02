#= =================== =#
#  init IMAS structures #
#= =================== =#

function init(ids::IMAS.IDS, time::Real)
    error("Function init() not defined for ids of type $(typeof(ids))")
end

"""
    init(equilibrium::IMAS.equilibrium, time::Real=0.0; B0, R0, ϵ, δ, κ, beta_n, ip, x_point::Union{Vector, NTuple{2}, Bool}=false)

Initialize equilibrium IDS based on some basic Miller geometry parameters
"""
function init(equilibrium::IMAS.equilibrium, time::Real=0.0;
        B0, R0, ϵ, δ, κ, beta_n, ip,
        x_point::Union{Vector,NTuple{2},Bool}=false)
    time_index = get_time_index(equilibrium.time_slice, time)
    eqt = equilibrium.time_slice[time_index]
    eqt.boundary.minor_radius = ϵ * R0
    eqt.boundary.geometric_axis.r = R0
    eqt.boundary.elongation = κ
    eqt.boundary.triangularity = δ
    set_field_time_array(equilibrium.vacuum_toroidal_field, :b0, time_index, B0)
    equilibrium.vacuum_toroidal_field.r0 = R0
    eqt.global_quantities.ip = ip
    eqt.global_quantities.beta_normal = beta_n
    if x_point === true
        x_point = (R0 * (1 - 1.1 * δ * ϵ), -R0 * 1.1 * κ * ϵ)
    end
    if isa(x_point, Union{Vector,Tuple})
        resize!(eqt.boundary.x_point, 1)
        eqt.boundary.x_point[1].r = x_point[1]
        eqt.boundary.x_point[1].z = x_point[2]
    end
    return equilibrium
end

"""
    init(radial_build::IMAS.radial_build; Bmax_OH=nothing, Bmax_TF=nothing, layers...)

Initialize radial_build IDS based on center stack layers (thicknesses) and maximum fields

NOTE: index and material follows from standard naming of layers
*  0 ...gap... : vacuum
*  1 OH: ohmic coil
*  2 TF: toroidal field coil
*  3 ...shield...: neutron shield
*  4 ...blanket...: neutron blanket
*  5 ...wall....: 
* -1 ...vessel...: 
"""
function init(radial_build::IMAS.radial_build; Bmax_OH=nothing, Bmax_TF=nothing, layers...)
    if Bmax_OH !== nothing
        radial_build.oh_b_field_max = Bmax_OH
    end
    if Bmax_TF !== nothing
        radial_build.tf_b_field_max = Bmax_TF
    end
    # assign layers
    resize!(radial_build.center_stack, length(layers))
    for (klayer, (layer_name, layer_thickness)) in enumerate(layers)
        radial_build.center_stack[klayer].thickness = layer_thickness
        radial_build.center_stack[klayer].name = replace(String(layer_name), "_" => " ")
        if occursin("gap", lowercase(radial_build.center_stack[klayer].name))
            radial_build.center_stack[klayer].index = 0
            radial_build.center_stack[klayer].material = "vacuum"
        elseif uppercase(radial_build.center_stack[klayer].name) == "OH"
            radial_build.center_stack[klayer].index = 1
        elseif uppercase(radial_build.center_stack[klayer].name) == "TF"
            radial_build.center_stack[klayer].index = 2
        elseif occursin("shield", lowercase(radial_build.center_stack[klayer].name))
            radial_build.center_stack[klayer].index = 3
        elseif occursin("blanket", lowercase(radial_build.center_stack[klayer].name))
            radial_build.center_stack[klayer].index = 4
        elseif occursin("wall", lowercase(radial_build.center_stack[klayer].name))
            radial_build.center_stack[klayer].index = 5
        elseif occursin("vessel", lowercase(radial_build.center_stack[klayer].name))
            radial_build.center_stack[klayer].index = -1
        end
    end
    return radial_build
end

"""
    init(radial_build::IMAS.radial_build, eqt::IMAS.equilibrium__time_slice; is_nuclear_facility=true)

Simple initialization of radial_build IDS based on equilibrium time_slice
"""
function init(radial_build::IMAS.radial_build, eqt::IMAS.equilibrium__time_slice; is_nuclear_facility=true)
    Bmax_OH = nothing
    Bmax_TF = eqt.profiles_1d.f[end] / eqt.boundary.geometric_axis.r

    rmin = eqt.boundary.geometric_axis.r - eqt.boundary.minor_radius
    rmax = eqt.boundary.geometric_axis.r + eqt.boundary.minor_radius

    if is_nuclear_facility
        nlayers = 6
        gap = (rmax - rmin) / 10.0
        rmin -= gap
        rmax += gap
        dr = rmin / nlayers
        init(radial_build,
            Bmax_OH=Bmax_OH,
            Bmax_TF=Bmax_TF,
            gap_TF=dr * 2.0,
            OH=dr,
            TF=dr,
            inner_shield=dr / 2.0,
            inner_blanket=dr,
            inner_wall=dr / 2.0,
            vacuum_vessel=rmax - rmin,
            outer_wall=dr / 2.0,
            outer_blanket=dr,
            outer_shield=dr / 2.0)

    else
        nlayers = 4.5
        gap = (rmax - rmin) / 10.0
        rmin -= gap
        rmax += gap
        dr = rmin / nlayers
        init(radial_build,
            Bmax_OH=Bmax_OH,
            Bmax_TF=Bmax_TF,
            gap_TF=dr * 2.0,
            OH=dr,
            TF=dr,
            inner_wall=dr / 2.0,
            vacuum_vessel=rmax - rmin,
            outer_wall=dr / 2.0)
    end

    return radial_build
end

#= ============ =#
#  AbstractActor #
#= ============ =#
abstract type AbstractActor end

"""
    Take a step with a given actor
"""
function Base.step(actor::AbstractActor)
    error("Function step() not defined for actor of type $(typeof(actor))")
end

"""
    store output data in IDS
    NOTE: actors should take dd and output dd, operating in place
"""
function finalize(actor::AbstractActor)
    error("Function finalize() not defined for actor of type $(typeof(actor))")
end

#= =========== =#
#  Equilibrium  #
#= =========== =#

abstract type EquilibriumActor <: AbstractActor end

include("actors/solovev_equilibrium_actor.jl")

export SolovevEquilibriumActor

#= ===== =#
#  Coils  #
#= ===== =#

abstract type CoilsActor <: AbstractActor end

include("actors/coils_actor.jl")

export PFcoilsOptActor

#= ============ =#
#  Radial Build  #
#= ============ =#

abstract type RadialBuildActor <: AbstractActor end

include("actors/radial_build_actor.jl")

