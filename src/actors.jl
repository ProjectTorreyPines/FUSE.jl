import DataStructures

#= =================== =#
#  init IMAS structures #
#= =================== =#
"""
    dummy `init`
"""
function init(ids::IMAS.IDS)
    error("Function init() not defined for IDS of type $(typeof(ids))")
end

function init_from_gasc(dd::IMAS.dd, filename, case, no_small_gaps=true; eq_kw = Dict(), cp_kw = Dict(), rb_kw = Dict(), verbose=false)
    gasc = read_GASC(filename, case + 1) # +1 to account for GASC (PYTHON) and FUSE (JULIA) numbering

    eq = dd.equilibrium
    cp = dd.core_profiles
    bd = dd.build

    R0 = gasc["INPUTS"]["radial build"]["majorRadius"]
    Z0 = 0.0
    ϵ = 1 / gasc["INPUTS"]["radial build"]["aspectRatio"]
    κ = gasc["OUTPUTS"]["plasma parameters"]["elongation"]
    δ = gasc["INPUTS"]["plasma parameters"]["triangularity"]
    B0 = gasc["INPUTS"]["conductors"]["magneticFieldOnAxis"]
    ip = gasc["INPUTS"]["plasma parameters"]["plasmaCurrent"] * 1E6
    βn = gasc["OUTPUTS"]["plasma parameters"]["betaN"]
    ejima = gasc["INPUTS"]["plasma parameters"]["ejimaCoeff"]
    x_point = true
    symmetric = false
    resolution = 129

    # equilibrium
    init(eq; B0, R0, Z0, ϵ, δ, κ, beta_n = βn, ip, x_point = x_point, eq_kw...)
    eqactor = SolovevEquilibriumActor(dd, symmetric = symmetric)
    step(eqactor, verbose = false)
    finalize(eqactor, resolution, (maximum([R0 * (1 - ϵ * 2), 0.0]), R0 * (1 + ϵ * 2)), (-R0 * ϵ * κ * 1.5, R0 * ϵ * κ * 1.5))

    # core profiles
    init(cp, ejima = ejima; cp_kw...)

    # build
    norm = gasc["OUTPUTS"]["radial build"]["innerPlasmaRadius"]

    radial_build = DataStructures.OrderedDict()
    radial_build["gap_OH"] = gasc["OUTPUTS"]["radial build"]["innerSolenoidRadius"]
    radial_build["OH"] = gasc["INPUTS"]["radial build"]["rbOH"] * norm

    radial_build["hfs_gap_TF"] = gasc["INPUTS"]["radial build"]["gapTFOH"] * norm
    radial_build["hfs_TF"] = gasc["INPUTS"]["radial build"]["rbTF"] * norm
    if no_small_gaps
        radial_build["hfs_TF"] += radial_build["hfs_gap_TF"]
        pop!(radial_build, "hfs_gap_TF")
    end

    radial_build["hfs_gap_shield"] = gasc["INPUTS"]["radial build"]["gapBlanketCoil"] * norm
    radial_build["hfs_shield"] = gasc["INPUTS"]["radial build"]["rbInnerShield"] * norm
    if no_small_gaps
        radial_build["hfs_shield"] += radial_build["hfs_gap_shield"]
        pop!(radial_build, "hfs_gap_shield")
    end
    radial_build["hfs_blanket"] = gasc["INPUTS"]["radial build"]["rbInnerBlanket"] * norm

    radial_build["hfs_wall"] = gasc["INPUTS"]["radial build"]["gapInnerBlanketWall"] * norm
    radial_build["plasma"] = (gasc["INPUTS"]["radial build"]["majorRadius"] - sum(values(radial_build))) * 2
    radial_build["lfs_wall"] = gasc["INPUTS"]["radial build"]["gapOuterBlanketWall"] * norm

    radial_build["lfs_blanket"] = gasc["INPUTS"]["radial build"]["rbOuterBlanket"] * norm
    radial_build["lfs_shield"] = gasc["INPUTS"]["radial build"]["rbOuterShield"] * norm
    radial_build["lfs_gap_shield"] = gasc["INPUTS"]["radial build"]["gapBlanketCoil"] * norm
    if no_small_gaps
        radial_build["lfs_shield"] += radial_build["lfs_gap_shield"]
        pop!(radial_build, "lfs_gap_shield")
    end

    radial_build["lfs_TF"] = radial_build["hfs_TF"]
    radial_build["lfs_gap_TF"] = gasc["INPUTS"]["radial build"]["gapTFOH"] * norm
    if no_small_gaps
        radial_build["lfs_TF"] += radial_build["lfs_gap_TF"]
        pop!(radial_build, "lfs_gap_TF")
    end

    radial_build["gap_cryostat"] = radial_build["gap_OH"] * 3
    init(bd, radial_build; verbose)

    # TF coils
    bd.tf.coils_n = 16

    # cross-section outlines
    build_cx(bd, dd.equilibrium.time_slice[], 3)

    return dd
end


#= ============ =#
#  AbstractActor #
#= ============ =#
abstract type AbstractActor end

"""
    Base.step(actor::AbstractActor)

Placeholder function to take a step with a given actor
"""
function Base.step(actor::AbstractActor)
    error("Function step() not defined for actor of type $(typeof(actor))")
end

"""
    finalize(actor::AbstractActor)

Dummy function to finalize actor and store its in a IDS
"""
function finalize(actor::AbstractActor)
    actor
end

#= =========== =#
#  Equilibrium  #
#= =========== =#
include("actors/equilibrium_actor.jl")

export SolovevEquilibriumActor

#= ===== =#
#  Coils  #
#= ===== =#
include("actors/coils_actor.jl")

export PFcoilsOptActor

#= ===== =#
#  Build  #
#= ===== =#
include("actors/build_actor.jl")

#= ========= =#
#  Transport  #
#= ========= =#

include("actors/transport_actor.jl")
