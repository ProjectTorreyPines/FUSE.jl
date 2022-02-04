#= =================== =#
#  init IMAS structures #
#= =================== =#
"""
    dummy `init`
"""
function init(ids::IMAS.IDS)
    error("Function init() not defined for IDS of type $(typeof(ids))")
end

function init_from_gasc(dd::IMAS.dd, filename, case; eq_kw = Dict(), cp_kw = Dict(), rb_kw = Dict())
    gasc = read_GASC(filename, case)

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
    norm = 0.570 / gasc["OUTPUTS"]["radial build"]["innerSolenoidRadius"]
    tmp = DataStructures.OrderedDict()
    tmp["gap_OH"] = gasc["OUTPUTS"]["radial build"]["innerSolenoidRadius"] * norm
    tmp["OH"] = (gasc["OUTPUTS"]["radial build"]["outerSolenoidRadius"] - gasc["OUTPUTS"]["radial build"]["innerSolenoidRadius"]) * norm
    tmp["hfs_TF"] = (gasc["OUTPUTS"]["radial build"]["outerRadiusTF"] - gasc["OUTPUTS"]["radial build"]["innerRadiusTF"]) * norm
    tmp["hfs_shield"] = gasc["OUTPUTS"]["radial build"]["rbInnerShield"] * norm
    tmp["hfs_blanket"] = gasc["OUTPUTS"]["radial build"]["rbInnerBlanket"] * norm
    tmp["hfs_wall"] = 0.1
    tmp["vacuum_vessel"] = (gasc["INPUTS"]["radial build"]["majorRadius"] - sum(values(tmp))) * 2
    tmp["lfs_wall"] = 0.1
    tmp["lfs_blanket"] = gasc["OUTPUTS"]["radial build"]["rbOuterBlanket"] * norm
    tmp["lfs_shield"] = gasc["OUTPUTS"]["radial build"]["rbOuterShield"] * norm
    tmp["lfs_TF"] = tmp["hfs_TF"] * norm
    tmp["gap_cryostat"] = tmp["gap_OH"] * 3 * norm
    init(bd, tmp)

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
