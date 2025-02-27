import EGGO

#= ========= =#
#  ActorEGGO  #
#= ========= =#
Base.@kwdef mutable struct FUSEparameters__ActorEGGO{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    #== actor parameters ==#
    model::Entry{Symbol} = Entry{Symbol}("-", "Neural network model to be used")
    #== display and debugging parameters ==#
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
    debug::Entry{Bool} = Entry{Bool}("-", "Print debug information withing EGGO solve"; default=false)
end

mutable struct ActorEGGO{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorEGGO{P}
    act::ParametersAllActors{P}
    green::Dict{Symbol,Any}
    basis_functions::Dict{Symbol,Any}
    basis_functions_1d::Dict{Any,Any}
    NNmodel::Dict{Any,Any}
end

"""
    ActorEGGO(dd::IMAS.dd, act::ParametersAllActors; kw...)

Runs the Fixed boundary equilibrium solver EGGO
"""
function ActorEGGO(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorEGGO(dd, act.ActorEGGO, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorEGGO(dd::IMAS.dd{D}, par::FUSEparameters__ActorEGGO{P}, act::ParametersAllActors{P}; kw...) where {D<:Real,P<:Real}
    logging_actor_init(ActorEGGO)
    par = par(kw...)

    model_name = :d3d_efit01
    green = EGGO.get_greens_function_tables(model_name)
    basis_functions = EGGO.get_basis_functions(model_name, green)
    NNmodel = EGGO.get_model(model_name)
    basis_functions_1d, _ = EGGO.get_basis_functions_1d(model_name)

    return ActorEGGO(dd, par, act, green, basis_functions, basis_functions_1d, NNmodel)
end

"""
    _step(actor::ActorEGGO)

Runs EGGO on the r_z boundary, equilibrium pressure and equilibrium j_tor
"""
function _step(actor::ActorEGGO{D,P}) where {D<:Real,P<:Real}
    dd = actor.dd

    eqt = dd.equilibrium.time_slice[]
    eqt1d = eqt.profiles_1d

    wall = Dict{Symbol,Vector{Float64}}()
    wall[:rlim], wall[:zlim] = IMAS.first_wall(dd.wall)

    Rb_target, Zb_target = eqt.boundary.outline.r, eqt.boundary.outline.z
    b0 = eqt.global_quantities.vacuum_toroidal_field.b0
    r0 = eqt.global_quantities.vacuum_toroidal_field.r0
    pend = eqt.profiles_1d.pressure[end]

    psi_norm = collect(range(0.0, 1.0, actor.green[:nw]))
    pp_target = IMAS.interp1d(eqt1d.psi_norm, eqt1d.dpressure_dpsi).(psi_norm) * 2π
    ffp_target = IMAS.interp1d(eqt1d.psi_norm, eqt1d.f_df_dpsi).(psi_norm) * 2π
    ecurrt_target = fill(0.0, 6)

    Jt, psirz, Ip, fcurrt =
        EGGO.predict_model(Rb_target, Zb_target, pp_target, ffp_target, ecurrt_target, actor.NNmodel, actor.green, actor.basis_functions, actor.basis_functions_1d)

    EGGO.get_surfaces(eqt, Matrix(transpose(psirz)), Ip, fcurrt, actor.green, wall, Rb_target, Zb_target, pp_target, ffp_target, ecurrt_target, b0, r0, pend)

    return actor
end
