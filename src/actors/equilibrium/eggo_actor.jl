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
    green::Dict
    basis_functions::Dict
    basis_functions_1d::Dict
    NNmodel::Dict
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
    basis_functions = EGGO.get_basis_functions(model_name)
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

    # prepare inputs
    wall = Dict{Symbol,Vector{Float64}}()
    wall[:rlim], wall[:zlim] = IMAS.first_wall(dd.wall)
    psi_norm = collect(range(0.0, 1.0, actor.green[:nw]))
    pp_target = IMAS.interp1d(eqt1d.psi_norm, eqt1d.dpressure_dpsi).(psi_norm) * 2π
    ffp_target = IMAS.interp1d(eqt1d.psi_norm, eqt1d.f_df_dpsi).(psi_norm) * 2π
    ecurrt_target = fill(0.0, 6)
    bound_mxh = IMAS.MXH(eqt.boundary.outline.r, eqt.boundary.outline.z, 4)
    pp_fit, ffp_fit = EGGO.fit_ppffp(pp_target, ffp_target, actor.basis_functions_1d)

    # make actual prediction
    Jt, psirz, Ip, fcurrt = EGGO.predict_model(bound_mxh, pp_fit, ffp_fit, ecurrt_target, actor.NNmodel, actor.green, actor.basis_functions)

    # average out EGGO solution with previous time slice(s)
    # until EGGO becomes a bit more robust
    n = 2 # number of time slices to average
    i = IMAS.index(eqt)
    if i > n
        d = 1.0
        for k in 1:n
            eqt0 = dd.equilibrium.time_slice[i-k]
            eqt02d = findfirst(:rectangular, eqt0.profiles_2d)
            if size(eqt02d.psi) == size(psirz)
                psirz .+= (eqt02d.psi ./ 2pi)
                d += 1
            end
        end
        psirz ./= d
    end

    # pp' and ff' that were actually used in EGGO
    pp = zero(actor.basis_functions_1d[:pp][1, :])
    for k in eachindex(pp_fit)
        pp .+= pp_fit[k] .* actor.basis_functions_1d[:pp][k, :]
    end
    ffp = zero(actor.basis_functions_1d[:ffp][1, :])
    for k in eachindex(ffp_fit)
        ffp .+= ffp_fit[k] .* actor.basis_functions_1d[:ffp][k, :]
    end

    # ff' and p' edge offsets
    b0 = eqt.global_quantities.vacuum_toroidal_field.b0
    r0 = eqt.global_quantities.vacuum_toroidal_field.r0
    pend = eqt.profiles_1d.pressure[end]

    # fill out eqt quantities
    EGGO.fill_eqt(eqt, psirz, actor.green, wall, pp, ffp, b0, r0, pend)
    return actor
end
