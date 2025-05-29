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
    use_vacuumfield_green:: Entry{Bool}("-", "Use Vacuum Fields green's function tables")
    #== display and debugging parameters ==#
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
    debug::Entry{Bool} = Entry{Bool}("-", "Print debug information withing EGGO solve"; default=false)
end

mutable struct ActorEGGO{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorEGGO{P}}
    act::ParametersAllActors{P}
    green::Dict
    basis_functions::Dict
    basis_functions_1d::Dict
    bf1d_itp::Dict
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
    par = OverrideParameters(par; kw...)

    model_name = :d3d_efit01efit02cake02
    green = EGGO.get_greens_function_tables(model_name)
    basis_functions = EGGO.get_basis_functions(model_name)
    NNmodel = EGGO.get_model(model_name)
    basis_functions_1d, bf1d_itp = EGGO.get_basis_functions_1d(model_name)
    return ActorEGGO(dd, par, act, green, basis_functions, basis_functions_1d,bf1d_itp, NNmodel)
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
    pp_fit, ffp_fit = EGGO.fit_ppffp(pp_target, ffp_target, actor.basis_functions_1d)
    coeff_max = 0.1
    pp_index = findfirst(x -> x > coeff_max || x<-1*coeff_max, pp_fit[2:end]) 
    ffp_index = findfirst(x -> x > coeff_max || x<-1*coeff_max, ffp_fit[2:end]) 
    pp_index = isnothing(pp_index) ? length(pp_fit) : pp_index
    ffp_index = isnothing(ffp_index) ? length(ffp_fit)  : ffp_index

    #Ip_pp = -IMAS.trapz(eqt1d.area,eqt1d.dpressure_dpsi*2π/ eqt1d.gm9)
    #Ip_ffp = -IMAS.trapz(eqt1d.area,eqt1d.dpressure_dpsi*eqt1d.gm1*2π/eqt1d.gm9)
    #@show(Ip_pp,Ip_ffp)
    #pp_index = 4
    #ffp_index = 6
    #pp_fit*=0.0
    #ffp_fit*=0.0
    #pp_fit[1:pp_index], ffp_fit[1:ffp_index] = EGGO.fit_ppffp(pp_target, ffp_target, actor.basis_functions_1d,pp_index,ffp_index)
    # make actual prediction
    Ip_target = eqt.global_quantities.ip
    coils = VacuumFields.IMAS_pf_active__coils(dd; green_model=:quad, zero_currents=false)
    if actor.use_vacuumfield_green
        green[:ggridfc] = VacuumFields.Green_table(green[:rgrid], green[:zgrid], coils)
    end
    Jt, psirz, Ip = EGGO.predict_model_from_boundary(eqt.boundary.outline.r, eqt.boundary.outline.z, pp_fit, ffp_fit, actor.NNmodel, actor.green, actor.basis_functions,coils,Ip_target)

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
    Ψaxis,Raxis,Zaxis,Ψbnd,ffp,pp = EGGO.get_ΨaxisΨbndffppp(psirz, actor.green, actor.basis_functions,actor.basis_functions_1d, actor.bf1d_itp, wall, pp_fit, ffp_fit,Ip_target)    

    # ff' and p' edge offsets
    b0 = eqt.global_quantities.vacuum_toroidal_field.b0
    r0 = eqt.global_quantities.vacuum_toroidal_field.r0
    pend = eqt.profiles_1d.pressure[end]
    # fill out eqt quantities
    EGGO.fill_eqt(eqt, psirz, actor.green, wall, pp, ffp, b0, r0, pend,Ψbnd,Ψaxis,Raxis,Zaxis)
    return actor
end

