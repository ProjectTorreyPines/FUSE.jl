import EGGO
import VacuumFields

#= ========= =#
#  ActorEGGO  #
#= ========= =#
@actor_parameters_struct ActorEGGO{T} begin
    #== actor parameters ==#
    model::Entry{Symbol} = Entry{Symbol}("-", "Neural network model to be used")
    use_vacuumfield_green::Entry{Bool} = Entry{Bool}("-", "Use Vacuum Fields green's function tables"; default=false)
    decimate_boundary::Entry{Int} = Entry{Int}("-", "Parameter to decimate number of boundary points"; default=4)
    timeslice_average::Entry{Int} = Entry{Int}("-", "Number of time slices to average"; default=0)
    #== display and debugging parameters ==#
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
    debug::Entry{Bool} = Entry{Bool}("-", "Print debug information withing EGGO solve"; default=false)
end



mutable struct ActorEGGO{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorEGGO{P}}
    act::ParametersAllActors{P}
    green::EGGO.GreenFunctionTables{Float64}
    basis_functions::EGGO.BasisFunctions{Float64}
    basis_functions_1d::EGGO.BasisFunctions1D{Float64}
    bf1d_itp::EGGO.BasisFunctions1Dinterp
    coils::Vector{<:VacuumFields.AbstractCoil}
    NNmodel::EGGO.NeuralNetModel{Float64}
end

"""
    ActorEGGO(dd::IMAS.dd, act::ParametersAllActors; kw...)

Runs the ML-based free boundary equilibrium solver EGGO
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
    coils = VacuumFields.MultiCoils(dd.pf_active)
    green.ggridfc_vf = VacuumFields.Green_table(green.rgrid, green.zgrid, coils)
    return ActorEGGO(dd, par, act, green, basis_functions, basis_functions_1d, bf1d_itp, coils, NNmodel)
end

function _step(actor::ActorEGGO{D,P}) where {D<:Real,P<:Real}
    dd = actor.dd
    par = actor.par

    eqt = dd.equilibrium.time_slice[]
    eqt1d = eqt.profiles_1d

    # prepare inputs
    rlim, zlim = IMAS.first_wall(dd.wall)
    wall = EGGO.Wall(rlim,zlim)
    psi_norm = range(0.0, 1.0, actor.green.nw)
    pp_target = IMAS.interp1d(eqt1d.psi_norm, eqt1d.dpressure_dpsi).(psi_norm) * 2π
    ffp_target = IMAS.interp1d(eqt1d.psi_norm, eqt1d.f_df_dpsi).(psi_norm) * 2π
    pp_fit, ffp_fit = EGGO.fit_ppffp(pp_target, ffp_target, actor.basis_functions_1d)

    # make actual prediction
    Ip_target = eqt.global_quantities.ip
    Rb_target = eqt.boundary.outline.r[1:par.decimate_boundary:end]
    Zb_target = eqt.boundary.outline.z[1:par.decimate_boundary:end]
    Rb_target[end] = Rb_target[1]
    Zb_target[end] = Zb_target[1]

    _, psirz, _ = EGGO.predict_model_from_boundary(
        Rb_target,
        Zb_target,
        pp_fit,
        ffp_fit,
        actor.NNmodel,
        actor.green,
        actor.basis_functions,
        actor.coils,
        Ip_target,
        par.use_vacuumfield_green
    )

    # average out EGGO solution with previous time slice(s)
    # until EGGO becomes a bit more robust
    n = par.timeslice_average # number of time slices to average
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
    Ψaxis, Raxis, Zaxis, Ψbnd, ffp, pp = EGGO.get_ΨaxisΨbndffppp(psirz, actor.green, actor.basis_functions, actor.basis_functions_1d, actor.bf1d_itp, wall, pp_fit, ffp_fit)#,Ip_target)    

    # ff' and p' edge offsets
    b0 = eqt.global_quantities.vacuum_toroidal_field.b0
    r0 = eqt.global_quantities.vacuum_toroidal_field.r0
    pend = eqt.profiles_1d.pressure[end]

    # fill out eqt quantities
    EGGO.fill_eqt(eqt, psirz, actor.green, wall, pp, ffp, b0, r0, pend, Ψbnd, Ψaxis, Raxis, Zaxis)

    eqt.global_quantities.free_boundary = 1 # free boundary

    return actor
end

