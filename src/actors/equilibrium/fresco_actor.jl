import FRESCO

#= =========== =#
#  ActorFRESCO  #
#= =========== =#
Base.@kwdef mutable struct FUSEparameters__ActorFRESCO{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    #== actor parameters ==#
    control::Switch{Symbol} = Switch{Symbol}([:vertical, :shape], "-", ""; default=:shape)
    relax::Entry{Float64} = Entry{Float64}("-", "Relaxation on the Picard iterations"; default=0.5)
    #tolerance::Entry{Float64} = Entry{Float64}("-", "Tolerance for terminating iterations"; default=1e-4)
    #== data flow parameters ==#
    # fixed_grid::Switch{Symbol} = Switch{Symbol}([:poloidal, :toroidal], "-", "Fix P and Jt on this rho grid"; default=:toroidal)
    #== display and debugging parameters ==#
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
    debug::Entry{Bool} = Entry{Bool}("-", "Print debug information withing FRESCO solve"; default=false)
    #== IMAS psi grid settings ==#
    nR::Entry{Int} = Entry{Int}("-", "Grid resolution along R"; default=129)
    nZ::Entry{Int} = Entry{Int}("-", "Grid resolution along Z"; default=129)
end

mutable struct ActorFRESCO{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorFRESCO{P}
    canvas::Union{Nothing,FRESCO.Canvas}
end

"""
    ActorFRESCO(dd::IMAS.dd, act::ParametersAllActors; kw...)

Runs the Fixed boundary equilibrium solver FRESCO
"""
function ActorFRESCO(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorFRESCO(dd, act.ActorFRESCO; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorFRESCO(dd::IMAS.dd{D}, par::FUSEparameters__ActorFRESCO{P}; kw...) where {D<:Real,P<:Real}
    logging_actor_init(ActorFRESCO)
    par = par(kw...)
    return ActorFRESCO(dd, par, nothing)
end

"""
    _step(actor::ActorFRESCO)

Runs FRESCO on the r_z boundary, equilibrium pressure and equilibrium j_tor
"""
function _step(actor::ActorFRESCO)
    dd = actor.dd
    par = actor.par
    eqt = dd.equilibrium.time_slice[]
    eqt1d = eqt.profiles_1d

    # Ip_target = eqt.global_quantities.ip
    #
    # if par.fixed_grid === :poloidal
    #     rho = sqrt.(eqt1d.psi_norm)
    #     rho[1] = 0.0
    #     P = (FRESCO.FE(rho, eqt1d.pressure), :poloidal)
    #     # don't allow current to change sign
    #     Jt = (FRESCO.FE(rho, [sign(j) == sign(Ip_target) ? j : 0.0 for j in eqt1d.j_tor]), :poloidal)
    #     Pbnd = eqt1d.pressure[end]
    # elseif par.fixed_grid === :toroidal
    #     rho = eqt1d.rho_tor_norm
    #     P = (FRESCO.FE(rho, eqt1d.pressure), :toroidal)
    #     # don't allow current to change sign
    #     Jt = (FRESCO.FE(rho, [sign(j) == sign(Ip_target) ? j : 0.0 for j in eqt1d.j_tor]), :toroidal)
    #     Pbnd = eqt1d.pressure[end]
    # end

    gpp = IMAS.interp1d(eqt1d.psi_norm, eqt1d.dpressure_dpsi, :cubic)
    gffp = IMAS.interp1d(eqt1d.psi_norm, eqt1d.f_df_dpsi, :cubic)
    profile = FRESCO.PprimeFFprime(x -> gpp(x), x -> gffp(x))

    actor.canvas = FRESCO.Canvas(dd, par.nR, par.nZ)

    FRESCO.solve!(actor.canvas, profile, 100, 3; relax=par.relax, par.debug, par.control)

    # using Plots
    # p1 = Plots.heatmap(Rs, Zs, psi', aspect_ratio=:equal, xrange=(0,12.5), yrange=(-9,9), size=(400,500))
    # Plots.plot!(p1, coils, label=nothing)
    # Plots.contour!(p1, Rs, Zs, psi', levels=[C.Ψbnd], color=:white)
    # display(p1)

    return actor
end

# finalize by converting FRESCO canvas to dd.equilibrium
function _finalize(actor::ActorFRESCO)
    canvas = actor.canvas
    dd = actor.dd
    eq = dd.equilibrium
    eqt = eq.time_slice[]
    eqt1d = eqt.profiles_1d
    eq2d = resize!(eqt.profiles_2d, 1)[1]

    Raxis, Zaxis, _ = FRESCO.find_axis(canvas)
    eqt.global_quantities.magnetic_axis.r = Raxis
    eqt.global_quantities.magnetic_axis.z = Zaxis
    eqt1d.psi = range(canvas.Ψaxis, canvas.Ψbnd, length(eqt1d.psi))
    # p, p', f, ff' don't change

    eq2d.grid_type.index = 1
    eq2d.grid.dim1 = canvas.Rs
    eq2d.grid.dim2 = canvas.Zs
    eq2d.psi = canvas.Ψ

    return actor
end
