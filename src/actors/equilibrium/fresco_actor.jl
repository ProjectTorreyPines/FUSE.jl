import FRESCO

#= =========== =#
#  ActorFRESCO  #
#= =========== =#
Base.@kwdef mutable struct FUSEparameters__ActorFRESCO{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    #== actor parameters ==#
    control::Switch{Symbol} = Switch{Symbol}([:vertical, :shape], "-", "Vertical control algorithm to be used"; default=:shape)
    number_of_iterations::Entry{Tuple{Int,Int}} = Entry{Tuple{Int,Int}}("-", "Number of outer and inner iterations"; default=(100, 3))
    relax::Entry{Float64} = Entry{Float64}("-", "Relaxation on the Picard iterations"; default=0.5)
    tolerance::Entry{Float64} = Entry{Float64}("-", "Tolerance for terminating iterations"; default=1e-4)
    nR::Entry{Int} = Entry{Int}("-", "Grid resolution along R"; default=129)
    nZ::Entry{Int} = Entry{Int}("-", "Grid resolution along Z"; default=129)
    reuse_Green_table::Entry{Bool} = Entry{Bool}("-", "Reuse Green function table"; default=false)
    #== display and debugging parameters ==#
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
    debug::Entry{Bool} = Entry{Bool}("-", "Print debug information withing FRESCO solve"; default=false)
end

mutable struct ActorFRESCO{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorFRESCO{P}
    act::ParametersAllActors{P}
    canvas::Union{Nothing,FRESCO.Canvas}
    profile::Union{Nothing,FRESCO.PressureJt}
end

"""
    ActorFRESCO(dd::IMAS.dd, act::ParametersAllActors; kw...)

Runs the Fixed boundary equilibrium solver FRESCO
"""
function ActorFRESCO(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorFRESCO(dd, act.ActorFRESCO, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorFRESCO(dd::IMAS.dd{D}, par::FUSEparameters__ActorFRESCO{P}, act::ParametersAllActors{P}; kw...) where {D<:Real,P<:Real}
    logging_actor_init(ActorFRESCO)
    par = par(kw...)
    return ActorFRESCO(dd, par, act, nothing, nothing)
end

"""
    _step(actor::ActorFRESCO)

Runs FRESCO on the r_z boundary, equilibrium pressure and equilibrium j_tor
"""
function _step(actor::ActorFRESCO{D,P}) where {D<:Real,P<:Real}
    dd = actor.dd
    par = actor.par
    act = actor.act

    # FRESCO requires wall information
    if IMAS.hasdata(dd.wall)
        fw_r, fw_z = IMAS.first_wall(dd.wall)
    elseif IMAS.hasdata(dd.pf_active) && findfirst(:rectangular, dd.equilibrium.time_slice[].profiles_2d) !== nothing
        fw_r, fw_z = IMAS.first_wall(dd.equilibrium.time_slice[], dd.pf_active)
    elseif IMAS.hasdata(dd.pf_active)
        fw_r, fw_z = IMAS.first_wall(dd.pf_active)
    else
        error("FRESCO needs at least pf_active.coils to be filled")
    end
    ΔR = maximum(fw_r) - minimum(fw_r)
    ΔZ = maximum(fw_z) - minimum(fw_z)
    Rs = range(max(0.01, minimum(fw_r) - ΔR / 100), maximum(fw_r) + ΔR / 100, par.nR)
    Zs = range(minimum(fw_z) - ΔZ / 100, maximum(fw_z) + ΔZ / 100, par.nZ)

    if actor.canvas !== nothing && par.reuse_Green_table
        Green_table = actor.canvas._Gvac
        Nr, Nz, _ = size(Green_table)
        @assert Nr == par.nR && Nz == par.nZ "nR and/or nZ have changed; set reuse_Green_table to false"
    else
        Green_table = D[;;;]
    end
    actor.canvas = FRESCO.Canvas(dd, Rs, Zs; load_pf_passive=false, Green_table,
                                 act.ActorPFactive.strike_points_weight, act.ActorPFactive.x_points_weight)

    actor.profile = FRESCO.PressureJt(dd)

    FRESCO.solve!(actor.canvas, actor.profile, par.number_of_iterations...; par.relax, par.debug, par.control, par.tolerance)

    return actor
end

# finalize by converting FRESCO canvas to dd.equilibrium
function _finalize(actor::ActorFRESCO)
    canvas = actor.canvas
    profile = actor.profile
    dd = actor.dd
    eq = dd.equilibrium
    eqt = eq.time_slice[]
    eqt1d = eqt.profiles_1d
    eq2d = resize!(eqt.profiles_2d, 1)[1]

    eqt.global_quantities.magnetic_axis.r = canvas.Raxis
    eqt.global_quantities.magnetic_axis.z = canvas.Zaxis
    eqt.global_quantities.psi_boundary = canvas.Ψbnd
    eqt.global_quantities.psi_axis = canvas.Ψaxis

    Npsi = length(eqt1d.psi)
    eqt1d.psi = range(canvas.Ψaxis, canvas.Ψbnd, Npsi)
    eqt1d.dpressure_dpsi = FRESCO.pprime.(Ref(canvas), Ref(profile), eqt1d.psi_norm)
    eqt1d.f_df_dpsi = FRESCO.ffprime.(Ref(canvas), Ref(profile), eqt1d.psi_norm)

    @ddtime(eq.vacuum_toroidal_field.b0 = eqt.global_quantities.vacuum_toroidal_field.b0)
    eq.vacuum_toroidal_field.r0 = eqt.global_quantities.vacuum_toroidal_field.r0

    fend = eqt.global_quantities.vacuum_toroidal_field.b0 * eqt.global_quantities.vacuum_toroidal_field.r0
    f2 = 2 * IMAS.cumtrapz(eqt1d.psi, eqt1d.f_df_dpsi)
    f2 .= f2 .- f2[end] .+ fend^2
    eqt1d.f = sign(fend) .* sqrt.(f2)

    pend = eqt1d.pressure[end]
    eqt1d.pressure = IMAS.cumtrapz(eqt1d.psi, eqt1d.dpressure_dpsi)
    eqt1d.pressure .+= pend .- eqt1d.pressure[end]

    eq2d.grid_type.index = 1
    eq2d.grid.dim1 = collect(range(canvas.Rs[1], canvas.Rs[end], Npsi))
    eq2d.grid.dim2 = collect(range(canvas.Zs[1], canvas.Zs[end], Npsi))
    FRESCO.update_interpolation!(canvas)
    eq2d.psi = [canvas._Ψitp(r, z) for r in eq2d.grid.dim1, z in eq2d.grid.dim2]

    return actor
end
