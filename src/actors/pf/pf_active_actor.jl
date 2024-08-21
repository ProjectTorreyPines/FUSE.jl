import MXHEquilibrium
import Optim
using LinearAlgebra
import VacuumFields
import VacuumFields: GS_IMAS_pf_active__coil

#= ============= =#
#  ActorPFactive  #
#= ============= =#
Base.@kwdef mutable struct FUSEparameters__ActorPFactive{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    green_model::Switch{Symbol} = Switch{Symbol}(options_green_model, "-", "Model used for the coils Green function calculations"; default=:quad)
    update_equilibrium::Entry{Bool} = Entry{Bool}("-", "Overwrite target equilibrium with the one that the coils can actually make"; default=false)
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
end

mutable struct ActorPFactive{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorPFactive{P}
    eqt_out::IMAS.equilibrium__time_slice{D}
    boundary_control_points::Vector{VacuumFields.FluxControlPoint{Float64}}
    flux_control_points::Vector{VacuumFields.FluxControlPoint{Float64}}
    saddle_control_points::Vector{VacuumFields.SaddleControlPoint{Float64}}
    λ_regularize::Float64
    cost::Float64
    setup_cache::Any
end

"""
    ActorPFactive(dd::IMAS.dd, act::ParametersAllActors; kw...)

Finds the optimal coil currents to match the equilibrium boundary shape

!!! note

    Manupulates data in `dd.pf_active`
"""
function ActorPFactive(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorPFactive(dd, act.ActorPFactive; kw...)
    finalize(step(actor))
    return actor
end

function ActorPFactive(dd::IMAS.dd, par::FUSEparameters__ActorPFactive; kw...)
    logging_actor_init(ActorPFactive)
    par = par(kw...)

    boundary_control_points, flux_control_points, saddle_control_points = default_control_points(dd.equilibrium.time_slice[], dd.pulse_schedule.position_control)

    return ActorPFactive(
        dd,
        par,
        dd.equilibrium.time_slice[],
        boundary_control_points,
        flux_control_points,
        saddle_control_points,
        -1.0,
        NaN,
        nothing
    )
end

"""
    _step(actor::ActorPFactive)

Find currents that satisfy boundary and flux/saddle constraints in a least-square sense
"""
function _step(actor::ActorPFactive{T}) where {T<:Real}
    dd = actor.dd

    # setup (with caching, used for coil optimization)
    if actor.setup_cache === nothing
        eqt = dd.equilibrium.time_slice[]
        actor.setup_cache = setup(actor, eqt)
    end
    fixed_eq, image_eq, fixed_coils, pinned_coils, optim_coils, ψbound = actor.setup_cache

    # determine a good regularization value
    if actor.λ_regularize < 0.0
        if fixed_eq === nothing
            actor.λ_regularize = 0.0
        else
            actor.λ_regularize = VacuumFields.optimal_λ_regularize(
                vcat(pinned_coils, optim_coils),
                fixed_eq,
                image_eq;
                flux_cps = vcat(actor.boundary_control_points, actor.flux_control_points),
                saddle_cps = actor.saddle_control_points,
                ψbound,
                fixed_coils
            )
        end
    end

    # Find coil currents
    _, actor.cost = VacuumFields.find_coil_currents!(
        vcat(pinned_coils, optim_coils),
        fixed_eq,
        image_eq;
        flux_cps = vcat(actor.boundary_control_points, actor.flux_control_points),
        saddle_cps = actor.saddle_control_points,
        ψbound,
        fixed_coils,
        actor.λ_regularize
    )

    return actor
end

"""
    _finalize(actor::ActorPFactive)

Update actor.eqt_out 2D equilibrium PSI based on coils currents
"""
function _finalize(actor::ActorPFactive{D,P}) where {D<:Real,P<:Real}
    dd = actor.dd
    par = actor.par

    # evaluate current limits
    pf_current_limits(dd.pf_active, dd.build)

    # evaluate eq_out
    eqt_in = dd.equilibrium.time_slice[]
    eqt2d_in = findfirst(:rectangular, eqt_in.profiles_2d)
    actor.eqt_out = eqt_out = deepcopy(eqt_in)
    setfield!(actor.eqt_out, :_parent, WeakRef(dd.equilibrium.time_slice))
    eqt2d_out = findfirst(:rectangular, eqt_out.profiles_2d)
    if !ismissing(eqt_in.global_quantities, :ip)
        # convert dd.pf_active to coils for VacuumFields calculation
        coils = VacuumFields.IMAS_pf_active__coils(dd; par.green_model, zero_currents=false)

        # convert equilibrium to MXHEquilibrium.jl format, since this is what VacuumFields uses
        EQfixed = IMAS2Equilibrium(eqt_in)

        # update ψ map
        scale_eq_domain_size = 1.5
        Rgrid = range(EQfixed.r[1] / scale_eq_domain_size, EQfixed.r[end] * scale_eq_domain_size; length=length(EQfixed.r))
        Zgrid = range(EQfixed.z[1] * scale_eq_domain_size, EQfixed.z[end] * scale_eq_domain_size; length=length(EQfixed.z))
        eqt2d_out.grid.dim1 = Rgrid
        eqt2d_out.grid.dim2 = Zgrid
        eqt2d_out.psi = collect(VacuumFields.fixed2free(EQfixed, coils, Rgrid, Zgrid)')
    end

    if par.do_plot
        display(plot(actor))
    end

    # update equilibrium psi2d and retrace flux surfaces
    if par.update_equilibrium
        eqt2d_in.psi = collect(VacuumFields.fixed2free(EQfixed, coils, eqt2d_in.grid.dim1, eqt2d_in.grid.dim2)')
        IMAS.flux_surfaces(eqt_in)
    end

    return actor
end

function default_control_points(eqt::IMAS.equilibrium__time_slice, pc::IMAS.pulse_schedule__position_control; saddle_weight::Float64=0.01, strike_weight::Float64=0.01)
    psib = eqt.global_quantities.psi_boundary
    if ismissing(eqt.global_quantities, :ip) # field nulls
        fixed_eq = nothing
        rb, zb = eqt.boundary.outline.r, eqt.boundary.outline.z
        boundary_control_points = VacuumFields.FluxControlPoints(rb, zb, psib)
    else # solutions with plasma
        fixed_eq = IMAS2Equilibrium(eqt)
        boundary_control_points = VacuumFields.boundary_control_points(fixed_eq, 0.999)
    end

    saddle_control_points = VacuumFields.SaddleControlPoint{Float64}[]
    if saddle_weight == 0.0
        # pass
    elseif !isempty(pc.x_point)
        # we favor taking the x-points from the pulse schedule, if available
        for x_point in pc.x_point
            r = @ddtime(x_point.r.reference)
            if r == 0.0 || isnan(r)
                continue
            end
            z = @ddtime(x_point.z.reference)
            push!(saddle_control_points, VacuumFields.SaddleControlPoint(r, z, saddle_weight))
        end
    else
        for x_point in eqt.boundary.x_point
            push!(saddle_control_points, VacuumFields.SaddleControlPoint(x_point.r, x_point.z, saddle_weight))
        end
    end

    psib = eqt.global_quantities.psi_boundary
    flux_control_points = VacuumFields.FluxControlPoint{Float64}[]
    if strike_weight == 0.0
        # pass
    elseif !isempty(pc.strike_point)
        # we favor taking the strike points from the pulse schedule, if available
        flux_control_points = VacuumFields.FluxControlPoint{Float64}[]
        for strike_point in pc.strike_point
            r = @ddtime(strike_point.r.reference)
            if r == 0.0 || isnan(r)
                continue
            end
            z = @ddtime(strike_point.z.reference)
            push!(flux_control_points, VacuumFields.FluxControlPoint(r, z, psib, strike_weight))
        end
    else
        for strike_point in eqt.boundary.strike_point
            push!(flux_control_points, VacuumFields.FluxControlPoint(strike_point.r, strike_point.z, psib, strike_weight))
        end
    end

    return boundary_control_points, flux_control_points, saddle_control_points
end

function setup(actor::ActorPFactive, eqt::IMAS.equilibrium__time_slice)
    if ismissing(eqt.global_quantities, :ip) # field nulls
        fixed_eq = nothing
        image_eq = nothing

    else # solutions with plasma
        fixed_eq = IMAS2Equilibrium(eqt)
        image_eq = VacuumFields.Image(fixed_eq)
    end

    # Get coils organized by their function and initialize them
    fixed_coils, pinned_coils, optim_coils = fixed_pinned_optim_coils(actor; zero_currents=true)

    return (fixed_eq=fixed_eq, image_eq=image_eq, fixed_coils=fixed_coils, pinned_coils=pinned_coils, optim_coils=optim_coils, ψbound=eqt.global_quantities.psi_boundary)
end

"""
    fixed_pinned_optim_coils(actor::ActorPFactive{D,P}; zero_currents::Bool) where {D<:Real,P<:Real}

Returns tuple of GS_IMAS_pf_active__coil structs organized by their function:

  - fixed: fixed position and current
  - pinned: coils with fixed position but current is optimized
  - optim: coils that have theri position and current optimized
"""
function fixed_pinned_optim_coils(actor::ActorPFactive{D,P}; zero_currents::Bool) where {D<:Real,P<:Real}
    dd = actor.dd
    par = actor.par

    fixed_coils = GS_IMAS_pf_active__coil{D,D}[]
    pinned_coils = GS_IMAS_pf_active__coil{D,D}[]
    optim_coils = GS_IMAS_pf_active__coil{D,D}[]
    for coil in dd.pf_active.coil
        if IMAS.is_ohmic_coil(coil)
            coil_tech = dd.build.oh.technology
        else
            coil_tech = dd.build.pf_active.technology
        end
        if zero_currents
            @ddtime(coil.current.data = 0.0)   # zero currents for all coils
        end
        pfcoil = GS_IMAS_pf_active__coil(coil, coil_tech, par.green_model)
        if :shaping ∉ [IMAS.index_2_name(coil.function)[f.index] for f in coil.function]
            push!(fixed_coils, pfcoil)
        elseif coil.identifier == "optim"
            push!(optim_coils, pfcoil)
        elseif coil.identifier == "fixed"
            push!(fixed_coils, pfcoil)
        else
            push!(pinned_coils, pfcoil)
        end
    end

    # push!(fixed_coils, popat!(optim_coils,1))
    # VacuumFields.imas(fixed_coils[end]).identifier = "fixed"
    # push!(fixed_coils, popat!(optim_coils,length(optim_coils)))
    # VacuumFields.imas(fixed_coils[end]).identifier = "fixed"

    return (fixed_coils=fixed_coils, pinned_coils=pinned_coils, optim_coils=optim_coils)
end