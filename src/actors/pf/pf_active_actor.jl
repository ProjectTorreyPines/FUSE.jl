import MXHEquilibrium
import Optim
using LinearAlgebra
import VacuumFields

#= ============= =#
#  ActorPFactive  #
#= ============= =#
Base.@kwdef mutable struct FUSEparameters__ActorPFactive{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    green_model::Switch{Symbol} = Switch{Symbol}(options_green_model, "-", "Model used for the coils Green function calculations"; default=:simple)
    update_equilibrium::Entry{Bool} = Entry{Bool}("-", "Overwrite target equilibrium with the one that the coils can actually make"; default=false)
    do_plot::Entry{Bool} = Entry{Bool}("-", "Plot"; default=false)
end

mutable struct ActorPFactive{D,P} <: ReactorAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorPFactive{P}
    eq_out::IMAS.equilibrium{D}
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
    return ActorPFactive(
        dd,
        par,
        dd.equilibrium,
        VacuumFields.FluxControlPoint{Float64}[],
        VacuumFields.FluxControlPoint{Float64}[],
        VacuumFields.SaddleControlPoint{Float64}[],
        0.0,
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

    if actor.setup_cache === nothing
        eqt = dd.equilibrium.time_slice[]
        actor.setup_cache = setup(actor, eqt)
    end
    fixed_eq, image_eq, fixed_coils, pinned_coils, optim_coils, ψbound = actor.setup_cache

    # determine a good regularization value
    if actor.λ_regularize < 0.0
        actor.λ_regularize = VacuumFields.optimal_λ_regularize(
            vcat(pinned_coils, optim_coils),
            fixed_eq,
            image_eq,
            vcat(actor.boundary_control_points, actor.flux_control_points),
            actor.saddle_control_points;
            ψbound,
            fixed_coils
        )
    end

    # Find coil currents
    _, actor.cost = VacuumFields.find_coil_currents!(
        vcat(pinned_coils, optim_coils),
        fixed_eq,
        image_eq,
        vcat(actor.boundary_control_points, actor.flux_control_points),
        actor.saddle_control_points;
        ψbound,
        fixed_coils,
        actor.λ_regularize
    )

    return actor
end

"""
    _finalize(actor::ActorPFactive)

Update actor.eq_out 2D equilibrium PSI based on coils currents
"""
function _finalize(actor::ActorPFactive{D,P}) where {D<:Real,P<:Real}
    dd = actor.dd
    par = actor.par

    # evaluate current limits
    pf_current_limits(dd.pf_active, dd.build)

    if par.update_equilibrium || par.do_plot
        actor.eq_out = IMAS.lazycopy(dd.equilibrium)

        eqt_in = dd.equilibrium.time_slice[]
        eqt2d_in = findfirst(:rectangular, eqt_in.profiles_2d)
        eqt_out = actor.eq_out.time_slice[dd.global_time]
        eqt2d_out = findfirst(:rectangular, eqt_out.profiles_2d)

        if !ismissing(eqt_in.global_quantities, :ip)
            # convert dd.pf_active to coils for VacuumFields calculation
            coils = IMAS_pf_active__coils(dd; par.green_model)

            # convert equilibrium to MXHEquilibrium.jl format, since this is what VacuumFields uses
            EQfixed = IMAS2Equilibrium(eqt_in)

            # update ψ map
            scale_eq_domain_size = 1.0
            Rgrid = range(EQfixed.r[1] / scale_eq_domain_size, EQfixed.r[end] * scale_eq_domain_size; length=length(EQfixed.r))
            Zgrid = range(EQfixed.z[1] * scale_eq_domain_size, EQfixed.z[end] * scale_eq_domain_size; length=length(EQfixed.z))
            eqt2d_out.grid.dim1 = Rgrid
            eqt2d_out.grid.dim2 = Zgrid
            eqt2d_out.psi = collect(VacuumFields.fixed2free(EQfixed, coils, Rgrid, Zgrid)')
        end

        if par.do_plot
            display(plot(actor; equilibrium=true))
        end

        # update equilibrium psi2d and retrace flux surfaces
        if par.update_equilibrium
            eqt2d_in.psi = eqt2d_out.psi
            IMAS.flux_surfaces(eqt_in)
            actor.eq_out = eq_in
        end
    end

    return actor
end

function setup(actor::ActorPFactive, eqt::IMAS.equilibrium__time_slice)
    # Get coils organized by their function and initialize them
    fixed_coils, pinned_coils, optim_coils = fixed_pinned_optim_coils(actor)

    if ismissing(eqt.global_quantities, :ip) # field nulls
        fixed_eq = nothing
        image_eq = nothing
        rb, zb = eqt.boundary.outline.r, eqt.boundary.outline.z
        actor.boundary_control_points = VacuumFields.FluxControlPoints(rb, zb, eqt.global_quantities.psi_boundary)

    else # solutions with plasma
        fixed_eq = IMAS2Equilibrium(eqt)
        image_eq = VacuumFields.Image(fixed_eq)
        actor.boundary_control_points = VacuumFields.boundary_control_points(fixed_eq, 0.999)
    end

    saddle_weight = length(actor.boundary_control_points) / length(eqt.boundary.x_point)
    actor.saddle_control_points = [VacuumFields.SaddleControlPoint(x_point.r, x_point.z, saddle_weight) for x_point in eqt.boundary.x_point]

    return (fixed_eq=fixed_eq, image_eq=image_eq, fixed_coils=fixed_coils, pinned_coils=pinned_coils, optim_coils=optim_coils, ψbound=eqt.global_quantities.psi_boundary)
end

"""
    fixed_pinned_optim_coils(actor::ActorPFactive{D,P}) where {D<:Real,P<:Real}

Returns tuple of GS_IMAS_pf_active__coil structs organized by their function:

  - fixed: fixed position and current
  - pinned: coils with fixed position but current is optimized
  - optim: coils that have theri position and current optimized
"""
function fixed_pinned_optim_coils(actor::ActorPFactive{D,P}) where {D<:Real,P<:Real}
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
        if coil.identifier == "pinned"
            push!(pinned_coils, GS_IMAS_pf_active__coil(coil, coil_tech, par.green_model))
        elseif coil.identifier == "optim"
            push!(optim_coils, GS_IMAS_pf_active__coil(coil, coil_tech, par.green_model))
        elseif coil.identifier == "fixed"
            push!(fixed_coils, GS_IMAS_pf_active__coil(coil, coil_tech, par.green_model))
        else
            push!(pinned_coils, GS_IMAS_pf_active__coil(coil, coil_tech, par.green_model))
            #     error("$(IMAS.location(coil)).identifier=`$(coil.identifier)` is not valid. Accepted values are [\"optim\", \"pinned\", \"fixed\"]")
        end
    end

    return fixed_coils, pinned_coils, optim_coils
end