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
    eq_in::IMAS.equilibrium{D}
    eq_out::IMAS.equilibrium{D}
    λ_regularize::Float64
    flux_control_points::Vector{VacuumFields.FluxControlPoint{Float64}}
    saddle_control_points::Vector{VacuumFields.SaddleControlPoint{Float64}}
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
    return ActorPFactive(dd, par, dd.equilibrium, dd.equilibrium, 0.0, VacuumFields.FluxControlPoint{Float64}[], VacuumFields.SaddleControlPoint{Float64}[])
end

"""
    _step(actor::ActorPFactive)

Find currents that satisfy boundary and flux/saddle constraints in a least-square sense
"""
function _step(actor::ActorPFactive{T}) where {T<:Real}
    eqt = actor.eq_in.time_slice[]

    # Get coils (as GS3_IMAS_pf_active__coil) organized by their function and initialize them
    fixed_coils, pinned_coils, optim_coils = fixed_pinned_optim_coils(actor)

    # Find coil currents
    actor.λ_regularize = VacuumFields.find_coil_currents!(
        eqt, vcat(pinned_coils, optim_coils), fixed_coils;
        actor.λ_regularize, actor.flux_control_points, actor.saddle_control_points)

    return actor
end

"""
    _finalize(actor::ActorPFactive)

Update actor.eq_out 2D equilibrium PSI based on coils currents
"""
function _finalize(actor::ActorPFactive{D,P}) where {D<:Real,P<:Real}
    dd = actor.dd
    par = actor.par

    if par.update_equilibrium || par.do_plot
        actor.eq_out = IMAS.lazycopy(dd.equilibrium)

        eqt_in = actor.eq_in.time_slice[]
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
            actor.eq_out = actor.eq_in
        end
    end

    # evaluate current limits
    pf_current_limits(dd.pf_active, dd.build)

    return actor
end

"""
    VacuumFields.find_coil_currents!(
        eqt::IMAS.equilibrium__time_slice,
        coils::Vector{GS3_IMAS_pf_active__coil{D,D}},
        fixed_coils::Vector{GS3_IMAS_pf_active__coil{D,D}};
        λ_regularize::Real,
        flux_control_points::Vector{<:VacuumFields.FluxControlPoint}=VacuumFields.FluxControlPoint{Float64}[],
        saddle_control_points::Vector{<:VacuumFields.SaddleControlPoint}=VacuumFields.SaddleControlPoint{Float64}[]
    ) where {D<:Real}

Find PF/OH coil currents given an equilibrium time slice
"""
function VacuumFields.find_coil_currents!(
    eqt::IMAS.equilibrium__time_slice,
    coils::Vector{GS3_IMAS_pf_active__coil{D,D}},
    fixed_coils::Vector{GS3_IMAS_pf_active__coil{D,D}};
    λ_regularize::Real,
    flux_control_points::Vector{<:VacuumFields.FluxControlPoint}=VacuumFields.FluxControlPoint{Float64}[],
    saddle_control_points::Vector{<:VacuumFields.SaddleControlPoint}=VacuumFields.SaddleControlPoint{Float64}[]
) where {D<:Real}

    if ismissing(eqt.global_quantities, :ip) # field nulls
        psib = eqt.global_quantities.psi_boundary
        rb, zb = eqt.boundary.outline.r, eqt.boundary.outline.z
        boundary_control_points = VacuumFields.FluxControlPoints(rb, zb, psib)

    else # solutions with plasma
        fixed_eq = IMAS2Equilibrium(eqt)
        _, psib = MXHEquilibrium.psi_limits(fixed_eq)
        boundary_control_points = VacuumFields.boundary_control_points(fixed_eq, 0.999)
    end

    VacuumFields.find_coil_currents!(coils, fixed_eq, vcat(boundary_control_points, flux_control_points), saddle_control_points; ψbound=psib, fixed_coils, λ_regularize)

    return λ_regularize
end

function fixed_pinned_optim_coils(actor::ActorPFactive{D,P}) where {D<:Real,P<:Real}
    return fixed_pinned_optim_coils(actor, :currents, GS3_IMAS_pf_active__coil)
end
