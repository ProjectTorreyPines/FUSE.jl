import MXHEquilibrium
import Optim
using LinearAlgebra
import VacuumFields
import VacuumFields: GS_IMAS_pf_active__coil

#= ============= =#
#  ActorPFactive  #
#= ============= =#
@actor_parameters_struct ActorPFactive{T} begin
    green_model::Switch{Symbol} = Switch{Symbol}(options_green_model, "-", "Model used for the coils Green function calculations"; default=:quad)
    update_equilibrium::Entry{Bool} = Entry{Bool}("-", "Overwrite target equilibrium with the one that the coils can actually make"; default=false)
    x_points_weight::Entry{Float64} = Entry{Float64}("-", "Weight givent to x-point constraints"; default=0.1)
    strike_points_weight::Entry{Float64} = Entry{Float64}("-", "Weight givent to strike-point constraints"; default=0.1)
    magnetic_probe_weight::Entry{Float64} = Entry{Float64}("-", "Weight givent to magnetic probes measurements"; default=1.0)
    flux_loop_weight::Entry{Float64} = Entry{Float64}("-", "Weight givent to flux loops measurements"; default=1.0)
    boundary_weight::Entry{Float64} = Entry{Float64}("-", "Weight givent to the boundary"; default=1.0)
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
end

mutable struct ActorPFactive{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorPFactive{P}}
    eqt_out::IMAS.equilibrium__time_slice{D}
    iso_control_points::Vector{VacuumFields.IsoControlPoint{D}}
    saddle_control_points::Vector{VacuumFields.SaddleControlPoint{D}}
    flux_control_points::Vector{VacuumFields.FluxControlPoint{D}}
    field_control_points::Vector{VacuumFields.FieldControlPoint{D}}
    λ_regularize::Float64
    cost::Float64
    setup_cache::Any
end

"""
    ActorPFactive(dd::IMAS.dd, act::ParametersAllActors; kw...)

Optimizes poloidal field coil currents to match target equilibrium constraints.

The actor solves for optimal PF coil currents that satisfy boundary shape, X-point positions,
and strike point locations in a least-squares sense. It uses Green's functions to model
electromagnetic coupling between coils and plasma, with optional regularization to prevent
unrealistic current solutions.

Key features:

  - Supports different Green's function calculation methods (quadrature, filament models)
  - Handles boundary shape control through iso-flux control points
  - Maintains X-point and strike point positions with configurable weights
  - Includes automatic regularization parameter selection for stable solutions
  - Can update target equilibrium with achievable coil configuration

Control points:

  - Boundary iso-flux surfaces for plasma shape control
  - X-point locations as saddle point constraints
  - Strike point locations as flux surface constraints
  - Magnetic axis position and flux value

!!! note

    Manipulates data in `dd.pf_active`
"""
function ActorPFactive(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorPFactive(dd, act.ActorPFactive; kw...)
    finalize(step(actor))
    return actor
end

function ActorPFactive(dd::IMAS.dd, par::FUSEparameters__ActorPFactive; kw...)
    logging_actor_init(ActorPFactive)
    par = OverrideParameters(par; kw...)

    eqt_control_points = VacuumFields.equilibrium_control_points(dd.equilibrium.time_slice[], dd.pulse_schedule.position_control;
        par.boundary_weight,
        par.x_points_weight,
        par.strike_points_weight)
    mag_control_points = VacuumFields.magnetic_control_points(dd.magnetics;
        par.flux_loop_weight,
        par.magnetic_probe_weight)

    return ActorPFactive(
        dd,
        par,
        dd.equilibrium.time_slice[],
        [eqt_control_points.iso_cps; mag_control_points.iso_cps],
        eqt_control_points.saddle_cps,
        [eqt_control_points.flux_cps; mag_control_points.flux_cps],
        mag_control_points.field_cps,
        -1.0,
        NaN,
        nothing)
end

"""
    _step(actor::ActorPFactive)

Solves the least-squares optimization for optimal PF coil currents.

The step function performs the core optimization by:

 1. Setting up coil system (fixed, pinned, and optimizable coils) with caching
 2. Computing optimal regularization parameter if not specified
 3. Solving the regularized least-squares problem for coil currents
 4. Evaluating cost function to assess solution quality

The optimization balances satisfaction of electromagnetic constraints with
regularization to prevent unrealistic current magnitudes.
"""
function _step(actor::ActorPFactive{T}) where {T<:Real}
    dd = actor.dd

    # setup coils information (with caching)
    if actor.setup_cache === nothing
        eqt = dd.equilibrium.time_slice[]
        if ismissing(eqt.global_quantities, :ip) # field nulls
            fixed_eq = nothing
            image_eq = nothing

        else # solutions with plasma
            fixed_eq = IMAS2Equilibrium(eqt)
            image_eq = VacuumFields.Image(fixed_eq)
        end

        # Get coils organized by their function and initialize them
        fixed_coils, pinned_coils, optim_coils = fixed_pinned_optim_coils(actor; zero_currents=true)

        actor.setup_cache =
            (fixed_eq=fixed_eq, image_eq=image_eq, fixed_coils=fixed_coils, pinned_coils=pinned_coils, optim_coils=optim_coils, ψbound=eqt.global_quantities.psi_boundary)
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
                iso_cps=actor.iso_control_points,
                saddle_cps=actor.saddle_control_points,
                flux_cps=actor.flux_control_points,
                field_cps=actor.field_control_points,
                ψbound,
                fixed_coils)
        end
    end

    # Find coil currents
    _, actor.cost = VacuumFields.find_coil_currents!(
        vcat(pinned_coils, optim_coils),
        fixed_eq,
        image_eq;
        iso_cps=actor.iso_control_points,
        saddle_cps=actor.saddle_control_points,
        flux_cps=actor.flux_control_points,
        field_cps=actor.field_control_points,
        ψbound,
        fixed_coils,
        actor.λ_regularize)

    return actor
end

"""
    _finalize(actor::ActorPFactive)

Generates output equilibrium with updated magnetic flux from optimized coil currents.

Updates the 2D equilibrium flux map (ψ) based on the computed coil currents, evaluates
current limits for the PF system, and optionally updates the target equilibrium if
requested. The output equilibrium shows what the coil system can actually achieve.
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
        fw = IMAS.first_wall(dd.wall)
        IMAS.flux_surfaces(eqt_in, fw.r, fw.z)
    end

    return actor
end

function update_eq_out(actor::ActorPFactive{D,P}) where {D<:Real,P<:Real}
    dd = actor.dd
    par = actor.par

    # evaluate eq_out
    eqt_in = dd.equilibrium.time_slice[]
    eqt2d_out = findfirst(:rectangular, actor.eqt_out.profiles_2d)
    coils = VacuumFields.IMAS_pf_active__coils(dd; par.green_model, zero_currents=false)
    EQfixed = IMAS2Equilibrium(eqt_in)
    Rgrid = eqt2d_out.grid.dim1
    Zgrid = eqt2d_out.grid.dim2
    eqt2d_out.psi = collect(VacuumFields.fixed2free(EQfixed, coils, Rgrid, Zgrid)')

    return actor
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