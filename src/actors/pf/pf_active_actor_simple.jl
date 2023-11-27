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
    return ActorPFactive(dd, par, dd.equilibrium, dd.equilibrium, 0.0)
end

"""
    _step(actor::ActorPFactive)

Optimize coil currents to produce equilibrium at current time
"""
function _step(actor::ActorPFactive{T}) where {T<:Real}
    dd = actor.dd

    actor.eq_in = dd.equilibrium
    actor.eq_out = IMAS.lazycopy(dd.equilibrium)

    # Get coils (as GS3_IMAS_pf_active__coil) organized by their function and initialize them
    fixed_coils, pinned_coils, optim_coils = fixed_pinned_optim_coils(actor)

    weight_strike = 1.0

    # run rail type optimizer
    λ_regularize = find_currents(
        actor.eq_in, dd, pinned_coils, optim_coils, fixed_coils;
        actor.λ_regularize,
        weight_strike)
    actor.λ_regularize = λ_regularize

    return actor
end

"""
    _finalize(actor::ActorPFactive)

Update actor.eq_out 2D equilibrium PSI based on coils positions and currents
"""
function _finalize(actor::ActorPFactive{D,P}) where {D<:Real,P<:Real}
    dd = actor.dd
    par = actor.par

    if par.update_equilibrium || par.do_plot

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
            eqt2d_out.psi = VacuumFields.fixed2free(EQfixed, coils, Rgrid, Zgrid)
        end
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

    return actor
end

function find_currents(
    eq::IMAS.equilibrium,
    dd::IMAS.dd{D},
    pinned_coils::Vector{GS3_IMAS_pf_active__coil{D,D}},
    optim_coils::Vector{GS3_IMAS_pf_active__coil{D,D}},
    fixed_coils::Vector{GS3_IMAS_pf_active__coil{D,D}};
    λ_regularize::Real,
    weight_strike::Real) where {D<:Real}

    bd = dd.build
    pc = dd.pulse_schedule.position_control

    eqt = eq.time_slice[]
    # field nulls
    if ismissing(eqt.global_quantities, :ip)
        # find ψp
        ψp_constant = eqt.global_quantities.psi_boundary
        rb, zb = IMAS.boundary(pc, 1)
        Bp_fac, ψp, Rp, Zp = VacuumFields.field_null_on_boundary(ψp_constant, rb, zb, fixed_coils)
        # solutions with plasma
    else
        fixed_eq = IMAS2Equilibrium(eqt)
        # private flux regions
        Rx = Float64[]
        Zx = Float64[]
        if weight_strike > 0.0
            private = IMAS.flux_surface(eqt, eqt.profiles_1d.psi[end], false)
            vessel = IMAS.get_build_layer(bd.layer; type=_plasma_)
            for (pr, pz) in private
                indexes, crossings = IMAS.intersection(vessel.outline.r, vessel.outline.z, pr, pz)
                for cr in crossings
                    push!(Rx, cr[1])
                    push!(Zx, cr[2])
                end
            end
            if isempty(Rx)
                @warn "weight_strike>0 but no strike point found"
            end
        end
        # find ψp
        Bp_fac, ψp, Rp, Zp = VacuumFields.ψp_on_fixed_eq_boundary(fixed_eq, fixed_coils; Rx, Zx, fraction_inside=1.0 - 1E-4)
        # weight more near the x-points
        h = (maximum(Zp) - minimum(Zp)) / 2.0
        o = (maximum(Zp) + minimum(Zp)) / 2.0
        weight = sqrt.(((Zp .- o) ./ h) .^ 2 .+ h) / h # asda
        # give each strike point the same weight as the lcfs
        weight[end-length(Rx)+1:end] .= length(Rp) / (1 + length(Rx)) * weight_strike
        if all(weight .== 1.0)
            weight = Float64[]
        end
    end

    coils = vcat(pinned_coils, optim_coils)
    currents = VacuumFields.currents_to_match_ψp(Bp_fac, ψp, Rp, Zp, coils; weights=weight, λ_regularize=0.0, return_cost=false)

    return λ_regularize
end

function fixed_pinned_optim_coils(actor::ActorPFactive{D,P}) where {D<:Real,P<:Real}
    return fixed_pinned_optim_coils(actor, :currents, GS3_IMAS_pf_active__coil)
end