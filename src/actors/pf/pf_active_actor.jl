import MXHEquilibrium
import Optim
using LinearAlgebra

#= =============== =#
#  ActorPFcoilsOpt  #
#= =============== =#
options_green_model = [
    :point => "one filament per coil",
    :simple => "like :point, but OH coils have three filaments",
    :corners => "like :simple, but PF coils have filaments at the four corners",
    :realistic => "possibly hundreds of filaments per coil (very slow!)"
]

options_optimization_scheme = [
    :none => "Do not optimize",
    :currents => "Find optimial coil currents but do not change coil positions",
    :rail => "Find optimial coil positions"
]

Base.@kwdef mutable struct FUSEparameters__ActorPFcoilsOpt{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    green_model::Switch{Symbol} = Switch{Symbol}(options_green_model, "-", "Model used for the coils Green function calculations"; default=:simple)
    symmetric::Entry{Bool} = Entry{Bool}("-", "Force PF coils location to be up-down symmetric"; default=true)
    weight_currents::Entry{T} = Entry{T}("-", "Weight of current limit constraint"; default=2.0)
    weight_strike::Entry{T} = Entry{T}("-", "Weight given to matching the strike-points"; default=0.1)
    weight_lcfs::Entry{T} = Entry{T}("-", "Weight given to matching last closed flux surface"; default=1.0)
    weight_null::Entry{T} = Entry{T}("-", "Weight given to get field null for plasma breakdown"; default=1.0)
    maxiter::Entry{Int} = Entry{Int}("-", "Maximum number of optimizer iterations"; default=1000)
    optimization_scheme::Switch{Symbol} = Switch{Symbol}(options_optimization_scheme, "-", "Type of PF coil optimization to carry out"; default=:rail)
    update_equilibrium::Entry{Bool} = Entry{Bool}("-", "Overwrite target equilibrium with the one that the coils can actually make"; default=false)
    do_plot::Entry{Bool} = Entry{Bool}("-", "Plot"; default=false)
    verbose::Entry{Bool} = Entry{Bool}("-", "Verbose"; default=false)
end

Base.@kwdef mutable struct PFcoilsOptTrace
    params::Vector{Vector{Float64}} = Vector{Float64}[]
    cost_lcfs::Vector{Float64} = Float64[]
    cost_currents::Vector{Float64} = Float64[]
    cost_oh::Vector{Float64} = Float64[]
    cost_1to1::Vector{Float64} = Float64[]
    cost_spacing::Vector{Float64} = Float64[]
    cost_total::Vector{Float64} = Float64[]
end

function Base.empty!(trace::PFcoilsOptTrace)
    for field in fieldnames(typeof(trace))
        empty!(getfield(trace, field))
    end
end

mutable struct ActorPFcoilsOpt{D,P} <: ReactorAbstractActor
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorPFcoilsOpt{P}
    eq_in::IMAS.equilibrium{D}
    eq_out::IMAS.equilibrium{D}
    λ_regularize::Float64
    trace::PFcoilsOptTrace
end

"""
    ActorPFcoilsOpt(dd::IMAS.dd, act::ParametersAllActors; kw...)

Finds the optimal coil currents and locations of the poloidal field coils
to match the equilibrium boundary shape and obtain a field-null region at plasma start-up.

!!! note

    Manupulates data in `dd.pf_active`
"""
function ActorPFcoilsOpt(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorPFcoilsOpt(dd, act.ActorPFcoilsOpt; kw...)

    par = actor.par

    if par.optimization_scheme == :none
        if par.do_plot
            plot(actor.eq_in; cx=true)
            plot!(dd.build)
            display(plot!(actor.pf_active))
        end

    else
        old_update_equilibrium = par.update_equilibrium
        par.update_equilibrium = false

        if par.optimization_scheme == :currents
            # find coil currents
            finalize(step(actor))

        elseif par.optimization_scheme == :rail
            # optimize coil location and currents
            finalize(step(actor))

            if par.do_plot
                display(plot(actor.trace, :cost; title="Evolution of cost"))
                display(plot(actor.trace, :params; title="Evolution of optimized parameters"))
            end
        end

        if par.do_plot
            # final time slice
            time_index = length(dd.equilibrium.time)
            display(plot(dd.pf_active, :currents; time=dd.equilibrium.time[time_index], title="Current limits at t=$(dd.equilibrium.time[time_index]) s"))
            display(plot(actor; equilibrium=true, time_index))
            # field null time slice
            if par.weight_null > 0.0
                display(plot(actor; equilibrium=true, rail=true, time_index=1))
            end
        end

        par.update_equilibrium = old_update_equilibrium
        finalize(actor)
    end

    return actor
end

function ActorPFcoilsOpt(dd::IMAS.dd, par::FUSEparameters__ActorPFcoilsOpt; kw...)
    logging_actor_init(ActorPFcoilsOpt)
    par = par(kw...)

    λ_regularize = 1E-3
    trace = PFcoilsOptTrace()

    # reset pf coil rails
    n_coils = [rail.coils_number for rail in dd.build.pf_active.rail]
    init_pf_active!(dd.pf_active, dd.build, n_coils)

    return ActorPFcoilsOpt(dd, par, dd.equilibrium, dd.equilibrium, λ_regularize, trace)
end

"""
    _step(actor::ActorPFcoilsOpt)

Optimize coil currents and positions to produce sets of equilibria while minimizing coil currents
"""
function _step(actor::ActorPFcoilsOpt{T}) where {T<:Real}
    dd = actor.dd
    par = actor.par

    # reset pf_active to work on a possibly updated equilibrium
    empty!(actor.trace)
    actor.eq_in = dd.equilibrium
    actor.eq_out = IMAS.lazycopy(dd.equilibrium)

    # Get coils (as GS_IMAS_pf_active__coil) organized by their function and initialize them
    fixed_coils, pinned_coils, optim_coils = fixed_pinned_optim_coils(actor, par.optimization_scheme)
    coils = vcat(pinned_coils, optim_coils, fixed_coils)
    for coil in coils
        coil.current_time = actor.eq_in.time
        coil.current_data = zeros(T, size(actor.eq_in.time))
    end

    # run rail type optimizer
    if par.optimization_scheme in (:rail, :currents)
        (λ_regularize, trace) = optimize_coils_rail(
            actor.eq_in, dd, pinned_coils, optim_coils, fixed_coils;
            par.symmetric,
            actor.λ_regularize,
            par.weight_lcfs,
            par.weight_null,
            par.weight_currents,
            par.weight_strike,
            par.maxiter,
            par.verbose)
    else
        error("Supported ActorPFcoilsOpt optimization_scheme are `:currents` or `:rail`")
    end
    actor.λ_regularize = λ_regularize
    actor.trace = trace

    # transfer the results to IMAS.pf_active
    for coil in coils
        transfer_info_GS_coil_to_IMAS(dd.build, coil)
    end

    return actor
end

"""
    _finalize(actor::ActorPFcoilsOpt)

Update actor.eq_out 2D equilibrium PSI based on coils positions and currents
"""
function _finalize(actor::ActorPFcoilsOpt{D,P}) where {D<:Real,P<:Real}
    dd = actor.dd
    par = actor.par

    update_equilibrium = par.update_equilibrium

    coils = GS_IMAS_pf_active__coil{D,D}[]
    for (k, coil) in enumerate(dd.pf_active.coil)
        if k <= dd.build.pf_active.rail[1].coils_number
            coil_tech = dd.build.oh.technology
        else
            coil_tech = dd.build.pf_active.technology
        end
        push!(coils, GS_IMAS_pf_active__coil(coil, coil_tech, par.green_model))
    end

    # update equilibrium
    for time_index in eachindex(actor.eq_in.time_slice)
        if ismissing(actor.eq_in.time_slice[time_index].global_quantities, :ip)
            continue
        end
        for coil in coils
            coil.time_index = time_index
        end

        # convert equilibrium to MXHEquilibrium.jl format, since this is what VacuumFields uses
        EQfixed = IMAS2Equilibrium(actor.eq_in.time_slice[time_index])

        # update ψ map
        scale_eq_domain_size = 1.0
        R = range(EQfixed.r[1] / scale_eq_domain_size, EQfixed.r[end] * scale_eq_domain_size; length=length(EQfixed.r))
        Z = range(EQfixed.z[1] * scale_eq_domain_size, EQfixed.z[end] * scale_eq_domain_size; length=length(EQfixed.z))
        ψ_f2f = transpose(VacuumFields.fixed2free(EQfixed, coils, R, Z))
        actor.eq_out.time_slice[time_index].profiles_2d[1].grid.dim1 = R
        actor.eq_out.time_slice[time_index].profiles_2d[1].grid.dim2 = Z
        if false # hack to force up-down symmetric equilibrium
            actor.eq_out.time_slice[time_index].profiles_2d[1].psi = (ψ_f2f .+ ψ_f2f[1:end, end:-1:1]) ./ 2.0
        else
            actor.eq_out.time_slice[time_index].profiles_2d[1].psi = copy(ψ_f2f)
        end
    end

    # update psi
    if update_equilibrium
        for time_index in eachindex(actor.eq_out.time_slice)
            if !ismissing(actor.eq_out.time_slice[time_index].global_quantities, :ip)
                psi1 = actor.eq_out.time_slice[time_index].profiles_2d[1].psi
                actor.eq_in.time_slice[time_index].profiles_2d[1].psi = psi1
                IMAS.flux_surfaces(actor.eq_in.time_slice[time_index])
            end
        end
    end

    return actor
end

function pack_rail(bd::IMAS.build, λ_regularize::Float64, symmetric::Bool)::Vector{Float64}
    distances = []
    for rail in bd.pf_active.rail
        if rail.name !== "OH"
            # not symmetric
            if !symmetric
                coil_distances = collect(range(-1.0, 1.0; length=rail.coils_number + 2))[2:end-1]
                # even symmetric
            elseif mod(rail.coils_number, 2) == 0
                coil_distances = collect(range(-1.0, 1.0; length=rail.coils_number + 2))[2+Int(rail.coils_number // 2):end-1]
                # odd symmetric
            else
                coil_distances = collect(range(-1.0, 1.0; length=rail.coils_number + 2))[2+Int((rail.coils_number - 1) // 2)+1:end-1]
            end
            append!(distances, coil_distances)
        end
    end
    oh_height_off = []
    for rail in bd.pf_active.rail
        if rail.name == "OH"
            push!(oh_height_off, 1.0)
            if !symmetric
                push!(oh_height_off, 0.0)
            end
        end
    end
    packed = vcat(distances, oh_height_off, log10(λ_regularize))

    return packed
end

function unpack_rail!(packed::Vector, optim_coils::Vector, symmetric::Bool, bd::IMAS.build)
    λ_regularize = packed[end]
    if symmetric
        n_oh_params = 1
    else
        n_oh_params = 2
    end
    if any(rail.name == "OH" for rail in bd.pf_active.rail)
        oh_height_off = packed[end-n_oh_params:end-1]
        distances = packed[1:end-n_oh_params]
    else
        oh_height_off = []
        distances = packed[1:end-1]
    end

    if length(optim_coils) != 0 # optim_coils have zero length in case of the `static` optimization
        kcoil = 0
        koptim = 0
        koh = 0

        for rail in bd.pf_active.rail
            if rail.name == "OH"
                # mirror OH size when it reaches maximum extent of the rail
                oh_height_off[1] = mirror_bound(oh_height_off[1], -1.0, 1.0)
                if !symmetric
                    offset = oh_height_off[2]
                else
                    offset = 0.0
                end
                z_oh, height_oh = size_oh_coils(rail.outline.z, rail.coils_cleareance, rail.coils_number, oh_height_off[1], offset)
                for k in 1:rail.coils_number
                    koptim += 1
                    koh += 1
                    optim_coils[koptim].z = z_oh[koh]
                    optim_coils[koptim].height = height_oh
                end
            else
                r_interp = IMAS.interp1d(rail.outline.distance, rail.outline.r)
                z_interp = IMAS.interp1d(rail.outline.distance, rail.outline.z)
                # not symmetric
                if !symmetric
                    dkcoil = rail.coils_number
                    coil_distances = distances[kcoil+1:kcoil+dkcoil]
                    # even symmetric
                elseif mod(rail.coils_number, 2) == 0
                    dkcoil = Int(rail.coils_number // 2)
                    coil_distances = distances[kcoil+1:kcoil+dkcoil]
                    coil_distances = vcat(-reverse(coil_distances), coil_distances)
                    # odd symmetric
                else
                    dkcoil = Int((rail.coils_number - 1) // 2)
                    coil_distances = distances[kcoil+1:kcoil+dkcoil]
                    coil_distances = vcat(-reverse(coil_distances), 0.0, coil_distances)
                end
                kcoil += dkcoil

                # mirror coil position when they reach the end of the rail
                coil_distances = mirror_bound.(coil_distances, -1.0, 1.0)

                # get coils r and z from distances
                r_coils = r_interp.(coil_distances)
                z_coils = z_interp.(coil_distances)

                # assign to optim coils
                for k in 1:rail.coils_number
                    koptim += 1
                    optim_coils[koptim].r = r_coils[k]
                    optim_coils[koptim].z = z_coils[k]
                end
            end
        end
    end

    return 10^λ_regularize
end

function optimize_coils_rail(
    eq::IMAS.equilibrium,
    dd::IMAS.dd{D},
    pinned_coils::Vector{GS_IMAS_pf_active__coil{D,D}},
    optim_coils::Vector{GS_IMAS_pf_active__coil{D,D}},
    fixed_coils::Vector{GS_IMAS_pf_active__coil{D,D}};
    symmetric::Bool,
    λ_regularize::Real,
    weight_lcfs::Real,
    weight_null::Real,
    weight_currents::Real,
    weight_strike::Real,
    maxiter::Integer,
    verbose::Bool) where {D<:Real}

    bd = dd.build
    pc = dd.pulse_schedule.position_control

    fixed_eqs = []
    weights = []
    for time_index in eachindex(eq.time_slice)
        eqt = eq.time_slice[time_index]
        # field nulls
        if ismissing(eqt.global_quantities, :ip)
            # find ψp
            ψp_constant = eqt.global_quantities.psi_boundary
            rb, zb = IMAS.boundary(pc, 1)
            Bp_fac, ψp, Rp, Zp = VacuumFields.field_null_on_boundary(ψp_constant, rb, zb, fixed_coils)
            push!(fixed_eqs, (Bp_fac, ψp, Rp, Zp))
            push!(weights, Float64[])
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
            Bp_fac, ψp, Rp, Zp = VacuumFields.ψp_on_fixed_eq_boundary(fixed_eq, fixed_coils; Rx, Zx, fraction_inside=1.0 - 1E-6)
            push!(fixed_eqs, (Bp_fac, ψp, Rp, Zp))
            # weight more near the x-points
            h = (maximum(Zp) - minimum(Zp)) / 2.0
            o = (maximum(Zp) + minimum(Zp)) / 2.0
            weight = sqrt.(((Zp .- o) ./ h) .^ 2 .+ h) / h
            # give each strike point the same weight as the lcfs
            weight[end-length(Rx)+1:end] .= length(Rp) / (1 + length(Rx)) * weight_strike
            if all(weight .== 1.0)
                weight = Float64[]
            end
            push!(weights, weight)
        end
    end

    packed = pack_rail(bd, λ_regularize, symmetric)
    trace = PFcoilsOptTrace()

    oh_indexes = [coil.pf_active__coil.name == "OH" for coil in vcat(pinned_coils, optim_coils)]

    packed_tmp = [packed]
    function placement_cost(packed; do_trace=false)
        packed_tmp[1] = packed

        index = findall(.>(1.0), abs.(packed[1:end-1]))
        if length(index) > 0
            cost_1to1 = sum(abs.(packed[index]) .- 1.0)
        else
            cost_1to1 = 0.0
        end

        λ_regularize = unpack_rail!(packed, optim_coils, symmetric, bd)
        coils = vcat(pinned_coils, optim_coils)

        all_cost_lcfs = []
        all_cost_currents = []
        all_cost_oh = []
        for (time_index, (fixed_eq, weight)) in enumerate(zip(fixed_eqs, weights))
            for coil in vcat(pinned_coils, optim_coils, fixed_coils)
                coil.time_index = time_index
            end
            currents, cost_lcfs0 = VacuumFields.currents_to_match_ψp(fixed_eq..., coils; weights=weight, λ_regularize, return_cost=true)
            current_densities = currents .* [coil.turns_with_sign / area(coil) for coil in coils]
            oh_max_current_densities = [coil_J_B_crit(bd.oh.max_b_field, coil.coil_tech)[1] for coil in coils[oh_indexes]]
            pf_max_current_densities = [coil_J_B_crit(coil_selfB(coil), coil.coil_tech)[1] for coil in coils[oh_indexes.==false]]
            max_current_densities = vcat(oh_max_current_densities, pf_max_current_densities)
            fraction_max_current_densities = abs.(current_densities ./ max_current_densities)
            #currents cost
            push!(all_cost_currents, norm(exp.(fraction_max_current_densities * weight_currents) / exp(1)) / length(currents))
            # boundary and oh costs
            if ismissing(eq.time_slice[time_index].global_quantities, :ip)
                push!(all_cost_lcfs, cost_lcfs0 * weight_null)
                push!(all_cost_oh, 0.0)
            else
                push!(all_cost_lcfs, cost_lcfs0 * weight_lcfs)
                # OH cost (to avoid OH shrinking too much)
                index = findfirst(rail -> rail.name === "OH", bd.pf_active.rail)
                if index !== nothing
                    oh_rail_length = diff(collect(extrema(bd.pf_active.rail[index].outline.z)))[1]
                    total_oh_coils_length = sum(coil.height for coil in coils[oh_indexes.==true])
                    cost_oh = oh_rail_length / total_oh_coils_length
                else
                    cost_oh = 0.0
                end
                push!(all_cost_oh, cost_oh)
            end
        end
        cost_lcfs = norm(all_cost_lcfs) / length(all_cost_lcfs)
        cost_currents = norm(all_cost_currents) / length(all_cost_currents)
        cost_oh = norm(all_cost_oh) / length(all_cost_oh)
        # Spacing between PF coils
        cost_spacing = 0.0
        if length(optim_coils) > 0
            for (k1, c1) in enumerate(optim_coils)
                for (k2, c2) in enumerate(optim_coils)
                    if k1 < k2 && c1.pf_active__coil.name != "OH" && c2.pf_active__coil.name != "OH"
                        cost_spacing += exp(sqrt((c1.width + c2.width)^2 + (c1.height + c2.height)^2) / sqrt((c1.r - c2.r)^2 + (c1.z - c2.z)^2))
                    end
                end
            end
            cost_spacing = cost_spacing / length(optim_coils)^2
        end

        cost_lcfs_2 = cost_lcfs^2 * 10000.0
        cost_currents_2 = cost_currents^2
        cost_oh_2 = cost_oh^2
        cost_1to1_2 = cost_1to1^2
        cost_spacing_2 = cost_spacing^2

        # total cost
        cost = sqrt(cost_lcfs_2 + cost_currents_2 + cost_oh_2 + cost_1to1_2 + cost_spacing_2)
        if do_trace
            push!(trace.params, packed)
            push!(trace.cost_lcfs, sqrt(cost_lcfs_2))
            push!(trace.cost_currents, sqrt(cost_currents_2))
            push!(trace.cost_oh, sqrt(cost_oh_2))
            push!(trace.cost_1to1, sqrt(cost_1to1_2))
            push!(trace.cost_spacing, sqrt(cost_spacing_2))
            push!(trace.cost_total, cost)
        end
        if isnan(cost)
            display(plot([p[end] for p in trace.params]))
            if isnan(cost_lcfs)
                error("optimize_coils_rail cost_lcfs is NaN")
            end
            if isnan(cost_currents)
                error("optimize_coils_rail cost_currents is NaN")
            end
            if isnan(cost_oh)
                error("optimize_coils_rail cost_oh is NaN")
            end
            if isnan(cost_1to1)
                error("optimize_coils_rail cost_1to1 is NaN")
            end
            if isnan(cost_spacing)
                error("optimize_coils_rail cost_spacing is NaN")
            end
        end
        return cost

    end

    function clb(x)
        placement_cost(packed_tmp[1]; do_trace=true)
        return false
    end

    if maxiter == 0
        placement_cost(packed)
        λ_regularize = unpack_rail!(packed, optim_coils, symmetric, bd)
    else
        res = Optim.optimize(
            placement_cost,
            packed,
            Optim.NelderMead(),
            Optim.Options(; time_limit=60 * 2, iterations=maxiter, callback=clb, g_tol=1E-4, show_trace=verbose);
            autodiff=:forward
        )
        if verbose
            println(res)
        end
        packed = Optim.minimizer(res)
        λ_regularize = unpack_rail!(packed, optim_coils, symmetric, bd)
    end

    return λ_regularize, trace
end

"""
    fixed_pinned_optim_coils(actor::ActorPFcoilsOpt, optimization_scheme::Symbol)

Returns tuple of GS_IMAS_pf_active__coil coils organized by their function:

  - fixed: fixed position and current
  - pinned: coils with fixed position but current is optimized
  - optim: coils that have theri position and current optimized
"""
function fixed_pinned_optim_coils(actor::ActorPFcoilsOpt{D,P}, optimization_scheme::Symbol) where {D<:Real,P<:Real}
    dd = actor.dd
    par = actor.par

    fixed_coils = GS_IMAS_pf_active__coil{D,D}[]
    pinned_coils = GS_IMAS_pf_active__coil{D,D}[]
    optim_coils = GS_IMAS_pf_active__coil{D,D}[]
    for (k, coil) in enumerate(dd.pf_active.coil)
        if k <= dd.build.pf_active.rail[1].coils_number
            coil_tech = dd.build.oh.technology
        else
            coil_tech = dd.build.pf_active.technology
        end
        if coil.identifier == "pinned"
            push!(pinned_coils, GS_IMAS_pf_active__coil(coil, coil_tech, par.green_model))
        elseif (coil.identifier == "optim") && (optimization_scheme == :currents)
            push!(pinned_coils, GS_IMAS_pf_active__coil(coil, coil_tech, par.green_model))
        elseif coil.identifier == "optim"
            push!(optim_coils, GS_IMAS_pf_active__coil(coil, coil_tech, par.green_model))
        elseif coil.identifier == "fixed"
            push!(fixed_coils, GS_IMAS_pf_active__coil(coil, coil_tech, par.green_model))
        else
            error("Accepted type of coil.identifier are only \"optim\", \"pinned\", or \"fixed\"")
        end
    end
    return fixed_coils, pinned_coils, optim_coils
end

#= ======== =#
#  plotting  #
#= ======== =#
"""
    plot_pfcoilsactor_cx(actor::ActorPFcoilsOpt; time_index=1, equilibrium=true, rail=true)

Plot ActorPFcoilsOpt optimization cross-section
"""
@recipe function plot_ActorPFcoilsOpt_cx(
    actor::ActorPFcoilsOpt{D,P};
    time_index=nothing,
    equilibrium=true,
    build=true,
    coils_flux=false,
    rail=false,
    plot_r_buffer=1.6) where {D<:Real, P<:Real}

    @assert typeof(time_index) <: Union{Nothing,Integer}
    @assert typeof(equilibrium) <: Bool
    @assert typeof(build) <: Bool
    @assert typeof(coils_flux) <: Bool
    @assert typeof(rail) <: Bool
    @assert typeof(plot_r_buffer) <: Real

    dd = actor.dd
    par = actor.par

    if time_index === nothing
        time_index = length(actor.eq_out.time_slice)
    end
    time = actor.eq_out.time_slice[time_index].time

    # if there is no equilibrium then treat this as a field_null plot
    field_null = false
    if length(actor.eq_out.time_slice[time_index].profiles_2d) == 0 || ismissing(actor.eq_out.time_slice[time_index].profiles_2d[1], :psi)
        coils_flux = equilibrium
        field_null = true
    end

    # when plotting coils_flux the build is not visible anyways
    if coils_flux
        build = false
    end

    # setup plotting area
    xlim = [0.0, maximum(dd.build.layer[end].outline.r)]
    ylim = [minimum(dd.build.layer[end].outline.z), maximum(dd.build.layer[end].outline.z)]
    xlim --> xlim * plot_r_buffer
    ylim --> ylim
    aspect_ratio --> :equal

    # plot build
    if build
        @series begin
            exclude_layers --> [:oh]
            dd.build
        end
    end

    # plot coils_flux
    if coils_flux
        ngrid = 129
        R = range(xlim[1], xlim[2]; length=ngrid)
        Z = range(ylim[1], ylim[2]; length=Int(ceil(ngrid * (ylim[2] - ylim[1]) / (xlim[2] - xlim[1]))))

        coils = GS_IMAS_pf_active__coil{D,D}[]
        for (k, coil) in enumerate(dd.pf_active.coil)
            if k <= dd.build.pf_active.rail[1].coils_number
                coil_tech = dd.build.oh.technology
            else
                coil_tech = dd.build.pf_active.technology
            end
            coil = GS_IMAS_pf_active__coil(coil, coil_tech, par.green_model)
            coil.time_index = time_index
            push!(coils, coil)
        end

        # ψ coil currents
        ψbound = actor.eq_out.time_slice[time_index].global_quantities.psi_boundary
        ψ = VacuumFields.coils_flux(2 * pi, coils, R, Z)

        ψmin = minimum(x -> isnan(x) ? Inf : x, ψ)
        ψmax = maximum(x -> isnan(x) ? -Inf : x, ψ)
        ψabsmax = maximum(x -> isnan(x) ? -Inf : x, abs.(ψ))

        if field_null
            clims = (-ψabsmax / 10 + ψbound, ψabsmax / 10 + ψbound)
        else
            clims = (ψmin, ψmax)
        end

        @series begin
            seriestype --> :contourf
            c --> :diverging
            colorbar_entry --> false
            levels --> range(clims[1], clims[2]; length=21)
            linewidth --> 0.0
            R, Z, transpose(ψ)
        end

        if field_null
            @series begin
                seriestype --> :contour
                colorbar_entry --> false
                levels --> [ψbound]
                linecolor --> :gray
                R, Z, transpose(ψ)
            end
        end

        @series begin
            wireframe --> true
            exclude_layers --> [:oh]
            dd.build
        end
    end

    # plot equilibrium
    if equilibrium
        if field_null
            pc = dd.pulse_schedule.position_control
            @series begin
                cx := true
                label --> "Field null region"
                seriescolor --> :red
                IMAS.boundary(pc, 1)
            end
        else
            @series begin
                cx := true
                label --> "Final"
                seriescolor --> :red
                actor.eq_out.time_slice[time_index]
            end
            @series begin
                cx := true
                label --> "Target"
                seriescolor --> :blue
                lcfs --> true
                linestyle --> :dash
                actor.eq_in.time_slice[time_index]
            end
        end
    end

    # plot pf_active coils
    @series begin
        time --> time
        dd.pf_active
    end

    # plot optimization rails
    if rail
        @series begin
            label --> (build ? "Coil opt. rail" : "")
            dd.build.pf_active.rail
        end
    end

end

"""
    plot_pfcoilsactor_trace(trace::PFcoilsOptTrace, what::Symbol=:cost; start_at::Int=1)

Plot ActorPFcoilsOpt optimization trace

Attributes:

  - what::Symbol=:cost or :currents or individual fields of the PFcoilsOptTrace structure
  - start_at=::Int=1 index of the first element of the trace to start plotting
"""
@recipe function plot_ActorPFcoilsOpt_trace(
    trace::PFcoilsOptTrace,
    what::Symbol=:cost;
    start_at=1)

    @assert typeof(start_at) <: Integer

    start_at = minimum([start_at, length(trace.cost_total)])
    x = start_at:length(trace.cost_total)
    legend --> :bottomleft
    if what == :cost
        if sum(trace.cost_lcfs[start_at:end]) > 0.0
            data = trace.cost_lcfs[start_at:end]
            index = data .> 0.0
            @series begin
                label --> "ψ"
                yscale --> :log10
                x[index], data[index]
            end
        end
        if sum(trace.cost_currents[start_at:end]) > 0.0
            data = trace.cost_currents[start_at:end]
            index = data .> 0.0
            @series begin
                label --> "currents"
                yscale --> :log10
                x[index], data[index]
            end
        end
        if sum(trace.cost_oh[start_at:end]) > 0.0
            data = trace.cost_oh[start_at:end]
            index = data .> 0.0
            @series begin
                label --> "oh"
                yscale --> :log10
                x[index], data[index]
            end
        end
        if sum(trace.cost_1to1[start_at:end]) > 0.0
            data = trace.cost_1to1[start_at:end]
            index = data .> 0.0
            @series begin
                label --> "1to1"
                yscale --> :log10
                x[index], data[index]
            end
        end
        if sum(trace.cost_spacing[start_at:end]) > 0.0
            data = trace.cost_spacing[start_at:end]
            index = data .> 0.0
            @series begin
                label --> "spacing"
                yscale --> :log10
                x[index], data[index]
            end
        end
        @series begin
            label --> "total"
            yscale --> :log10
            linestyle --> :dash
            color --> :black
            # ylim --> [minimum(trace.cost_total[start_at:end]) / 10,maximum(trace.cost_total[start_at:end])]
            x, trace.cost_total[start_at:end]
        end

    elseif what == :params
        nparams = length(getfield(trace, what)[1]) - 1

        for k in 1:nparams
            @series begin
                label --> "#$k"
                x, [getfield(trace, what)[i][k] for i in eachindex(trace.cost_total)][start_at:end]
            end
        end

    else
        @series begin
            if occursin("cost_", String(what))
                yscale --> :log10
            end
            label --> String(what)
            x, getfield(trace, what)[start_at:end]
        end
    end
end
