import MXHEquilibrium
import Optim
using LinearAlgebra

#= =============== =#
#  ActorPFcoilsOpt  #
#= =============== =#
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
    weight_currents::Entry{T} = Entry{T}("-", "Weight of current limit constraint"; default=2.0) # current limit tolerance
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

mutable struct ActorPFcoilsOpt{D,P} <: ReactorAbstractActor{D,P}
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
            plot!(dd.build; legend=false)
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
            display(plot(dd.pf_active, :currents; time0=dd.equilibrium.time[time_index], title="PF currents at at t=$(dd.equilibrium.time[time_index]) s"))
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
    init_pf_active!(dd.pf_active, dd.build, dd.equilibrium.time_slice[], n_coils)

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
    fixed_coils, pinned_coils, optim_coils = fixed_pinned_optim_coils(actor)
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
            par.verbose,
            do_trace=par.do_plot)
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

    coils = GS_IMAS_pf_active__coil{D,D}[]
    for coil in dd.pf_active.coil
        if IMAS.is_ohmic_coil(coil)
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
        R = range(EQfixed.r[1] / scale_eq_domain_size, EQfixed.r[end] * scale_eq_domain_size, length(EQfixed.r))
        Z = range(EQfixed.z[1] * scale_eq_domain_size, EQfixed.z[end] * scale_eq_domain_size, length(EQfixed.z))
        eqt2d_out = findfirst(:rectangular, actor.eq_out.time_slice[time_index].profiles_2d)
        eqt2d_out.grid.dim1 = R
        eqt2d_out.grid.dim2 = Z
        eqt2d_out.psi = collect(VacuumFields.fixed2free(EQfixed, coils, R, Z)')
    end

    # update psi
    if par.update_equilibrium
        for time_index in eachindex(actor.eq_out.time_slice)
            if !ismissing(actor.eq_out.time_slice[time_index].global_quantities, :ip)
                eqt2d_out = findfirst(:rectangular, actor.eq_out.time_slice[time_index].profiles_2d)
                eqt2d_in = findfirst(:rectangular, actor.eq_in.time_slice[time_index].profiles_2d)
                eqt2d_in.psi = eqt2d_out.psi
                IMAS.flux_surfaces(actor.eq_in.time_slice[time_index])
            end
        end
    end

    return actor
end

function pack_rail(bd::IMAS.build, λ_regularize::Float64, symmetric::Bool)
    distances = Float64[]
    lbounds = Float64[]
    ubounds = Float64[]
    for rail in bd.pf_active.rail
        if rail.name !== "OH"
            # not symmetric
            if !symmetric
                coil_distances = collect(range(-1.0, 1.0, rail.coils_number + 2))[2:end-1]
                # even symmetric
            elseif mod(rail.coils_number, 2) == 0
                coil_distances = collect(range(-1.0, 1.0, rail.coils_number + 2))[2+Int(rail.coils_number // 2):end-1]
                # odd symmetric
            else
                coil_distances = collect(range(-1.0, 1.0, rail.coils_number + 2))[2+Int((rail.coils_number - 1) // 2)+1:end-1]
            end
            append!(distances, coil_distances)
            append!(lbounds, coil_distances .* 0.0 .- 1.0)
            append!(ubounds, coil_distances .* 0.0 .+ 1.0)
        end
    end
    oh_height_off = Float64[]
    for rail in bd.pf_active.rail
        if rail.name == "OH"
            push!(oh_height_off, 1.0)
            push!(lbounds, 1.0 - 1.0 / rail.coils_number / 2.0)
            push!(ubounds, 1.0)
            if !symmetric
                push!(oh_height_off, 0.0)
                push!(lbounds, -1.0 / rail.coils_number / 2.0)
                push!(ubounds, 1.0 / rail.coils_number / 2.0)
            end
        end
    end
    packed = vcat(distances, oh_height_off, log10(λ_regularize))
    push!(lbounds, -25.0)
    push!(ubounds, -5.0)

    return packed, (lbounds, ubounds)
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
        oh_height_off = Float64[]
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
    verbose::Bool,
    do_trace::Bool) where {D<:Real}

    bd = dd.build
    pc = dd.pulse_schedule.position_control

    eq_cp_psibs = []
    for time_index in eachindex(eq.time_slice)
        eqt = eq.time_slice[time_index]
        # field nulls
        if ismissing(eqt.global_quantities, :ip)
            psib = eqt.global_quantities.psi_boundary
            rb, zb = IMAS.boundary(pc, 1)
            flux_cps = VacuumFields.FluxControlPoints(rb, zb, psib)
            push!(eq_cp_psibs, (nothing, flux_cps, psib))
            # solutions with plasma
        else
            fixed_eq = IMAS2Equilibrium(eqt)
            _, psib = MXHEquilibrium.psi_limits(fixed_eq)

            # private flux regions
            Rx = Float64[]
            Zx = Float64[]
            if weight_strike > 0.0
                private, _ = IMAS.flux_surface(eqt, eqt.profiles_1d.psi[end], :open)
                vessel = IMAS.get_build_layer(bd.layer; type=_plasma_)
                for (pr, pz) in private
                    _, crossings = IMAS.intersection(vessel.outline.r, vessel.outline.z, pr, pz)
                    for cr in crossings
                        push!(Rx, cr[1])
                        push!(Zx, cr[2])
                    end
                end
                if isempty(Rx)
                    @warn "weight_strike>0 but no strike point found"
                end
            end

            bnd_cps = VacuumFields.boundary_control_points(fixed_eq, 0.999)
            strike_cps = VacuumFields.FluxControlPoints(Rx, Zx, psib)

            # weight more near the x-points
            Zmax, Zmin = -Inf, Inf
            for cp in bnd_cps
                cp.Z > Zmax && (Zmax = cp.Z)
                cp.Z < Zmin && (Zmin = cp.Z)
            end
            h = 0.5 * (Zmax - Zmin)
            o = 0.5 * (Zmax + Zmin)
            for cp in bnd_cps
                cp.weight = sqrt(((cp.Z - o) / h) ^ 2 + h) / h
            end

            # give each strike point the same weight as the lcfs
            for cp in strike_cps
                cp.weight = length(bnd_cps) / (1 + length(strike_cps)) * weight_strike
            end
            flux_cps = vcat(bnd_cps, strike_cps)

            push!(eq_cp_psibs, (fixed_eq, flux_cps, psib))
        end
    end

    packed, bounds = pack_rail(bd, λ_regularize, symmetric)
    trace = PFcoilsOptTrace()

    oh_indexes = [IMAS.is_ohmic_coil(coil.pf_active__coil) for coil in vcat(pinned_coils, optim_coils)]

    # MAKE RUN WITH eq_cps

    function placement_cost(packed::Vector{Float64}; do_trace::Bool=false)
        λ_regularize = unpack_rail!(packed, optim_coils, symmetric, bd)
        coils = vcat(pinned_coils, optim_coils)

        all_cost_lcfs = 0.0
        all_const_currents = Float64[]
        for (time_index, (fixed_eq, flux_cps, psib)) in enumerate(eq_cp_psibs)
            for coil in vcat(pinned_coils, optim_coils, fixed_coils)
                coil.time_index = time_index
            end

            # find best currents to match boundary
            currents, cost_lcfs0 = VacuumFields.find_coil_currents!(coils, fixed_eq, flux_cps; ψbound=psib, fixed_coils, λ_regularize, return_cost=true)

            # currents cost
            current_densities = currents .* [coil.turns_with_sign / area(coil) for coil in coils]
            oh_max_current_densities = [coil_J_B_crit(bd.oh.max_b_field, coil.coil_tech)[1] for coil in coils[oh_indexes]] ./ weight_currents
            pf_max_current_densities = [coil_J_B_crit(coil_selfB(coil), coil.coil_tech)[1] for coil in coils[oh_indexes.==false]] ./ weight_currents
            max_current_densities = vcat(oh_max_current_densities, pf_max_current_densities)
            append!(all_const_currents, log10.(abs.(current_densities) ./ (1.0 .+ max_current_densities)))

            # boundary costs
            if ismissing(eq.time_slice[time_index].global_quantities, :ip)
                all_cost_lcfs += cost_lcfs0 * weight_null
            else
                all_cost_lcfs += cost_lcfs0 * weight_lcfs
            end
        end

        # Spacing between PF coils
        const_spacing = 0.0
        if length(optim_coils) > 0
            for (k1, c1) in enumerate(optim_coils)
                for (k2, c2) in enumerate(optim_coils)
                    if k1 < k2
                        d = sqrt((c1.r - c2.r)^2 + (c1.z - c2.z)^2)
                        s = sqrt((c1.width + c2.width)^2 + (c1.height + c2.height)^2)
                        if IMAS.is_ohmic_coil(c1.pf_active__coil) && IMAS.is_ohmic_coil(c2.pf_active__coil)
                        else
                            const_spacing = max(const_spacing, s - d)
                        end
                    end
                end
            end
        end

        cost = all_cost_lcfs
        if do_trace
            push!(trace.params, packed)
            push!(trace.cost_lcfs, all_cost_lcfs)
            push!(trace.cost_oh, 0.0)
            push!(trace.cost_currents, maximum(all_const_currents))
            push!(trace.cost_1to1, 0.0)
            push!(trace.cost_spacing, maximum(const_spacing))
            push!(trace.cost_total, cost)
        end

        return cost, vcat(all_const_currents, const_spacing), [0.0]
    end

    if maxiter == 0
        placement_cost(packed)
        λ_regularize = unpack_rail!(packed, optim_coils, symmetric, bd)

    else
        x_optimum = deepcopy(packed)
        function logger(st, x_optimum, x_tol=1E-3)
            α = 0.1
            x_optimum .= x_optimum .* (1.0 - α) .+ Metaheuristics.minimizer(st) .* α
            if norm(Metaheuristics.minimizer(st) .- x_optimum) < x_tol
                st.termination_status_code = Metaheuristics.OBJECTIVE_DIFFERENCE_LIMIT
                st.stop = true
            end
        end

        options = Metaheuristics.Options(; seed=1, iterations=1000)
        algorithm = Metaheuristics.ECA(; N=length(packed), options)
        res = Metaheuristics.optimize(x -> placement_cost(x; do_trace), bounds, algorithm; logger=st -> logger(st, x_optimum))
        packed = Metaheuristics.minimizer(res)
        if verbose
            println(res)
        end
        λ_regularize = unpack_rail!(packed, optim_coils, symmetric, bd)
    end

    return λ_regularize, trace
end

function fixed_pinned_optim_coils(actor::ActorPFcoilsOpt{D,P}) where {D<:Real,P<:Real}
    return fixed_pinned_optim_coils(actor, actor.par.optimization_scheme, GS_IMAS_pf_active__coil)
end