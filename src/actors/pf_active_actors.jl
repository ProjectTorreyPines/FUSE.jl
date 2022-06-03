using Equilibrium
import Optim
import VacuumFields
using LinearAlgebra

#= =============== =#
#  ActorPFcoilsOpt  #
#= =============== =#
Base.@kwdef mutable struct PFcoilsOptTrace
    params::Vector{Vector{Real}} = Vector{Real}[]
    cost_ψ::Vector{Real} = Real[]
    cost_currents::Vector{Real} = Real[]
    cost_total::Vector{Real} = Real[]
end

mutable struct ActorPFcoilsOpt <: AbstractActor
    eq_in::IMAS.equilibrium
    eq_out::IMAS.equilibrium
    pf_active::IMAS.pf_active
    bd::IMAS.build
    symmetric::Bool
    λ_regularize::Real
    trace::PFcoilsOptTrace
    green_model::Symbol
end

function ParametersActor(::Type{Val{:ActorPFcoilsOpt}})
    par = ParametersActor(nothing)
    options = [
        :point => "one filament per coil",
        :simple => "like :point, but OH coils have three filaments",
        :corners => "like :simple, but PF coils have filaments at the four corners",
        :realistic => "hundreds of filaments per coil (very slow!)",
    ]
    par.green_model = Switch(options, "", "Model used for the Greens function calculation"; default=:simple)
    par.symmetric = Entry(Bool, "", "PF coils location should be up-down symmetric"; default=true)
    par.λ_currents = Entry(Real, "", "Weight of current limit constraint"; default=0.5)
    par.λ_strike = Entry(Real, "", "Weight given to matching the strike-points"; default=0.0)
    par.λ_ψ = Entry(Real, "", "Weight given to matching last closed flux surface"; default=1.0)
    par.λ_null = Entry(Real, "", "Weight given to get field null for plasma breakdown"; default=1E-3)
    options = [
        :none => "Do not optimize",
        :currents => "Find optimial coil currents but do not change coil positions",
        :rail => "Find optimial coil positions"
    ]
    par.optimization_scheme = Switch(options, "", "Optimization to carry out"; default=:currents)
    par.update_equilibrium = Entry(Bool, "", "Overwrite target equilibrium with the one that the coils can actually make"; default=false)
    par.do_plot = Entry(Bool, "", "plot"; default=false)
    par.verbose = Entry(Bool, "", "verbose"; default=false)
    return par
end

"""
    ActorPFcoilsOpt(dd::IMAS.dd, act::ParametersActor; kw...)

This actor fiends the coil currents and locations of the poloidal field coils for both the field-null during start-up and to match the equilibrium boundary shape.

!!! note 
    Stores data in `dd.pf_active`
"""
function ActorPFcoilsOpt(dd::IMAS.dd, act::ParametersActor; kw...)
    par = act.ActorPFcoilsOpt(kw...)
    actor = ActorPFcoilsOpt(dd; par.green_model, par.symmetric)

    if par.optimization_scheme == :none
        if par.do_plot
            plot(actor.eq_in)
            plot!(actor.bd)
            display(plot!(actor.pf_active))
        end

    else
        if par.optimization_scheme == :currents
            # find coil currents
            step(actor; par.λ_ψ, par.λ_null, par.λ_currents, par.λ_strike, par.verbose, maxiter=1000, optimization_scheme=:currents)
            finalize(actor)

        elseif par.optimization_scheme == :rail
            # optimize coil location and currents
            step(actor; par.λ_ψ, par.λ_null, par.λ_currents, par.λ_strike, par.verbose, maxiter=1000, optimization_scheme=:rail)
            finalize(actor)

            if par.do_plot
                display(plot(actor.trace, :cost))
                display(plot(actor.trace, :params))
            end
        end

        if par.do_plot
            # field null time slice
            display(plot(actor.pf_active, :currents, time=dd.equilibrium.time[1]))
            display(plot(actor, equilibrium=true, rail=true, time_index=1))
            # final time slice
            display(plot(actor.pf_active, :currents, time=dd.equilibrium.time[end]))
            display(plot(actor, equilibrium=true, time_index=length(dd.equilibrium.time)))
        end

        finalize(actor; update_eq_in=par.update_equilibrium)
    end

    return actor
end

function ActorPFcoilsOpt(dd::IMAS.dd; kw...)
    return ActorPFcoilsOpt(dd.equilibrium, dd.build, dd.pf_active; kw...)
end

function ActorPFcoilsOpt(
    eq_in::IMAS.equilibrium,
    bd::IMAS.build,
    pf::IMAS.pf_active;
    λ_regularize=1E-3,
    green_model=:simple,
    symmetric=false)

    # reset pf coil rails
    n_coils = [rail.coils_number for rail in bd.pf_active.rail]
    init_pf_active(pf, bd, n_coils)

    # basic constructors
    eq_out = deepcopy(eq_in)

    # constructor
    pfactor = ActorPFcoilsOpt(eq_in, eq_out, pf, bd, symmetric, λ_regularize, PFcoilsOptTrace(), green_model)

    return pfactor
end

"""
    step(pfactor::ActorPFcoilsOpt;
        symmetric=pfactor.symmetric,
        λ_regularize=pfactor.λ_regularize,
        λ_ψ=1.0,
        λ_null=1E-3,
        λ_currents=0.5,
        λ_strike=0.0,
        maxiter=10000,
        optimization_scheme=:rail,
        verbose=false)

Optimize coil currents and positions to produce sets of equilibria while minimizing coil currents
"""
function step(pfactor::ActorPFcoilsOpt;
    symmetric=pfactor.symmetric,
    λ_regularize=pfactor.λ_regularize,
    λ_ψ=1.0,
    λ_null=1E-3,
    λ_currents=0.5,
    λ_strike=0.0,
    maxiter=10000,
    optimization_scheme=:rail,
    verbose=false)

    fixed_coils, pinned_coils, optim_coils = fixed_pinned_optim_coils(pfactor, optimization_scheme)
    coils = vcat(pinned_coils, optim_coils, fixed_coils)
    for coil in coils
        coil.current_time = pfactor.eq_in.time
        coil.current_data = zeros(size(pfactor.eq_in.time))
    end

    bd = pfactor.bd
    # run rail type optimizer
    if optimization_scheme in [:rail, :currents]
        (λ_regularize, trace) = optimize_coils_rail(pfactor.eq_in; pinned_coils, optim_coils, fixed_coils, symmetric, λ_regularize, λ_ψ, λ_null, λ_currents, λ_strike, bd, maxiter, verbose)
    else
        error("Supported ActorPFcoilsOpt optimization_scheme are `:currents` or `:rail`")
    end
    pfactor.λ_regularize = λ_regularize
    pfactor.trace = trace

    # transfer the results to IMAS.pf_active
    for coil in coils
        transfer_info_GS_coil_to_IMAS(bd, coil)
    end

    return pfactor
end

"""
    finalize(pfactor::ActorPFcoilsOpt; scale_eq_domain_size = 1.0)

Update pfactor.eq_out 2D equilibrium PSI based on coils positions and currents
"""
function finalize(pfactor::ActorPFcoilsOpt; scale_eq_domain_size=1.0, update_eq_in=false)
    coils = GS_IMAS_pf_active__coil[]
    for (k, coil) in enumerate(pfactor.pf_active.coil)
        if k <= pfactor.bd.pf_active.rail[1].coils_number
            coil_tech = pfactor.bd.oh.technology
        else
            coil_tech = pfactor.bd.pf_active.technology
        end
        push!(coils, GS_IMAS_pf_active__coil(coil, coil_tech, pfactor.green_model))
    end

    # update equilibrium
    for time_index in 1:length(pfactor.eq_in.time_slice)
        if ismissing(pfactor.eq_in.time_slice[time_index].global_quantities, :ip)
            continue
        end
        for coil in coils
            coil.time_index = time_index
        end

        # convert equilibrium to Equilibrium.jl format, since this is what VacuumFields uses
        EQfixed = IMAS2Equilibrium(pfactor.eq_in.time_slice[time_index])

        # # update ψ map
        R = range(EQfixed.r[1] / scale_eq_domain_size, EQfixed.r[end] * scale_eq_domain_size, length=length(EQfixed.r))
        Z = range(EQfixed.z[1] * scale_eq_domain_size, EQfixed.z[end] * scale_eq_domain_size, length=length(EQfixed.z))
        ψ_f2f = VacuumFields.fixed2free(EQfixed, coils, R, Z)
        pfactor.eq_out.time_slice[time_index].profiles_2d[1].grid.dim1 = R
        pfactor.eq_out.time_slice[time_index].profiles_2d[1].grid.dim2 = Z
        pfactor.eq_out.time_slice[time_index].profiles_2d[1].psi = transpose(ψ_f2f)
    end

    # update psi
    if update_eq_in
        for time_index in 1:length(pfactor.eq_out.time_slice)
            if !ismissing(pfactor.eq_out.time_slice[time_index].global_quantities, :ip)
                psi1 = pfactor.eq_out.time_slice[time_index].profiles_2d[1].psi
                pfactor.eq_in.time_slice[time_index].profiles_2d[1].psi = psi1
                IMAS.flux_surfaces(pfactor.eq_in.time_slice[time_index])
            end
        end
    end
end

#= ==================================== =#
#  IMAS.pf_active__coil to VacuumFields  #
#= ==================================== =#
mutable struct GS_IMAS_pf_active__coil <: VacuumFields.AbstractCoil
    pf_active__coil::IMAS.pf_active__coil
    r::Real
    z::Real
    width::Real
    height::Real
    turns_with_sign::Real
    spacing::Real
    coil_tech::Union{IMAS.build__oh__technology,IMAS.build__pf_active__technology}
    current_data::Vector{<:Real}
    current_time::Vector{<:Real}
    time_index::Int
    green_model::Symbol
end

function GS_IMAS_pf_active__coil(
    pf_active__coil::IMAS.pf_active__coil,
    coil_tech::Union{IMAS.build__oh__technology,IMAS.build__pf_active__technology},
    green_model::Symbol)
    return GS_IMAS_pf_active__coil(pf_active__coil,
        pf_active__coil.element[1].geometry.rectangle.r,
        pf_active__coil.element[1].geometry.rectangle.z,
        pf_active__coil.element[1].geometry.rectangle.width,
        pf_active__coil.element[1].geometry.rectangle.height,
        pf_active__coil.element[1].turns_with_sign,
        get_spacing_from_turns(pf_active__coil),
        coil_tech,
        pf_active__coil.current.data,
        pf_active__coil.current.time,
        1,
        green_model)
end

function Base.getproperty(coil::GS_IMAS_pf_active__coil, field::Symbol)
    if field == :current
        return getfield(coil, :current_data)[coil.time_index]
    else
        return getfield(coil, field)
    end
end

function Base.setproperty!(coil::GS_IMAS_pf_active__coil, field::Symbol, value)
    if field == :current
        getfield(coil, :current_data)[coil.time_index] = value
    else
        setfield!(coil, field, value)
    end
    if field in [:width, :height, :spacing]
        s = sign(getfield(coil, :turns_with_sign))
        turns = Int(ceil(coil.width .* coil.height ./ coil.spacing .^ 2))
        setfield!(coil, :turns_with_sign, s * turns)
    end
end

function transfer_info_GS_coil_to_IMAS(bd::IMAS.build, coil::GS_IMAS_pf_active__coil)
    pf_active__coil = coil.pf_active__coil
    pf_active__coil.element[1].geometry.rectangle.r = coil.r
    pf_active__coil.element[1].geometry.rectangle.z = coil.z
    pf_active__coil.element[1].geometry.rectangle.width = coil.width
    pf_active__coil.element[1].geometry.rectangle.height = coil.height
    pf_active__coil.element[1].turns_with_sign = coil.turns_with_sign
    pf_active__coil.b_field_max = range(0.1, 30, step=0.1)
    pf_active__coil.temperature = [-1, coil.coil_tech.temperature]
    pf_active__coil.current_limit_max = [abs(coil_J_B_crit(b, coil.coil_tech)[1] * area(coil) / coil.turns_with_sign) for b in pf_active__coil.b_field_max, t in pf_active__coil.temperature]
    pf_active__coil.b_field_max_timed.time = coil.current_time
    if pf_active__coil.name == "OH"
        pf_active__coil.b_field_max_timed.data = [bd.oh.max_b_field for time_index in 1:length(coil.current_time)]
    else
        pf_active__coil.b_field_max_timed.data = [coil_selfB(coil, time_index) for time_index in 1:length(coil.current_time)]
    end
    pf_active__coil.current.time = coil.current_time
    pf_active__coil.current.data = coil.current_data
end

function set_turns_from_spacing!(coil::GS_IMAS_pf_active__coil)
    pf_active__coil = getfield(coil, :pf_active__coil)
    return set_turns_from_spacing!(pf_active__coil, coil.spacing)
end

function set_turns_from_spacing!(pf_active__coil::IMAS.pf_active__coil, spacing::Real)
    s = sign(pf_active__coil.element[1].turns_with_sign)
    set_turns_from_spacing!(pf_active__coil, spacing, s)
end

function set_turns_from_spacing!(pf_active__coil::IMAS.pf_active__coil, spacing::Real, s::Int)
    pf_active__coil.element[1].turns_with_sign = s * Int(ceil(IMAS.area(pf_active__coil) / spacing^2))
end

function get_spacing_from_turns(coil::GS_IMAS_pf_active__coil)
    pf_active__coil = getfield(coil, :pf_active__coil)
    return get_spacing_from_turns(pf_active__coil)
end

function get_spacing_from_turns(pf_active__coil::IMAS.pf_active__coil)
    return sqrt((pf_active__coil.element[1].geometry.rectangle.width * pf_active__coil.element[1].geometry.rectangle.height) / abs(pf_active__coil.element[1].turns_with_sign))
end

function area(coil::GS_IMAS_pf_active__coil)
    return IMAS.area(coil.pf_active__coil)
end

"""
    VacuumFields.Green(coil::GS_IMAS_pf_active__coil, R::Real, Z::Real)

Calculates coil green function at given R and Z coordinate
"""
function VacuumFields.Green(coil::GS_IMAS_pf_active__coil, R::Real, Z::Real)
    if coil.green_model == :point # fastest
        return VacuumFields.Green(coil.r, coil.z, R, Z, coil.turns_with_sign)

    elseif coil.green_model in [:corners, :simple] # medium
        if coil.pf_active__coil.name == "OH"
            n = 3
            z_filaments = range(coil.z - (coil.height - coil.width / 2.0) / 2.0, coil.z + (coil.height - coil.width / 2.0) / 2.0, length=n)
            green = []
            for z in z_filaments
                push!(green, VacuumFields.Green(coil.r, z, R, Z, coil.turns_with_sign / n))
            end
            return sum(green)

        elseif coil.green_model == :corners
            return VacuumFields.Green(VacuumFields.ParallelogramCoil(coil.r, coil.z, coil.width / 2.0, coil.height / 2.0, 0.0, 90.0, nothing), R, Z, coil.turns_with_sign / 4)

        elseif coil.green_model == :simple
            return VacuumFields.Green(coil.r, coil.z, R, Z, coil.turns_with_sign)
        end

    elseif coil.green_model == :realistic # high-fidelity
        return VacuumFields.Green(VacuumFields.ParallelogramCoil(coil.r, coil.z, coil.width, coil.height, 0.0, 90.0, coil.spacing), R, Z)

    else
        error("GS_IMAS_pf_active__coil coil.green_model can only be (in order of accuracy) :realistic, :corners, :simple, and :point")
    end
end

# step
function pack_rail(bd::IMAS.build, λ_regularize::Float64, symmetric::Bool)::Vector{Float64}
    distances = []
    for rail in bd.pf_active.rail
        if rail.name !== "OH"
            # not symmetric
            if !symmetric
                coil_distances = collect(range(-1.0, 1.0, length=rail.coils_number + 2))[2:end-1]
                # even symmetric
            elseif mod(rail.coils_number, 2) == 0
                coil_distances = collect(range(-1.0, 1.0, length=rail.coils_number + 2))[2+Int(rail.coils_number // 2):end-1]
                # odd symmetric
            else
                coil_distances = collect(range(-1.0, 1.0, length=rail.coils_number + 2))[2+Int((rail.coils_number - 1) // 2)+1:end-1]
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

"""
    coil_selfB(coil::GS_IMAS_pf_active__coil, time_index::Int=0)
PF coil self-induced magnetic field
NOTE: infinite wire approximation
"""
function coil_selfB(coil::GS_IMAS_pf_active__coil, time_index::Int=0)
    if time_index == 0
        time_index = coil.time_index
    end
    b = abs.(constants.μ_0 * coil.current_data[time_index] * coil.turns_with_sign / (2pi * min(coil.width, coil.height)))
    if b < 0.1
        return 0.1
    else
        return b
    end
end

function optimize_coils_rail(
    eq::IMAS.equilibrium;
    pinned_coils::Vector{GS_IMAS_pf_active__coil},
    optim_coils::Vector{GS_IMAS_pf_active__coil},
    fixed_coils::Vector{GS_IMAS_pf_active__coil},
    symmetric::Bool,
    λ_regularize::Real,
    λ_ψ::Real,
    λ_null::Real,
    λ_currents::Real,
    λ_strike::Real,
    bd::IMAS.build,
    maxiter::Int,
    verbose::Bool)

    R0 = eq.vacuum_toroidal_field.r0
    fixed_eqs = []
    weights = []
    for time_index in 1:length(eq.time_slice)
        eqt = eq.time_slice[time_index]
        # field nulls
        if ismissing(eqt.global_quantities, :ip)
            # find ψp
            Bp_fac, ψp, Rp, Zp = VacuumFields.field_null_on_boundary(eqt.global_quantities.psi_boundary,
                eqt.boundary.outline.r,
                eqt.boundary.outline.z,
                fixed_coils)
            push!(fixed_eqs, (Bp_fac, ψp, Rp, Zp))
            push!(weights, Float64[])
            # solutions with plasma
        else
            fixed_eq = IMAS2Equilibrium(eqt)
            # private flux regions
            Rx = []
            Zx = []
            if λ_strike > 0.0
                private = IMAS.flux_surface(eqt, eqt.profiles_1d.psi[end], false)
                vessel = IMAS.get_build(bd, type=_plasma_)
                for (pr, pz) in private
                    pvx, pvy = IMAS.intersection(vessel.outline.r, vessel.outline.z, pr, pz; as_list_of_points=false)
                    append!(Rx, pvx)
                    append!(Zx, pvy)
                end
            end
            # find ψp
            Bp_fac, ψp, Rp, Zp = VacuumFields.ψp_on_fixed_eq_boundary(fixed_eq, fixed_coils; Rx, Zx)
            push!(fixed_eqs, (Bp_fac, ψp, Rp, Zp))
            # give each strike point the same weight as the lcfs
            weight = Rp .* 0.0 .+ 1.0
            weight[end-length(Rx)+1:end] .= length(Rp) / (1 + length(Rx)) * λ_strike
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
            cost_1to1 = sum(abs.(packed[index]) .- 1.0) * 10
        else
            cost_1to1 = 0.0
        end

        λ_regularize = unpack_rail!(packed, optim_coils, symmetric, bd)
        coils = vcat(pinned_coils, optim_coils)

        all_cost_ψ = []
        all_cost_currents = []
        all_cost_oh = []
        for (time_index, (fixed_eq, weight)) in enumerate(zip(fixed_eqs, weights))
            for coil in vcat(pinned_coils, optim_coils, fixed_coils)
                coil.time_index = time_index
            end
            currents, cost_ψ0 = VacuumFields.currents_to_match_ψp(fixed_eq..., coils; weights=weight, λ_regularize, return_cost=true)
            current_densities = currents .* [coil.turns_with_sign / area(coil) for coil in coils]
            oh_max_current_densities = [coil_J_B_crit(bd.oh.max_b_field, coil.coil_tech)[1] for coil in coils[oh_indexes]]
            pf_max_current_densities = [coil_J_B_crit(coil_selfB(coil), coil.coil_tech)[1] for coil in coils[oh_indexes.==false]]
            max_current_densities = vcat(oh_max_current_densities, pf_max_current_densities)
            fraction_max_current_densities = abs.(current_densities ./ max_current_densities)
            #currents cost
            push!(all_cost_currents, norm(exp.(fraction_max_current_densities / λ_currents) / exp(1)) / length(currents))
            # boundary cost
            if ismissing(eq.time_slice[time_index].global_quantities, :ip)
                push!(all_cost_ψ, cost_ψ0 / λ_null)
                push!(all_cost_oh, 0.0)
            else
                #OH cost
                oh_current_densities = current_densities[oh_indexes]
                avg_oh = Statistics.mean(oh_current_densities)
                cost_oh = norm(oh_current_densities .- avg_oh) / avg_oh
                push!(all_cost_ψ, cost_ψ0 / λ_ψ)
                push!(all_cost_oh, cost_oh)
            end
        end
        cost_ψ = norm(all_cost_ψ) / length(all_cost_ψ)
        cost_currents = norm(all_cost_currents) / length(all_cost_currents)
        cost_oh = norm(all_cost_oh) / length(all_cost_oh)
        #spacing
        cost_spacing = 0.0
        for (k1, c1) in enumerate(optim_coils)
            for (k2, c2) in enumerate(optim_coils)
                if k1 < k2
                    cost_spacing += 1.0 / sqrt((c1.r - c2.r)^2 + (c1.z - c2.z)^2)
                end
            end
        end
        cost_spacing = cost_spacing / R0
        # total cost
        cost = sqrt(cost_ψ^2 + cost_currents^2 + 0.1 * cost_oh^2 + cost_1to1^2 + 10 * cost_spacing^2)
        if do_trace
            push!(trace.params, packed)
            push!(trace.cost_ψ, cost_ψ)
            push!(trace.cost_currents, cost_currents)
            push!(trace.cost_total, cost)
        end
        if isnan(cost)
            error("optimize_coils_rail cost is NaN")
        end
        return cost

    end

    function clb(x)
        placement_cost(packed_tmp[1]; do_trace=true)
        false
    end

    if maxiter == 0
        placement_cost(packed)
        λ_regularize = unpack_rail!(packed, optim_coils, symmetric, bd)
    else
        res = Optim.optimize(placement_cost, packed, Optim.NelderMead(), Optim.Options(time_limit=60 * 2, iterations=maxiter, callback=clb, g_tol=1E-4); autodiff=:forward)
        if verbose
            println(res)
        end
        packed = Optim.minimizer(res)
        λ_regularize = unpack_rail!(packed, optim_coils, symmetric, bd)
    end

    return λ_regularize, trace
end

"""
    fixed_pinned_optim_coils(pfactor, optimization_scheme)

Returns tuple of GS_IMAS_pf_active__coil coils organized by their function:
- fixed: fixed position and current
- pinned: coils with fixed position but current is optimized
- optim: coils that have theri position and current optimized
"""
function fixed_pinned_optim_coils(pfactor, optimization_scheme)
    fixed_coils = GS_IMAS_pf_active__coil[]
    pinned_coils = GS_IMAS_pf_active__coil[]
    optim_coils = GS_IMAS_pf_active__coil[]
    for (k, coil) in enumerate(pfactor.pf_active.coil)
        if k <= pfactor.bd.pf_active.rail[1].coils_number
            coil_tech = pfactor.bd.oh.technology
        else
            coil_tech = pfactor.bd.pf_active.technology
        end
        if coil.identifier == "pinned"
            push!(pinned_coils, GS_IMAS_pf_active__coil(coil, coil_tech, pfactor.green_model))
        elseif (coil.identifier == "optim") && (optimization_scheme == :currents)
            push!(pinned_coils, GS_IMAS_pf_active__coil(coil, coil_tech, pfactor.green_model))
        elseif coil.identifier == "optim"
            push!(optim_coils, GS_IMAS_pf_active__coil(coil, coil_tech, pfactor.green_model))
        elseif coil.identifier == "fixed"
            push!(fixed_coils, GS_IMAS_pf_active__coil(coil, coil_tech, pfactor.green_model))
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
    plot_pfcoilsactor_cx(pfactor::ActorPFcoilsOpt; time_index=1, equilibrium=true, rail=true)

Plot ActorPFcoilsOpt optimization cross-section
"""
@recipe function plot_ActorPFcoilsOpt_cx(pfactor::ActorPFcoilsOpt; time_index=nothing, equilibrium=true, build=true, coils_flux=false, rail=false, plot_r_buffer=1.6)

    if time_index === nothing
        time_index = length(pfactor.eq_out.time_slice)
    end
    time = pfactor.eq_out.time_slice[time_index].time

    # if there is no equilibrium then treat this as a field_null plot
    field_null = false
    if length(pfactor.eq_out.time_slice[time_index].profiles_2d) == 0 || ismissing(pfactor.eq_out.time_slice[time_index].profiles_2d[1], :psi)
        coils_flux = equilibrium
        field_null = true
    end

    # when plotting coils_flux the build is not visible anyways
    if coils_flux
        build = false
    end

    # setup plotting area
    xlim = [0.0, maximum(pfactor.bd.layer[end].outline.r)]
    ylim = [minimum(pfactor.bd.layer[end].outline.z), maximum(pfactor.bd.layer[end].outline.z)]
    xlim --> xlim * plot_r_buffer
    ylim --> ylim
    aspect_ratio --> :equal

    # plot build
    if build
        @series begin
            exclude_layers --> [:oh]
            pfactor.bd
        end
    end

    # plot coils_flux
    if coils_flux
        ngrid = 129
        R = range(xlim[1], xlim[2], length=ngrid)
        Z = range(ylim[1], ylim[2], length=Int(ceil(ngrid * (ylim[2] - ylim[1]) / (xlim[2] - xlim[1]))))

        coils = GS_IMAS_pf_active__coil[]
        for (k, coil) in enumerate(pfactor.pf_active.coil)
            if k <= pfactor.bd.pf_active.rail[1].coils_number
                coil_tech = pfactor.bd.oh.technology
            else
                coil_tech = pfactor.bd.pf_active.technology
            end
            coil = GS_IMAS_pf_active__coil(coil, coil_tech, pfactor.green_model)
            coil.time_index = time_index
            push!(coils, coil)
        end

        # ψ coil currents
        ψbound = pfactor.eq_out.time_slice[time_index].global_quantities.psi_boundary
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
            levels --> range(clims[1], clims[2], length=21)
            linewidth --> 0.0
            R, Z, transpose(ψ)
        end

        if field_null
            @series begin
                seriestype --> :contour
                colorbar_entry --> false
                levels --> [ψbound]
                linecolor --> :black
                R, Z, transpose(ψ)
            end
        end

        @series begin
            wireframe --> true
            exclude_layers --> [:oh]
            pfactor.bd
        end
    end

    # plot equilibrium
    if equilibrium
        if field_null
            @series begin
                label --> "Field null region"
                seriescolor --> :red
                pfactor.eq_out.time_slice[time_index]
            end
        else
            @series begin
                label --> "Final"
                seriescolor --> :red
                pfactor.eq_out.time_slice[time_index]
            end
            @series begin
                label --> "Target"
                seriescolor --> :blue
                lcfs --> true
                linestyle --> :dash
                pfactor.eq_in.time_slice[time_index]
            end
        end
    end

    # plot pf_active coils
    @series begin
        time --> time
        pfactor.pf_active
    end

    # plot optimization rails
    if rail
        @series begin
            label --> (build ? "Coil opt. rail" : "")
            pfactor.bd.pf_active.rail
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
@recipe function plot_ActorPFcoilsOpt_trace(trace::PFcoilsOptTrace, what::Symbol=:cost; start_at=1)
    start_at = minimum([start_at, length(trace.cost_total)])
    x = start_at:length(trace.cost_total)
    legend --> :bottomleft
    if what == :cost
        if sum(trace.cost_ψ[start_at:end]) > 0.0
            @series begin
                label --> "ψ"
                yscale --> :log10
                x, trace.cost_ψ[start_at:end]
            end
        end
        if sum(trace.cost_currents[start_at:end]) > 0.0
            @series begin
                label --> "currents"
                yscale --> :log10
                x, trace.cost_currents[start_at:end]
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
                x, [getfield(trace, what)[i][k] for i in 1:length(trace.cost_total)][start_at:end]
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
