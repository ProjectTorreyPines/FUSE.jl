import NeoclassicalTransport
import GACODE

#= ================================================= =#
#  Split of the pedestal pressure into Te_ped/Ti_ped  #
#= ================================================= =#
# EPED-like models only constrain the total pedestal pressure: how that pressure is shared between
# electrons and ions at the pedestal top is selected by the `Ti_ped_model` parameter of the pedestal actors:
#
#   :ratio       Ti_ped = r * Te_ped                                      (r from `T_ratio_pedestal`)
#   :increment   Ti_ped - Ti_sep = r * (Te_ped - Te_sep)
#   :heat_flux   pedestal power balance, ΔT = Q Δr / (e n χ S), with χi/χe = `chi_i_over_chi_e`
#   :fingerprint pedestal power balance with χi = χi_neo + f χe, where f is the turbulent χi/χe
#                of the dominant pedestal mode (`fingerprint_mode`)
#
# For EPED-like models the electron χ is not a free parameter: it is the one that makes the
# separatrix + pedestal temperature increments match the predicted pedestal pressure.

"""
    PEDESTAL_FINGERPRINTS

Approximate transport "fingerprints" of the pedestal modes: turbulent χi/χe and expected De/χe.
From Kotschenreuther et al., Nucl. Fusion 59, 096001 (2019), Table 1:

  - MHD-like (KBM):  χi/χe = 1,     De/χe = 2/3
  - MTM:             χi/χe ≈ 1/10,  De/χe ≈ 1/10
  - ETG:             χi/χe ≈ 1/10,  De/χe ≈ 1/20
  - ITG/TEM:         χe/χi ≈ 1/4–1 (χi/χe = 2 is used), De/(χi+χe) ≈ -0.1 ± 0.33 (no De/χe check)

The χi/χe ratios are applied to the turbulent part of χi only, the neoclassical χi is added separately.
"""
const PEDESTAL_FINGERPRINTS = Dict{Symbol,@NamedTuple{chi_i_over_chi_e::Float64, De_over_chi_e::Float64}}(
    :MHD => (chi_i_over_chi_e=1.0, De_over_chi_e=2.0 / 3.0),
    :MTM => (chi_i_over_chi_e=0.1, De_over_chi_e=0.1),
    :ETG => (chi_i_over_chi_e=0.1, De_over_chi_e=0.05),
    :ITG_TEM => (chi_i_over_chi_e=2.0, De_over_chi_e=NaN)
)

"""
    pedestal_ti_split(dd::IMAS.DD, par, rho_top::Real, Ti_over_Te::Real; warn::Bool=true)

Collects what is needed to split the pedestal pressure between electrons and ions according to `par.Ti_ped_model`.

`rho_top` is the pedestal top in rho_tor_norm and `Ti_over_Te` the ratio from `T_ratio_pedestal`.
For `:heat_flux` and `:fingerprint` the power crossing the pedestal top, densities and flux-surface area
at mid-pedestal, the pedestal radial extent (in minor radius) and, if requested, the neoclassical ion χ are evaluated.

Returns a NamedTuple to be used with `pedestal_temperatures` (pressure closure) or `pedestal_ti` (Te_ped given).
"""
function pedestal_ti_split(dd::IMAS.DD, par, rho_top::Real, Ti_over_Te::Real; warn::Bool=true)
    model = par.Ti_ped_model
    Te_sep = Float64(par.Te_sep)
    Ti_sep = Float64(separatrix_ti(par, Ti_over_Te))

    f = NaN
    chi_i_neo = 0.0
    ae = NaN
    ai = NaN
    De = NaN
    De_over_chi_e_expected = NaN
    Qe = NaN
    Qi = NaN

    if model ∈ (:heat_flux, :fingerprint)
        if model == :heat_flux
            f = Float64(par.chi_i_over_chi_e)
            @assert f > 0.0 "Ti_ped_model=:heat_flux requires chi_i_over_chi_e > 0"
        elseif par.fingerprint_mode == :custom
            f = Float64(par.chi_i_over_chi_e)
        else
            fp = PEDESTAL_FINGERPRINTS[par.fingerprint_mode]
            f = fp.chi_i_over_chi_e
            De_over_chi_e_expected = fp.De_over_chi_e
        end

        cp1d = dd.core_profiles.profiles_1d[]
        eqt = dd.equilibrium.time_slice[]
        rho = cp1d.grid.rho_tor_norm
        rho_mid = (1.0 + rho_top) / 2.0

        # pedestal radial extent in minor radius and flux-surface area at mid-pedestal
        eq1d = eqt.profiles_1d
        rmin = IMAS.interp1d(eq1d.rho_tor_norm, 0.5 .* (eq1d.r_outboard .- eq1d.r_inboard))
        dr = rmin(1.0) - rmin(rho_top)
        surface = IMAS.interp1d(rho, cp1d.grid.surface)(rho_mid)
        ne_mid = IMAS.interp1d(rho, cp1d.electrons.density_thermal)(rho_mid)
        ni_mid = IMAS.interp1d(rho, sum(ion.density_thermal for ion in cp1d.ion if IMAS.hasdata(ion, :density_thermal)))(rho_mid)

        # power crossing the pedestal top
        total_source = IMAS.total_sources(dd.core_sources, cp1d; time0=dd.global_time, fields=[:power_inside, :total_ion_power_inside, :particles_inside])
        Qe_raw = IMAS.interp1d(rho, total_source.electrons.power_inside)(rho_top)
        Qi_raw = IMAS.interp1d(rho, total_source.total_ion_power_inside)(rho_top)
        Qtot = Qe_raw + Qi_raw
        @assert Qtot > 0.0 "Ti_ped_model=:$(model) requires positive power crossing the pedestal (Qe + Qi = $Qtot W)"
        # NOTE: each channel is floored to 1% of the total, otherwise its pedestal temperature increment vanishes
        Qe = max(Qe_raw, 0.01 * Qtot)
        Qi = max(Qi_raw, 0.01 * Qtot)
        ae = Qe * dr / (IMAS.mks.e * ne_mid * surface)
        ai = Qi * dr / (IMAS.mks.e * ni_mid * surface)

        if model == :fingerprint && par.chi_i_neo_model == :changhinton
            chi_i_neo = chi_i_neo_changhinton(eqt, cp1d, rho_mid; warn)
        end
        @assert f > 0.0 || chi_i_neo > 0.0 "Ti_ped_model=:fingerprint requires either a turbulent χi/χe > 0 or a neoclassical χi (chi_i_neo_model=:changhinton)"

        # particle diffusivity at mid-pedestal (only used as a fingerprint consistency check)
        Γe = IMAS.interp1d(rho, total_source.electrons.particles_inside)(rho_top)
        Δne = IMAS.interp1d(rho, cp1d.electrons.density_thermal)(rho_top) - cp1d.electrons.density_thermal[end]
        if Γe > 0.0 && Δne > 0.0
            De = Γe * dr / (Δne * surface)
        end
    end

    return (; model, Te_sep, Ti_sep, Ti_over_Te=Float64(Ti_over_Te), rho_top=Float64(rho_top), f, chi_i_neo, ae, ai, Qe, Qi, De, De_over_chi_e_expected)
end

"""
    chi_i_neo_changhinton(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, rho::Real; warn::Bool=true)

Neoclassical ion heat diffusivity [m²/s] from the Chang-Hinton model at `rho` (on the current profiles),
defined as `q_i = e nᵢ χᵢ dTᵢ/dr` with the same flux convention (power = flux * surface) used by ActorNeoclassical.

For hydrogenic main ions the Chang-Hinton flux is linear in dTᵢ/dr, so χᵢ does not depend on the gradient.
If χᵢ cannot be evaluated (e.g. flat Tᵢ) it returns 0.0, with a warning if `warn=true`.
"""
function chi_i_neo_changhinton(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, rho::Real; warn::Bool=true)
    rho_cp = cp1d.grid.rho_tor_norm
    i_cp = IMAS.argmin_abs(rho_cp, rho)
    i_eq = IMAS.argmin_abs(eqt.profiles_1d.rho_tor_norm, rho)

    sol = NeoclassicalTransport.changhinton(eqt, cp1d, rho, 1)

    ne = cp1d.electrons.density_thermal[i_cp]
    Te = cp1d.electrons.temperature[i_cp]
    rhos = GACODE.rho_s(cp1d, eqt)[i_cp]
    a = eqt.boundary.minor_radius
    q_i = GACODE.gyrobohm_energy_flux(ne, Te, rhos, a) * sol.ENERGY_FLUX_i * GACODE.volume_prime_miller_correction(eqt)[i_eq]

    eq1d = eqt.profiles_1d
    rmin = IMAS.interp1d(eq1d.rho_tor_norm, 0.5 .* (eq1d.r_outboard .- eq1d.r_inboard)).(rho_cp)
    dTidr = -IMAS.gradient(rmin, cp1d.ion[1].temperature)[i_cp]

    # NOTE: Chang-Hinton is formulated with ni = ne
    chi = q_i / (IMAS.mks.e * ne * dTidr)
    if !isfinite(chi) || chi < 0.0
        warn && @warn "Chang-Hinton neoclassical χi at rho=$(round(rho; digits=3)) is not valid (χi=$chi, dTi/dr=$dTidr): χi_neo=0 assumed"
        return 0.0
    end
    return chi
end

"""
    pedestal_chi_e(Pstar::Real, ne::Real, ni::Real, ae::Real, ai::Real, chi_i_neo::Real, f::Real)

Electron χ such that the pedestal temperature increments `ΔTe = ae/χe` and `ΔTi = ai/(χi_neo + f χe)`
add up to the pressure increment `Pstar = ne ΔTe + ni ΔTi` [eV m⁻³], i.e. the positive root of

    Pstar f χe² + (Pstar χi_neo - ne ae f - ni ai) χe - ne ae χi_neo = 0

Returns `NaN` if there is no positive solution.
"""
function pedestal_chi_e(Pstar::Real, ne::Real, ni::Real, ae::Real, ai::Real, chi_i_neo::Real, f::Real)
    if Pstar <= 0.0
        return NaN
    end
    A = Pstar * f
    B = Pstar * chi_i_neo - ne * ae * f - ni * ai
    C = -ne * ae * chi_i_neo
    if A > 0.0
        # A > 0 and C <= 0: there is always one and only one positive root
        sq = sqrt(B^2 - 4.0 * A * C)
        return B < 0.0 ? (-B + sq) / (2.0 * A) : 2.0 * C / (-B - sq)
    elseif B > 0.0
        # f = 0: only neoclassical ion transport, solution exists if ions alone do not exceed Pstar
        return -C / B
    else
        return NaN
    end
end

"""
    pedestal_temperatures(split::NamedTuple, p_ped::Real, ne_ped::Real, ni_ped::Real)

Pedestal-top Te and Ti [eV] such that `e (ne_ped Te_ped + ni_ped Ti_ped) = p_ped` [Pa], split according to `split.model`.

Returns `(Te_ped, Ti_ped, chi_e, chi_i)`, with the diffusivities only defined for `:heat_flux` and `:fingerprint`.
"""
function pedestal_temperatures(split::NamedTuple, p_ped::Real, ne_ped::Real, ni_ped::Real)
    p = p_ped / IMAS.mks.e # eV m⁻³
    Pstar = p - ne_ped * split.Te_sep - ni_ped * split.Ti_sep
    chi_e = NaN
    chi_i = NaN

    if split.model == :ratio
        Te_ped = p / (ne_ped + ni_ped * split.Ti_over_Te)
        Ti_ped = Te_ped * split.Ti_over_Te

    elseif split.model == :increment
        ΔTe = Pstar / (ne_ped + ni_ped * split.Ti_over_Te)
        Te_ped = split.Te_sep + ΔTe
        Ti_ped = split.Ti_sep + split.Ti_over_Te * ΔTe

    elseif split.model ∈ (:heat_flux, :fingerprint)
        chi_e = pedestal_chi_e(Pstar, ne_ped, ni_ped, split.ae, split.ai, split.chi_i_neo, split.f)
        chi_i = split.chi_i_neo + split.f * chi_e
        Te_ped = split.Te_sep + split.ae / chi_e
        Ti_ped = split.Ti_sep + split.ai / chi_i

    else
        error("Ti_ped_model can only be :ratio, :increment, :heat_flux, or :fingerprint")
    end

    if Pstar <= 0.0 || !isfinite(Te_ped) || !isfinite(Ti_ped)
        # pedestal pressure does not exceed the separatrix temperatures at the pedestal-top density
        # (or the ions alone would exceed it): keep the separatrix Ti/Te ratio
        if split.model != :ratio
            @warn "Ti_ped_model=:$(split.model): no valid pedestal Te/Ti split (pressure above separatrix = $(round(Pstar * IMAS.mks.e)) Pa), separatrix Ti/Te assumed"
        end
        ratio = split.Ti_sep / split.Te_sep
        Te_ped = p / (ne_ped + ni_ped * ratio)
        Ti_ped = Te_ped * ratio
        chi_e = NaN
        chi_i = NaN
    end

    return (; Te_ped, Ti_ped, chi_e, chi_i)
end

"""
    pedestal_ti(split::NamedTuple, Te_ped::Real)

Pedestal-top Ti [eV] for a given pedestal-top Te [eV] (for models that set Te_ped directly, like WPED).

For `:heat_flux` and `:fingerprint` the electron χ is the one that carries `split.Qe` at the given `Te_ped - Te_sep`.

Returns `(Ti_ped, chi_e, chi_i)`
"""
function pedestal_ti(split::NamedTuple, Te_ped::Real)
    ΔTe = Te_ped - split.Te_sep
    if split.model == :ratio
        return (Ti_ped=Te_ped * split.Ti_over_Te, chi_e=NaN, chi_i=NaN)
    elseif ΔTe <= 0.0
        # no pedestal: keep the separatrix Ti/Te ratio
        return (Ti_ped=Te_ped * split.Ti_sep / split.Te_sep, chi_e=NaN, chi_i=NaN)
    elseif split.model == :increment
        return (Ti_ped=split.Ti_sep + split.Ti_over_Te * ΔTe, chi_e=NaN, chi_i=NaN)
    elseif split.model ∈ (:heat_flux, :fingerprint)
        chi_e = split.ae / ΔTe
        chi_i = split.chi_i_neo + split.f * chi_e
        return (Ti_ped=split.Ti_sep + split.ai / chi_i, chi_e=chi_e, chi_i=chi_i)
    else
        error("Ti_ped_model can only be :ratio, :increment, :heat_flux, or :fingerprint")
    end
end

"""
    pedestal_ti_split_result(split::NamedTuple, temperatures::NamedTuple, Ti_ped_history::Vector{Float64}=[temperatures.Ti_ped])

Merges the split inputs with the resulting pedestal temperatures, the inferred De/χe (fingerprint consistency check),
and the trial Ti_ped of each evaluation of the self-consistent exchange/χi_neo solve (see `Ti_ped_self_consistent`)
"""
function pedestal_ti_split_result(split::NamedTuple, temperatures::NamedTuple, Ti_ped_history::Vector{Float64}=[Float64(temperatures.Ti_ped)])
    De_over_chi_e = split.De / temperatures.chi_e
    return merge(split, temperatures, (; De_over_chi_e, Ti_ped_history=copy(Ti_ped_history)))
end

"""
    solve_pedestal_ti(residual::Function, T_lo::Real, T_hi::Real; atol::Real=0.5, maxiter::Int=60)

Bisection for the root of `residual(T)` in `[T_lo, T_hi]`, where `residual(T)` returns `Ti_new(T) - T` [eV].
Returns `(T, converged, history)`, with `history` the trial T of each evaluation;
`converged=false` if the residual has the same sign at both ends of the bracket.
"""
function solve_pedestal_ti(residual::Function, T_lo::Real, T_hi::Real; atol::Real=0.5, maxiter::Int=60)
    history = Float64[]
    r(T) = (push!(history, T); residual(T))
    a, b = Float64(T_lo), Float64(T_hi)
    ra, rb = r(a), r(b)
    if ra == 0.0
        return (T=a, converged=true, history)
    elseif rb == 0.0
        return (T=b, converged=true, history)
    elseif sign(ra) == sign(rb)
        return (T=abs(ra) < abs(rb) ? a : b, converged=false, history)
    end
    for _ in 1:maxiter
        m = (a + b) / 2.0
        if b - a < 2.0 * atol
            return (T=m, converged=true, history)
        end
        rm = r(m)
        if rm == 0.0
            return (T=m, converged=true, history)
        elseif sign(rm) == sign(ra)
            a, ra = m, rm
        else
            b, rb = m, rm
        end
    end
    return (T=(a + b) / 2.0, converged=false, history)
end

"""
    save_temperatures(cp1d::IMAS.core_profiles__profiles_1d)

Copy of the electron and ion temperature profiles, to be used with `restore_temperatures!`
"""
function save_temperatures(cp1d::IMAS.core_profiles__profiles_1d)
    Te = copy(cp1d.electrons.temperature)
    Ti = [ismissing(ion, :temperature) ? nothing : copy(ion.temperature) for ion in cp1d.ion]
    return (; Te, Ti)
end

"""
    restore_temperatures!(cp1d::IMAS.core_profiles__profiles_1d, saved::NamedTuple)

Restores the temperature profiles saved with `save_temperatures`
"""
function restore_temperatures!(cp1d::IMAS.core_profiles__profiles_1d, saved::NamedTuple)
    cp1d.electrons.temperature = copy(saved.Te)
    for (ion, Ti) in zip(cp1d.ion, saved.Ti)
        if Ti !== nothing
            ion.temperature = copy(Ti)
        end
    end
    return cp1d
end
