using FUSE
using Test

# ============================================================================
#  Following GSLite's equilibrium (evolve_equilibrium=false)
#
#  With FUSE_EVOLVE_EQ=0, receive! commits only GSLite's map, axis and boundary,
#  and follow_gslite_equilibrium! recomputes the rest of the slice on that map.
# ============================================================================

# Analytic map: a plasma whose ψ rises outward from near (1.70, 0), and a divertor coil below the
# wall at (1.40, -1.45). Between them sits an X-point. Softened so the grid never meets the coil.
const _C_COIL = 0.43
_psi_xpoint(r, z) = ((r - 1.70) / 0.62)^2 + (z / 1.05)^2 - 1.0 + _C_COIL * log((r - 1.40)^2 + (z + 1.45)^2 + 0.05^2)

# its X-point, by Newton on ∇ψ = 0
function _xpoint(; r=1.5, z=-1.0)
    for _ in 1:50
        d2 = (r - 1.40)^2 + (z + 1.45)^2 + 0.05^2
        g = [2(r - 1.70) / 0.62^2 + 2_C_COIL * (r - 1.40) / d2, 2z / 1.05^2 + 2_C_COIL * (z + 1.45) / d2]
        hrz = -4_C_COIL * (r - 1.40) * (z + 1.45) / d2^2
        H = [2/0.62^2+_C_COIL*(2/d2-4(r - 1.40)^2/d2^2) hrz; hrz 2/1.05^2+_C_COIL*(2/d2-4(z + 1.45)^2/d2^2)]
        r, z = [r, z] .- H \ g
    end
    return r, z
end

"""
The dd as receive! leaves it on an exchange: GSLite's map on its 33x33 grid, with the axis, boundary
psi and boundary outline found on it, in a slice whose 1-D profiles came from an earlier exchange.
core_profiles and pulse_schedule are FUSE's.
"""
function _followed_dd(; p_scale=1.0, j_scale=1.0, j_peaking=1.3)
    dd = IMAS.dd()
    resize!(dd.wall.description_2d, 1)
    limiter = resize!(dd.wall.description_2d[1].limiter.unit, 1)[1]
    θ = range(0, 2π, length=97)
    limiter.outline.r = 1.70 .+ 0.72 .* cos.(θ)
    limiter.outline.z = 1.25 .* sin.(θ)
    fw_r, fw_z = IMAS.first_wall(dd.wall)

    dd.global_time = 0.0
    dd.equilibrium.time = [0.0]
    eqt = resize!(dd.equilibrium.time_slice, 1)[1]
    eqt.time = 0.0
    eqt.global_quantities.vacuum_toroidal_field.r0 = 1.6955
    eqt.global_quantities.vacuum_toroidal_field.b0 = 2.0
    p2d = resize!(eqt.profiles_2d, 1)[1]
    p2d.grid_type.index = 1

    # the slice of an earlier exchange
    p2d.grid.dim1 = collect(range(0.84, 2.54, length=129))
    p2d.grid.dim2 = collect(range(-1.6, 1.6, length=129))
    p2d.psi = [_psi_xpoint(r, z) for r in p2d.grid.dim1, z in p2d.grid.dim2]
    r, z, itp = IMAS.ψ_interpolant(p2d)
    ar = IMAS.find_magnetic_axis(r, z, itp, 1.0)
    pb = IMAS.find_psi_boundary(r, z, p2d.psi, itp(ar.RA, ar.ZA), :increasing, ar.RA, ar.ZA, fw_r, fw_z;
        raise_error_on_not_open=false, raise_error_on_not_closed=false)
    eqt1d = eqt.profiles_1d
    eqt1d.psi = collect(range(itp(ar.RA, ar.ZA), pb.last_closed, length=101))
    psin = range(0.0, 1.0, 101)
    eqt1d.pressure = 2.0e4 .* (1.0 .- psin) .^ 2
    eqt1d.dpressure_dpsi = IMAS.gradient(eqt1d.psi, eqt1d.pressure)
    eqt1d.f_df_dpsi = -0.5 .* (1.0 .- psin)
    eqt1d.f = fill(2.0 * 1.6955, 101)
    IMAS.flux_surfaces(eqt, fw_r, fw_z)
    ip = eqt.global_quantities.ip

    # receive!'s commit of GSLite's map
    R = collect(range(0.84, 2.54, length=33))
    Z = collect(range(-1.6, 1.6, length=33))
    p2d.grid.dim1 = R
    p2d.grid.dim2 = Z
    p2d.psi = [_psi_xpoint(r, z) for r in R, z in Z]
    empty!(p2d, :phi)
    rgrid, zgrid = range(R[1], R[end], 33), range(Z[1], Z[end], 33)
    itp = FUSE.Interpolations.cubic_spline_interpolation((rgrid, zgrid), p2d.psi; extrapolation_bc=FUSE.Interpolations.Line())
    ar = IMAS.find_magnetic_axis(rgrid, zgrid, itp, 1.0)
    pb = IMAS.find_psi_boundary(rgrid, zgrid, p2d.psi, itp(ar.RA, ar.ZA), :increasing, ar.RA, ar.ZA, fw_r, fw_z;
        raise_error_on_not_open=false, raise_error_on_not_closed=false)
    surfaces = IMAS.trace_simple_surfaces([itp(ar.RA, ar.ZA), pb.last_closed], rgrid, zgrid, p2d.psi, itp, ar.RA, ar.ZA, fw_r, fw_z)
    eqt.boundary.outline.r = surfaces[end].r
    eqt.boundary.outline.z = surfaces[end].z
    eqt.global_quantities.magnetic_axis.r = ar.RA
    eqt.global_quantities.magnetic_axis.z = ar.ZA
    eqt.global_quantities.psi_axis = itp(ar.RA, ar.ZA)
    eqt.global_quantities.psi_boundary = pb.last_closed

    # FUSE's own profiles, and GSLite's Ip
    dd.core_profiles.time = [0.0]
    cp1d = resize!(dd.core_profiles.profiles_1d, 1)[1]
    cp1d.time = 0.0
    cp1d.grid.rho_tor_norm = collect(range(0.0, 1.0, 101))
    ρ = cp1d.grid.rho_tor_norm
    cp1d.pressure = p_scale .* 3.0e4 .* (1.0 .- ρ .^ 2) .^ 1.5
    cp1d.j_tor = j_scale .* 5.6e6 .* (1.0 .- ρ .^ 2) .^ j_peaking
    dd.pulse_schedule.flux_control.time = [-1.0, 1.0]
    dd.pulse_schedule.flux_control.i_plasma.reference = [1.1 * ip, 1.1 * ip]
    return dd
end

@testset "eq-follow traces GSLite's map on a finer grid" begin
    dd = _followed_dd()
    eqt = dd.equilibrium.time_slice[]
    xr, xz = _xpoint()
    to_xpoint(outline) = minimum(hypot.(outline.r .- xr, outline.z .- xz))
    @test to_xpoint(eqt.boundary.outline) > 0.03  # receive!'s boundary, on GSLite's grid
    # No fallback: on GSLite's grid flux_surfaces' first trial surface is smaller than a cell and
    # comes back empty, so the slice got a `nothing` boundary psi (the 2026-11-07 runs, from 4.07 s).
    @test_logs min_level = Logging.Warn FUSE.follow_gslite_equilibrium!(dd)
    @test size(eqt.profiles_2d[1].psi) == (FUSE.GSLITE_FOLLOW_GRID, FUSE.GSLITE_FOLLOW_GRID)
    @test to_xpoint(eqt.boundary.outline) < 0.01
end

@testset "eq-follow carries GSLite's Ip" begin
    dd = _followed_dd()
    FUSE.follow_gslite_equilibrium!(dd)
    @test dd.equilibrium.time_slice[].global_quantities.ip ≈ IMAS.get_from(dd, Val(:ip), :pulse_schedule) rtol = 1e-3
end

@testset "eq-follow sends FUSE's own betap and li" begin
    function followed(; kw...)
        dd = _followed_dd(; kw...)
        FUSE.follow_gslite_equilibrium!(dd)
        return dd.equilibrium.time_slice[].global_quantities
    end
    ref = followed()
    # betap = <p> / (μ0 Ip² / 2L²) and li = L² ∫Bp² dV / (V μ0² Ip²), with FUSE's pressure, current and Ip
    @test followed(; p_scale=2.0).beta_pol ≈ 2 * ref.beta_pol rtol = 1e-6
    twice = followed(; j_scale=2.0)
    @test twice.li_1 ≈ ref.li_1 rtol = 1e-6
    @test twice.beta_pol ≈ ref.beta_pol / 4 rtol = 1e-6
    @test followed(; j_peaking=3.0).li_1 > ref.li_1 + 0.1
end

@testset "eq-follow with FUSE's profiles equal to the map's gives IMAS's betap and li" begin
    dd = _followed_dd()
    FUSE.follow_gslite_equilibrium!(dd)
    eqt = dd.equilibrium.time_slice[]
    fw_r, fw_z = IMAS.first_wall(dd.wall)
    # FUSE's current = the map's enclosed current, and FUSE's pressure = the slice's
    i_map = [FUSE.trapz(s.ll, s.Bp) / IMAS.mks.μ_0 for s in IMAS.trace_surfaces(eqt, fw_r, fw_z)]
    cp1d = dd.core_profiles.profiles_1d[]
    ρ = eqt.profiles_1d.rho_tor_norm
    cp1d.j_tor = IMAS.gradient(cp1d.grid.area, IMAS.interp1d(ρ, i_map, :cubic).(cp1d.grid.rho_tor_norm))
    cp1d.pressure = IMAS.interp1d(ρ, eqt.profiles_1d.pressure, :cubic).(cp1d.grid.rho_tor_norm)
    FUSE.follow_gslite_equilibrium!(dd)
    betap, li = eqt.global_quantities.beta_pol, eqt.global_quantities.li_1
    IMAS.flux_surfaces(eqt, fw_r, fw_z)  # IMAS's own, from the map's Bp
    # this map isn't a Grad-Shafranov solution, so its enclosed current isn't the slice's Ip: compare
    # li Ip² and betap Ip², which are ∫Bp² dV and <p> up to the same constants
    gq = eqt.global_quantities
    @test betap * i_map[end]^2 ≈ gq.beta_pol * gq.ip^2 rtol = 1e-2
    @test li * i_map[end]^2 ≈ gq.li_1 * gq.ip^2 rtol = 1e-2
end

@testset "eq-follow keeps the slice when it can't follow" begin
    dd = _followed_dd()
    eqt = dd.equilibrium.time_slice[]
    before = deepcopy(eqt)
    empty!(dd.pulse_schedule.flux_control.i_plasma, :reference)
    @test_logs (:warn, r"^eq-follow: failed") FUSE.follow_gslite_equilibrium!(dd)
    @test eqt.profiles_2d[1].psi == before.profiles_2d[1].psi
    @test eqt.boundary.outline.z == before.boundary.outline.z
end

@testset "eq-follow logs the line the image build checks" begin
    dd = _followed_dd()
    @test_logs (:info, r"^eq-follow: betap=\S+ li=\S+ ") FUSE.follow_gslite_equilibrium!(dd; verbose=true)
end
