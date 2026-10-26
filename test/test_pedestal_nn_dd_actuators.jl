#= ============================================================ =#
#  build_fuse29_actuators(source=:dd) on a synthetic dd: coil     #
#  currents as stored (per turn = PTDATA), the convention guard,  #
#  the beam-geometry torque fallback, gas_injection valves and    #
#  the loader repair of omas-divided F-coils. Needs the ONNX      #
#  bundle (FUSE_PEDESTAL_NN_DIR); skipped otherwise.              #
#    cd FUSE && julia --project test/test_pedestal_nn_dd_actuators.jl
#= ============================================================ =#
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using Test
using FUSE
import IMAS

have_bundle = withenv(FUSE.PEDESTAL_NN_AUTOFETCH_ENV => "0") do
    FUSE.pedestal_nn_dir_complete(FUSE.resolve_pedestal_nn_dir())
end
if !have_bundle
    @warn "fuse29 bundle not on disk; skipping test_pedestal_nn_dd_actuators"
else
    nn = FUSE.load_pedestal_nn()
    t = collect(0.0:0.01:1.0)

    function synthetic_dd(; divided::Bool=false) # divided = F-coils as the omas D3D mapping ≥ v0.96.0 exports them
        dd = IMAS.dd()
        dd.global_time = 0.5
        dd.equilibrium.time = [0.0, 1.0]
        dd.equilibrium.vacuum_toroidal_field.b0 = [-2.0, -2.0]
        # coils: F6A (one 55-turn element) at PTDATA -3000 A, F1A (58 turns), ECOILA (48 single-turn elements) at -7000 A,
        # F3A (58 turns) at PTDATA -200 A: below the 500 A threshold on its own, converted with the group
        for (name, turns, nel, ptdata) in (("F6A", 55.0, 1, -3000.0), ("F1A", 58.0, 1, -1500.0), ("ECOILA", 1.0, 48, -7000.0), ("F3A", 58.0, 1, -200.0))
            coil = resize!(dd.pf_active.coil, length(dd.pf_active.coil) + 1)[end]
            coil.name = name
            for el in resize!(coil.element, nel)
                el.turns_with_sign = turns
            end
            coil.current.time = t
            coil.current.data = fill(divided ? ptdata / turns : ptdata, length(t))
        end
        # one co-injected 80 keV D beam, 2 MW, R_tang = 1.15 m
        unit = resize!(dd.nbi.unit, 1)[1]
        unit.name = "30L"
        unit.power_launched.time = t
        unit.power_launched.data = fill(2.0e6, length(t))
        unit.energy.time = t
        unit.energy.data = fill(80.0e3, length(t))
        unit.species.a = 2.0
        bg = resize!(unit.beamlets_group, 1)[1]
        bg.tangency_radius = 1.15
        bg.direction = 1
        # gas valve A at 10 Torr·L/s
        valve = resize!(dd.gas_injection.valve, 1)[1]
        valve.name = "GASA"
        valve.flow_rate.time = t
        valve.flow_rate.data = fill(10.0 * FUSE.PA_M3_PER_TORR_L, length(t))
        return dd
    end

    @testset "fuse29 :dd actuators read the coil currents as stored" begin
        dd = synthetic_dd()
        u, live = FUSE.build_fuse29_actuators(nn, dd; source=:dd, time=0.5)
        idx(name) = findfirst(==(name), nn.actuator_names)
        @test live[idx("f6a")] && isapprox(u[idx("f6a")], -3000.0; rtol=1e-5)
        @test live[idx("f1a")] && isapprox(u[idx("f1a")], -1500.0; rtol=1e-5)
        @test live[idx("f3a")] && isapprox(u[idx("f3a")], -200.0; rtol=1e-5)
        @test live[idx("ecoila")] && isapprox(u[idx("ecoila")], -7000.0; rtol=1e-5)
        @test live[idx("bt")] && u[idx("bt")] == -2.0f0
        @test live[idx("gasa_cal")] && isapprox(u[idx("gasa_cal")], 10.0; rtol=1e-4)
        @test !live[idx("gasb_cal")]
        # torque from geometry: P R / v, v = sqrt(2 E e / m_D)
        v = sqrt(2 * 80.0e3 * IMAS.mks.e / (2.0 * IMAS.mks.m_u))
        @test live[idx("tinj")] && isapprox(u[idx("tinj")], 2.0e6 * 1.15 / v; rtol=1e-5)
        @test !live[idx("pinj")] # no pulse_schedule in this dd
        @test isempty(FUSE.fuse29_check_units(nn, u))
        @test getfield(dd, :_aux)[:fuse29_coil_convention_checked]
    end

    @testset "convention guard" begin
        dd = synthetic_dd(; divided=true) # F-coils at PTDATA / turns
        aux = getfield(dd, :_aux)
        msg = @test_logs (:warn, r"looks divided by turns_with_sign") FUSE.fuse29_coil_convention_check(nn, dd.pf_active.coil, aux)
        @test occursin("divided", msg)
        @test FUSE.fuse29_coil_convention_check(nn, dd.pf_active.coil, aux) === nothing # once per dd
        dd2 = synthetic_dd()
        for coil in dd2.pf_active.coil
            coil.current.data .*= coil.element[1].turns_with_sign # amp-turns stored as current
        end
        @test occursin("multiplied", @test_logs (:warn, r"multiplied") FUSE.fuse29_coil_convention_check(nn, dd2.pf_active.coil, getfield(dd2, :_aux)))
        dd3 = synthetic_dd()
        @test FUSE.fuse29_coil_convention_check(nn, dd3.pf_active.coil, getfield(dd3, :_aux)) === nothing
    end

    @testset "loader repair of omas-divided F-coils" begin
        dd = synthetic_dd(; divided=true)
        @test_logs (:warn, r"divided by turns_with_sign") FUSE.repair_pf_active_omas_turns!(dd.pf_active)
        @test isapprox(dd.pf_active.coil[1].current.data[1], -3000.0; rtol=1e-6)
        @test dd.pf_active.coil[3].current.data[1] == -7000.0 # 48 single-turn elements: untouched
        @test isapprox(dd.pf_active.coil[4].current.data[1], -200.0; rtol=1e-6) # low-current F3A repaired with the group
        ddp = synthetic_dd()
        FUSE.repair_pf_active_omas_turns!(ddp.pf_active) # correct dd: no change, no warning
        @test isapprox(ddp.pf_active.coil[1].current.data[1], -3000.0; rtol=1e-6)
        @test isapprox(ddp.pf_active.coil[4].current.data[1], -200.0; rtol=1e-6)
    end

    @testset "model torque used when core_sources are current" begin
        dd = synthetic_dd()
        src = resize!(dd.core_sources.source, 1)[1]
        src.identifier.index = 2
        gq = resize!(src.global_quantities, 1)[1]
        gq.time = 0.5
        gq.torque_tor = 1.234
        u, live = FUSE.build_fuse29_actuators(nn, dd; source=:dd, time=0.5)
        @test isapprox(u[findfirst(==("tinj"), nn.actuator_names)], 1.234; rtol=1e-6)
        u2, _ = FUSE.build_fuse29_actuators(nn, dd; source=:dd, time=0.9) # 0.4 s away: fall back to geometry
        @test !isapprox(u2[findfirst(==("tinj"), nn.actuator_names)], 1.234; rtol=1e-3)
    end
end
