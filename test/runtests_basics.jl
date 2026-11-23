using FUSE
using Test

@testset "act_migrations" begin
    act = FUSE.ParametersActors()
    schema_keys = Set(string(k) for k in keys(act))
    for old_name in keys(FUSE.ACT_ACTOR_MIGRATIONS)
        @test old_name ∉ schema_keys
    end
end

@testset "xpoints_from_pulse_schedule" begin
    dd = FUSE.IMAS.dd()
    pc = dd.pulse_schedule.position_control
    @test FUSE.xpoints_from_pulse_schedule(pc, 0.0) == :none
    pc.time = [0.0, 1.0]
    resize!(pc.x_point, 2)
    pc.x_point[1].r.reference = [1.5, 1.5]
    pc.x_point[1].z.reference = [2.0, 2.1]
    pc.x_point[2].r.reference = [1.5, 1.5]
    pc.x_point[2].z.reference = [-2.0, -2.1]
    @test FUSE.xpoints_from_pulse_schedule(pc, 0.0) == :double
    @test FUSE.xpoints_from_pulse_schedule(pc, 3.0) == :lower
    @test FUSE.xpoints_from_pulse_schedule(pc, -3.0) == :upper
    # an x-point that is not active at the current time
    pc.x_point[1].z.reference = [NaN, 2.1]
    @test FUSE.xpoints_from_pulse_schedule(pc, 0.0) == :lower
    resize!(pc.x_point, 1)
    pc.x_point[1].r.reference = [1.5, 1.5]
    pc.x_point[1].z.reference = [2.0, 2.1]
    @test FUSE.xpoints_from_pulse_schedule(pc, 0.0) == :upper
end

@testset "GC" begin
    act = FUSE.ParametersActors()
    @test parent(act.ActorNeoclassical) === act
    GC.gc()
    @test parent(act.ActorNeoclassical) === act
end

@testset "apply_act_overrides!" begin
    act = FUSE.ParametersActors()
    act.ActorFRESCO.nR = 65
    act.ActorFRESCO.nZ = 65
    overrides = FUSE.unset_parameters!(FUSE.ParametersActors())
    overrides.ActorFRESCO.nR = 129 # a default value is forced too
    FUSE.apply_act_overrides!(act, overrides)
    @test act.ActorFRESCO.nR == 129
    @test act.ActorFRESCO.nZ == 65
end

@testset "optimization_parameters" begin
    ini = FUSE.ParametersInits()

    ini.core_profiles.zeff = 2.0 ↔ [1.2, 2.5]
    @test ini.core_profiles.zeff == 2.0
    @test typeof(ini.core_profiles.zeff) <: Float64

    ini.equilibrium.ngrid = 200 ↔ [129, 257]
    @test ini.equilibrium.ngrid == 200
    @test typeof(ini.equilibrium.ngrid) <: Int

    ini.build.symmetric = true ↔ (false, true)
    @test ini.build.symmetric == true
    @test typeof(ini.build.symmetric) <: Bool
end

@testset "evolve_densities :zeff" begin
    dd = IMAS.json2imas(joinpath(pkgdir(FUSE), "sample", "D3D_standard_Lmode.json"); show_warnings=false)
    cp1d = dd.core_profiles.profiles_1d[]

    # fast electrons are not among the ions that get :zeff: they are registered as :fixed
    cp1d.electrons.density_fast = zero(cp1d.electrons.density_thermal)
    evolve_densities = FUSE.setup_density_evolution_electron_flux_match_fixed_zeff(cp1d)
    @test !haskey(evolve_densities, :electrons_fast)
    @test FUSE.check_evolve_densities(cp1d, evolve_densities) === nothing
    @test evolve_densities[:electrons_fast] == :fixed

    # a fast ion can be :fixed too
    evolve_densities[:D_fast] = :fixed
    @test FUSE.check_evolve_densities(cp1d, evolve_densities) === nothing

    # every thermal ion must be :zeff
    evolve_densities[:C] = :fixed
    @test_throws AssertionError FUSE.check_evolve_densities(cp1d, evolve_densities)
end

@testset "costing" begin
    value = FUSE.IMAS.Measurements.value

    # CPI table: every month is read (the last year may be incomplete)
    CPI = FUSE.load_inflation_rate()
    months = [:Jan, :Feb, :Mar, :Apr, :May, :Jun, :Jul, :Aug, :Sep, :Oct, :Nov, :Dec]
    @test all(!any(ismissing, CPI[1:end-1, month]) for month in months)
    # annual averages published by the U.S. Bureau of Labor Statistics
    @test CPI[findfirst(==(2023), CPI.Year), "Year Avg"] ≈ 304.702 atol = 1e-3
    @test CPI[findfirst(==(2024), CPI.Year), "Year Avg"] ≈ 313.689 atol = 1e-3

    # ARIES: land and buildings of a plant that is given by its powers alone
    dd = IMAS.dd()
    bop = dd.balance_of_plant
    bop.time = [0.0]
    bop.power_plant.power_cycle_type = "rankine"
    bop.power_plant.total_heat_supplied = [2.5e9]
    bop.power_plant.power_electric_generated = [1.0e9]
    bop.power_electric_net = [8.0e8]
    act = FUSE.ParametersActors()
    act.ActorCosting.model = :ARIES
    FUSE.ActorCosting(dd, act)
    systems = dd.costing.cost_direct_capital.system
    facility = systems[findfirst(sys -> sys.name == "facility", systems)]
    cost(name) = value(facility.subsystem[findfirst(sub -> sub.name == name, facility.subsystem)].cost)
    par = act.ActorCostingARIES
    da = FUSE.DollarAdjust{FUSE.IMAS.Measurements.Measurement{Float64}}(act.ActorCosting.future_inflation_rate, act.ActorCosting.construction_start_year, missing, missing)
    @test cost("land") ≈ value(FUSE.cost_direct_capital_ARIES(Val(:land), par.land_space, 8.0e8, da))
    @test cost("buildings") ≈ value(FUSE.cost_direct_capital_ARIES(Val(:buildings), par.land_space, par.building_volume, 1.0e9, 2.5e9, 8.0e8, da))

    # yearly cost of blanket replacement: the capital cost of the blanket over its lifetime
    @test value(FUSE.cost_operations_ARIES(Val(:blanket_replacement), 350.0, 7.0, da)) == 50.0
    @test value(FUSE.cost_operations_ARIES(Val(:blanket_replacement), nothing, 7.0, da)) == 0.0
end

@testset "ActorPowerNeeds EU_DEMO" begin
    dd = IMAS.dd()
    act = FUSE.ParametersActors()
    FUSE.ActorPowerNeeds(dd, act; model=:EU_DEMO)
    needs = Dict(sys.name => sys.power for sys in dd.balance_of_plant.power_electric_plant_operation.system)
    @test needs["cryostat"] == [30e6]
    @test needs["tritium_handling"] == [15e6]
end

@testset "equalize_symmetric_areas!" begin
    # coils 2 and 3 mirror each other; coil 1 sits next to coil 2, within a coil size of the mirror image of coil 3
    # coil 4 has no mirror; coils 5 and 6 mirror each other
    r = [5.1, 5.0, 5.0, 8.0, 1.0, 1.0]
    z = [1.0, 1.0, -1.0, 0.5, 2.0, -2.0]
    areas = [1.0, 2.0, 3.0, 0.5, 4.0, 6.0]
    @test FUSE.equalize_symmetric_areas!(areas, r, z) == [1.0, 3.0, 3.0, 0.5, 6.0, 6.0]
end

@testset "Checkpoint" begin
    a_old = a = [1]
    b_old = b = [2]
    dd_old = dd = IMAS.dd()
    ini_old = ini = FUSE.ParametersInits()

    # checkpoint empty
    empty!(FUSE.checkpoint)

    # checkin
    @checkin :test1 a b dd

    a[1] = 100

    @test a[1] == 100
    @test a === a_old
    @test b[1] == 2
    @test b === b_old
    @test dd === dd_old

    # add to the same checking
    @checkin :test1 ini

    # partial checkout
    @checkout :test1 a b

    @test a[1] == 1
    @test a !== a_old
    @test b[1] == 2
    @test b !== b_old

    # partial checkout
    @checkout :test1 dd ini

    @test dd !== dd_old
    @test ini !== ini_old

    # checkpoint empty
    empty!(FUSE.checkpoint)
    @assert isempty(FUSE.checkpoint)
end