using FUSE
using Test

@testset "act_migrations" begin
    act = FUSE.ParametersActors()
    schema_keys = Set(string(k) for k in keys(act))
    for old_name in keys(FUSE.ACT_ACTOR_MIGRATIONS)
        @test old_name ∉ schema_keys
    end
end

@testset "GC" begin
    act = FUSE.ParametersActors()
    @test parent(act.ActorNeoclassical) === act
    GC.gc()
    @test parent(act.ActorNeoclassical) === act
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

    # yearly cost of blanket replacement: the cost of the blanket (in dollars of the current year) over its lifetime
    da.year_assessed = FUSE.Dates.year(FUSE.Dates.now())
    expected = value(FUSE.future_dollars(350.0 / 7.0, da))
    @test value(FUSE.cost_operations_ARIES(Val(:blanket_replacement), 350.0, 7.0, da)) ≈ expected
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