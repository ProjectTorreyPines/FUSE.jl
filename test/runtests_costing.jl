using Test
using FUSE
using Measurements: value, uncertainty

@testset "costing_year_uncertainty" begin
    cpi = FUSE.load_inflation_rate()
    cpi_year_min = minimum(cpi.Year)
    cpi_year_max = maximum(cpi.Year)

    # Regression: uncertain historical year must work and propagate uncertainty
    da_hist = FUSE.DollarAdjust{FUSE.Measurement{Float64}}(0.025 ± 0.0, 2015.0 ± 5.0, 2016, missing)
    val_hist = FUSE.future_dollars(100.0, da_hist)
    @test val_hist isa FUSE.Measurement{Float64}
    @test !isnan(value(val_hist))
    @test uncertainty(val_hist) > 0.0
    cpi_hist = FUSE.load_inflation_rate()
    cpi_2016 = cpi_hist[findfirst(==(2016), cpi_hist.Year), "Year Avg"]
    cpi_grad_2015 = FUSE.cpi_year_average_gradient(cpi_hist, 2015.0)
    expected_hist_unc = abs((100.0 / cpi_2016) * cpi_grad_2015 * 5.0)
    @test isapprox(uncertainty(val_hist), expected_hist_unc; rtol=1e-12)

    # Zero uncertainty should match scalar path
    da_scalar = FUSE.DollarAdjust{Float64}(0.025, 2015.0, 2016, missing)
    val_scalar = FUSE.future_dollars(100.0, da_scalar)
    da_zero = FUSE.DollarAdjust{FUSE.Measurement{Float64}}(0.025 ± 0.0, 2015.0 ± 0.0, 2016, missing)
    val_zero = FUSE.future_dollars(100.0, da_zero)
    @test isapprox(value(val_zero), val_scalar; rtol=1e-12)

    # Hard errors for uncertain ranges outside available historical CPI interval
    da_too_early = FUSE.DollarAdjust{FUSE.Measurement{Float64}}(0.025 ± 0.0, (cpi_year_min - 1.0) ± 0.5, 2016, missing)
    @test_throws ErrorException FUSE.future_dollars(100.0, da_too_early)

    da_overlap = FUSE.DollarAdjust{FUSE.Measurement{Float64}}(0.025 ± 0.0, cpi_year_max ± 2.0, 2016, missing)
    @test_throws ErrorException FUSE.future_dollars(100.0, da_overlap)

    # Fully future uncertain year range remains supported
    da_future = FUSE.DollarAdjust{FUSE.Measurement{Float64}}(0.025 ± 0.0, (cpi_year_max + 5.0) ± 1.0, 2016, missing)
    val_future = FUSE.future_dollars(100.0, da_future)
    @test val_future isa FUSE.Measurement{Float64}
end
