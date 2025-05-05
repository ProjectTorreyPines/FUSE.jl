import Test

#= ========== =#
#  test_cases  #
#= ========== =#

import OrderedCollections

struct TestCases <: AbstractDict{String,Any}
    data::OrderedCollections.OrderedDict{String,Any}
end

TestCases() = TestCases(OrderedCollections.OrderedDict{String,Any}())

function Base.getindex(tc::TestCases, key::AbstractString)
    if key ∉ keys(tc.data)
        error("Test case \"$(key)\" not found. Valid test cases are: $(sort!(collect(keys(tc.data))))")
    end
    return getindex(tc.data, key)
end

function Base.setindex!(tc::TestCases, value, key::AbstractString)
    return tc.data[key] = value
end

Base.length(tc::TestCases) = length(tc.data)

function Base.iterate(tc::TestCases, state=1)
    return iterate(tc.data, state)
end

Base.haskey(tc::TestCases, key::AbstractString) = haskey(tc.data, key)

Base.keys(tc::TestCases) = keys(tc.data)

Base.get(tc::TestCases, key::AbstractString, default) = get(tc.data, key, default)

function Base.delete!(tc::TestCases, key::AbstractString)
    return delete!(tc.data, key)
end

function Base.show(io::IO, ::MIME"text/plain", test_cases::TestCases)
    n = maximum(length, keys(test_cases))
    for testname in sort!(collect(keys(test_cases)))
        (args, kw) = test_cases[testname]
        kw_str = join(["$k=$(repr(v))" for (k, v) in kw], ", ")
        printstyled(io, "$(rpad(testname, n))"; bold=true)
        println(io, "   ini, act = FUSE.case_parameters($(join(repr.(args), ", "))$(isempty(kw_str) ? "" : "; ")$kw_str)")
    end
end

const test_cases = TestCases()

test_cases["ITER_ods"] = ([:ITER], Dict(:init_from => :ods))
test_cases["ITER_scalars"] = ([:ITER], Dict(:init_from => :scalars))
test_cases["D3D_Hmode"] = ([:D3D, :H_mode], Dict())
test_cases["D3D_Lmode"] = ([:D3D, :L_mode], Dict())
test_cases["D3D"] = ([:D3D, :default], Dict())
test_cases["FPP"] = ([:FPP], Dict())
test_cases["CAT"] = ([:CAT], Dict())
test_cases["JET_HDB5"] = ([:HDB5], Dict(:tokamak => :JET, :case => 500))
test_cases["ARC"] = ([:ARC], Dict())
test_cases["SPARC"] = ([:SPARC], Dict(:init_from => :ods))
test_cases["KDEMO"] = ([:KDEMO], Dict())
test_cases["KDEMO_compact"] = ([:KDEMO_compact], Dict())
test_cases["DTT"] = ([:DTT], Dict())
test_cases["EXCITE"] = ([:EXCITE], Dict())
test_cases["MANTA"] = ([:MANTA], Dict())
test_cases["UNIT"] = ([:UNIT], Dict())

"""
    test(testname::String, dd::IMAS.DD)

Convenience function used to run test cases as done by CI
"""
function test(testname::String, dd::IMAS.DD)
    if testname ∉ keys(test_cases)
        error("`$testname` not recognized. Valid tests are: $(collect(keys(test_cases)))")
    end

    args, kw = test_cases[testname]
    dd, ini, act = test(dd, args...; kw...)

    return (dd=dd, ini=ini, act=act)
end

function test(dd::IMAS.DD, args...; kw...)
    ini, act = case_parameters(args...; kw...)

    ini_act_tests_customizations!(ini, act)

    Test.@testset "init" begin
        init(dd, ini, act)
    end

    Test.@testset "sol" begin
        IMAS.sol(dd)
    end

    Test.@testset "whole_facility" begin
        ActorWholeFacility(dd, act)
    end

    return (dd=dd, ini=ini, act=act)
end

function test_ini_act_save_load(args...; kw...)
    ini, act = case_parameters(args...; kw...)

    Test.@testset "ini_dict" begin
        dict2ini(ini2dict(ini))
    end

    Test.@testset "act_dict" begin
        dict2act(act2dict(act))
    end

    Test.@testset "ini_yaml" begin
        ini.general.dd = missing # general.dd cannot be serialized
        ini_str = SimulationParameters.par2ystr(ini; skip_defaults=true, show_info=false)
        ini2 = SimulationParameters.ystr2par(ini_str, ParametersInits())
    end

    Test.@testset "act_yaml" begin
        act_str = SimulationParameters.par2ystr(act; skip_defaults=true, show_info=false)
        act2 = SimulationParameters.ystr2par(act_str, ParametersActors())
    end

    Test.@testset "ini_json" begin
        ini.general.dd = missing # general.dd cannot be serialized
        ini_str = SimulationParameters.par2jstr(ini)
        ini2 = SimulationParameters.jstr2par(ini_str, ParametersInits())
    end

    Test.@testset "act_json" begin
        act_str = SimulationParameters.par2jstr(act)
        act2 = SimulationParameters.jstr2par(act_str, ParametersActors())
    end

    Test.@testset "ini_hdf5" begin
        tmpdir = mktempdir()
        ini.general.dd = missing # general.dd cannot be serialized
        ini_str = SimulationParameters.par2hdf(ini, joinpath(tmpdir, "ini.h5"))
        ini2 = SimulationParameters.hdf2par(joinpath(tmpdir, "ini.h5"), ParametersInits())
    end

    Test.@testset "act_hdf5" begin
        tmpdir = mktempdir()
        act_str = SimulationParameters.par2hdf(act, joinpath(tmpdir, "act.h5"))
        act2 = SimulationParameters.hdf2par(joinpath(tmpdir, "act.h5"), ParametersActors())
    end

    return (ini=ini, act=act)
end

function ini_act_tests_customizations!(ini::ParametersAllInits, act::ParametersAllActors)
    # speedup the tests
    act.ActorStationaryPlasma.max_iterations = 2
    # use full model for ActorThermalPlant if environmental variable `FUSE_WITH_EXTENSIONS` is set
    if get(ENV, "FUSE_WITH_EXTENSIONS", "false") == "true"
        act.ActorThermalPlant.model = :network
    end
    return (ini=ini, act=act)
end
