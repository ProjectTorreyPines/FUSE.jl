#= ========== =#
#  test_cases  #
#= ========== =#

include("_toksys.jl")

import OrderedCollections

struct TestCases <: AbstractDict{String,Any}
    data::OrderedCollections.OrderedDict{String,Any}
end

TestCases() = TestCases(OrderedCollections.OrderedDict{String,Any}())

function Base.getindex(tc::TestCases, key::AbstractString)
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
test_cases["D3D_Hmode"] = ([:D3D], Dict(:scenario => :H_mode, :use_scenario_sources => true))
test_cases["D3D_Lmode"] = ([:D3D], Dict(:scenario => :L_mode, :use_scenario_sources => false))
test_cases["D3D"] = ([:D3D], Dict(:scenario => :default))
test_cases["FPP"] = ([:FPP], Dict())
test_cases["CAT"] = ([:CAT], Dict())
test_cases["JET_HDB5"] = ([:HDB5], Dict(:tokamak => :JET, :case => 500))
test_cases["ARC"] = ([:ARC], Dict())
test_cases["SPARC"] = ([:SPARC], Dict())
test_cases["KDEMO"] = ([:KDEMO], Dict())
test_cases["DTT"] = ([:DTT], Dict())
test_cases["EXCITE"] = ([:EXCITE], Dict())
test_cases["MANTA"] = ([:MANTA], Dict())

"""
    test(testname::String, dd::IMAS.DD; sol::Bool=false, whole_facility::Bool=false)

Convenience function used to run test cases as done by CI
"""
function test(testname::String, dd::IMAS.DD; sol::Bool=false, whole_facility::Bool=false)
    if testname ∉ keys(test_cases)
        error("`$testname` not recognized. Valid tests are: $(collect(keys(test_cases)))")
    end

    args, kw = test_cases[testname]
    ini, act = case_parameters(args...; kw...)

    init(dd, ini, act)

    if sol
        IMAS.sol(dd)
    end

    if whole_facility
        FUSE.ActorWholeFacility(dd, act)
    end

    return act
end

#= ===== =#
#  cases  #
#= ===== =#
# NOTE only called once at precompile time, kernel needs to be restarted to include new file in cases
for filename in readdir(joinpath(@__DIR__))
    if !startswith(splitpath(filename)[end], "_") && endswith(filename, ".jl")
        include(filename)
    end
end

"""
    case_parameters(case::Symbol, args...; kw...)

return `ini` and `act` for a use-case

NOTE: if case starts with `test__` then the regression test cases are loaded
"""
function case_parameters(case::Symbol, args...; kw...)
    case_str = string(case)
    if startswith(case_str, "test__")
        tmp, kw = test_cases[split(case_str, "__"; limit=2)[end]]
        case = tmp[1]
        if length(tmp) == 1
            args = []
        else
            args = tmp[2:end]
        end
    end

    if length(methods(case_parameters, (Type{Val{case}},))) == 0
        error("case `$case` does not exist.\nPossible options are:\n\n$(join(["$method" for method in methods(case_parameters)],"\n"))")
    end

    ini, act = case_parameters(Val{case}, args...; kw...)

    if startswith(case_str, "test__")
        act.ActorStationaryPlasma.max_iter = 1
    end

    return ini, act
end

"""
    case_parameters_creation_from_ini(ini_case::AbstractParameters)

Useful function to create a case_parameters function for a case from an ini
"""
function case_parameters_creation_from_ini(ini_case::ParametersAllInits)
    ini_base = typeof(ini_case)()
    for field in keys(ini_base)
        code_string = ""
        for sub_field in keys(getproperty(ini_base, field))
            value_case = getproperty(getproperty(ini_case, field), sub_field, missing)
            value_base = getproperty(getproperty(ini_base, field), sub_field, missing)
            if value_case !== value_base
                if typeof(value_case) <: Symbol
                    code_string = "ini.$field.$sub_field = :$value_case"
                elseif typeof(value_case) <: String
                    code_string = "ini.$field.$sub_field = \"$value_case\""
                elseif typeof(value_case) <: OrderedCollections.OrderedDict
                    code_string = "ini.$field.$sub_field = OrderedCollections.OrderedDict(\n    " * join([":$k => $v," for (k, v) in value_case], "\n    ") * "\n)"
                else
                    code_string = "ini.$field.$sub_field = $value_case"
                end
                println(code_string)
            end
        end
        if !isempty(code_string)
            println()
        end
    end
end
