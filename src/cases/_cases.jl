#= ===== =#
#  cases  #
#= ===== =#
# NOTE only called once at precompile time, kernel needs to be restarted to include new file in cases
for filename in readdir(joinpath(@__DIR__))
    if !startswith(splitpath(filename)[end], "_") && endswith(filename, ".jl")
        include(filename)
    end
end

const use_cases = Dict()
use_cases["ITER_ods"] = ([:ITER], Dict(:init_from => :ods))
use_cases["ITER_scalars"] = ([:ITER], Dict(:init_from => :scalars))
use_cases["D3D_Hmode"] = ([:D3D, :H_mode], Dict())
use_cases["D3D_Lmode"] = ([:D3D, :L_mode], Dict())
use_cases["D3D"] = ([:D3D, :default], Dict())
use_cases["FPP"] = ([:FPP], Dict())
use_cases["CAT"] = ([:CAT], Dict())
use_cases["JET_HDB5"] = ([:HDB5], Dict(:tokamak => :JET, :database_case => 500))
use_cases["ARC"] = ([:ARC], Dict())
use_cases["SPARC"] = ([:SPARC], Dict(:init_from => :ods))
use_cases["KDEMO"] = ([:KDEMO], Dict())
use_cases["KDEMO_compact"] = ([:KDEMO_compact], Dict())
use_cases["DTT"] = ([:DTT], Dict())
use_cases["EXCITE"] = ([:EXCITE], Dict())
use_cases["MANTA"] = ([:MANTA], Dict())
use_cases["UNIT"] = ([:UNIT], Dict())

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

    if length(methods(case_parameters, (Type{Val{case}}, map(typeof, args)...))) == 0
        error("case `$case` does not exist.\nPossible options are:\n\n$(join(["$method" for method in methods(case_parameters)],"\n"))")
    end

    ini, act = case_parameters(Val{case}, args...; kw...)

    if startswith(case_str, "test__")
        ini_act_tests_customizations!(ini, act)
    end

    consistent_ini_act!(ini, act)

    set_new_base!(ini)
    set_new_base!(act)

    return ini, act
end
