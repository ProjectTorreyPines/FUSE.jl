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
        ini_act_tests_customizations!(ini, act)
    end

    set_new_base!(ini)
    set_new_base!(act)

    return ini, act
end
