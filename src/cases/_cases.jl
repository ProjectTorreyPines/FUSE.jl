
import OrderedCollections

const test_cases = OrderedCollections.OrderedDict{String,Any}()

test_cases["ITER_ods"] = ([:ITER], Dict(:init_from => :ods))
test_cases["ITER_scalars"] = ([:ITER], Dict(:init_from => :scalars))
test_cases["D3D_Hmode"] = ([:D3D], Dict(:scenario => :H_mode, :use_ods_sources => true))
test_cases["D3D_Lmode"] = ([:D3D], Dict(:scenario => :L_mode, :use_ods_sources => false))
test_cases["D3D"] = ([:D3D], Dict(:scenario => :default))
test_cases["FPP"] = ([:FPP], Dict())
test_cases["CAT"] = ([:CAT], Dict())
test_cases["JET_HDB5"] = ([:HDB5], Dict(:tokamak => :JET, :case => 500))
test_cases["ARC"] = ([:ARC], Dict())
test_cases["SPARC"] = ([:SPARC], Dict())
test_cases["KDEMO"] = ([:KDEMO], Dict())
test_cases["DTT"] = ([:DTT], Dict())
test_cases["EXCITE"] = ([:EXCITE], Dict())
test_cases["EXCITE"] = ([:EXCITE], Dict())
test_cases["STEP"] = ([:STEP], Dict(:init_from => :scalars, :pf_from=>:ods))
#= ===== =#
#  cases  #
#= ===== =#
# NOTE only called once at precompile time, kernel needs to be restarted to include new file in cases
for filename in readdir(joinpath(@__DIR__))
    if !startswith(splitpath(filename)[end], "_") && endswith(filename, ".jl")
        include(filename)
    end
end

function case_parameters(case::Symbol; kw...)
    if length(methods(case_parameters, (Type{Val{case}},))) == 0
        throw(InexistentParameterException([case]))
    end
    return case_parameters(Val{case}; kw...)
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
