"""
    case_parameters(::Type{Val{:FPP}}; kw...)::Tuple{ParametersAllInits,ParametersAllActors}

Disambiguates between different FPP versions
"""
function case_parameters(design::Type{Val{:FPP}}; version::Symbol, kw...)::Tuple{ParametersAllInits,ParametersAllActors}
    if version in [:v1, :v1_demount]
        ini, act = case_parameters(:FPPv1; version, kw...)
    elseif version in [:v2]
        ini, act = case_parameters(:FPPv2)
    end
    return ini, act
end
