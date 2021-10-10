struct StopIteration <: Exception end

function no_Dual(x)
    if typeof(x) <: ForwardDiff.Dual
        x = x.value
        return no_Dual(x)
    else
        return x
    end
end