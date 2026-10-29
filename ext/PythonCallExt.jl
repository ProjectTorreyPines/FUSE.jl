module PythonCallExt

using FUSE
import PythonCall
import PythonCall: Py, pyimport, pyconvert, pyisinstance, pybuiltins

"""
    _solwidth_import(solwidth_path::AbstractString)

Import the `solwidth` package (and its flat `solwidth_data` module), adding `solwidth_path` to `sys.path` if given
"""
function _solwidth_import(solwidth_path::AbstractString)
    if !isempty(solwidth_path)
        path = abspath(expanduser(solwidth_path))
        sys = pyimport("sys")
        if !(path in pyconvert(Vector{String}, sys.path))
            sys.path.insert(0, path)
        end
    end
    solwidth = pyimport("solwidth")
    # importable once solwidth/__init__ has put the package directory on sys.path
    solwidth_data = pyimport("solwidth_data")
    return solwidth, solwidth_data
end

function _from_py(value::Py)
    if pyisinstance(value, pybuiltins.bool)
        return pyconvert(Bool, value)
    elseif pyisinstance(value, pybuiltins.str)
        return pyconvert(String, value)
    end
    try
        return pyconvert(Float64, pybuiltins.float(value))
    catch
        return pyconvert(String, pybuiltins.str(value))
    end
end

function FUSE._solwidth_solve(model::AbstractString, pars::AbstractDict, solwidth_path::AbstractString)
    solwidth, solwidth_data = _solwidth_import(solwidth_path)

    kw = Dict{Symbol,Any}(Symbol(k) => v for (k, v) in pars)
    py_pars = solwidth_data.SolverParameters(; kw...)
    py_result = solwidth.solve_one(solwidth.load_model(model), py_pars)

    result = Dict{String,Any}()
    for (k, v) in pyimport("dataclasses").asdict(py_result).items()
        result[pyconvert(String, k)] = _from_py(v)
    end
    # solve_one attaches these outside the dataclass fields
    for extra in ("status", "nfev")
        if PythonCall.pyhasattr(py_result, extra)
            result[extra] = _from_py(PythonCall.pygetattr(py_result, extra))
        end
    end
    return result
end

end
