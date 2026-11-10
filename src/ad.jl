# ============================================================ #
# Forward-mode automatic differentiation through whole actors  #
# ============================================================ #
# An actor that receives a `dd` whose element type is a ForwardDiff.Dual carries the derivatives
# of everything it computes with respect to the seeded inputs. Three rules are followed throughout:
#
#  1. Float64 runs are untouched: every Dual-specific branch is selected by the element type.
#  2. Inner implicit solves (flux matching, pedestal and equilibrium solves) are converged on the
#     value parts with the same algorithm, iterates and stopping decisions as a Float64 run. The
#     derivatives of the solution then follow from the implicit function theorem applied to the
#     converged residual, not from differentiating the iterations.
#  3. Control flow (convergence tests, mode switches, event times, grid indices) only ever looks at
#     value parts.

"""
    strip_partials(x)

Returns `x` with its derivative part set to zero (identity for non-Dual numbers)
"""
strip_partials(x::Real) = x
strip_partials(x::D) where {D<:ForwardDiff.Dual} = D(ForwardDiff.value(x))

"""
    dual_values(x)

Value part of a Dual number or array of Dual numbers (identity otherwise)
"""
dual_values(x::Real) = x
dual_values(x::ForwardDiff.Dual) = ForwardDiff.value(x)
dual_values(x::AbstractArray{<:Real}) = x
dual_values(x::AbstractArray{<:ForwardDiff.Dual}) = ForwardDiff.value.(x)

"""
    partials_matrix(x::AbstractVector{<:ForwardDiff.Dual})

Matrix of size (length(x), N) with the N partial derivatives of each element of `x`
"""
function partials_matrix(x::AbstractVector{D}) where {T,V,N,D<:ForwardDiff.Dual{T,V,N}}
    out = Matrix{V}(undef, length(x), N)
    for (i, xi) in enumerate(x)
        for k in 1:N
            out[i, k] = ForwardDiff.partials(xi, k)
        end
    end
    return out
end

"""
    with_partials(::Type{D}, values::AbstractVector{<:Real}, partials::AbstractMatrix{<:Real})

Vector of Dual numbers of type `D` built from values and a (length(values), N) matrix of partials
"""
function with_partials(::Type{D}, values::AbstractVector{<:Real}, partials::AbstractMatrix{<:Real}) where {T,V,N,D<:ForwardDiff.Dual{T,V,N}}
    @assert size(partials) == (length(values), N)
    return D[D(V(values[i]), ForwardDiff.Partials{N,V}(ntuple(k -> V(partials[i, k]), N))) for i in eachindex(values)]
end

"""
    rethrow_if_dual_type_error(::Type{D}, e)

Several actors catch any exception and carry on with a previous or default state. With Dual numbers
that would turn a missing method or a failed conversion into a silently wrong derivative, so those
errors are always rethrown when the element type is a Dual. No effect for other element types.
"""
rethrow_if_dual_type_error(::Type{<:Real}, e) = nothing
function rethrow_if_dual_type_error(::Type{<:ForwardDiff.Dual}, e)
    if e isa MethodError || e isa InexactError || e isa TypeError || e isa ForwardDiff.DualMismatchError
        rethrow(e)
    elseif e isa ErrorException && occursin("Dual", e.msg)
        rethrow(e)
    end
    return nothing
end
