function no_Dual(x)
    if typeof(x) <: ForwardDiff.Dual
        x = x.value
        return no_Dual(x)
    else
        return x
    end
end

function unwrap(v, inplace=false)
    unwrapped = inplace ? v : copy(v)
    for i in 2:length(v)
        while (unwrapped[i] - unwrapped[i - 1] >= pi)
            unwrapped[i] -= 2pi
        end
        while (unwrapped[i] - unwrapped[i - 1] <= -pi)
            unwrapped[i] += 2pi
        end
    end
    return unwrapped
end

function atan_eq(r, z, r0, z0)
    if r[1] == r[end] && z[1] == z[end]
        r = r[1:end - 1]
        z = z[1:end - 1]
    end
    θ = unwrap(atan.(z .- z0, r .- r0))
    if θ[2] < θ[1]
        r = reverse(r)
        z = reverse(z)
        θ = reverse(θ)
    end
    return r, z, θ
end

function two_curves_same_θ(r1, z1, r2, z2, scheme=:cubic)
    r0 = (sum(r1) / length(r1) + sum(r2) / length(r2)) / 2.0
    z0 = (sum(z1) / length(z1) + sum(z2) / length(z2)) / 2.0
    r1, z1, θ1 = atan_eq(r1, z1, r0, z0)
    r2, z2, θ2 = atan_eq(r2, z2, r0, z0)
    if length(θ2) > length(θ1)
        r1 = IMAS.interp(vcat(θ1 .- 2 * π, θ1, θ1 .+ 2 * π), vcat(r1, r1, r1), scheme=scheme).(θ2)
        z1 = IMAS.interp(vcat(θ1 .- 2 * π, θ1, θ1 .+ 2 * π), vcat(z1, z1, z1), scheme=scheme).(θ2)
        θ = θ2
    else
        r2 = IMAS.interp(vcat(θ2 .- 2 * π, θ2, θ2 .+ 2 * π), vcat(r2, r2, r2), scheme=scheme).(θ1)
        z2 = IMAS.interp(vcat(θ2 .- 2 * π, θ2, θ2 .+ 2 * π), vcat(z2, z2, z2), scheme=scheme).(θ1)
        θ = θ1
    end
    return r1, z1, r2, z2, θ
end
