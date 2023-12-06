"""
    step(t::T)::T where T

Unitary step triggered at t=0
"""
function step(t::T)::T where {T}
    return t >= 0.0
end

"""
    step(t::T, t_start::Float64)::T where T

Unitary step triggered at t=t_start
"""
function step(t::T, t_start::Float64)::T where {T}
    return step(t - t_start)
end

"""
    pulse(t::T)::T where T

Unitary pulse with width of 1, starting at t=0
"""
function pulse(t::T)::T where {T}
    a = step(t)
    b = step(-t + 1)
    return a * b
end

"""
    pulse(t::T, t_start::Float64, Δt::Float64)::T where T

Unitary pulse with given width Δt, starting at t=t_start
"""
function pulse(t::T, t_start::Float64, Δt::Float64)::T where {T}
    return pulse((t - t_start) / Δt)
end

"""
    ramp(t::T)::T where T

Unitary ramp from t=0 to t=1
"""
function ramp(t::T)::T where {T}
    a = t * (t < 1) * (t > 0)
    b = t >= 1
    return a + b
end

"""
    ramp(t::T, ramp_fraction::Float64)::T where {T}

Unitary ramp

The `ramp_fraction` defines the fraction of ramp with respect to 1.0 and must be between [0.0,1.0]

NOTE: This function is designed as is to be able to switch between `ramp(t, ramp_fraction)` and `trap(t, ramp_fraction)`.
"""
function ramp(t::T, ramp_fraction::Float64)::T where {T}
    @assert 0 <= ramp_fraction <= 1 "ramp `ramp_fraction` must be between [0.0,1.0]"
    k = 1.0 / ramp_fraction
    if ramp_fraction == 0.0
        return pulse(t)
    else
        return ramp(t * k)
    end
end

"""
    ramp(t::T, t_start::Float64, Δt::Float64)::T where T

Unitary ramp of duration Δt, starting at t=t_start
"""
function ramp(t::T, t_start::Float64, Δt::Float64)::T where {T}
    return ramp((t - t_start) / Δt)
end

"""
    trap(t::T, ramp_fraction::Float64)::T where T

Unitary trapezoid

The `ramp_fraction` defines the fraction of ramp with respect to flattop and must be between [0.0,0.5]
"""
function trap(t::T, ramp_fraction::Float64)::T where {T}
    @assert 0 <= ramp_fraction <= 0.5 "trap `ramp_fraction` must be between [0.0,0.5]"
    k = 1.0 / ramp_fraction
    if ramp_fraction == 0.0
        return pulse(t)
    else
        a = ramp(t * k) * (t < 0.5)
        b = ramp(-t * k + k) * (t >= 0.5)
        return a + b
    end
end

"""
    trap(t::T, t_start::Float64, Δt::Float64, ramp_fraction::Float64)::T where T

Unitary trapezoid of duration Δt, starting at t=t_start

The `ramp_fraction` defines the fraction of ramp with respect to flattop and must be between [0.0,0.5]
"""
function trap(t::T, t_start::Float64, Δt::Float64, ramp_fraction::Float64)::T where {T}
    return trap((t - t_start) / Δt, ramp_fraction)
end

"""
    gaus(t::T)::T where {T}

Unitary gaussian
"""
function gaus(t::T)::T where {T}
    return exp(-t^2 / 2.0)
end

"""
    gaus(t::T, t_start::Float64, Δt::Float64)::T where {T}

Unitary gaussian centered at t_start and with standard deviation Δt
"""
function gaus(t::T, t_start::Float64, Δt::Float64)::T where {T}
    return gaus((t - t_start) / Δt)
end

"""
    beta(t::T, mode::Float64)::T where {T}

Unitary beta distribution

The `mode` [-1.0, 1.0] defines how skewed the distribution is
"""
function beta(t::T, mode::Float64)::T where {T}
    @assert -1 <= mode <= 1 "beta `mode` must be between [-1.0,1.0]"
    if t < 0 || t > 1
        return 0  # Outside the support
    end

    α = 2.0  # keeping alpha constant

    # mode expressed from 0 to 1
    mode = (mode / 2.0) + 0.5

    # Special conditions where the mode is at the boundaries
    if mode == 0.0
        return t == 0.0 ? 1.0 : 0.0  # spike at 0
    elseif mode == 1.0
        return t == 1.0 ? 1.0 : 0.0  # spike at 1
    end

    if mode > 0.5
        t = -t + 1.0
        mode = -(mode - 0.5) + 0.5
    end

    # For other cases, we find the corresponding β from the mode
    β = ((α - 1) / mode) - (α - 2)

    # Calculate the value of the beta distribution at its mode
    peak_value = (mode^(α - 1) * (1 - mode)^(β - 1))

    # normalize
    return (t^(α - 1) * (1 - t)^(β - 1)) / peak_value
end

"""
    beta(t::T, t_start::Float64, Δt::Float64, mode::Float64)::T where {T}

Unitary beta distribution of duration Δt, starting at t=t_start

The `mode` [-1.0, 1.0] defines how skewed the distribution is
"""
function beta(t::T, t_start::Float64, Δt::Float64, mode::Float64)::T where {T}
    return beta((t - t_start) / Δt, mode)
end

# ====================== #
# TrajectoryPerturbation #
# ====================== #
mutable struct TrajectoryPerturbation
    t_start::Float64
    t_end::Float64
    nominal::Function
    upper::Function
    lower::Function
    t_perturb::Vector{Float64}
    v_perturb::Vector{Float64}
    correction::Any
end

"""
    TrajectoryPerturbation(t_start, t_end, nominal::Function, upper::Function, lower::Function, tt::Vector{Float64}, v::Vector{Float64}; method::Symbol=:linear)

Perturbation of time trace with vector of time `t` and values `v`

The perturbations are done as deviations from a nominal case and are all expressed in normalized way (between 0.0 and 1.0):

  - in time: between `t_start` and `t_end`

  - in value: between `upper` and `lower`
"""
function TrajectoryPerturbation(t_start, t_end, nominal::Function, upper::Function, lower::Function, t::AbstractVector{Float64}, v::AbstractVector{Float64}; method::Symbol=:linear)
    @assert length(t) == length(v)
    t_perturb = [0.0, t..., 1.0, 1E6] .* (t_end - t_start) .+ t_start
    nn = nominal.(t_perturb)
    uu = upper.(t_perturb)
    ll = lower.(t_perturb)

    vv = [0.0, v..., 0.0, 0.0]
    v_perturb = vv .* (uu .- nn) .* (vv .> 0) .- vv .* (ll .- nn) .* (vv .< 0)

    correction = IMAS.interp1d(t_perturb, v_perturb, method)

    return TrajectoryPerturbation(t_start, t_end, nominal, upper, lower, t_perturb, v_perturb, correction)
end

function TrajectoryPerturbation(t_start, t_end, nominal, upper, lower, n; method)
    t = range(0.0, 1.0, n + 2)[2:end-1]
    v = t .* 0.0
    return TrajectoryPerturbation(t_start, t_end, nominal, upper, lower, t, v; method)
end

function (tp::TrajectoryPerturbation)(t)
    out = tp.nominal.(t) + tp.correction.(t) .* (t .> tp.t_start) .* (t .< tp.t_end)
    out = max.(tp.lower.(t), min.(tp.upper.(t), out))
    return out
end

@recipe function plot_TrajectoryPerturbation(tp::TrajectoryPerturbation)
    @series begin
        label := "perturbed"
        lw := 2.0
        t -> tp(t)
    end
    @series begin
        primary := false
        seriestype := :scatter
        tp.t_perturb, tp.nominal.(tp.t_perturb) .+ tp.v_perturb
    end
    @series begin
        primary := false
        label := "upper/lower"
        ribbon := t -> (tp.upper(t) .- tp.lower(t)) / 2.0
        linewidth := 0.0
        t -> (tp.upper(t) .+ tp.lower(t)) / 2.0
    end
    @series begin
        label := "nominal"
        lw := 2.0
        xlim --> (tp.t_start, tp.t_end)
        linestyle := :dash
        tp.nominal
    end
end