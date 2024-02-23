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
    gaus(t::T, order::Float64=1.0)::T where {T}

Unitary gaussian
"""
function gaus(t::T, order::Float64=1.0)::T where {T}
    return exp(-(t^2 / 2.0)^order)
end

"""
    gaus(t::T, t_start::Float64, Δt::Float64, order::Float64=1.0)::T where {T}

Unitary gaussian centered at t_start and with standard deviation Δt
"""
function gaus(t::T, t_start::Float64, Δt::Float64, order::Float64=1.0)::T where {T}
    return gaus((t - t_start) / Δt, order)
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
