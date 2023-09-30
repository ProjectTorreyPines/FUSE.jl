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
    ramp(t::T, t_start::Float64, Δt::Float64)::T where T

Unitary ramp over the range of width Δt, starting at t=t_start
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
    @assert 0 <= ramp_fraction <= 0.5 "trap ramp_fraction must be between [0.0,0.5]"
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

Unitary trapezoid, with time shifted by t_start and scaled by Δt

The `ramp_fraction` defines the fraction of ramp with respect to flattop and must be between [0.0,0.5]
"""
function trap(t::T, t_start::Float64, Δt::Float64, ramp_fraction::Float64)::T where {T}
    return trap((t - t_start) / Δt, ramp_fraction)
end
