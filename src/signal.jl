"""
    step(t::T)::T where T

Unitary step triggered at t=0
"""
function step(t::T)::T where T
    return @. t >= 0.0
end

"""
    step(t::T, t_start::Float64)::T where T

Unitary step triggered at t=t_start
"""
function step(t::T, t_start::Float64)::T where T
    return step(t - t_start)
end

"""
    pulse(t::T)::T where T

Unitary pulse with width of 1, starting at t=0
"""
function pulse(t::T)::T where T
    a = step(t)
    b = step(-t + 1)
    return @. a * b
end

"""
    pulse(t::T, t_start::Float64, Δt::Float64)::T where T

Unitary pulse with given width Δt, starting at t=t_start
"""
function pulse(t::T, t_start::Float64, Δt::Float64)::T where T
    return pulse((t - t_start) / Δt)
end

"""
    ramp(t::T)::T where T

Unitary ramp from t=0 to t=1
"""
function ramp(t::T)::T where T
    a = @. t * (t < 1) * (t > 0)
    b = @. t >= 1
    return @. a + b
end

"""
    ramp(t::T, t_start::Float64, Δt::Float64)::T where T

Unitary ramp over the range of width Δt, starting at t=t_start
"""
function ramp(t::T, t_start::Float64, Δt::Float64)::T where T
    return ramp((t - t_start) / Δt)
end

"""
    trap(t::T, flattop_fraction::Float64)::T where T

Unitary trapezoid. The `flattop_fraction` defines the fraction of flat top within [0,1]
"""
function trap(t::T, flattop_fraction::Float64)::T where T
    @assert 0 <= flattop_fraction <= 1 "trap flattop_fraction must be between [0,1]"
    k = 1 / ((1 - flattop_fraction) / 2)
    if flattop_fraction == 1.0
        return pulse(t)
    else
        a = @. ramp(t * k) * (t < 0.5)
        b = @. ramp(-t * k + k) * (t >= 0.5)
        return @. a + b
    end
end

"""
    trap(t::T, t_start::Float64, Δt::Float64, flattop_fraction::Float64)::T where T

Unitary trapezoid, shifted by t_start and scaled by Δt. The `flattop_fraction` defines the fraction of flat top within [0,1]
"""
function trap(t::T, t_start::Float64, Δt::Float64, flattop_fraction::Float64)::T where T
    return trap((t - t_start) / Δt, flattop_fraction)
end
