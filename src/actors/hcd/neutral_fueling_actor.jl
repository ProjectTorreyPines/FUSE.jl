#= ====== =#
#  PELLET  #
#= ====== =#
Base.@kwdef mutable struct FUSEparameters__ActorNeutralFueling{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    τp::Entry{Float64} = Entry{Float64}("-", "Particle confinement time"; default=1.0)
    T_wall::Entry{Float64} = Entry{Float64}("-", "Wall temperature (eV)"; default=10.0)
    model::Switch{Symbol} = Switch{Symbol}([:neucg, :none], "-", "Neutral gas fueling model"; default=:neucg)
end

mutable struct ActorNeutralFueling{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorNeutralFueling{P}
    function ActorNeutralFueling(dd::IMAS.dd{D}, par::FUSEparameters__ActorNeutralFueling{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorNeutralFueling)
        par = par(kw...)
        return new{D,P}(dd, par)
    end
end

""" 
    g0(x::Float64)::Float64

    Interpolation of G function (Equation 12 in Burrel et al.) pulled from onetwo. Used where QuadGK calculation is slow
"""
function g0(x::Float64)::Float64
    g0f1(x) = (9.9995058e-1 + x*(7.7628443e+2 + x*(-4.3649686e+3 + x*6.1480022e+4))) /
              (1.0 + x*7.8621464e+2)
    
    g0f2(x) = (9.9703418e-1 + x*(7.7783588e+1 + x*(3.9012546e+2 + x*(-8.9205431e+2 + x*1.0037943e+3)))) /
              (1.0 + x*(8.4398466e+1 + x*7.1649574e+2))
    
    g0f3(x) = (9.7552510e-1 + x*(7.0647154 + x*(-4.0953920 + x*(9.0556774e-1 - x*5.8788928e-2)))) /
              (1.0 + x*(1.1344206e+1 + x*1.5956867e+1))
    
    g0f4(x) = (8.4513992e-1 + x*(-2.2811875e-1 + x*(2.5926818e-2 + x*(-1.4549910e-3 + x*3.3570582e-5)))) /
              (1.0 + x*(2.0520190 + x*(6.1145865e-1 + x*1.2572764e-1)))
    
    g0f5(x) = (1.9937501e-3 + x*(-3.2073160e-4 + x*(1.7925104e-5 - x*3.4571807e-7))) /
              (1.0 + x*(-3.3316230e-1 + x*3.6461690e-2))
    
    g0f6(x) = (-2.4903487e-7 + x*7.3163546e-9) /
              (1.0 + x*(-3.2915104e-1 + x*(4.3949080e-2 + x*(-2.9908526e-3 + x*(1.0457349e-4 - x*1.5316628e-6)))))

    # Conditional logic
    if x > 8.0
        if x > 15.0
            if x > 23.0
                return 0.0  # 23 < x
            else
                return g0f6(x)  # 15 < x < 23
            end
        else
            return g0f5(x)  # 8 < x < 15
        end
    elseif x > 1.0
        return g0f4(x)  # 1 < x < 8
    elseif x > 0.1
        return g0f3(x)  # 0.1 < x < 1
    elseif x > 0.01
        return g0f2(x)  # 0.01 < x < 0.1
    else
        return g0f1(x)  # 0.0 < x < 0.01
    end
end


""" 
    cxr(x::Vector{Float64})

This function calculates the charge exchange rate for hydrogen atoms interacting with protons in units of cm**3/s.
x is in units of keV for 1.0e-3 .le. x .le. 100, the the formula is taken from the paper by r.l. freeman and e.m. jones
clm-r 137 culham laboratory 1974 for x.lt.1.0e-3, a rate coefficient derived from an analytic average
over the approximate cross section, sigma=0.6937e-14*(1.0-0.155*LOG10 (e/1ev))**2 is used.  this cross section is
an approximation to that given by Riviere, Nuclear Fusion 11,363(1971).
"""
function cxr(x::Vector{Float64})
    tc = log.(x) .+ 6.9077553
    dum = 0.2530205e-4 .- tc .* 0.8230751e-6
    tc  = -0.1841757e2 .+ tc .* (0.528295 .- tc .* (0.2200477 .- tc.* (0.9750192e-1 .- tc .*
           (0.1749183e-1 .- tc .* (0.4954298e-3 .+
           tc .* (0.2174910e-3 .- tc.*dum))))))
   return  exp.(tc)
end

"""     
    eir(x::Vector{Float64})

This function calculates the ionization rate of atomic hydrogen by ele
impact in units of cm**3/s.  x is in units of keV
the formula is taken from the paper by r.l. freeman and e.m. jones
clm-r 137 culham laboratory 1974
"""
function eir(x::Vector{Float64})
    ta = log.(x) .+ 6.9077553
    ta = (-0.3173850e2 .+
          ta .* (0.1143818e2 .- ta .* (0.3833998e1 .- ta .* (0.7046692 .- ta
          .* (0.7431486e-1 .- ta .* (0.4153749e-2 .-ta .* 0.9486967e-4))))))
    return exp.(ta)   
end

"""
    nuslv1(sn1::Vector{Float64}, k11::Matrix{Float64}, nr::Float64)

This function solves the 1 species integral equation with the
Gauss-Seidel iteration technique.      
"""
function nuslv1(sn1::Vector{Float64}, k11::Matrix{Float64}, nr::Int)

    nxmax = 300

    f1 = zeros(nr)
    pn1 = zeros(nr)
    maxit = 40
    tol = 1.0e-4

    # Set up factor with diagonal term, and initial guess
    for i in 1:nr
        f1[i] = 1.0 / (1.0 -k11[i,i])
        pn1[i] = sn1[i]
    end

    # Now iterate until the solution converges
    for it in 1:maxit
        del = 0.0
        for i in 1:nr
            tn1 = sn1[i]
            
            for j in 1:nr
                if j == i
                    continue
                end
                tn1 += k11[j,i] * pn1[j]
            end

            tn1 *= f1[i]
            del1 = 0.0
            if pn1[i] != 0.0
                del1 = abs((tn1 - pn1[i]) / pn1[i])
            elseif tn1 != 0.0
                del1 = 1.0
            end
            del = max(del, del1)
            pn1[i] = tn1
        end
        if del < tol
            return pn1
        end
    end
    return error("Could not solve for neutral density")
end

"""
    neucg(dd::IMAS.dd, par::FUSEparameters__ActorNeutralFueling)

    This function calculates K matrix of equation 11 in Burrel 1977
"""
function get_At(r::Vector{Float64}, n::Int, r0::Float64,odelr::Float64, ra::Vector{Float64}, nra::Int, a1::Vector{Float64})
    nram1 = nra - 1
    
    function a1f(i::Int, rdum::Float64)
        return a1[i] + (a1[i+1] - a1[i]) * (rdum - ra[i]) * odelr
    end

    function x(rdummy::Float64)
        return sqrt(rdummy^2 - r0^2)
    end
    
    r0sq = r0^2
    a1int = 0.0
    singloid = r0 * odelr
    i = min(Int(round(singloid)) + 1, nram1)
    a11 = a1f(i, r0)
    x1 = 0.0
    at1 = zeros(n)
    for j in 1:n
        while r[j] > ra[i+1]
            i += 1
            a12 = a1[i]
            x2 = x(ra[i])
            a1int += (a11 + a12) * (x2 - x1)
            a11 = a12
            x1 = x2
            if i > nra
                break
            end
        end
        
        
        a12 = a1f(i, r[j])
        x2 = x(r[j])
        a1int += (a11 + a12) * (x2 - x1)
        at1[j] = 0.5 * a1int
        a11 = a12
        x1 = x2
    end
    return at1
end

"""
    neucg(dd::IMAS.dd, par::FUSEparameters__ActorNeutralFueling)

This function calculates K matrix of equation 9 in Burrel 1977
"""
function get_K_matrix(atrr01::Matrix{Float64},
                      vth::Vector{Float64},
                      r::Vector{Float64},
                      nr::Int,
                      theta::Vector{Float64},
                      nt::Int,
                      sint::Vector{Float64})
                  
    k11 = zeros(nr, nr)
    wt = theta[2] - theta[1]
    G0(x) = 2.0/sqrt(π)*QuadGK.quadgk(η-> η.^0 .* exp.(-η.^2 .- x ./ η), 0., Inf)[1]

    for i in 1:nr-1
        for j in (i + 1):nr-1
            if i == 1
                apm = atrr01[j, 1]
                apovp = apm / vth[j]
                g0ap = g0(apovp)
    
                k11[j, 1] = π * g0ap / r[j]
                apovp = apm / vth[1]
                g0ap = g0(apovp)
                k11[1, j] = π * g0ap / r[j]
            else
                sum1p, sum1r, sum1wp, sum1wr= 0.0, 0.0, 0.0, 0.0
                pgt2 = r[j] * r[j]
                
                for k in 1:nt
                    r0 = r[i] * sint[k]
                    wtx = wt / sqrt(pgt2 - r0 * r0)
                    
                    atplt = atrr01[i, k]
                    atpgt = atrr01[j, k]
                    ap = atpgt + atplt
                    am = atpgt - atplt
                    
                    # r = rholt=r(i), p=rhogt=r(j)
                    apovp = abs(ap / vth[j])
                    amovp = abs(am / vth[j])
 
                    g0ap = g0(apovp)
                    g0am = g0(amovp)

                    sum1p += wtx * (g0ap + g0am)

                    
                    # r = rhogt=r(j), p=rholt=r(i)
                    apovp = abs(ap / vth[i])
                    amovp = abs(am / vth[i])
                    g0ap = g0(apovp)
                    g0am = g0(amovp)

                    sum1r += wtx * (g0ap + g0am)                    
                end
                
                k11[j, i] = sum1p
                k11[i, j] = sum1r
            end
    
        end
    end
    return k11
end

"""
    neucg(dd::IMAS.dd, par::FUSEparameters__ActorNeutralFueling)

This model calculates neutral density source fueling based on the neucg 
model. K. Burrell,  Journal of Computational Physics 27.1 (1978): 88-102.
"""
function neucg(dd::IMAS.dd, par::FUSEparameters__ActorNeutralFueling)

    τp = par.τp
    T_wall = par.T_wall
    
    surface = dd.core_profiles.profiles_1d[].grid.surface
    volume = dd.core_profiles.profiles_1d[].grid.volume

    Te = dd.core_profiles.profiles_1d[].electrons.temperature/1e3
    rho = dd.core_profiles.profiles_1d[].grid.rho_tor_norm
    ne = dd.core_profiles.profiles_1d[].electrons.density/1e6
    ni  = dd.core_profiles.profiles_1d[].ion[1].density/1e6
    mi = dd.core_profiles.profiles_1d[].ion[1].element[1].a
    kappa = dd.equilibrium.time_slice[].profiles_1d.elongation[end]
    rmin = 0.5.*(dd.equilibrium.time_slice[].profiles_1d.r_outboard-
                 dd.equilibrium.time_slice[].profiles_1d.r_inboard)
    amin = 0.5.*(dd.equilibrium.time_slice[].profiles_1d.r_outboard[end]-
                 dd.equilibrium.time_slice[].profiles_1d.r_inboard[end])
    Rmaj = 0.5.*(dd.equilibrium.time_slice[].profiles_1d.r_outboard[end]+
                 dd.equilibrium.time_slice[].profiles_1d.r_inboard[end])
    raneut = sqrt.(kappa).*amin
    rhoa = dd.equilibrium.time_slice[].profiles_1d.rho_tor[end]
    ratef = raneut / rhoa
    eionr = ratef .* eir(Te)
    eirate = eir(dd.core_profiles.profiles_1d[].electrons.temperature/1e3).*ne

    vth = 1e2 .* sqrt.(2*1.6e-19*(Te*1e3) ./ (mi*1.6e-27))
 
    a1 = ne .* eionr .+ ni .* cxr(Te/mi)
    b1 = ni .* cxr(Te/mi) 
    vwall  = 100.0 * sqrt.(2.0*1.6e-19*T_wall ./ (mi*1.67e-27))
    
    flux1 = IMAS.trapz(surface, ni)/τp
    swall = 2.0*flux1/vwall

    theta = collect(0:0.05*(pi/2):pi/2)
    
    r = rho*1e2*rhoa
    ra = deepcopy(r)
    nr = length(r)
    nra = length(ra)
    nt = length(theta)
    
    odelr = 1.0/(r[2] - r[1])
    
    sint = sin.(theta)
    
    k11 = zeros(nr, nr)
    atrr01 = zeros(nr, nt)
    atar01 = zeros(nr, nt)

    for i in 1:nr
        k11[i, i] = 0.0
        
        for k in 1:nt
            r0 = r[i] * sint[k]
            atrr01[i:nr, k] .= get_At(r[i:nr], nr - i + 1, r0, odelr,ra,nra, a1)
            atar01[i, k] = atrr01[nr, k]
        end
    
        if i == 1
            for kk in 2:nt
                for ii in 1:nr
                atrr01[ii, kk] = atrr01[ii, 1]
                end
                atar01[1, kk] = atrr01[nr, kk]
            end
        end
    end

    
    G(x,n) = 2.0/sqrt(π)*QuadGK.quadgk(η-> η.^n .* exp.(-η.^2 .- x ./ η), 0., Inf)[1]
    G1h = zeros(nr, nt)
    Ap = zeros(nr, nt)
    Am = zeros(nr,nt)
    h = zeros(nr)
    for (ir,rr) in enumerate(rho)
        for (it,tt) in enumerate(theta)
            Ap[ir,it] = abs.(atrr01[ir,it] + atar01[ir,it])
            Am[ir,it] = abs.(atrr01[ir,it] - atar01[ir,it])
            G1h[ir,it] = G(Ap[ir,it]/vwall,1)+G(Am[ir,it]/vwall,1)
        end
    end


    for (ir,rr) in enumerate(rho)
        h[ir] = swall*IMAS.trapz(theta,G1h[ir,:])
    end

    k11 = get_K_matrix(atrr01,vth,r,nr,theta,nt,sint)
    k11 = b1./vth./sqrt(π) .* k11 

    pn11 = nuslv1(h,k11, nr)

    cm3_to_m3 = 1e6
    fn11 = pn11 * cm3_to_m3
    sione = fn11.*eirate

    correction = (trapz(volume,ni)/τp)/trapz(volume,sione/cm3_to_m3)
 
    return sione*correction,fn11*correction
    
end

"""
    ActorNeutralFueling(dd::IMAS.dd, act::ParametersAllActors; kw...)

Estimates the neutral fueling deposition
"""
function ActorNeutralFueling(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorNeutralFueling(dd, act.ActorNeutralFueling; kw...)
    step(actor)
    finalize(actor)
    return actor
end


function _step(actor::ActorNeutralFueling)
    dd = actor.dd
    par = actor.par

    if par.model == :neucg
        Sneut, nneut = neucg(dd, par)
    else
        return actor
    end

    cp1d = dd.core_profiles.profiles_1d[]
    resize!(cp1d.neutral, 1)[1]
    actor.dd.core_profiles.profiles_1d[].neutral[1].density = nneut

    source = resize!(dd.core_sources.source, :gas_puff, "identifier.name" => "gas"; wipe=false)
    IMAS.new_source(source, source.identifier.index, "gas", cp1d.grid.rho_tor_norm, cp1d.grid.volume, cp1d.grid.area; electrons_particles=Sneut)


    return actor
end
