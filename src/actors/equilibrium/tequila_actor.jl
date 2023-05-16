import TEQUILA

#= =========== =#
#  ActorTEQUILA  #
#= =========== =#
Base.@kwdef mutable struct FUSEparameters__ActorTEQUILA{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    free_boundary::Entry{Bool} = Entry(Bool, "-", "Convert fixed boundary equilibrium to free boundary one"; default=true)
    number_of_radial_grid_points::Entry{Int} = Entry(Int, "-", "Number of TEQUILA radial grid points"; default=20)
    number_of_fourier_modes::Entry{Int} = Entry(Int, "-", "Number of modes for Fourier decomposition"; default=20)
    number_of_MXH_harmonics::Entry{Int} = Entry(Int, "-", "Number of Fourier harmonics in MXH representation of flux surfaces"; default=10)
    number_of_iterations::Entry{Int} = Entry(Int, "-", "Number of TEQUILA iterations"; default=5)
    psi_norm_boundary_cutoff::Entry{Float64} = Entry(Float64, "-", "Cutoff psi_norm for determining boundary"; default=0.999)
    do_plot::Entry{Bool} = Entry(Bool, "-", "Plot before and after actor"; default=false)
end

mutable struct ActorTEQUILA <: PlasmaAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorTEQUILA
    shot::Union{Nothing, TEQUILA.Shot}
    psib::Real
end

"""
    ActorTEQUILA(dd::IMAS.dd, act::ParametersAllActors; kw...)

Runs the Fixed boundary equilibrium solver TEQUILA
"""
function ActorTEQUILA(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorTEQUILA(kw...)
    actor = ActorTEQUILA(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function ActorTEQUILA(dd::IMAS.dd, par::FUSEparameters__ActorTEQUILA; kw...)
    logging_actor_init(ActorTEQUILA)
    par = par(kw...)
    ActorTEQUILA(dd, par, nothing, 0.0)
end

"""
    prepare(dd::IMAS.dd, :ActorTEQUILA, act::ParametersAllActors; kw...)

Prepare dd to run ActorTEQUILA
* Copy pressure from core_profiles to equilibrium
* Copy j_parallel from core_profiles to equilibrium
"""
function prepare(dd::IMAS.dd, ::Type{Val{:ActorTEQUILA}}, act::ParametersAllActors; kw...)
    eq1d = dd.equilibrium.time_slice[].profiles_1d
    cp1d = dd.core_profiles.profiles_1d[]
    eq1d.j_tor = IMAS.interp1d(cp1d.grid.psi_norm, cp1d.j_tor).(eq1d.psi_norm)
    eq1d.pressure = IMAS.interp1d(cp1d.grid.psi_norm, cp1d.pressure).(eq1d.psi_norm)
end

"""
    step(actor::ActorTEQUILA)

Runs TEQUILA on the r_z boundary, equilibrium pressure and equilibrium j_tor
"""
function _step(actor::ActorTEQUILA)
    dd = actor.dd
    eqt = dd.equilibrium.time_slice[]
    eq1d = eqt.profiles_1d
    par = actor.par

    psin = eq1d.psi_norm
    actor.psib = IMAS.interp1d(psin, eq1d.psi)(par.psi_norm_boundary_cutoff)
    r_bound, z_bound, _ = IMAS.flux_surface(eqt, actor.psib, true)

    mxh = IMAS.MXH(r_bound, z_bound, par.number_of_MXH_harmonics; optimize_fit=true)

    rho_pol = sqrt.(psin)
    rho_pol[1] = 0.0
    dp_dψ = TEQUILA.FE(rho_pol, eq1d.dpressure_dpsi)
    Jt_R = TEQUILA.FE(rho_pol, eq1d.j_tor .* eq1d.gm9)
    pbnd = eq1d.pressure[end]
    fbnd = eq1d.f[end]

    # TEQUILA shot
    shot = TEQUILA.Shot(par.number_of_radial_grid_points, par.number_of_fourier_modes, mxh;
                        dp_dψ, Jt_R, pbnd, fbnd)
    actor.shot = TEQUILA.solve(shot, par.number_of_iterations)

    if par.do_plot
        p=plot(dd.equilibrium;cx=true, label="before")
        for (idx,s) in enumerate(eachcol(actor.shot.surfaces[:,2:end]))
            if idx == length(actor.shot.surfaces[1,2:end])
                plot!(p,IMAS.MXH(s), color=:red,label="after TEQUILA")
            else
                plot!(p,IMAS.MXH(s), color=:red)
            end
        end
        display(p)
    end

    return actor
end

# finalize by converting TEQUILA shot to dd.equilibrium
function _finalize(actor::ActorTEQUILA)
    try
        tequila2imas(actor.shot, actor.dd.equilibrium;
                     psib=actor.psib, free_boundary=actor.par.free_boundary)
    catch e
        display(plot(actor.shot))
        rethrow(e)
    end
    return actor
end

function tequila2imas(shot::TEQUILA.Shot, eq::IMAS.equilibrium; psib=0.0, free_boundary=false)

    eqt = eq.time_slice[]
    eq1d = eqt.profiles_1d

    R0 = shot.R0fe(1.0)
    Z0 = shot.Z0fe(1.0)
    eq.vacuum_toroidal_field.r0 = R0
    @ddtime(eq.vacuum_toroidal_field.b0 = shot.fbnd / R0)
    eqt.boundary.geometric_axis.r = R0
    eqt.boundary.geometric_axis.z = Z0

    Rax = shot.R0fe(0.0)
    Zax = shot.Z0fe(0.0)
    eqt.global_quantities.magnetic_axis.r = Rax
    eqt.global_quantities.magnetic_axis.z = Zax

    n_grid = 10 * length(shot.ρ)

    psit = shot.C[2:2:end, 1]
    psii = range(psit[1], psit[end], n_grid)
    eq1d.psi = psii .+ psib
    eq1d.pressure = MXHEquilibrium.pressure.(Ref(shot), psii)
    eq1d.dpressure_dpsi = MXHEquilibrium.pressure_gradient.(Ref(shot), psii)
    eq1d.f = TEQUILA.Fpol.(Ref(shot), psii)
    eq1d.f_df_dpsi = TEQUILA.dFpol_dψ.(Ref(shot), psii)

    resize!(eqt.profiles_2d, 1)
    eq2d = eqt.profiles_2d[1]

    R0 = shot.surfaces[1, end]
    Z0 = shot.surfaces[2, end]
    ϵ  = shot.surfaces[3, end]
    κ  = shot.surfaces[4, end]
    a = min(1.2 * R0 * ϵ, R0) # 20% bigger than plasma, but a no bigger than R0
    b = κ * a
    Rgrid = range(R0 - a, R0 + a, n_grid)
    Zgrid = range(Z0 - b, Z0 + b, n_grid)

    psirz = zeros(n_grid, n_grid)

    eq2d.grid_type.index = 1
    eq2d.grid.dim1 = Rgrid
    eq2d.grid.dim2 = Zgrid

    if free_boundary
        # convert from fixed to free boundary equilibrium
        bnd = TEQUILA.MXH(shot.surfaces[:,end])
        Rb, _ = bnd()
        psirz .= VacuumFields.fixed2free(shot, Int(ceil(length(Rb)/2)), Rgrid, Zgrid; ψbound=psib)
    else
        for (i, r) in enumerate(Rgrid)
            for (j, z) in enumerate(Zgrid)
                psirz[i, j] = shot(r, z)
            end
        end
    end
    eq2d.psi = psirz .+ psib

    IMAS.flux_surfaces(eqt)
end
