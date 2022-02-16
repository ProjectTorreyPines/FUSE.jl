#= ==================== =#
#  init equilibrium IDS  #
#= ==================== =#

"""
    function init_equilibrium(
        eq::IMAS.equilibrium;
        B0::Real,
        R0::Real,
        Z0::Real,
        ϵ::Real,
        κ::Real,
        δ::Real,
        βn::Real,
        ip::Real,
        x_point::Union{Vector,NTuple{2},Bool} = false)

Initialize equilibrium IDS based on some basic Miller geometry parameters
"""
function init_equilibrium(
    eq::IMAS.equilibrium;
    B0::Real,
    R0::Real,
    Z0::Real,
    ϵ::Real,
    κ::Real,
    δ::Real,
    βn::Real,
    ip::Real,
    x_point::Union{Vector,NTuple{2},Bool} = false)

    eqt = resize!(eq.time_slice)
    eqt.boundary.minor_radius = ϵ * R0
    eqt.boundary.geometric_axis.r = R0
    eqt.boundary.geometric_axis.z = Z0
    eqt.boundary.elongation = κ
    eqt.boundary.triangularity = δ
    eqt.global_quantities.ip = ip
    eqt.global_quantities.beta_normal = βn
    if x_point === true
        x_point = (R0 * (1 - 1.1 * δ * ϵ), -R0 * 1.1 * κ * ϵ)
    end
    if isa(x_point, Union{Vector,Tuple})
        resize!(eqt.boundary.x_point, 1)
        eqt.boundary.x_point[1].r = x_point[1]
        eqt.boundary.x_point[1].z = x_point[2]
    end
    eq.vacuum_toroidal_field.r0 = R0
    @ddtime eq.vacuum_toroidal_field.b0 = B0

    return eq
end

function init_equilibrium(dd::IMAS.dd, gasc::GASC)
    gasc = gasc.solution

    R0 = gasc["INPUTS"]["radial build"]["majorRadius"]
    Z0 = 0.0
    ϵ = 1 / gasc["INPUTS"]["radial build"]["aspectRatio"]
    κ = gasc["OUTPUTS"]["plasma parameters"]["elongation"]
    δ = gasc["INPUTS"]["plasma parameters"]["triangularity"]
    B0 = gasc["INPUTS"]["conductors"]["magneticFieldOnAxis"]
    ip = gasc["INPUTS"]["plasma parameters"]["plasmaCurrent"] * 1E6
    βn = gasc["OUTPUTS"]["plasma parameters"]["betaN"]
    x_point = true
    symmetric = false
    resolution = 129

    # initialize dd
    init_equilibrium(dd.equilibrium; B0, R0, Z0, ϵ, κ, δ, βn = βn, ip, x_point = x_point)

    # equilibrium solver
    eqactor = SolovevEquilibriumActor(dd, symmetric = symmetric)
    step(eqactor, verbose = false)
    finalize(eqactor, resolution, (maximum([R0 * (1 - ϵ * 2), 0.0]), R0 * (1 + ϵ * 2)), (-R0 * ϵ * κ * 1.5, R0 * ϵ * κ * 1.5))

    return dd
end

function init_equilibrium(dd::IMAS.dd, par::Parameters)
    if par.general.init_from == :scalars
        # init equilibrium
        init_equilibrium(
            dd.equilibrium;
            B0 = par.equilibrium.B0,
            R0 = par.equilibrium.R0,
            Z0 = par.equilibrium.Z0,
            ϵ = par.equilibrium.ϵ,
            κ = par.equilibrium.κ,
            δ = par.equilibrium.δ,
            βn = par.equilibrium.βn,
            ip = par.equilibrium.ip,
            x_point = par.equilibrium.x_point)

        # equilibrium
        eqactor = SolovevEquilibriumActor(dd, symmetric = par.equilibrium.symmetric)
        step(eqactor, verbose = false)
        finalize(eqactor, par.equilibrium.ngrid)

    elseif par.general.init_from == :ods
        dd1 = IMAS.json2imas(par.ods.filename)
        dd1.equilibrium.time = [t for t in dd1.equilibrium.time]
        dd.global_time = dd1.global_time = dd1.equilibrium.time[end]
        IMAS.flux_surfaces(dd1.equilibrium)
        dd.equilibrium = dd1.equilibrium

    elseif par.general.init_from == :gasc
        gasc = GASC(par.gasc.filename, par.gasc.case)
        init_equilibrium(dd, gasc)
    end

    # field null surface
    if par.equilibrium.field_null_surface > 0.0
        pushfirst!(dd.equilibrium.time_slice, field_null_surface(dd.equilibrium.time_slice[], par.equilibrium.field_null_surface))
        pushfirst!(dd.equilibrium.vacuum_toroidal_field.b0, @ddtime(dd.equilibrium.vacuum_toroidal_field.b0))
        pushfirst!(dd.equilibrium.time, -1.0)
        dd.equilibrium.time_slice[1].time = -1.0
    end

    return dd
end

"""
    field_null_surface(eqt, scale = 0.25, abs_psi_boundary = 0.1)

Return field null surface by scaling an existing equilibrium time_slice
"""
function field_null_surface(eqt::IMAS.equilibrium__time_slice, scale::Real = 0.25, abs_psi_boundary::Real = 0.1)
    eqb = IMAS.equilibrium__time_slice()
    eqb.global_quantities.psi_boundary = sign(eqt.profiles_1d.psi[1] - eqt.profiles_1d.psi[end]) * abs_psi_boundary
    eqb.boundary.outline.r, eqb.boundary.outline.z, _ = IMAS.flux_surface(eqt, eqt.profiles_1d.psi[1] * (1 - scale) + eqt.profiles_1d.psi[end] * scale)
    eqb.boundary.outline.r .-= minimum(eqb.boundary.outline.r) .- minimum(IMAS.flux_surface(eqt, eqt.profiles_1d.psi[end])[1])
    eqb.profiles_1d.psi = [eqb.global_quantities.psi_boundary]
    eqb.profiles_1d.f = [eqt.profiles_1d.f[end]]
    return eqb
end