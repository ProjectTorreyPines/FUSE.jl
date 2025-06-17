#= ==================== =#
#  init equilibrium IDS  #
#= ==================== =#
"""
    init_equilibrium!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())

Initialize `dd.equilibrium` starting from `ini` and `act` parameters
"""
function init_equilibrium!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())
    TimerOutputs.reset_timer!("init_equilibrium")
    TimerOutputs.@timeit timer "init_equilibrium" begin
        init_from = ini.general.init_from

        if init_from == :ods
            if !isempty(dd1.equilibrium.time_slice) && ini.equilibrium.boundary_from == :ods
                dd.equilibrium = deepcopy(dd1.equilibrium)
                eqt = dd.equilibrium.time_slice[]
                fw = IMAS.first_wall(dd.wall)
                if findfirst(:rectangular, eqt.profiles_2d) !== nothing
                    IMAS.flux_surfaces(eqt, fw.r, fw.z)
                end
            else
                init_from = :scalars
            end
        end

        # the pressure and j_tor to be used by equilibrium solver need to be set in dd.core_profiles
        # while the equilibrium boundary neds to be set in init_pulse_schedule
        if isempty(dd.core_profiles.profiles_1d)
            cp1d = resize!(dd.core_profiles.profiles_1d)
            if init_from == :ods
                # take p and j from input equilibrium ods
                cp1d.grid.rho_tor_norm = eqt.profiles_1d.rho_tor_norm
                cp1d.grid.psi = eqt.profiles_1d.psi
                cp1d.j_tor = eqt.profiles_1d.j_tor
                cp1d.pressure = eqt.profiles_1d.pressure
            else
                # guess pressure and j_tor from input current and peak pressure
                rhon = range(0.0, 1.0, 129)
                psin = rhon .^ 2
                cp1d.grid.rho_tor_norm = rhon
                cp1d.grid.psi = psin
                cp1d.j_tor = ini.equilibrium.ip .* (1.0 .- psin) ./ @ddtime(dd.pulse_schedule.position_control.geometric_axis.r.reference)

                # estimate ini.equilibrium.pressure_core
                if ismissing(ini.equilibrium, :pressure_core)
                    dd0 = IMAS.dd()
                    init_core_profiles!(dd0, ini, act, dd1)
                    ini.equilibrium.pressure_core = dd0.core_profiles.profiles_1d[].pressure[1]
                end

                cp1d.pressure = ini.equilibrium.pressure_core .* (1.0 .- psin).^2
            end
        end

        # solve equilibrium
        if !(init_from == :ods && ini.equilibrium.boundary_from == :ods) || findfirst(:rectangular, eqt.profiles_2d) === nothing
            act_copy = deepcopy(act)
            act_copy.ActorCHEASE.rescale_eq_to_ip = true
            ActorEquilibrium(dd, act_copy; ip_from=:pulse_schedule)
        end

        # field null surface
        if ini.equilibrium.field_null_surface > 0.0
            field_null_surface!(dd.pulse_schedule.position_control, dd.equilibrium, ini.equilibrium.field_null_surface)
        end

        return dd
    end
end

"""
    field_null_surface!(pc::IMAS.pulse_schedule__position_control, eq::IMAS.equilibrium, scale::Real=0.5, ψp_constant::Real=0.1)

Setup field null surface as pulse_schedule.position_control.boundary_outline and insert equilibrium time slice at time=-Inf
"""
function field_null_surface!(pc::IMAS.pulse_schedule__position_control, eq::IMAS.equilibrium, scale::Real=0.75, ψp_constant::Real=0.1)
    eqt = eq.time_slice[]

    # get coordinates for flux-null boundary at t=-Inf
    mxh = IMAS.MXH(eqt.boundary.outline.r, eqt.boundary.outline.z, 0)
    mxh.κ = 1.0
    mxh.c0 = 0.0
    mxh.ϵ *= scale
    pr, pz = mxh(length(pc.boundary_outline); adaptive=false)

    # insert flux-null boundary at t=-Inf
    if !isempty(pc.boundary_outline) && pc.time[1] == -Inf
        for (k, (r, z)) in enumerate(zip(pr, pz))
            pc.boundary_outline[k].r.reference[1] = r
            pc.boundary_outline[k].z.reference[1] = z
        end
    else
        pushfirst!(pc.time, -Inf)
        for (k, (r, z)) in enumerate(zip(pr, pz))
            pushfirst!(pc.boundary_outline[k].r.reference, r)
            pushfirst!(pc.boundary_outline[k].z.reference, z)
        end
    end

    # insert an equilibrium time slice at t=-Inf
    if !isempty(eq.time_slice) && eq.time[1] == -Inf
        eqb = empty!(eq.time_slice[1])
    else
        eqb = IMAS.equilibrium__time_slice()
        pushfirst!(eq.time_slice, eqb)
        pushfirst!(eq.time, -Inf)
        pushfirst!(eq.vacuum_toroidal_field.b0, 0.0)
    end
    eqb.time = -Inf

    # set B0 and psi_boundary for equilibrium time slice at t=-Inf
    eq.vacuum_toroidal_field.b0[1] = @ddtime(eq.vacuum_toroidal_field.b0)
    eqb.global_quantities.ip = 0.0
    eqb.global_quantities.psi_boundary = ψp_constant
    eqb.profiles_1d.psi = [ψp_constant]
    eqb.boundary.outline.r = pr
    eqb.boundary.outline.z = pz
    eqb.global_quantities.magnetic_axis.r = mxh.R0
    eqb.global_quantities.magnetic_axis.z = mxh.Z0

    return nothing
end
