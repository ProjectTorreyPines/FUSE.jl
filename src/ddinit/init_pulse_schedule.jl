#= ======================= =#
#  init pulse_schedule IDS  #
#= ======================= =#
"""
    init_pulse_schedule!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())

Initialize `dd.pulse_schedule` starting from `ini` and `act` parameters
"""
function init_pulse_schedule!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())
    TimerOutputs.reset_timer!("init_pulse_schedule")
    TimerOutputs.@timeit timer "init_pulse_schedule" begin
        init_from = ini.general.init_from

        ps = dd.pulse_schedule
        ps1 = dd1.pulse_schedule

        # ip
        if init_from == :ods && IMAS.hasdata(ps1.flux_control)
            ps.flux_control = deepcopy(ps1.flux_control)
        else
            time, data = get_time_dependent(ini.equilibrium, :ip)
            ps.flux_control.time = time
            ps.flux_control.i_plasma.reference = data
        end

        # B0
        if init_from == :ods && IMAS.hasdata(ps1.tf)
            ps.tf = deepcopy(ps1.tf)
        else
            # R0 should not be time dependent for definition of B0
            if !isempty(dd.build.layer)
                plasma = IMAS.get_build_layer(dd.build.layer; type=_plasma_)
                R0 = (plasma.start_radius + plasma.end_radius) / 2.0
            elseif typeof(getfield(ini.equilibrium, :R0).value) <: Union{Function,TimeData}
                error("`ini.equilibrium.R0` should not be time dependent")
            else
                mxhb = MXHboundary(ini, dd1)
                R0 = mxhb.mxh.R0
            end

            # B0
            time, data = get_time_dependent(ini.equilibrium, :B0)
            ps.tf.time = time
            ps.tf.b_field_tor_vacuum.reference = data
            ps.tf.r0 = R0
        end

        # position_control
        if init_from == :ods && IMAS.hasdata(ps1.position_control)
            ps.position_control = deepcopy(ps1.position_control)
        else
            time, _ = get_time_dependent(ini.equilibrium, [:R0, :Z0, :ϵ, :κ, :δ, :ζ, :tilt, :𝚶, :xpoints, :MXH_params, :rz_points])
            if !ismissing(ini.rampup, :ends_at)
                time = filter(t -> t > ini.rampup.ends_at, time)
                pushfirst!(time, ini.rampup.ends_at)
                if !ismissing(ini.rampup, :diverted_at)
                    pushfirst!(time, ini.rampup.diverted_at)
                end
                pushfirst!(time, 0.0)
                pushfirst!(time, -Inf)
            end
            ini_time_simulation_start = ini.time.simulation_start
            dd1_time_backup = dd1.global_time
            for (k, time0) in enumerate(time)
                if !ismissing(ini.time, :pulse_shedule_time_basis) && time0 < ini.time.pulse_shedule_time_basis[1]
                    ini.time.simulation_start = ini.time.pulse_shedule_time_basis[1]
                    dd1.global_time = ini.time.pulse_shedule_time_basis[1]
                elseif ismissing(ini.time, :pulse_shedule_time_basis) && time0 < ini_time_simulation_start
                    # This is necessary because equilibrium quantities may not be defined at < simulation_start as it happens for example sometimes when starting from ODS
                    ini.time.simulation_start = ini_time_simulation_start
                    dd1.global_time = ini_time_simulation_start
                else
                    ini.time.simulation_start = time0
                    dd1.global_time = time0
                end

                # get MXHboundary representation
                mxhb = MXHboundary(ini, dd1)

                if ismissing(ini.rampup, :ends_at)
                    init_pulse_schedule_postion_control(ps.position_control, mxhb, time0)
                else
                    if isempty(dd.wall)
                        wr = wall_radii(mxhb.mxh.R0, mxhb.mxh.minor_radius, ini.build.plasma_gap)
                    else
                        wr = wall_radii(IMAS.first_wall(dd.wall)..., mxhb.mxh.Z0)
                    end
                    mxh_bore, mxh_lim2div = limited_to_diverted(0.75, mxhb, wr.r_hfs, wr.r_lfs, ini.rampup.side)
                    if time0 <= 0.0
                        init_pulse_schedule_postion_control(ps.position_control, mxh_bore, time0)
                    elseif time0 == ini.rampup.diverted_at
                        init_pulse_schedule_postion_control(ps.position_control, mxh_lim2div, ini.rampup.diverted_at)
                    else
                        init_pulse_schedule_postion_control(ps.position_control, mxhb, time0)
                    end
                end
                if k == length(time) - 1 && time[k+1] == Inf
                    init_pulse_schedule_postion_control(ps.position_control, mxhb, Inf)
                    break
                end
            end
            ini.time.simulation_start = ini_time_simulation_start
            dd1.global_time = dd1_time_backup
        end

        # density & zeff
        if init_from == :ods && IMAS.hasdata(ps1.density_control)
            ps.density_control = deepcopy(ps1.density_control)
        else
            time, data = get_time_dependent(ini.core_profiles, [:zeff, :ne_value])
            dd.pulse_schedule.density_control.time = time
            dd.pulse_schedule.density_control.zeff.reference = data.zeff
            dd.pulse_schedule.density_control.zeff_pedestal.reference = data.zeff

            if ini.core_profiles.ne_setting == :greenwald_fraction_ped
                dd.pulse_schedule.density_control.n_e_pedestal_greenwald_fraction.reference = data.ne_value
            elseif ini.core_profiles.ne_setting == :ne_ped
                dd.pulse_schedule.density_control.n_e_pedestal.reference = data.ne_value
            elseif ini.core_profiles.ne_setting == :ne_line
                dd.pulse_schedule.density_control.n_e_line.reference = data.ne_value
            elseif ini.core_profiles.ne_setting == :greenwald_fraction
                dd.pulse_schedule.density_control.n_e_greenwald_fraction.reference = data.ne_value
            end
        end

        if ini.core_profiles.ne_setting in (:ne_ped, :greenwald_fraction_ped)
            @assert act.ActorPedestal.density_match == :ne_ped "ini.core_profiles.ne_setting=:$(ini.core_profiles.ne_setting) requires act.ActorPedestal.density_match=:ne_ped"
        else
            @assert act.ActorPedestal.density_match == :ne_line "ini.core_profiles.ne_setting=:$(ini.core_profiles.ne_setting) requires act.ActorPedestal.density_match=:ne_line"
        end

        # EC
        if init_from == :ods && IMAS.hasdata(ps1.ec)
            ps.ec = deepcopy(ps1.ec)
        else
            resize!(ps.ec.beam, length(ini.ec_launcher))
            ps.ec.time, powers_launched = get_time_dependent(ini.ec_launcher, :power_launched)
            for k in eachindex(ini.ec_launcher)
                ps.ec.beam[k].power_launched.reference = powers_launched[k]
            end
        end

        # IC
        if init_from == :ods && IMAS.hasdata(ps1.ic)
            ps.ic = deepcopy(ps1.ic)
        else
            resize!(ps.ic.antenna, length(ini.ic_antenna))
            ps.ic.time, powers_launched = get_time_dependent(ini.ic_antenna, :power_launched)
            for k in eachindex(ini.ic_antenna)
                ps.ic.antenna[k].power.reference = powers_launched[k]
            end
        end

        # LH
        if init_from == :ods && IMAS.hasdata(ps1.lh)
            ps.lh = deepcopy(ps1.lh)
        else
            resize!(ps.lh.antenna, length(ini.lh_antenna))
            ps.lh.time, powers_launched = get_time_dependent(ini.lh_antenna, :power_launched)
            for k in eachindex(ini.lh_antenna)
                ps.lh.antenna[k].power.reference = powers_launched[k]
            end
        end

        # NB
        if init_from == :ods && IMAS.hasdata(ps1.nbi)
            ps.nbi = deepcopy(ps1.nbi)
        else
            resize!(ps.nbi.unit, length(ini.nb_unit))
            ps.nbi.time, powers_launched = get_time_dependent(ini.nb_unit, :power_launched)
            energies = [nb_unit.beam_energy for nb_unit in ini.nb_unit]
            for k in eachindex(ini.nb_unit)
                ps.nbi.unit[k].power.reference = powers_launched[k]
                ps.nbi.unit[k].energy.reference = fill(energies[k], size(ps.nbi.time))
            end
        end

        # PL
        if init_from == :ods && IMAS.hasdata(ps1.pellet)
            ps.pellet = deepcopy(ps1.pellet)
        else
            resize!(ps.pellet.launcher, length(ini.pellet_launcher))
            ps.pellet.time, frequencies = get_time_dependent(ini.pellet_launcher, :frequency)
            for k in eachindex(ini.pellet_launcher)
                ps.pellet.launcher[k].frequency.reference = frequencies[k]
            end
        end

        return dd
    end
end

function get_time_dependent(par::AbstractParameters, field::Symbol)
    value = getfield(par, field).value

    # if it is a time dependent quantity
    if typeof(value) <: TimeData
        time = value.time
        data = value.data
    elseif typeof(value) <: Function
        time = time_range = collect(SimulationParameters.time_range(par))
        data = value.(time_range)
        if !(eltype(data) <: Number)
            data, mapping = SimulationParameters.encode_array(data)
            data = [mapping[Int(d)] for d in data]
        end

        # if it is a constant
    else
        time = Float64[-Inf, SimulationParameters.global_time(par), Inf]
        data = [value, value, value]
    end

    return time, data
end

function get_time_dependent(pars::ParametersVector, field::Symbol)
    all_times = Float64[]

    for par in pars
        time, data = get_time_dependent(par, field)
        append!(all_times, time)
        all_times = sort!(unique(all_times))
    end

    return all_times, [get_time_dependent(par, field, all_times) for par in pars]
end

function get_time_dependent(par::AbstractParameters, fields::Vector{Symbol})
    all_times = Float64[]

    for field in fields
        time, data = get_time_dependent(par, field)
        append!(all_times, time)
        all_times = sort!(unique(all_times))
    end

    return all_times, NamedTuple{Tuple(fields)}([get_time_dependent(par, field, all_times) for field in fields])
end

function get_time_dependent(par::AbstractParameters, field::Symbol, all_times::Vector{Float64})
    value = getfield(par, field).value

    if typeof(value) <: TimeData && value.time == all_times
        all_data = value.data
    elseif typeof(value) <: Union{Function,TimeData}
        all_data = data = value.(all_times)
        if !(eltype(data) <: Number)
            data, mapping = SimulationParameters.encode_array(data)
            all_data = [mapping[d] for d in Int.(data)]
        end
    else
        all_data = fill(value, size(all_times))
    end

    return all_data
end

"""
    init_pulse_schedule_postion_control(pc::IMAS.pulse_schedule__position_control, mxhb::MXHboundary, time0::Float64)

Initialize pulse_schedule.postion_control based on MXH boundary and number of x_points
"""
function init_pulse_schedule_postion_control(pc::IMAS.pulse_schedule__position_control, mxhb::MXHboundary, time0::Float64)
    # X-point information
    # NOTE: always fill for both X-points and set things to NaN if no X-point
    # NOTE: upper X-point always in first slot, lower X-point in second slot
    resize!(pc.x_point, 2; wipe=false)
    rxu = rxl = zxu = zxl = NaN
    for k in eachindex(mxhb.ZX)
        if mxhb.ZX[k] > mxhb.mxh.Z0
            rxu = mxhb.RX[k]
            zxu = mxhb.ZX[k]
        else
            rxl = mxhb.RX[k]
            zxl = mxhb.ZX[k]
        end
    end
    IMAS.set_time_array(pc.x_point[1].r, :reference, time0, rxu)
    IMAS.set_time_array(pc.x_point[1].z, :reference, time0, zxu)
    IMAS.set_time_array(pc.x_point[2].r, :reference, time0, rxl)
    IMAS.set_time_array(pc.x_point[2].z, :reference, time0, zxl)

    # from MXHboundary to pr,pz
    pr, pz = IMAS.resample_plasma_boundary(mxhb.r_boundary, mxhb.z_boundary; n_points=100)
    if isempty(mxhb.RX)
        # boundary without x-points
        mxh = mxhb.mxh
    else
        # boundary with x-points parametrized with MXH
        mxh = IMAS.MXH(pr, pz, 2)
    end

    # scalars from mxh
    IMAS.set_time_array(pc.minor_radius, :reference, time0, mxh.minor_radius)
    IMAS.set_time_array(pc.geometric_axis.r, :reference, time0, mxh.R0)
    IMAS.set_time_array(pc.geometric_axis.z, :reference, time0, mxh.Z0)
    IMAS.set_time_array(pc.elongation, :reference, time0, mxh.elongation)
    IMAS.set_time_array(pc.tilt, :reference, time0, mxh.tilt)
    IMAS.set_time_array(pc.triangularity, :reference, time0, mxh.triangularity)
    IMAS.set_time_array(pc.squareness, :reference, time0, mxh.squareness)
    IMAS.set_time_array(pc.ovality, :reference, time0, mxh.ovality)
    IMAS.set_time_array(pc.twist, :reference, time0, mxh.twist)

    # boundary
    resize!(pc.boundary_outline, length(pr); wipe=false)
    for (k, (r, z)) in enumerate(zip(pr, pz))
        IMAS.set_time_array(pc.boundary_outline[k].r, :reference, time0, r)
        IMAS.set_time_array(pc.boundary_outline[k].z, :reference, time0, z)
    end

    return pc
end

"""
    limited_to_diverted(
        initial_minor_radius_fraction::Float64,
        mxhb_diverted::MXHboundary,
        r_inner_wall::Float64,
        r_outer_wall::Float64,
        wall_side::Symbol)

Generates starting circular boundary and transition limited-to-diverted mxh boundaries
"""
function limited_to_diverted(
    initial_minor_radius_fraction::Float64,
    mxhb_diverted::MXHboundary,
    r_inner_wall::Float64,
    r_outer_wall::Float64,
    wall_side::Symbol)

    @assert 0.0 < initial_minor_radius_fraction <= 1.0
    @assert wall_side in (:hfs, :lfs)

    # final diverted shape
    mxh_diverted = mxhb_diverted.mxh
    r_diverted, z_diverted = mxh_diverted()
    δ_diverted = sin(mxh_diverted.s[1])
    if δ_diverted > 0.0
        index_z_height = argmax(r_diverted)
    else
        index_z_height = argmin(r_diverted)
    end

    # starting circular shape
    t = range(0, -2pi, length(r_diverted))
    if wall_side == :lfs
        a = (r_outer_wall - minimum(r_diverted)) / 2.0
    else
        a = (maximum(r_diverted) - r_inner_wall) / 2.0
    end
    initial_minor_radius = a * initial_minor_radius_fraction
    r_bore = initial_minor_radius * cos.(t)
    z_bore = initial_minor_radius * sin.(t)
    if wall_side == :lfs
        r_bore .+= r_outer_wall .- maximum(r_bore)
        z_bore .+= z_diverted[index_z_height]
    else
        r_bore .+= r_inner_wall .- minimum(r_bore)
        z_bore .+= z_diverted[index_z_height]
    end
    mxh_bore = IMAS.MXH(r_bore, z_bore, 0)
    mxhb_bore = MXHboundary(mxh_bore; upper_x_point=false, lower_x_point=false)
    mxhb_bore.RX = deepcopy(mxhb_diverted.RX)
    mxhb_bore.ZX = deepcopy(mxhb_diverted.ZX)

    # diverted to limited shape
    r = [r_bore; r_diverted]
    z = [z_bore; z_diverted]
    hull = IMAS.convex_hull(r, z; closed_polygon=true)
    r = [x for (x, y) in hull]
    z = [y for (x, y) in hull]
    mxh_lim2div = IMAS.MXH(r, z, length(mxh_diverted.c))
    mxhb_lim2div = MXHboundary(mxh_lim2div; upper_x_point=false, lower_x_point=false)
    mxhb_lim2div.RX = deepcopy(mxhb_diverted.RX)
    mxhb_lim2div.ZX = deepcopy(mxhb_diverted.ZX)

    return (mxhb_bore=mxhb_bore, mxhb_lim2div=mxhb_lim2div)
end
