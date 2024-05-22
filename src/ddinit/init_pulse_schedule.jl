#= ======================= =#
#  init pulse_schedule IDS  #
#= ======================= =#
"""
    init_pulse_schedule!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd(); simplify_time_traces::Float64=0.1)

Initialize `dd.pulse_schedule` starting from `ini` and `act` parameters
"""
function init_pulse_schedule!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd(); simplify_time_traces::Float64=0.1)
    TimerOutputs.reset_timer!("init_pulse_schedule")
    TimerOutputs.@timeit timer "init_pulse_schedule" begin
        init_from = ini.general.init_from

        if init_from == :ods
            if IMAS.hasdata(dd1.pulse_schedule, :time) && length(dd1.pulse_schedule.time) > 0
                dd.pulse_schedule = deepcopy(dd1.pulse_schedule)
            else
                init_from = :scalars
            end
        end

        if init_from == :scalars
            # ip
            time, data = get_time_dependent(ini.equilibrium, :ip; simplify_time_traces)
            dd.pulse_schedule.flux_control.time = time
            dd.pulse_schedule.flux_control.i_plasma.reference = data

            # R0 should not be time dependent for definition of B0
            if !isempty(dd.build.layer)
                plasma = IMAS.get_build_layer(dd.build.layer; type=_plasma_)
                R0 = (plasma.start_radius + plasma.end_radius) / 2.0
            elseif typeof(getfield(ini.equilibrium, :R0).value) <: Function
                error("`ini.equilibrium.R0` should not be time dependent")
            else
                mxh = IMAS.MXH(ini, dd1)
                R0 = mxh.R0
            end

            # B0
            time, data = get_time_dependent(ini.equilibrium, :B0; simplify_time_traces)
            dd.pulse_schedule.tf.time = time
            dd.pulse_schedule.tf.b_field_tor_vacuum.reference = data
            dd.pulse_schedule.tf.r0 = R0

            # initialize position_control from mxh
            if ini.equilibrium.boundary_from == :scalars
                time, data = get_time_dependent(ini.equilibrium, [:ip, :R0, :Z0, :ϵ, :κ, :δ, :ζ, :𝚶, :xpoints]; simplify_time_traces)
                for (k, time0) in enumerate(time)
                    R0 = data.R0[k]
                    Z0 = data.Z0[k]
                    ϵ = data.ϵ[k]
                    κ = data.κ[k]
                    δ = data.δ[k]
                    ζ = data.ζ[k]
                    𝚶 = data.𝚶[k]
                    κ = ini_equilibrium_elongation_true(κ, ϵ)
                    mxh = IMAS.MXH(R0, Z0, ϵ, κ, 0.0, [𝚶, 0.0], [asin(δ), -ζ])
                    nx = n_xpoints(data.xpoints[k])
                    mxhb = fitMXHboundary(mxh, nx)
                    init_pulse_schedule_postion_control(dd.pulse_schedule.position_control, mxhb, time0)
                    if k == length(time) - 1 && time[k+1] == Inf
                        init_pulse_schedule_postion_control(dd.pulse_schedule.position_control, mxhb, Inf)
                        break
                    end
                end

            else
                # can handle rampup but not time-dependence after that
                nx = n_xpoints(ini.equilibrium.xpoints)
                mxh = IMAS.MXH(ini, dd1)
                ini.equilibrium(mxh)
                mxhb = fitMXHboundary(mxh, nx)

                if !ismissing(ini.rampup, :ends_at)
                    wall = wall_radii(ini.build.layers, mxh.R0)
                    mxh_bore, mxh_lim2div = limited_to_diverted(.75, mxhb, wall.r_hfs, wall.r_lfs, ini.rampup.side)
                    init_pulse_schedule_postion_control(dd.pulse_schedule.position_control, mxh_bore, -Inf)
                    init_pulse_schedule_postion_control(dd.pulse_schedule.position_control, mxh_bore, 0.0)
                    init_pulse_schedule_postion_control(dd.pulse_schedule.position_control, mxh_lim2div, ini.rampup.diverted_at)
                    init_pulse_schedule_postion_control(dd.pulse_schedule.position_control, mxhb, ini.rampup.ends_at)
                    init_pulse_schedule_postion_control(dd.pulse_schedule.position_control, mxhb, Inf)
                else
                    init_pulse_schedule_postion_control(dd.pulse_schedule.position_control, mxhb, -Inf)
                    init_pulse_schedule_postion_control(dd.pulse_schedule.position_control, mxhb, global_time(ini))
                    init_pulse_schedule_postion_control(dd.pulse_schedule.position_control, mxhb, Inf)
                end
            end

            # density & zeff
            time, data = get_time_dependent(ini.core_profiles, [:zeff, :greenwald_fraction_ped, :ne_ped, :greenwald_fraction]; simplify_time_traces)
            dd.pulse_schedule.density_control.time = time
            dd.pulse_schedule.density_control.zeff.reference = data.zeff
            dd.pulse_schedule.density_control.zeff_pedestal.reference = data.zeff
            if !(typeof(data.greenwald_fraction_ped) <: Vector{Missing})
                dd.pulse_schedule.density_control.n_e_pedestal_greenwald_fraction.reference = data.greenwald_fraction_ped
            end
            if !(typeof(data.ne_ped) <: Vector{Missing})
                dd.pulse_schedule.density_control.n_e_pedestal.reference = data.ne_ped
            end
            if !(typeof(data.greenwald_fraction) <: Vector{Missing})
                dd.pulse_schedule.density_control.n_e_greenwald_fraction.reference = data.greenwald_fraction
            end

            # NB
            resize!(dd.pulse_schedule.nbi.unit, length(ini.nb_unit))
            for (k, ini_nbu) in enumerate(ini.nb_unit)
                time, data = get_time_dependent(ini_nbu, [:power_launched, :rho_0, :width]; simplify_time_traces)
                dd.pulse_schedule.nbi.time = time
                dd.pulse_schedule.nbi.unit[k].power.reference = data.power_launched
                dd.pulse_schedule.nbi.unit[k].deposition_rho_tor_norm.reference = data.rho_0
                dd.pulse_schedule.nbi.unit[k].deposition_rho_tor_norm_width.reference = data.width
            end

            # EC
            resize!(dd.pulse_schedule.ec.beam, length(ini.ec_launcher))
            for (k, ini_ecb) in enumerate(ini.ec_launcher)
                time, data = get_time_dependent(ini_ecb, [:power_launched, :rho_0, :width]; simplify_time_traces)
                dd.pulse_schedule.ec.time = time
                dd.pulse_schedule.ec.beam[k].power_launched.reference = data.power_launched
                dd.pulse_schedule.ec.beam[k].deposition_rho_tor_norm.reference = data.rho_0
                dd.pulse_schedule.ec.beam[k].deposition_rho_tor_norm_width.reference = data.width
            end

            # IC
            resize!(dd.pulse_schedule.ic.antenna, length(ini.ic_antenna))
            for (k, ini_ica) in enumerate(ini.ic_antenna)
                time, data = get_time_dependent(ini_ica, [:power_launched, :rho_0, :width]; simplify_time_traces)
                dd.pulse_schedule.ic.time = time
                dd.pulse_schedule.ic.antenna[k].power.reference = data.power_launched
                dd.pulse_schedule.ic.antenna[k].deposition_rho_tor_norm.reference = data.rho_0
                dd.pulse_schedule.ic.antenna[k].deposition_rho_tor_norm_width.reference = data.width
            end

            # LH
            resize!(dd.pulse_schedule.lh.antenna, length(ini.lh_antenna))
            for (k, ini_lha) in enumerate(ini.lh_antenna)
                time, data = get_time_dependent(ini_lha, [:power_launched, :rho_0, :width]; simplify_time_traces)
                dd.pulse_schedule.lh.time = time
                dd.pulse_schedule.lh.antenna[k].power.reference = data.power_launched
                dd.pulse_schedule.lh.antenna[k].deposition_rho_tor_norm.reference = data.rho_0
                dd.pulse_schedule.lh.antenna[k].deposition_rho_tor_norm_width.reference = data.width
            end

            # PL
            resize!(dd.pulse_schedule.pellet.launcher, length(ini.pellet_launcher))
            for (k, ini_peln) in enumerate(ini.pellet_launcher)
                time, data = get_time_dependent(ini_peln, [:frequency, :rho_0, :width]; simplify_time_traces)
                dd.pulse_schedule.pellet.time = time
                dd.pulse_schedule.pellet.launcher[k].frequency.reference = data.frequency
                dd.pulse_schedule.pellet.launcher[k].deposition_rho_tor_norm.reference = data.rho_0
                dd.pulse_schedule.pellet.launcher[k].deposition_rho_tor_norm_width.reference = data.width
            end
        end

        return dd
    end
end

function get_time_dependent(par::AbstractParameters, field::Symbol; simplify_time_traces::Float64)
    @assert 0.0 <= simplify_time_traces <= 1.0 "get_time_dependent() simplify_time_traces must be between [0,1]"

    value = getfield(par, field).value

    # if it is a time dependent quantity
    if typeof(value) <: Function
        time = collect(SimulationParameters.time_range(par))
        data = value.(time)
        if !(eltype(data) <: Number)
            data, mapping = SimulationParameters.encode_array(data)
            if simplify_time_traces != 0.0
                time, data = IMAS.simplify_2d_path(time, Float64.(data), simplify_time_traces)
            end
            data = [mapping[Int(d)] for d in data]
        else
            if simplify_time_traces != 0.0
                time, data = IMAS.simplify_2d_path(time, data, simplify_time_traces)
            end
        end

        # if it is a constant
    else
        if simplify_time_traces != 0.0 || isempty(SimulationParameters.time_range(par))
            time = Float64[-Inf, SimulationParameters.global_time(par), Inf]
            data = [value, value, value]
        else
            time = SimulationParameters.time_range(par)
            data = fill(value, length(time))
        end
    end

    return time, data
end

function get_time_dependent(par::AbstractParameters, fields::Vector{Symbol}; simplify_time_traces::Float64)
    all_times = Float64[]

    for field in fields
        time, data = get_time_dependent(par, field; simplify_time_traces)
        append!(all_times, time)
    end

    all_times = sort!(unique(all_times))

    return all_times, NamedTuple{Tuple(fields)}([get_time_dependent(par, field, all_times) for field in fields])
end

function get_time_dependent(par::AbstractParameters, field::Symbol, all_times::Vector{Float64})
    value = getfield(par, field).value

    if typeof(value) <: Function
        time = collect(SimulationParameters.time_range(par))
        data = value.(time)
        if !(eltype(data) <: Number)
            data, mapping = SimulationParameters.encode_array(data)
            all_data = [mapping[Int(d)] for d in IMAS.interp1d(time, Float64.(data), :constant).(all_times)]
        else
            all_data = IMAS.interp1d(time, data, :constant).(all_times)
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
    # MXHboundary adds x-points
    pr = mxhb.r_boundary
    pz = mxhb.z_boundary

    # X-point information
    # NOTE: always fill for both X-points and set things to NaN if no X-point
    # NOTE: upper X-point always in first slot, lower X-point in second slot
    resize!(pc.x_point, 2; wipe=false)
    rxu = rxl = zxu = zxl = NaN
    if length(mxhb.RX) == 0
        # pass
    elseif length(mxhb.RX) == 1
        if mxhb.ZX[1] > mxhb.mxh.Z0
            rxu = mxhb.RX[1]
            zxu = mxhb.ZX[1]
        else
            rxl = mxhb.RX[1]
            zxl = mxhb.ZX[1]
        end
    elseif length(mxhb.RX) == 2
        if mxhb.ZX[1] > mxhb.mxh.Z0
            rxu = mxhb.RX[1]
            zxu = mxhb.ZX[1]
        else
            rxl = mxhb.RX[1]
            zxl = mxhb.ZX[1]
        end
        if mxhb.ZX[2] > mxhb.mxh.Z0
            rxu = mxhb.RX[2]
            zxu = mxhb.ZX[2]
        else
            rxl = mxhb.RX[2]
            zxl = mxhb.ZX[2]
        end
    else
        error("cannot handle more than two X-points")
    end
    IMAS.set_time_array(pc.x_point[1].r, :reference, time0, rxu)
    IMAS.set_time_array(pc.x_point[1].z, :reference, time0, zxu)
    IMAS.set_time_array(pc.x_point[2].r, :reference, time0, rxl)
    IMAS.set_time_array(pc.x_point[2].z, :reference, time0, zxl)

    # boundary with x-points parametrized with MXH
    mxh = IMAS.MXH(pr, pz, 2)
    pr, pz = IMAS.resample_plasma_boundary(pr, pz; n_points=100)
    IMAS.reorder_flux_surface!(pr, pz, argmax(pz))

    # scalars
    IMAS.set_time_array(pc.minor_radius, :reference, time0, mxh.ϵ * mxh.R0)
    IMAS.set_time_array(pc.geometric_axis.r, :reference, time0, mxh.R0)
    IMAS.set_time_array(pc.geometric_axis.z, :reference, time0, mxh.Z0)
    IMAS.set_time_array(pc.elongation, :reference, time0, mxh.κ)
    IMAS.set_time_array(pc.triangularity, :reference, time0, sin(mxh.s[1]))
    IMAS.set_time_array(pc.squareness, :reference, time0, -mxh.s[2])
    IMAS.set_time_array(pc.ovality, :reference, time0, mxh.c[1])

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
    mxhb_bore = FUSE.MXHboundary(mxh_bore; upper_x_point=false, lower_x_point=false)

    # diverted to limited shape
    r = [r_bore; r_diverted]
    z = [z_bore; z_diverted]
    hull = convex_hull(r, z; closed_polygon=true)
    r = [x for (x, y) in hull]
    z = [y for (x, y) in hull]
    mxh_lim2div = IMAS.MXH(r, z, length(mxh_diverted.c))
    mxhb_lim2div = FUSE.MXHboundary(mxh_lim2div; upper_x_point=any(mxhb_diverted.ZX .> mxhb_diverted.mxh.Z0), lower_x_point=any(mxhb_diverted.ZX .< mxhb_diverted.mxh.Z0))

    return (mxhb_bore=mxhb_bore, mxhb_lim2div=mxhb_lim2div)
end
