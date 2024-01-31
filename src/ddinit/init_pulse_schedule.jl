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
            if !ismissing(dd1.pulse_schedule, :time) && length(dd1.pulse_schedule.time) > 0
                dd.pulse_schedule = dd1.pulse_schedule
            else
                init_from = :scalars
            end
        end

        if init_from == :scalars
            time, data = get_time_dependent(ini.equilibrium, :ip; simplify_time_traces)
            dd.pulse_schedule.flux_control.i_plasma.reference.time = time
            dd.pulse_schedule.flux_control.i_plasma.reference.data = data

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

            time, data = get_time_dependent(ini.equilibrium, :B0; simplify_time_traces)
            dd.pulse_schedule.tf.b_field_tor_vacuum_r.reference.time = time
            dd.pulse_schedule.tf.b_field_tor_vacuum_r.reference.data = data .* R0

            # initialize position_control from mxh
            if ini.equilibrium.boundary_from == :scalars
                all_times = Float64[]
                for shape_parameter in [:B0, :ip, :R0, :Z0, :ϵ, :κ, :δ, :ζ, :𝚶, :xpoints]
                    time, data = get_time_dependent(ini.equilibrium, shape_parameter; simplify_time_traces)
                    append!(all_times, time)
                end
                all_times = sort!(unique(all_times))

                simulation_start_bkp = ini.time.simulation_start
                try
                    for (k, time) in enumerate(all_times)
                        ini.time.simulation_start = time
                        R0 = ini.equilibrium.R0
                        Z0 = ini.equilibrium.Z0
                        ϵ = ini.equilibrium.ϵ
                        κ = ini.equilibrium.κ
                        δ = ini.equilibrium.δ
                        ζ = ini.equilibrium.ζ
                        𝚶 = ini.equilibrium.𝚶
                        κ = ini_equilibrium_elongation_true(κ, ϵ)
                        mxh = IMAS.MXH(R0, Z0, ϵ, κ, 0.0, [𝚶, 0.0], [asin(δ), -ζ])
                        nx = n_xpoints(ini.equilibrium.xpoints)
                        mxhb = fitMXHboundary(mxh, nx)
                        init_pulse_schedule_postion_control(dd.pulse_schedule.position_control, mxhb, time)
                        if k == length(all_times) - 1 && all_times[k+1] == Inf
                            init_pulse_schedule_postion_control(dd.pulse_schedule.position_control, mxhb, Inf)
                            break
                        end
                    end
                finally
                    ini.time.simulation_start = simulation_start_bkp
                end

            else
                # NOT SETUP FOR TIME DEPENDENCE YET
                nx = n_xpoints(ini.equilibrium.xpoints) # XPOINTS NOT SETUP FOR TIME DEPENDENCE YET
                mxh = IMAS.MXH(ini, dd1)
                ini.equilibrium(mxh)
                mxhb = fitMXHboundary(mxh, nx)

                init_pulse_schedule_postion_control(dd.pulse_schedule.position_control, mxhb, -Inf)
                init_pulse_schedule_postion_control(dd.pulse_schedule.position_control, mxhb, ini.time.simulation_start)
                init_pulse_schedule_postion_control(dd.pulse_schedule.position_control, mxhb, Inf)
            end

            # NB
            resize!(dd.pulse_schedule.nbi.unit, length(ini.nb_unit))
            for (k, ini_nbu) in enumerate(ini.nb_unit)
                time, data = get_time_dependent(ini_nbu, :power_launched; simplify_time_traces)
                dd.pulse_schedule.nbi.unit[k].power.reference.time = time
                dd.pulse_schedule.nbi.unit[k].power.reference.data = data
            end

            # EC
            resize!(dd.pulse_schedule.ec.launcher, length(ini.ec_launcher))
            for (k, ini_ecb) in enumerate(ini.ec_launcher)
                time, data = get_time_dependent(ini_ecb, :power_launched; simplify_time_traces)
                dd.pulse_schedule.ec.launcher[k].power.reference.time = time
                dd.pulse_schedule.ec.launcher[k].power.reference.data = data
            end

            # IC
            resize!(dd.pulse_schedule.ic.antenna, length(ini.ic_antenna))
            for (k, ini_ica) in enumerate(ini.ic_antenna)
                time, data = get_time_dependent(ini_ica, :power_launched; simplify_time_traces)
                dd.pulse_schedule.ic.antenna[k].power.reference.time = time
                dd.pulse_schedule.ic.antenna[k].power.reference.data = data
            end

            # LH
            resize!(dd.pulse_schedule.lh.antenna, length(ini.lh_antenna))
            for (k, ini_lha) in enumerate(ini.lh_antenna)
                time, data = get_time_dependent(ini_lha, :power_launched; simplify_time_traces)
                dd.pulse_schedule.lh.antenna[k].power.reference.time = time
                dd.pulse_schedule.lh.antenna[k].power.reference.data = data
            end

        end

        return dd
    end
end

function get_time_dependent(par::AbstractParameters, field::Symbol; simplify_time_traces::Float64)
    @assert 0.0 <= simplify_time_traces <= 1.0 "get_time_dependent() simplify_time_traces must be between [0,1]"

    value = getfield(par, field).value

    if typeof(value) <: Function
        time = collect(SimulationParameters.top(par).time.pulse_shedule_time_basis)
        data = value.(time)
        if !(eltype(data) <: Number)
            data = Float64.(SimulationParameters.encode_array(data)[1])
        end
        if simplify_time_traces != 0.0
            time, data = IMAS.simplify_2d_path(time, data, simplify_time_traces)
        end
    else
        time = Float64[-Inf, SimulationParameters.top(par).time.simulation_start, Inf]
        data = [value, value, value]
    end
    return time, data
end

"""
    init_pulse_schedule_postion_control(
        pc::IMAS.pulse_schedule__position_control,
        mxhb::MXHboundary,
        time0::Float64)

Initialize pulse_schedule.postion_control based on MXH boundary and number of x_points
"""
function init_pulse_schedule_postion_control(
    pc::IMAS.pulse_schedule__position_control,
    mxhb::MXHboundary,
    time0::Float64)

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
    IMAS.set_time_array(pc.x_point[1].r.reference, :data, time0, rxu)
    IMAS.set_time_array(pc.x_point[1].z.reference, :data, time0, zxu)
    IMAS.set_time_array(pc.x_point[2].r.reference, :data, time0, rxl)
    IMAS.set_time_array(pc.x_point[2].z.reference, :data, time0, zxl)

    # boundary with x-points parametrized with MXH
    mxh = IMAS.MXH(pr, pz, 2)
    pr, pz = IMAS.resample_plasma_boundary(pr, pz; n_points=100)
    IMAS.reorder_flux_surface!(pr, pz, argmax(pz))

    # scalars
    IMAS.set_time_array(pc.minor_radius.reference, :data, time0, mxh.ϵ * mxh.R0)
    IMAS.set_time_array(pc.geometric_axis.r.reference, :data, time0, mxh.R0)
    IMAS.set_time_array(pc.geometric_axis.z.reference, :data, time0, mxh.Z0)
    IMAS.set_time_array(pc.elongation.reference, :data, time0, mxh.κ)
    IMAS.set_time_array(pc.triangularity.reference, :data, time0, sin(mxh.s[1]))
    IMAS.set_time_array(pc.squareness.reference, :data, time0, -mxh.s[2])

    # boundary
    resize!(pc.boundary_outline, length(pr); wipe=false)
    for (k, (r, z)) in enumerate(zip(pr, pz))
        IMAS.set_time_array(pc.boundary_outline[k].r.reference, :data, time0, r)
        IMAS.set_time_array(pc.boundary_outline[k].z.reference, :data, time0, z)
    end

    return pc
end