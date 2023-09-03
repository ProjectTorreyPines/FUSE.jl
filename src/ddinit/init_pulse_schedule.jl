#= ======================= =#
#  init pulse_schedule IDS  #
#= ======================= =#
"""
    init_pulse_schedule(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors; simplify_time_traces::Float64=0.1)

Initialize `dd.pulse_schedule` starting from `ini` and `act` parameters
"""
function init_pulse_schedule(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors; simplify_time_traces::Float64=0.1)
    TimerOutputs.reset_timer!("init_pulse_schedule")
    TimerOutputs.@timeit timer "init_pulse_schedule" begin
        init_from = ini.general.init_from
        dd.global_time = ini.time.simulation_start

        if init_from == :ods
            dd1 = IMAS.json2imas(ini.ods.filename)
            dd1.global_time = ini.time.simulation_start
            if !ismissing(dd1.pulse_schedule, :time) && length(dd1.pulse_schedule.time) > 0
                dd.pulse_schedule = dd1.pulse_schedule
            else
                init_from = :scalars
            end

            if ismissing(ini.equilibrium, :ip)
                ini.equilibrium.ip = dd1.equilibrium.time_slice[].global_quantities.ip
            end
            if ismissing(ini.equilibrium, :R0)
                ini.equilibrium.R0 = dd1.equilibrium.vacuum_toroidal_field.r0
            end
            if ismissing(ini.equilibrium, :B0)
                ini.equilibrium.B0 = @ddtime dd1.equilibrium.vacuum_toroidal_field.b0
            end

        else
            dd1 = IMAS.dd()
        end

        nx = n_xpoints(ini, dd1) # x-points are not time dependent yet

        if init_from == :scalars
            time, data = get_time_dependent(ini.equilibrium, :ip, simplify_time_traces)
            dd.pulse_schedule.flux_control.i_plasma.reference.time = time
            dd.pulse_schedule.flux_control.i_plasma.reference.data = data

            # R0 should not be time dependent for definition of B0
            if !isempty(dd.build.layer)
                plasma = IMAS.get_build_layer(dd.build.layer, type=_plasma_)
                R0 = (plasma.R_start + plasma.R_end) / 2.0
            elseif typeof(getfield(ini.equilibrium, :R0).value) <: Function
                error("`ini.equilibrium.R0` should not be time dependent")
            else
                mxh = IMAS.MXH(ini, dd1)
                R0 = mxh.R0
            end

            time, data = get_time_dependent(ini.equilibrium, :B0, simplify_time_traces)
            dd.pulse_schedule.tf.b_field_tor_vacuum_r.reference.time = time
            dd.pulse_schedule.tf.b_field_tor_vacuum_r.reference.data = data .* R0

            # initialize position_control from mxh
            if ini.equilibrium.boundary_from == :scalars
                shape_parameters = Dict{Symbol,Tuple{Vector{Float64},Vector{Float64}}}()
                all_times = Float64[]
                for shape_parameter in [:R0, :Z0, :ϵ, :κ, :δ, :ζ, :𝚶]
                    time, data = get_time_dependent(ini.equilibrium, shape_parameter, simplify_time_traces)
                    shape_parameters[shape_parameter] = time, data
                    append!(all_times, time)
                end
                all_times = sort!(unique(all_times))

                for (k, time) in enumerate(all_times)
                    R0 = IMAS.interp1d(shape_parameters[:R0][1], shape_parameters[:R0][2]).(time)
                    Z0 = IMAS.interp1d(shape_parameters[:Z0][1], shape_parameters[:Z0][2]).(time)
                    ϵ = IMAS.interp1d(shape_parameters[:ϵ][1], shape_parameters[:ϵ][2]).(time)
                    κ = IMAS.interp1d(shape_parameters[:κ][1], shape_parameters[:κ][2]).(time)
                    δ = IMAS.interp1d(shape_parameters[:δ][1], shape_parameters[:δ][2]).(time)
                    ζ = IMAS.interp1d(shape_parameters[:ζ][1], shape_parameters[:ζ][2]).(time)
                    𝚶 = IMAS.interp1d(shape_parameters[:𝚶][1], shape_parameters[:𝚶][2]).(time)
                    κ = ini_equilibrium_elongation_true(κ, ϵ)
                    mxh = IMAS.MXH(R0, Z0, ϵ, κ, 0.0, [𝚶, 0.0], [asin(δ), -ζ])
                    mxhb = fitMXHboundary(mxh, nx)
                    init_pulse_schedule_postion_control(dd.pulse_schedule.position_control, mxhb, time)
                    if k == length(all_times) - 1 && all_times[k+1] == Inf
                        init_pulse_schedule_postion_control(dd.pulse_schedule.position_control, mxhb, Inf)
                        break
                    end
                end

            else
                # NOT SETUP FOR TIME DEPENDENCE YET
                mxh = IMAS.MXH(ini, dd1)
                ini.equilibrium(mxh)
                mxhb = fitMXHboundary(mxh, nx)

                init_pulse_schedule_postion_control(dd.pulse_schedule.position_control, mxhb, ini.time.simulation_start)
                init_pulse_schedule_postion_control(dd.pulse_schedule.position_control, mxhb, Inf)
            end
        end

        return dd
    end
end

function get_time_dependent(par::AbstractParameters, field::Symbol, simplify_time_traces::Float64=0.0)
    @assert 0.0 <= simplify_time_traces <= 1.0 "get_time_dependent() simplify_time_traces must be between [0,1]"

    value = getfield(par, field).value

    if typeof(value) <: Function
        time = collect(SimulationParameters.top(par).time.pulse_shedule_time_basis)
        data = value.(time)
        if simplify_time_traces != 0.0
            time, data = IMAS.simplify_2d_path(time, data, simplify_time_traces)
        end
    else
        time = Float64[SimulationParameters.top(par).time.simulation_start, Inf]
        data = [value, value]
    end
    return time, data
end

"""
    init_pulse_schedule_postion_control(
        pc::IMAS.pulse_schedule__position_control,
        mxhb::FUSE.MXHboundary,
        time0::Float64)

Initialize pulse_schedule.postion_control based on MXH boundary and number of x_points
"""
function init_pulse_schedule_postion_control(
    pc::IMAS.pulse_schedule__position_control,
    mxhb::FUSE.MXHboundary,
    time0::Float64)

    # MXHboundary adds x-points
    pr = mxhb.r_boundary
    pz = mxhb.z_boundary

    # x-point information
    resize!(pc.x_point, length(mxhb.RX); wipe=false)
    for (k, (rx, zx)) in enumerate(zip(mxhb.RX, mxhb.ZX))
        IMAS.set_time_array(pc.x_point[k].r.reference, :data, time0, rx)
        IMAS.set_time_array(pc.x_point[k].z.reference, :data, time0, zx)
    end

    # boundary with x-points parametrized with MXH
    mxh = IMAS.MXH(pr, pz, 2)
    pr, pz = mxh(100; adaptive=false)

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