"""
    ip_control(dd::IMAS.dd, δt::Float64)

Initializes (or returns existing) plasma current controller
"""
function ip_control(dd::IMAS.dd, δt::Float64)
    index = findfirst(controller.name == "ip" for controller in dd.controllers.linear_controller)
    if index !== nothing
        ctrl_ip = dd.controllers.linear_controller[index]
    else
        ctrl_ip = resize!(dd.controllers.linear_controller, "name" => "ip")
        # guess PID gains based on resistivity
        cp1d = dd.core_profiles.profiles_1d[]
        conductivity_parallel = cp1d.conductivity_parallel
        f = (k, x) -> 1.0 / conductivity_parallel[k]
        η_avg = trapz(cp1d.grid.area, f) / cp1d.grid.area[end]
        P = η_avg * 5.0
        I = η_avg * 0.5
        D = 0.0
        # instantiate PID controller
        IMAS.pid_controller(ctrl_ip, P, I, D)
        if IMAS.fxp_request_service(ctrl_ip)
            @info("Running Ip controller via FXP")
        end
        # initial guess for Vloop
        Vloop0 = IMAS.get_from(dd, Val{:vloop}, :core_profiles)
        ip_control(ctrl_ip, dd; time0=dd.global_time - δt)
        ctrl_ip.inputs.data[1, 1] = Vloop0 / I / δt * 1.5
        ctrl_ip.outputs.data[1, 1] = Vloop0
    end
    return ctrl_ip
end

"""
    ip_control(ctrl_ip::IMAS.controllers__linear_controller, dd::IMAS.dd)

Activates controller for the plasma current
"""
function ip_control(ctrl_ip::IMAS.controllers__linear_controller, dd::IMAS.dd; time0::Float64=dd.global_time)
    ip_target_now = IMAS.get_from(dd, Val{:ip}, :pulse_schedule)
    ip_value_now = IMAS.get_from(dd, Val{:ip}, :core_profiles) # eventually this should be Ip from Rogowski coil

    # the output of the controller is stored in ctrl_ip.outputs.data
    ctrl_ip(ip_target_now, ip_value_now, time0)

    # eventually we should also update pf_active
    # ...

    return ctrl_ip
end
