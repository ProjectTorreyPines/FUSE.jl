"""
    controller(dd::IMAS.dd, ctrl_ip, ::Type{Val{:ip}})

Controller for the plasma current
"""
function controller(dd::IMAS.dd, ::Type{Val{:ip}})
    ip_target_now = IMAS.get_from(dd, Val{:ip}, :pulse_schedule)
    ip_value_now = IMAS.get_from(dd, Val{:ip}, :core_profiles) # eventually this should be Ip from Rogowski coil

    # the output of the controller is stored in ctrl_ip.outputs.data
    ctrl_ip = IMAS.controller(dd.controllers, "ip")
    vloop_now = ctrl_ip(ip_target_now, ip_value_now, dd.global_time)

    # eventually we should also update pf_active
end
