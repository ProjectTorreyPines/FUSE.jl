"""
    ip_control(ctrl_ip::IMAS.controllers__linear_controller, dd::IMAS.dd)

Controller for the plasma current
"""
function ip_control(ctrl_ip::IMAS.controllers__linear_controller, dd::IMAS.dd; time0::Float64=dd.global_time)
    ip_target_now = IMAS.get_from(dd, Val{:ip}, :pulse_schedule)
    ip_value_now = IMAS.get_from(dd, Val{:ip}, :core_profiles) # eventually this should be Ip from Rogowski coil

    # the output of the controller is stored in ctrl_ip.outputs.data
    ctrl_ip(ip_target_now, ip_value_now, time0)

    # eventually we should also update pf_active
end
