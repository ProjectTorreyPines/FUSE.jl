"""
    controller(dd::IMAS.dd, ctrl_ip, ::Type{Val{:ip}})

Controller for the plasma current
"""
function controller(dd::IMAS.dd, ctrl_ip::T, ::Type{Val{:ip}}) where {T<:Union{IMAS.controllers__linear_controller,IMAS.controllers__nonlinear_controller}}
    ip_target_now = IMAS.get_from(dd, Val{:ip}, :pulse_schedule)
    ip_value_now = IMAS.get_from(dd, Val{:ip}, :core_profiles) # eventually this should be Ip from Rogowski coil

    vloop_now = ctrl_ip(ip_target_now, ip_value_now, dd.global_time)

    @ddtime(dd.pulse_schedule.flux_control.loop_voltage.reference = vloop_now)
    # eventually we should also update the pf_active
end