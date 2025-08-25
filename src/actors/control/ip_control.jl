Base.@kwdef mutable struct FUSEparameters__ActorControllerIp{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    algorithm::Entry{Symbol} = Entry{Symbol}("-", "Algorithm to dispatch on"; default=:PID)
    P::Entry{T} = Entry{T}("-", "Proportional gain")
    I::Entry{T} = Entry{T}("-", "Integral gain")
    D::Entry{T} = Entry{T}("-", "Derivative gain")
    Vloop_initial::Entry{T} = Entry{T}("V", "Initial loop voltage")
end

mutable struct ActorControllerIp{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorControllerIp{P}}
    controller::Union{Nothing,IMAS.controllers__linear_controller{<:D},IMAS.controllers__nonlinear_controller{<:D}}
end

function ActorControllerIp(dd::IMAS.dd{D}, par::FUSEparameters__ActorControllerIp{P}; kw...) where {D<:Real,P<:Real}
    logging_actor_init(ActorControllerIp)
    par = OverrideParameters(par; kw...)
    return ActorControllerIp{D,P}(dd, par, nothing)
end

"""
    ActorControllerIp(dd::IMAS.dd, act::ParametersAllActors; kw...)

Controls the loop voltage Vloop to obtain the target plasma cuurent Ip pulse_schedule

Algorithm names that have "PID" in their name will be considered linear and operate on `dd.controllers.linear_controller` otherwise `dd.controllers.nonlinear_controller`

!!! note

    External controllers can be defined by dispatching on the controllers__linear_controller functors (note, yes "PID" in the control_algorithm)

        function (controller::IMAS.controllers__linear_controller{T})(control_name::Val{:ip}, control_algorithm::Val{:my_PID}, setpoint::T, value::T, time0::Float64) where {T<:Real}
            ...
        end

    or controllers__nonlinear_controller functors (note, no "PID" in the control_algorithm)

        function (controller::IMAS.controllers__nonlinear_controller{T})(control_name::Val{:ip}, control_algorithm::Val{:my_MPC}, setpoint::T, value::T, time0::Float64) where {T<:Real}
            ...
        end
"""
function ActorControllerIp(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorControllerIp(dd, act.ActorControllerIp; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorControllerIp; time0::Float64=dd.global_time)
    dd = actor.dd
    par = actor.par

    if actor.controller === nothing

        if contains(string(par.algorithm), "PID")
            controller_category = dd.controllers.linear_controller
        else
            controller_category = dd.controllers.nonlinear_controller
        end
        index = findfirst(controller.name == "ip" for controller in controller_category)

        if index !== nothing
            # if controller already exists in dd, just reuse it
            controller = controller_category[index]

        else
            # if controller does not exist as part of dd, then initialize one
            controller = resize!(controller_category, "name" => "ip")
            controller.input_names = ["Ip error"]
            controller.output_names = ["Vloop"]

            # initialize controller
            if contains(string(par.algorithm), "PID")
                IMAS.pid_controller(controller, par.P, par.I, par.D)
            else
                controller(Val(:ip), Val(par.algorithm))
            end

            # initial guess for Vloop if none provided
            if ismissing(par, :Vloop_initial)
                Ω = IMAS.plasma_lumped_resistance(dd)
                Ip0 = IMAS.get_from(dd, Val(:ip), :pulse_schedule)
                Vloop0 = Ip0 * Ω
            else
                Vloop0 = par.Vloop_initial
            end
            controller.inputs.time = [time0]
            controller.inputs.data = [0.0;;]
            controller.outputs.time = [time0]
            controller.outputs.data = [Vloop0;;]
        end

        actor.controller = controller

    else
        # run the controller
        ip_target_now = IMAS.get_from(dd, Val(:ip), :pulse_schedule; time0)
        ip_value_now = IMAS.get_from(dd, Val(:ip), :core_profiles) # eventually this should be Ip from a synthetic Rogowski coil
        actor.controller(Val(:ip), Val(par.algorithm), ip_target_now, ip_value_now, time0)
    end

    return actor
end
