import EPEDNN

#= ============= =#
#  ActorPedestal  #
#= ============= =#
Base.@kwdef mutable struct FUSEparameters__ActorPedestal{T<:Real} <: ParametersActorPlasma{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    #== actor parameters ==#
    rho_nml::Entry{T} = Entry{T}("-", "Defines rho at which the no man's land region starts")
    rho_ped::Entry{T} = Entry{T}("-", "Defines rho at which the pedestal region starts") # rho_nml < rho_ped
    density_match::Switch{Symbol} = Switch{Symbol}([:ne_line, :ne_ped], "-", "Matching density based on ne_ped or line averaged density"; default=:ne_ped)
    model::Switch{Symbol} = Switch{Symbol}([:EPED, :WPED], "-", "Pedestal model to use"; default=:EPED)
    #== data flow parameters ==#
    ip_from::Switch{Symbol} = switch_get_from(:ip)
    βn_from::Switch{Symbol} = switch_get_from(:βn)
    ne_from::Switch{Symbol} = switch_get_from(:ne_ped)
    zeff_ped_from::Switch{Symbol} = switch_get_from(:zeff_ped)
    #== display and debugging parameters ==#
    warn_nn_train_bounds::Entry{Bool} = Entry{Bool}("-", "EPED-NN raises warnings if querying cases that are certainly outside of the training range"; default=false)
end

mutable struct ActorPedestal{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorPedestal{P}
    ped_actor::Union{Nothing,ActorEPED{D,P}, ActorWPED{D,P}}
end

"""
    ActorPedestal(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evaluates the pedestal boundary condition (height and width)
"""
function ActorPedestal(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorPedestal(dd, act.ActorPedestal, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function blend_core_edge(actor::ActorEPED)
    return
end

function blend_core_edge(actor::ActorWPED)
    return
end


function ActorPedestal(dd::IMAS.dd, par::FUSEparameters__ActorPedestal, act::ParametersAllActors; kw...)
    logging_actor_init(ActorPedestal)
    par = par(kw...)
    if par.model == :EPED
        ped_actor = ActorEPED(dd, act.ActorEPED; ne_ped_from=par.ne_from, par.zeff_ped_from, par.βn_from, par.ip_from, par.rho_nml, par.rho_ped)
    elseif par.model == :WPED
        ped_actor = ActorWPED(dd, act.ActorWPED)
    end

    return ActorPedestal(dd, par, ped_actor)
end

"""
    _step(actor::ActorPedestal)

Runs pedestal actor to evaluate pedestal width and height
"""
function _step(actor::ActorPedestal{D,P}) where {D<:Real,P<:Real}

    dd = actor.dd
    par = actor.par

    eq = dd.equilibrium
    eqt = eq.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]

    if par.density_match == :ne_ped
        _finalize(_step(actor.ped_actor))

    elseif par.density_match == :ne_line
        error("not fully implemented yet")
        ne_line_wanted = @ddtime(dd.pulse_schedule.density_control.n_e_line.reference)
        ne_ped_over_ne_sep = cp1d.electrons.density_thermal[end] / neped
        if par.model == :EPED
            function cost_ne_ped_from_nel(neped, ne_line_wanted)
                actor.ped_actor.inputs.neped = neped / 1e19
                cp1d.electrons.density_thermal[end] = ne_ped_over_ne_sep * neped
                sol = actor.epedmod(actor.ped_actor.inputs; par.only_powerlaw, par.warn_nn_train_bounds)

                new_density = IMAS.blend_core_edge_Hmode(cp1d.electrons.density_thermal, cp1d.grid.rho_tor_norm, neped, max(sol.width.GH.H, 0.01), par.rho_nml, par.rho_ped)
                cp1d.electrons.density_thermal = new_density
                error = ((IMAS.line_average_density_middle_plasma(eqt, cp1d) - ne_line_wanted) / ne_line_wanted)^2
                @show error, IMAS.line_average_density_middle_plasma(eqt, cp1d)
                return error
            end
            res = Optim.optimize(x -> cost_ne_ped_from_nel(x, ne_line_wanted), ne_line_wanted / 10, ne_line_wanted * 10, Optim.GoldenSection(); rel_tol=1E-3)
            actor.ped_actor.inputs.neped = res.minimizer / 1e19
            println("ne_ped found == $(res.minimizer)")
        else
            error("not implemented yet for model == $(par.model)")
        end
        ped_actor = ActorEPED(dd, act.ActorEPED; ne_ped_from=par.ne_from, par.zeff_ped_from, par.βn_from, par.ip_from, rho_nml=0.9, rho_ped=0.95)
    end

    return actor
end

"""
    _finalize(actor::ActorPedestal)

Writes results to dd.summary.local.pedestal and possibly updates core_profiles
"""
function _finalize(actor::ActorPedestal)
    return actor
end

