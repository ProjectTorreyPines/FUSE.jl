 import PAM


# #= ====== =#
# #  PELLET  #
# #= ====== =#
 Base.@kwdef mutable struct FUSEparameters__ActorPAM{T<:Real} <: ParametersActor{T}
     _parent::WeakRef = WeakRef(nothing)
     _name::Symbol = :not_set
     _time::Float64 = NaN
     drift_model::Switch{Symbol} = Switch{Symbol}([:HPI2, :Parks, :none], "-", "drift model"; default=:none)
     time_from:: Switch{Symbol} = Switch{Symbol}([:pulse_schedule, :pellet_time, :none], "-", "initialize time for the pellet calculations"; default=:none)
     time_step::Entry{Float64} = Entry{Float64}("-", "Time step [s]"; default=0.00001)
     time_end::Entry{Float64} = Entry{Float64}("-", "Time end [s]"; default=0.0008)
     Bt_dependance::Entry{Bool} = Entry{Bool}("-", "Enable Bt dependance"; default=false)
     density_update::Entry{Bool} = Entry{Bool}("-", "Update plasma density"; default=false)
 end

 mutable struct ActorPAM{D,P} <: SingleAbstractActor{D,P}
     dd::IMAS.dd{D}
     par::FUSEparameters__ActorPAM{P}
     outputs::Union{PAM.Pellet1,Vector{<:PAM.Pellet1}}
 

    function ActorPAM(dd::IMAS.dd{D}, par::FUSEparameters__ActorPAM{P}; kw...) where {D<:Real,P<:Real}
         logging_actor_init(ActorPAM)
         par = par(kw...)
         return new{D,P}(dd, par, PAM.Pellet1[])
    end
 end

 """
     ActorPAM(dd::IMAS.dd, act::ParametersAllActors; kw...)

 Estimates the Pellet particle direction, ablation rate, density source deposition

 !!! note

     Reads data in `dd.pellet`, 'dd.equilibrium', and stores data in ...
 """
 function ActorPAM(dd::IMAS.dd, act::ParametersAllActors; kw...)
     par = act.ActorPAM(kw...)
     actor = ActorPAM(dd, par)
     step(actor)
     finalize(actor)
     return actor
 end

 function _step(actor::ActorPAM)
    dd = actor.dd
    par = actor.par

    if par.time_from == :none
        t0 = 0
    end

    inputs=(
         t_start = t0,
         t_finish = par.time_end,
         time_step = par.time_step, 
         drift_model=par.drift_model,
         BtDependance=par.Bt_dependance,
         plasma_update=par.density_update
      )
      
     
      
      
      
 
     actor.outputs = PAM.run_PAM(dd, inputs)
     

     return actor
 end
 """
    _finalize(actor::ActorPAM)

 Update plasma sources in dd
 """
 function _finalize(actor::ActorPAM)
    dd = actor.dd
    cs = dd.core_sources
    eqt = dd.equilibrium.time_slice[]
    output = actor.outputs
    rho = eqt.profiles_1d.rho_tor_norm
    volume = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.volume).(rho)
    area = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.area).(rho)   
    
   
   


    for (ps, pll) in zip(dd.pulse_schedule.pellet.launcher, dd.pellets.launcher)
         electrons_particles = vec(sum(output.density_source; dims=1))
        
         source = resize!(cs.source, :pellet, "identifier.name" => pll.name;  wipe=false)

         IMAS.new_source(
            source,
            source.identifier.index,
            pll.name,
            rho,
            volume,
            area;
            electrons_particles,
        )
    end

     return actor

     
 end    


