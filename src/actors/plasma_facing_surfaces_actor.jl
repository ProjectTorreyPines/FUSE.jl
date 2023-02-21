#= ============== =#
#  ActorPlasmaFacingSurfaces #
#= ============== =#

# should be the content below into a module to avoid unwanted redefinition of methods?
# can add auto export of all FUSEParameters and Actor structs if needed.
# Can inject directly struct def into module with "@newactor module ActorName"
# will permit addition of methods/parameters into local scope without constraints on naming

pfc_materials() = [:tungsten, :SiC, :graphite]

@actor_params struct DivertorParameters
    inner = DivertorLegParameters()
    outer = DivertorLegParameters()
    dome = DivertorDomeParameters()
end

@actor_params struct DivertorDomeParameters
    material = Entry(:tungsten, "dome surface material"; options=pfc_materials())
    𝓁_xpt = Entry(0.0, "distance between top of the dome and xpt")
    α_xpt = Entry(0.0, "fraction of distance between top of the dome and xpt")
end

@actor_params struct DivertorLegParameters
    target = DivertorTarget()
    baffle_cfr = CFRDivertorBaffle()
    baffle_pfr = PFRDivertorBaffle()
    l_leg = Entry(1.0, "length of the divertor leg"; units="m")
    d_heat_shield = Entry(0.25, "length of the divertor leg"; units="m")
    method = Entry(:manual, "method to define the divertor leg length"; options=[:manual, :d_heat_shield])
end

@actor_params struct CFRDivertorBaffle
    material = Entry(:tungsten, "baffle surface material"; options=pfc_materials())
    method = Entry(:xpt_equator, "method to define the baffle design point"; options=[:xpt_equator, :manual])
    𝓁_baffle = Entry(0.0, "poloidal distance from target boundary point to baffle design point"; units="m")
    r_baffle = Entry(0.0, "radial distance from separatrix to baffle design point"; units="m")
end

@actor_params struct PFRDivertorBaffle
    material = Entry(:tungsten, "baffle surface material"; options=pfc_materials())
    method = Entry(:none, "method to define the baffle design point"; options=[:none, :manual])
    𝓁_baffle = Entry(0.0, "poloidal distance from target boundary point to baffle design point"; units="m")
    r_baffle = Entry(0.0, "radial distance from separatrix to baffle design point"; units="m")
end

@actor_params mutable struct TargetBoundary
    λ_plasma = Entry(0.003, "plasma width"; units="m")
    β_plasma = Entry(5.0, "plasma width decay factor")
    Λ_buffer = Entry(0.05, "buffer between plasma and wall"; units="m")
end

@actor_params struct MainChamberWallBoundary
    λ_plasma = Entry(0.003, "plasma width"; units="m")
    β_plasma = Entry(5.0, "plasma width decay factor")
    Λ_buffer = Entry(0.05, "buffer between plasma and wall"; units="m")
end

@actor_params struct MainChamberWallDesign
    material = Entry(:tungsten, "baffle surface material"; options=pfc_materials())
    mid_plane = MainChamberWallBoundary()
    method = Entry(:plasma_width, "method to draw main chamber wall"; options=[:plasma_width, :conformal])
end


@actor_params struct DivertorTarget
    material = Entry(:tungsten, "target surface material"; options=pfc_materials())
    cfr = TargetBoundary() # cfr target boundary point 
    pfr = TargetBoundary() # pfr target boundary point
    Λ_peak = Entry(0.0, "position of the plasma peak with respect to separatrix "; units="m") # 
    type = Entry(:flat, "type of target"; options=[:flat])
    θ_target = Entry(60.0, "inclination of the target with respect to the separatrix normal")
end

@actor_params struct DivertorsDesign
    lower = DivertorParameters() #lower divertor parameters
    upper = DivertorParameters() #upper divertor parameters
end

@actor_params struct MainChamberWallsDesign
    outer = MainChamberWallDesign()
    inner = MainChamberWallDesign()
end

@actor_params struct ActorPlasmaFacingSurfaces
    divertors = DivertorsDesign()
    main_chamber_walls = MainChamberWallsDesign()
end

"""
    ActorDivertors(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evaluates divertor loading and deposited power

!!! note 
    Stores data in `dd.divertors`
"""
@new_actor ActorPlasmaFacingSurfaces

function _step(actor::ActorPlasmaFacingSurfaces)
    dd = actor.dd
    #empty!(dd.divertors)
end
