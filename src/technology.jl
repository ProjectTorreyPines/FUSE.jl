"""
Material properties
"""
Base.@kwdef struct MaterialProperties
    yield_strength::Float64 = NaN
    young_modulus::Float64 = NaN
    poisson_ratio::Float64 = NaN
end

const stainless_steel = MaterialProperties(;
    yield_strength=800E6, # Pa
    young_modulus=193.103448275E9, # Pa
    poisson_ratio=0.33
)

const pure_copper = MaterialProperties(;
    yield_strength=70E6, # Pa
    young_modulus=110E9, # Pa
    poisson_ratio=0.34
)

"""
    coil_technology(technology::Symbol, coil_type::Symbol)

Return coil parameters from technology $(supported_coils_techs) and coil type [:oh, :tf, :pf_active]"
"""
function coil_technology(coil_tech::Union{IMAS.build__pf_active__technology,IMAS.build__oh__technology,IMAS.build__tf__technology}, technology::Symbol, coil_type::Symbol)
    if coil_type ∉ (:oh, :tf, :pf_active)
        error("Supported coil type are [:oh, :tf, :pf_active]")
    end
    
    FusionMaterials.is_supported_material(technology, IMAS._tf_)

    if technology == :Copper
        coil_tech.material = "Copper"
        coil_tech.temperature = 293.0
        coil_tech.fraction_steel = 0.0
        coil_tech.ratio_SC_to_copper = 0.0
        coil_tech.fraction_void = 0.2

    elseif technology ∈ (:Nb3Sn, :NbTi, :ITER_Nb3Sn, :KDEMO_Nb3Sn, :ReBCO)
        if technology == :Nb3Sn
            coil_tech.temperature = 4.2
            coil_tech.material = "Nb3Sn"
            coil_tech.fraction_void = 0.1
        elseif technology == :NbTi 
            coil_tech.temperature =  4.2
            coil_tech.material = "NbTi"
            coil_tech.fraction_void = 0.2 # from Supercond. Sci. Technol. 36 (2023) 075009
        elseif technology == :ITER_Nb3Sn
            coil_tech.temperature = 4.2
            coil_tech.material = "ITER_Nb3Sn"
            coil_tech.fraction_void = 0.1 
        elseif technology == :KDEMO_Nb3Sn
            coil_tech.temperature = 4.2 
            coil_tech.material = "KDEMO_Nb3Sn"
            if coil_type == :tf
                coil_tech.fraction_void = 0.26 # from NF 55 (2015) 053027, Table 2
            end
        else
            coil_tech.temperature = 4.2
            coil_tech.material = "ReBCO"
        end
        coil_tech.fraction_steel = 0.5
        coil_tech.ratio_SC_to_copper = 1.0
        coil_tech.fraction_void = 0.1
    end

    if technology == :ITER_Nb3Sn
        if coil_type == :oh
            coil_tech.thermal_strain = -0.64
            coil_tech.JxB_strain = -0.05
            coil_tech.fraction_steel = 0.46
        elseif coil_type == :tf
            coil_tech.thermal_strain = -0.69
            coil_tech.JxB_strain = -0.13
            coil_tech.fraction_steel = 0.55
        elseif coil_type == :pf_active
            coil_tech.thermal_strain = -0.64
            coil_tech.JxB_strain = -0.05
            coil_tech.fraction_steel = 0.46
        end
    end

    coil_tech.thermal_strain = 0.0
    coil_tech.JxB_strain = 0.0

    return coil_tech
end

"""
    coil_J_B_crit(Bext, coil_tech::Union{IMAS.build__pf_active__technology,IMAS.build__oh__technology,IMAS.build__tf__technology})

Returns critical current density and magnetic field given an external magnetic field and coil technology
"""
function coil_J_B_crit(Bext, coil_tech::Union{IMAS.build__pf_active__technology,IMAS.build__oh__technology,IMAS.build__tf__technology})
    if coil_tech.material == "Copper"
        mat = FusionMaterials.Material(:Copper)
    elseif coil_tech.material == "Nb3Sn"
        mat = FusionMaterials.Material(:Nb3Sn; coil_tech, Bext)
    elseif coil_tech.material == "NbTi"
        mat = FusionMaterials.Material(:NbTi; coil_tech, Bext)
    elseif coil_tech.material == "KDEMO_Nb3Sn"
        mat = FusionMaterials.Material(:KDEMO_Nb3Sn; coil_tech, Bext)
    elseif coil_tech.material == "ITER_Nb3Sn"
        mat = FusionMaterials.Material(:ITER_Nb3Sn; coil_tech, Bext)
    elseif coil_tech.material == "ReBCO"
        mat = FusionMaterials.Material(:ReBCO; coil_tech, Bext)
    end
    Jcrit, Bcrit = mat.critical_current_density, mat.critical_magnetic_field

    return Jcrit, Bcrit
end

function GAMBL_blanket(bm::IMAS.blanket__module)
    layers = resize!(bm.layer, 3)

    n = 1
    layers[n].name = "First wall"
    layers[n].material = "Tungsten"
    layers[n].thickness = 0.02

    n = n + 1
    layers[n].name = "Breeder"
    layers[n].material = "lithium-lead"
    layers[n].thickness = 0.5

    n = n + 1
    layers[n].name = "Shield"
    layers[n].material = "Tungsten"
    layers[n].thickness = 0.05

    return bm
end
