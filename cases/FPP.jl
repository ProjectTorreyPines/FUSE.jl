"""
To generate a JSON file from a GASC run:

```python
    filename = "path_to_the_gasc_output.pkl"
    casename = os.path.splitext(os.path.basename(filename))[1]
    json = OMFITjson(casename + ".json", objects_encode=False)
    json.update(OMFITpickle(filename))
```
"""
function case_parameters(::Type{Val{:FPP}}; version::Symbol, init_from::Symbol)::Tuple{ParametersInit, ParametersActor}
    if version == :v1
        filename = "FPPv1.0_aspectRatio3.5_PBpR35.json"
        case = 0
    elseif version == :v1_demount
        filename = "FPPv1.0_aspectRatio3.5_PBpR35_demount.json"
        case = 0
    end

    gasc = GASC(joinpath(dirname(abspath(@__FILE__)), "..", "sample", filename), case)
    ini, act = case_parameters(gasc)
    ini.general.casename = "FPP_$(version)_$(init_from)"
    ini.general.init_from = init_from

    if init_from == :ods
        ini.ods.filename = joinpath(dirname(abspath(@__FILE__)), "..", "sample", "FPPv1.0_demount_eq.json")
        act.ActorCXbuild.rebuild_wall = true # false to use wall from ODS
        act.ActorHFSsizing.fixed_aspect_ratio = true
    end

    ini.core_profiles.rot_core = 0.0
    ini.core_profiles.bulk = :DT

    # ini.equilibrium.δ *= -1 ### for negative triangularity

    ini.core_profiles.zeff = 1.1 ↔ [1.1, 2.5]
    ini.core_profiles.greenwald_fraction = 0.9 ↔ [0.8, 0.95]
    ini.ec_launchers.power_launched = 45e6 ↔ [30e6, 100e6]

    ini.tf.shape = :princeton_D_scaled
    ini.tf.n_coils = 16

    ini.pf_active.n_oh_coils = 6
    ini.pf_active.n_pf_coils_inside = 0
    ini.pf_active.n_pf_coils_outside = 8

    ini.material.shield = "Tungsten"
    ini.material.blanket = "lithium-lead"

    act.ActorPFcoilsOpt.symmetric = true

    return set_new_base!(ini), set_new_base!(act)
end
