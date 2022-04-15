"""
To generate a JSON file from a GASC run:

```python
    filename = "path_to_the_gasc_output.pkl"
    casename = os.path.splitext(os.path.basename(filename))[1]
    json = OMFITjson(casename + ".json", objects_encode=False)
    json.update(OMFITpickle(filename))
```
"""

function case_parameters(::Type{Val{:FPP}}; version::Symbol=:v1, init_from::Symbol)
    if version == :v1
        filename = "FPPv1.0_aspectRatio3.5_PBpR35.json"
        case = 0
    elseif version == :v1_demount
        filename = "FPPv1.0_aspectRatio3.5_PBpR35_demount.json"
        case = 0
    end

    gasc = GASC(joinpath(dirname(abspath(@__FILE__)), "..", "sample", filename), case)
    ini, act = case_parameters(gasc)

    ini.general.casename = "FPP_$(init_from)"
    ini.general.init_from = init_from

    ini.gasc.no_small_gaps = true

    if init_from == :ods
        @assert version == :v0 "Can only init_from=:ods with version=:v0"
        ini.ods.filename = joinpath(dirname(abspath(@__FILE__)), "..", "sample", "fpp_gasc_59_step.json")
        ini.build.blanket = 0.9
        ini.build.shield = 0.5
        ini.build.vessel = 0.125
    else
        ini.core_profiles.rot_core = 0.0
        ini.core_profiles.bulk = :DT
    end

#    ini.oh.flattop_duration = 1000

    ini.tf.shape = :princeton_D_scaled
    ini.tf.n_coils = 16

    ini.pf_active.n_oh_coils = 6
    ini.pf_active.n_pf_coils_inside = 0
    ini.pf_active.n_pf_coils_outside = 4

    ini.material.shield = "Tungsten"
    ini.material.blanket = "FLiBe"

    act.PFcoilsOptActor.symmetric = true

    return set_new_base!(ini), set_new_base!(act)
end
