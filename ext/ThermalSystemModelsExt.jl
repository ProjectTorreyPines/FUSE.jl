module ThermalSystemModelsExt

# import FUSE namespace
import FUSE
for n in names(FUSE; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(FUSE), :eval, :include)
        @eval import FUSE: $n
    end
end

include(joinpath("..","src", "actors", "balance_plant", "thermal_system_models_ext.jl"))

end
