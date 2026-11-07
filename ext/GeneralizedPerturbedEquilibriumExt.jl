module GeneralizedPerturbedEquilibriumExt

# import FUSE namespace
import FUSE
for n in names(FUSE; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(FUSE), :eval, :include)
        @eval import FUSE: $n
    end
end

include(joinpath("..", "src", "actors", "stability", "gpec_actor_ext.jl"))

end
