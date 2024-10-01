using FUSE
using Test

@testset "init_expressions" begin
    all_expr_fields = intersect(FUSE.init_expressions()...)
    all_expr_fields_from_file = FUSE.load_init_expressions()
    @assert isempty(symdiff(all_expr_fields, all_expr_fields_from_file)) "expressions after init() have changed! run `make init_expressions`"
end
