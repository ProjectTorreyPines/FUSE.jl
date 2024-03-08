using FUSE
using Test

@testset "check_init_expressions" begin
    data_fields_all = union(FUSE.init_expressions()...)
    data_fields_from_file = FUSE.load_init_expressions()
    @assert isempty(symdiff(data_fields_from_file, data_fields_all)) "data fields after init have changed! run `make init_expressions`"
end
