using FUSE
using Test

@testset "check_init_primary_quanties" begin
    data_fields_all = union(FUSE.init_primary_quanties()...)
    data_fields_from_file = FUSE.load_init_primary_quanties()
    @assert isempty(symdiff(data_fields_from_file, data_fields_all)) "data fields after init have changed! run `make init_primary_quanties`"
end
