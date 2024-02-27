using FUSE
using Test

@testset "check_data_fields_after_init" begin
    data_fields_all = union(FUSE.write_data_fields_after_init()...)
    data_fields_from_file = FUSE.load_data_fields_after_init()
    @assert isempty(symdiff(data_fields_from_file, data_fields_all)) "datafields after init has changed"
end