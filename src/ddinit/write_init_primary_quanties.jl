"""
    init_primary_quanties(; save::Bool=false)

Writes the location of dd that have data in them (not expressions) running init on FUSE.test_cases
"""
function init_primary_quanties(; save::Bool=false)
    all_data_fields = Vector{OrderedCollections.OrderedSet{String}}()
    for (testname, (args, kw)) in FUSE.test_cases
        @info ("init primary quantities of $testname")
        ini, act = FUSE.case_parameters(args...; kw...)
        dd = IMAS.dd()
        init(dd, ini, act; purge_derived_quanties=false)

        data_fields, expr_fields = IMAS.data_and_expression_ulocations(dd)

        push!(all_data_fields, data_fields)
    end
    if save
        json_str = IMAS.IMASDD.JSON.json(sort!(union(all_data_fields...)), 4)
        open(joinpath(__FUSE__, "src", "ddinit", "data_fields_init_from_scalars.json"), "w") do file
            return write(file, json_str)
        end
    end
    return all_data_fields
end


"""
    function load_init_primary_quanties()

Returns the union of initializing from scalars for all cases as a OrderedSet
"""
function load_init_primary_quanties()
    json_str = read(joinpath(__FUSE__, "src", "ddinit", "data_fields_init_from_scalars.json"), String)
    return OrderedCollections.OrderedSet(IMAS.IMASDD.JSON.parse(json_str))
end

