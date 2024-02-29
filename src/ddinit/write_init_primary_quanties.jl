"""
    init_primary_quanties()

Writes the union of initializing from scalars for all cases
"""
function init_primary_quanties(; save_file=false)
    all_data_fields = Vector{OrderedCollections.OrderedSet{String}}()
    for (testname, (args, kw)) in FUSE.test_cases
        if :init_from ∈ keys(kw) && kw[:init_from] == :scalars
            ini, act = FUSE.case_parameters(args...; kw...)
            dd = IMAS.dd()
            init(dd, ini, act)

            data_fields, expr_fields = IMAS.data_and_expression_ulocations(dd)

            push!(all_data_fields, data_fields)
        end
    end
    if save_file
        json_str = IMAS.IMASDD.JSON.json(union(all_data_fields...), 4)
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

