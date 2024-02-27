"""
    write_data_fields_after_init()

Writes the union of initializing from scalars for all cases
"""
function write_data_fields_after_init()
    all_data_fields = Vector{OrderedCollections.OrderedSet{String}}()
    for (testname, (args, kw)) in FUSE.test_cases
        if :init_from ∈ keys(kw) && kw[:init_from] == :scalars
            ini, act = FUSE.case_parameters(args...; kw...)
            dd = IMAS.dd()
            init(dd, ini, act)
            
            data_fields,expr_fields = IMAS.data_and_expression_ulocations(dd)
            
            push!(all_data_fields, data_fields)
        end
    end
    json_str = JSON.json(union(all_data_fields...))
    open(joinpath(__FUSE__,"src","ddinit","data_fields_init_from_scalars.json"), "w") do file
       write(file, json_str)
    end
    return all_data_fields
end


"""
    function load_data_fields_after_init()

Returns the union of initializing from scalars for all cases as a OrderedSet
"""
function load_data_fields_after_init()
    json_str = read(joinpath(__FUSE__,"src","ddinit","data_fields_init_from_scalars.json"), String)
    return OrderedCollections.OrderedSet(JSON.parse(json_str))
end

