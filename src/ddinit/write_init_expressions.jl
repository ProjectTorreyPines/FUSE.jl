import JSON

"""
    init_expressions(; save::Bool=false)

Returns the set of ulocations that always have expressions after running init, for all FUSE.test_cases
"""
function init_expressions(; save::Bool=false)
    all_expr_fields = Vector{OrderedCollections.OrderedSet{String}}()
    for (testname, (args, kw)) in test_cases
        @info ("init primary quantities of $testname")
        ini, act = case_parameters(args...; kw...)
        dd = IMAS.dd()
        init(dd, ini, act; restore_expressions=false)

        data_fields, expr_fields = IMAS.data_and_expression_ulocations(dd)

        push!(all_expr_fields, expr_fields)
    end
    if save
        json_str = JSON.json(sort!(intersect(all_expr_fields...)), 4)
        open(joinpath(__FUSE__, "src", "ddinit", "init_expressions.json"), "w") do file
            return write(file, json_str)
        end
    end
    return all_expr_fields
end

"""
    load_init_expression_quantities()

Returns the intersection of initializing for all cases as a OrderedSet
"""
function load_init_expressions()
    json_str = read(joinpath(__FUSE__, "src", "ddinit", "init_expressions.json"), String)
    return OrderedCollections.OrderedSet(JSON.parse(json_str))
end

"""
    restore_init_expressions!(dd::IMAS.dd; verbose::Bool)

Turn into expressions Make fields that FUSE expects to be expressions after init
"""
function restore_init_expressions!(dd::IMAS.dd; verbose::Bool)
    restored_expressions = Set()
    for ulocation in load_init_expressions()
        restored = IMAS.selective_delete!(dd, IMAS.i2p(ulocation))
        if restored
            push!(restored_expressions, ulocation)
        end
    end
    if verbose && !isempty(restored_expressions)
        @info "The following ODS locations have been turned into expressions:"
        for ulocation in sort!(collect(restored_expressions))
            @info "-  $ulocation"
        end
    end
end