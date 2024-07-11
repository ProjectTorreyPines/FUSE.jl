using Pkg
Pkg.activate(@__DIR__)
import OrderedCollections
import Memoize
import JSON

IMASDDjl_path = joinpath(@__DIR__, "..", "..", "IMASDD", "src")

const imas_version = "3_40_0"
# JSON files come from https://github.com/gafusion/omas
if "OMAS_ROOT" ∉ keys(ENV) && "OMFIT_ROOT" ∈ keys(ENV) && isdir(joinpath(ENV["OMFIT_ROOT"], "omas", "omas"))
    ENV["OMAS_ROOT"] = joinpath(ENV["OMFIT_ROOT"], "omas")
end
if "OMAS_ROOT" ∈ keys(ENV)
    const omas_imas_structure_folder = joinpath(ENV["OMAS_ROOT"], "omas", "imas_structures")
    run(`sh -c "rm -rf $(@__DIR__)/data_structures"`)
    run(`sh -c "mkdir $(@__DIR__)/data_structures"`)
    run(`sh -c "cp -rf $(omas_imas_structure_folder)/$(imas_version)/*.json $(@__DIR__)/data_structures"`)
    run(`sh -c "rm $(@__DIR__)/data_structures/_*"`)
    println("Copying latest IMAS JSON (v$imas_version) from " * ENV["OMAS_ROOT"])
end

#= ====== =#
#  Header  #
#= ====== =#
include(joinpath(IMASDDjl_path, "data_header.jl"))

#= ====== =#
#  Header  #
#= ====== =#

const type_translator = Dict{String,DataType}("STR" => String, "INT" => Int, "FLT" => Float64, "CPX" => Complex{Float64})

"""
    imas2jl_data_type(imas_data_type::String)

Convert IMAS data type to Julia data type
"""
Memoize.@memoize function imas2jl_data_type(imas_data_type::String)
    # find data type and dimension
    (tp, dim) = split(imas_data_type, "_")
    dim = parse(Int, replace(dim, "D" => ""))
    # convert to Julia type
    if tp ∈ keys(type_translator)
        if dim == 0
            jldata_type = type_translator[tp]
            if jldata_type <: Integer
                jlzero = 0
            elseif jldata_type <: String
                jlzero = ""
            else
                jlzero = 0.0
            end
        else
            jldata_type = Array{type_translator[tp],dim}
            jlzero = Array{type_translator[tp]}(undef, (0 for k in 1:dim)...)
        end
    else
        throw(ArgumentError("`$tp` has not been mapped to Julia data type"))
    end
    return jldata_type, jlzero
end

"""
    imas_dd_ids_filenames(extras::Bool)::Vector{String}

Return list of JSON data structure filenames
"""
Memoize.@memoize function imas_dd_ids_filenames(extras::Bool)::Vector{String}
    if extras
        dirs = ("data_structures", "data_structures_extra")
    else
        dirs = ("data_structures",)
    end
    filenames = String[]
    for dir in dirs
        directory = joinpath(@__DIR__, dir)
        dir_filenames = readdir(directory)
        append!(filenames, [joinpath(directory, filename) for filename in dir_filenames if !startswith(filename, '_') && endswith(filename, ".json")])
    end
    return filenames
end

"""
    imas_dd_ids_names(extras::Bool=true)::Vector{String}

Return list of IDS names, possibly including extra structures
"""
function imas_dd_ids_names(extras::Bool=true)::Vector{String}
    return unique!([replace(basename(filename), ".json" => "") for filename in imas_dd_ids_filenames(extras)])
end

"""
    function imas_dd_ids(ids_name::String, extras::Bool=true)

Read the IMAS data structures in OMAS JSON format, possibly including extra structures
"""
Memoize.@memoize function imas_dd_ids(ids_name::String, extras::Bool=true)
    tmp = Dict{String,Any}()
    for filename in imas_dd_ids_filenames(extras)
        if basename(filename) == "$ids_name.json"
            merge!(tmp, JSON.parsefile(filename))
        end
    end
    return tmp
end

#= =================== =#
#  IMAS data structure  #
#= =================== =#
ids_names = imas_dd_ids_names()
ids_extras = String["balance_of_plant", "blanket", "build", "costing", "neutronics", "solid_mechanics", "requirements", "stability"]
append!(ids_names, ids_extras)
ids_names = unique(ids_names)

#= ================================== =#
#  Function sanitize extra structures  #
#= ================================== =#
function sanitize_extra_structures(ids_names, ids_extras)
    filename = joinpath(@__DIR__, "data_structures", "equilibrium.json")
    default = JSON.parsefile(filename)

    for ids_name in ids_names
        filename = joinpath(@__DIR__, "data_structures_extra", "$ids_name.json")
        if !isfile(filename)
            continue
        end
        extra = JSON.parsefile(filename; dicttype=OrderedCollections.OrderedDict)

        sort!(extra)

        if ids_name in ids_extras
            for (key, value) in default
                if any([match(r"equilibrium\." * item * r"\b", key) !== nothing for item in ("time",)])
                    extra[replace(key, r"^equilibrium\." => "$ids_name.")] = value = deepcopy(value)
                    value["full_path"] = replace(value["full_path"], r"^equilibrium/" => "$ids_name/")
                end
            end
        end

        open(filename, "w") do io
            return JSON.print(io, extra, 1)
        end
    end
end

#= ================================================= =#
#  Implementing selected IMAS DD as Julia structures  #
#= ================================================= =#
"""
    imas_julia_struct(desired_structure::Vector{String})

Translate the IMAS `desired_structure` entries into Julia structs
Note that `desired_structure` are fully qualified IMAS locations
"""
function imas_julia_struct(desired_structure::Vector{String})
    struct_commands = String[]

    all_info = Dict{String,Info}()
    all_info["global_time"] = Info(Tuple([]), "s", "FLT_0D", "Generic global time", true)
    all_info["dd"] = Info(Tuple([]), "-", "STRUCTURE", "Top level structure", false)

    branches = []
    timedep_structures = String[]
    push!(branches, "")

    # generate hierarchy of dicts
    ddjson = Dict() # this holds the information in json format
    ddjson_noextra = Dict() # this holds the information in json format, strictly from IMAS DD
    ddict = Dict() # this holds the information in hierarchical format
    old_path_1 = ""
    for sel in sort!(desired_structure)
        # ignore error fields
        if any(endswith(sel, x) for x in ("error_upper", "error_lower", "error_index"))
            continue
        end

        # split IMAS path in 
        path = split(sel, ".")

        # load imas data structure as needed
        if old_path_1 != path[1]
            merge!(ddjson, imas_dd_ids(string(path[1])))
            merge!(ddjson_noextra, imas_dd_ids(string(path[1]), false))
            old_path_1 = path[1]
        end

        # no obsolescent
        if ("lifecycle_status" in keys(ddjson[sel])) && (ddjson[sel]["lifecycle_status"] == "obsolescent")
            continue
        end

        # detect if this is a time dependendent array of structure
        if ("full_path" in keys(ddjson[sel])) && (endswith(ddjson[sel]["full_path"], "(itime)"))
            push!(timedep_structures, sel * "[:]")
            @assert sel * "[:].time" ∈ keys(ddjson) "Missing: $(sel * "[:].time") in JSON"
        end

        h = ddict
        for (k, item) in enumerate(path)
            if (k > 1) & (k == length(path))
                if ddjson[sel]["data_type"] ∈ ("STRUCTURE", "STRUCT_ARRAY")
                    continue
                end
                jldata_type, jlzero = imas2jl_data_type(ddjson[sel]["data_type"])
                if path[end] == "time" # leave `time` leaves as Float64
                    my_dtype = ":: $jldata_type"
                else
                    if contains("$jldata_type", "{")
                        my_dtype = replace(":: $jldata_type", "Float64" => "<:T")
                    else
                        my_dtype = replace(":: $jldata_type", "Float64" => "T")
                    end
                end
                my_zero = replace(repr(jlzero), "Float64" => "T")
                h[item] = (my_dtype, my_zero)
            end
            if item ∉ keys(h)
                h[item] = Dict()
                push!(branches, path[1:k])
            end
            h = h[item]
        end
    end

    matching_structs = OrderedCollections.OrderedDict{String,Vector{String}}()

    # generate source code of Julia structures
    # here we loop over each julia structure
    is_structarray = false
    for branch in reverse(branches)

        is_ggd = false
        if any(contains(path, "grid_ggd") for path in branch)
            is_ggd = true
            @show branch
        end

        h = ddict
        for item in branch
            h = h[item]
        end

        is_structarray = length(branch) > 0 && occursin("[", branch[end])
        if is_structarray
            if is_ggd
                struct_type = "IDSvectorRawElement{T}"
            elseif join(branch, ".") ∈ timedep_structures
                struct_type = "IDSvectorTimeElement{T}"
            else
                struct_type = "IDSvectorStaticElement{T}"
            end
        elseif length(branch) == 0
            struct_type = "DD{T}"
        elseif length(branch) == 1
            struct_type = "IDStop{T}"
        elseif is_ggd
            struct_type = "IDSraw{T}"
        else
            struct_type = "IDS{T}"
        end

        sep = "__"
        ulocation = join(branch, ".")
        struct_name_ = replace(join(branch, sep), "[:]" => "_")
        struct_name = replace(struct_name_, r"_$" => "")
        fields = String[]
        txt_parent = String[]
        inits = String[]
        for item in sort!(collect(keys(h)))

            item_ulocation = strip("$ulocation.$item", '.')
            info = get(ddjson, replace(item_ulocation, r"\[\:\]$" => ""), Dict{String,Any}())
            coordinates = Tuple(get(info, "coordinates", String[]))
            units = get(info, "units", "-")
            if occursin("[:]", item)
                data_type = get(info, "data_type", "STRUCT_ARRAY")
            else
                data_type = get(info, "data_type", "STRUCTURE")
            end
            data_type = replace(data_type, r"_TYPE$" => "")
            documentation = get(info, "documentation", "")
            extra = replace(item_ulocation, r"\[\:\]$" => "") ∉ keys(ddjson_noextra)
            all_info[item_ulocation] = Info(coordinates, units, data_type, documentation, extra)

            if typeof(h[item]) <: Tuple && typeof(h[item][1]) <: String # leaf
                push!(fields, "    var\"$(item)\" $(h[item][1])")
                push!(inits, h[item][2])

            else # branch
                if length(struct_name_) == 0 # top level
                    push!(fields, "    var\"$(item)\" :: $(item){T}")
                    push!(inits, "$(item){T}()")
                elseif occursin("[:]", item) # arrays of structs
                    item = replace(item, "[:]" => "")
                    push!(fields, "    var\"$(item)\" :: IDSvector{$(struct_name_)$(sep)$(item){T}}")
                    push!(inits, "IDSvector{$(struct_name_)$(sep)$(item){T}}()")
                else # structs
                    push!(fields, "    var\"$(item)\" :: $(struct_name_)$(sep)$(item){T}")
                    push!(inits, "$(struct_name_)$(sep)$(item){T}()")
                end
                push!(txt_parent, "    setfield!(ids.$(item), :_parent, WeakRef(ids))")
            end
        end

        txt_parent = join(txt_parent, "\n")
        if length(txt_parent) > 0
            txt_parent = "\n$(txt_parent)"
        end
        if length(struct_name) == 0
            struct_name = "dd"
            push!(fields, "    global_time :: Float64")
            push!(inits, "0.0")
            push!(fields, "    _aux :: Dict")
            push!(inits, "Dict()")
        end

        txt_inits = join(inits, ", ")
        txt_fields = join(fields, "\n")

        if txt_fields ∉ keys(matching_structs)
            matching_structs[txt_fields] = String[]
        end
        push!(matching_structs[txt_fields], struct_name)

        txt = fields
        push!(txt, "    _filled::Set{Symbol}")
        push!(txt, "    _frozen::Bool")
        push!(txt, "    _in_expression::Vector{Symbol}")
        push!(txt, "    _ref :: Union{Nothing,$struct_name}")
        push!(txt, "    _parent :: WeakRef")

        txt_constructor = """
        function $(struct_name){T}() where T
            ids = $struct_name{T}($(join(map(x -> split(x, "=")[1], split(txt_inits, ", ")), ", ")), Set{Symbol}(), false, Symbol[], nothing, WeakRef(nothing))$(txt_parent)
            return ids
        end
        """

        txt = join(txt, "\n")
        txt = """
        mutable struct $(struct_name){T} <: $(struct_type)
        $(txt)
        end

        $(rstrip(txt_constructor))

        $(struct_name)() = $(struct_name){Float64}()
        """

        push!(struct_commands, txt)

        # #uncomment to debug data structure
        # println(txt)
        # eval(Meta.parse(txt))
    end

    push!(struct_commands, "const _all_info = Dict{String,Info}(")
    for ulocation in sort!(collect(keys(all_info)))
        info = all_info[ulocation]
        push!(struct_commands, "$(repr(ulocation)) => $(repr(info)),")
    end
    push!(struct_commands, ")")

    return struct_commands
end

function field_name(field)
    return ":" * strip(split(field, "::")[1])[5:end-1]
end

sanitize_extra_structures(ids_names, ids_extras)

desired_structure = String[]
for ids_name in ids_names
    append!(desired_structure, keys(imas_dd_ids(ids_name)))
end

struct_commands = imas_julia_struct(desired_structure)

filename = abspath(joinpath(IMASDDjl_path, "dd.jl"))
open(filename, "w") do io
    return println(io, join(struct_commands, "\n"))
end
println("$filename generated!")
