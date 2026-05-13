Base.@kwdef mutable struct study_database_item
    name::Union{Nothing,String} = nothing # item name (will be used to name the group in HDF5)
    dd::Union{Nothing,IMAS.dd} = nothing
    ini::Union{Nothing,ParametersAllInits} = nothing
    act::Union{Nothing,ParametersAllActors} = nothing
    log::Union{Nothing,String} = nothing
    timer::Union{Nothing,String} = nothing
    error::Union{Nothing,String} = nothing
end

mutable struct study_database
    df::DataFrames.DataFrame
    items::Vector{study_database_item}
end


"""
    filter(pred::Function, db::study_database)

Return a new `study_database` containing only the rows of `db.df` for which
`pred(row)` is true, and the corresponding entries from `db.items`.

`pred` is called with a `DataFrameRow` (from `eachrow`). This function does not
mutate the original `study_database`.
"""
function Base.filter(pred::Function, db::study_database; kwargs...)
    @assert nrow(db.df) == length(db.items) "Mismatch between df rows and items length"
    @assert !(:__rownum__ in names(db.df)) "Temporary column collision"
    # make a shallow copy of the DataFrame and apply the mutating filter!
    df2 = copy(db.df)
    df2[!, :__rownum__] = 1:nrow(df2)
    DataFrames.filter!(pred, df2; kwargs...)
    inds = df2[!, :__rownum__]
    DataFrames.select!(df2, DataFrames.Not(:__rownum__))
    return study_database(df2, db.items[inds])
end

"""
    filter!(pred::Function, db::study_database)

In-place filter: keep only rows where `pred(row)` is true and drop the rest
from both `db.df` and `db.items`. Returns the modified `db`.
"""
function Base.filter!(pred::Function, db::study_database; kwargs...)
    @assert nrow(db.df) == length(db.items) "Mismatch between df rows and items length"
    @assert !(:__rownum__ in names(db.df)) "Temporary column collision"
    db.df[!, :__rownum__] = 1:nrow(db.df)
    DataFrames.filter!(pred, db.df; kwargs...)
    inds = db.df[!, :__rownum__]
    DataFrames.select!(db.df, DataFrames.Not(:__rownum__))
    db.items = db.items[inds]
    return db
end

"""
    filter(mask::AbstractVector{Bool}, db::study_database)

Return a new `study_database` selecting rows where `mask` is true. `mask` must be
an `AbstractVector{Bool}` (for example `Vector{Bool}` or `BitVector`). If your mask
contains `missing` values, coalesce them before calling this function.
"""
function Base.filter(mask::AbstractVector{Bool}, db::study_database)
    @assert nrow(db.df) == length(db.items) "Mismatch between df rows and items length"
    @assert length(mask) == nrow(db.df) "Mask length must equal number of rows in df"
    return study_database(db.df[mask, :], db.items[mask])
end

"""
    filter!(mask::AbstractVector{Bool}, db::study_database)

In-place mask filter: keep only rows where `mask` is true. `mask` must be an
`AbstractVector{Bool}` (e.g. `BitVector`). Returns the modified `db`.
"""
function Base.filter!(mask::AbstractVector{Bool}, db::study_database)
    @assert nrow(db.df) == length(db.items) "Mismatch between df rows and items length"
    @assert length(mask) == nrow(db.df) "Mask length must equal number of rows in df"
    db.df = db.df[mask, :]
    db.items = db.items[mask]
    return db
end

"""
    subset(db::study_database, args...; kwargs...)

DataFrames-style `subset` dispatch for `study_database` implemented the same
way as the `filter` wrapper: add temporary row numbers, delegate to
`DataFrames.subset`, then pick items and remove the temporary column.
"""
function DataFrames.subset(db::study_database, args...; kwargs...)
    @assert nrow(db.df) == length(db.items) "Mismatch between df rows and items length"
    @assert !(:__rownum__ in names(db.df)) "Temporary column collision"
    db.df[!, :__rownum__] = 1:nrow(db.df)
    df2 = DataFrames.subset(db.df, args...; kwargs...)
    inds = df2[!, :__rownum__]
    DataFrames.select!(df2, DataFrames.Not(:__rownum__))
    items_out = db.items[inds]
    return study_database(df2, items_out)
end

"""
    subset!(db::study_database, args...; kwargs...)

In-place `subset!` for `study_database`. Mutates `db.df` using
`DataFrames.subset!` and keeps `db.items` aligned by recording row numbers in
a temporary `:__rownum__` column.
"""
function DataFrames.subset!(db::study_database, args...; kwargs...)
    @assert nrow(db.df) == length(db.items) "Mismatch between df rows and items length"
    @assert !(:__rownum__ in names(db.df)) "Temporary column collision"
    db.df[!, :__rownum__] = 1:nrow(db.df)
    DataFrames.subset!(db.df, args...; kwargs...)
    inds = db.df[!, :__rownum__]
    DataFrames.select!(db.df, DataFrames.Not(:__rownum__))
    db.items = db.items[inds]
    return db
end


"""
    getindex(db::study_database, rows)

Row-only indexing for `study_database`. Supported forms:
- `db[i]`, `db[1:3]`, `db[[1,4,7]]`, `db[mask::Vector{Bool}]`.

These return a new `study_database` containing the selected rows and the
corresponding `items`. Two-argument indexing `db[rows, cols]` delegates to the
inner `DataFrame` and returns whatever `DataFrame` indexing would return.
"""
function Base.getindex(db::study_database, i::Integer)
    @assert nrow(db.df) == length(db.items) "Mismatch between df rows and items length"
    inds = i:i
    return study_database(db.df[inds, :], db.items[collect(inds)])
end

function Base.getindex(db::study_database, r::AbstractRange{<:Integer})
    @assert nrow(db.df) == length(db.items) "Mismatch between df rows and items length"
    return study_database(db.df[r, :], db.items[collect(r)])
end

function Base.getindex(db::study_database, idx::AbstractVector{<:Integer})
    @assert nrow(db.df) == length(db.items) "Mismatch between df rows and items length"
    return study_database(db.df[idx, :], db.items[idx])
end

function Base.getindex(db::study_database, mask::AbstractVector{Bool})
    @assert nrow(db.df) == length(db.items) "Mismatch between df rows and items length"
    @assert length(mask) == nrow(db.df) "Mask length must equal number of rows in df"
    return study_database(db.df[mask, :], db.items[mask])
end

# Delegate two-argument indexing to the inner DataFrame (returns DataFrame result)
function Base.getindex(db::study_database, rows, cols)
    return db.df[rows, cols]
end


function save_study_database(
    savedir::AbstractString,
    parent_group::AbstractString,
    dd::Union{Nothing,IMAS.dd},
    ini::Union{Nothing,ParametersAllInits},
    act::Union{Nothing,ParametersAllActors},
    log_io::IO; 
    error_info::Any=nothing,
    timer::Bool=true,
    varinfo::Bool=false,
    freeze::Bool=false,
    overwrite_groups::Bool=false,
    verbose::Bool=false,
    database_name::Union{Nothing,AbstractString}=nothing,
    kw...
)

    savedir = abspath(savedir)
    if !isdir(savedir)
        mkpath(savedir)
    end

    parent_group = IMAS.norm_hdf5_path(parent_group)

    # allow caller to supply a specific HDF5 filename; otherwise default to pid-based file
    if database_name !== nothing
        name = string(database_name)
        @assert endswith(lowercase(name), ".h5") "database_name must end with .h5"
		h5_filename =joinpath(savedir, name)
    else
        h5_filename = joinpath(savedir, "pid$(getpid())_output.h5")
    end

    function check_and_create_group(fid::HDF5.File, target_group::AbstractString)
        if haskey(fid, target_group)
            if target_group == "/"
                gparent = fid
            else
                if !overwrite_groups
                    error("Target group '$target_group' already exists in file '$(fid.filename)'. " *
                          "\n       Set `overwrite_groups`=true to replace the existing group.")
                else
                    verbose && @warn "Target group '$target_group' already exists. Overwriting it..."
                    HDF5.delete_object(fid, target_group)
                    gparent = HDF5.create_group(fid, target_group)
                end
            end
        else
            gparent = HDF5.create_group(fid, target_group)
        end
        attr = HDF5.attrs(gparent)
        attr["date_time"] = Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS")
        return gparent
    end

    function check_and_write(fid::HDF5.File, target_group::AbstractString, data)
        if haskey(fid, target_group)
            if target_group != "/"
                if !overwrite_groups
                    error("Target group '$target_group' already exists in file '$(fid.filename)'. " *
                          "\n       Set `overwrite_groups`=true to replace the existing group.")
                else
                    verbose && @warn "Target group '$target_group' already exists. Overwriting it..."
                    HDF5.delete_object(fid, target_group)
                end
            end
        end
        HDF5.write(fid, target_group, data)
        attr = HDF5.attrs(fid[target_group])
        return attr["date_time"] = Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS")
    end

    mode = isfile(h5_filename) ? "r+" : "w"

    HDF5.h5open(h5_filename, mode) do fid

        attr = HDF5.attrs(fid)
        attr["FUSE_version"] = string(pkgversion(FUSE))
        attr["date_time"] = Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS")
        attr["original_file_abs_path"] = abspath(fid.filename)
        attr["original_file_rel_path"] = relpath(fid.filename)

        if !haskey(fid, parent_group)
            HDF5.create_group(fid, parent_group)
        end
        attr = HDF5.attrs(fid[parent_group])
        attr["FUSE_version"] = string(pkgversion(FUSE))
        attr["date_time"] = Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS")
        attr["original_file_abs_path"] = abspath(fid.filename)
        attr["original_file_rel_path"] = relpath(fid.filename)

        # Write error information into the HDF5 file (instead of separate txt file)
        if error_info !== nothing
            error_str = ""
            if typeof(error_info) <: Exception
                io = IOBuffer()
                showerror(io, error_info, catch_backtrace())
                error_str = String(take!(io))
            else
                error_str = string(error_info)
            end
            check_and_write(fid, parent_group * "/error.txt", error_str)
        end

        if ini !== nothing
            gparent = check_and_create_group(fid, parent_group * "/ini.h5")
            SimulationParameters.par2hdf!(ini, gparent)
        end

        if dd !== nothing
            IMAS.imas2hdf(dd, h5_filename; mode="a", freeze, target_group=parent_group * "/dd.h5", overwrite=overwrite_groups, show_warnings=verbose)
        end

        if act !== nothing
            gparent = check_and_create_group(fid, parent_group * "/act.h5")
            SimulationParameters.par2hdf!(act, gparent)
        end

        # save timer output
        if timer
            check_and_write(fid, parent_group * "/timer.txt", string(FUSE.timer))
        end

        # save memory trace
        if parse(Bool, get(ENV, "FUSE_MEMTRACE", "false"))
            memtrace_string = String[]
            for (date, txt, kb) in FUSE.memtrace.data
                push!(memtrace_string, "$date $kb \"$txt\"")
            end
            check_and_write(fid, parent_group * "/memtrace.txt", memtrace_string)
        end

        # save vars usage
        if varinfo
            varinfo_string = string(FUSE.varinfo(FUSE; all=true, imported=true, recursive=true, sortby=:size, minsize=1024))
            check_and_write(fid, parent_group * "/varinfo.txt", varinfo_string)
        end

        # save log
        flush(log_io)
        seekstart(log_io)
        log_str = read(log_io, String)
        if !isempty(log_str)
            check_and_write(fid, parent_group * "/log.txt", log_str)
        end
    end

    return savedir
end

function save_study_database(h5path::AbstractString, db::study_database; kw...)
    @assert nrow(db.df) == length(db.items) "Mismatch between df rows and items length"
    @assert endswith(lowercase(h5path), ".h5") "The first argument must be a file path ending with .h5"

    db_path = abspath(h5path)
    savedir = dirname(db_path)
    if !isdir(savedir)
        mkpath(savedir)
    end

    database_name = basename(db_path)

    # Check and validate :gparent column, create/update if necessary
    if !hasproperty(db.df, :gparent)
        # If :gparent column doesn't exist, create it from item names
        parent_groups = String[]
        for (i, item) in pairs(db.items)
            if !isnothing(item.name) && !isempty(item.name)
                push!(parent_groups, IMAS.norm_hdf5_path(item.name))
            else
                Lpad = length(string(length(db.items)))
                auto_item_name = "/item" * lpad(i, Lpad, "0")
                push!(parent_groups, IMAS.norm_hdf5_path(auto_item_name))
                @warn "Item at index $i has an empty name; assigning autogenerated name \"$auto_item_name\""
            end
        end
        # Add the new column to the DataFrame
        db.df[!, :gparent] = parent_groups
    else
        # :gparent exists, process each value
        parent_groups = db.df[!, :gparent]

        # Process each :gparent value individually
        updated_parent_groups = String[]
        Lpad = length(string(length(db.items)))

        for (i, pg) in enumerate(parent_groups)
            if ismissing(pg) || !isa(pg, AbstractString) || isempty(pg)
                # Generate auto name for empty/missing values
                item = db.items[i]
                if !isnothing(item.name) && !isempty(item.name)
                    push!(updated_parent_groups, IMAS.norm_hdf5_path(item.name))
                else
                    auto_item_name = "/item" * lpad(i, Lpad, "0")
                    push!(updated_parent_groups, IMAS.norm_hdf5_path(auto_item_name))
                    @warn "Item at index $i has empty :gparent and name; assigning autogenerated name \"$auto_item_name\""
                end
            else
                # Just normalize existing non-empty values
                push!(updated_parent_groups, IMAS.norm_hdf5_path(pg))
            end
        end

        parent_groups = updated_parent_groups
        db.df[!, :gparent] = parent_groups
    end
    
    for i in eachindex(parent_groups)
        pg = parent_groups[i]
        item = db.items[i]

        log_buf = IOBuffer(item.log === nothing ? "" : item.log)
        try
            save_study_database(savedir, pg, item.dd, item.ini, item.act, log_buf;
                                error_info=item.error, database_name=database_name, kw...)
        finally
            close(log_buf)
        end
    end

    # Merge/append extract.csv into the given HDF5 file and write a local copy
    merged_df = copy(db.df)
    HDF5.h5open(db_path, "r+") do fid
        if haskey(fid, "/extract.csv")
            existing_df = coalesce.(CSV.read(IOBuffer(fid["/extract.csv"][]), DataFrame), NaN)
            merged_df = vcat(existing_df, merged_df; cols=:union)
            HDF5.delete_object(fid, "/extract.csv")
        end

        unique!(merged_df)
        sort!(merged_df, "gparent")

        io_buffer = IOBuffer()
        CSV.write(io_buffer, merged_df)
        csv_text = String(take!(io_buffer))
        HDF5.write(fid, "extract.csv", csv_text)
        attr = HDF5.attrs(fid["/extract.csv"])
        attr["date_time"] = Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS")
        attr["FUSE_version"] = string(pkgversion(FUSE))
    end

    return db_path
end


"""
    load_study_database(filename::AbstractString; kw...)

Loads all data in a combined HDF5 database. Additional keyword arguments are forwarded.

### Example:
```julia
    data = FUSE.load_study_database("database.h5")
```
"""
function load_study_database(filename::AbstractString; kw...)
    @assert HDF5.ishdf5(filename) "\"$filename\" is not the HDF5 format"

    HDF5.h5open(filename, "r") do H5_fid
        df = coalesce.(CSV.read(IOBuffer(H5_fid["/extract.csv"][]), DataFrame), NaN)
        return load_study_database(filename, df[!, :gparent]; kw...)
    end
end

"""
    load_study_database(filename::AbstractString, conditions::Function; kw...)

Loads a combined HDF5 database using a filtering condition. The condition is applied to the
extract.csv contents to select parent group paths for loading.

### Example:
```julia
    data = FUSE.load_study_database("database.h5", x -> x.status=="fail"; kw...)
    data = FUSE.load_study_database("database.h5", x -> x.R0>2 && x."<zeff>">1.5; kw...)
```
"""
function load_study_database(filename::AbstractString, conditions::Function; kw...)
    @assert HDF5.ishdf5(filename) "\"$filename\" is not the HDF5 format"

    HDF5.h5open(filename, "r") do H5_fid
        df = coalesce.(CSV.read(IOBuffer(H5_fid["/extract.csv"][]), DataFrame), NaN)
        parent_groups = filter(conditions, df)[!, :gparent]
        return load_study_database(filename, parent_groups; kw...)
    end
end

"""
    load_study_database(filename::AbstractString, parent_group::AbstractString, kw...)

Loads a combined HDF5 database for a single parent group path

### Example:
```julia
    data = FUSE.load_study_database("database.h5", "/case01"; kw...)
```
"""
function load_study_database(filename::AbstractString, parent_group::AbstractString, kw...)
    return load_study_database(filename, [parent_group]; kw...)
end

"""
    load_study_database(filename::AbstractString, parent_groups::Vector{<:AbstractString}; pattern::Regex=r"", kw...)

Loads a combined HDF5 database for the specified parent group paths.

### Example:
```julia
    data = FUSE.load_study_database("database.h5", ["/case01", "/case02"]; pattern=r"dd.h5")
```
"""
function load_study_database(filename::AbstractString, parent_groups::Vector{<:AbstractString}; pattern::Regex=r"", kw...)
    @assert HDF5.ishdf5(filename) "\"$filename\" is not the HDF5 format"

    parent_groups = IMAS.norm_hdf5_path.(parent_groups)

    H5_fid = HDF5.h5open(filename, "r")

    # Load dataframe (from extract)
    df = coalesce.(CSV.read(IOBuffer(H5_fid["/extract.csv"][]), DataFrame), NaN)
    df = DataFrames.subset(df, :gparent => ByRow(x -> x in parent_groups))

    Nparents = length(parent_groups)

    # Prepare a study_database with Nparents empty items and fill them in-place
    items = [study_database_item() for _ in 1:Nparents]

    for (k, gparent) in pairs(parent_groups)
        filterd_keys = filter(x -> occursin(pattern, x), keys(H5_fid[gparent]))
        items[k].name = lstrip(gparent, '/')
        for key in filterd_keys
            h5path = gparent * "/" * key
            if key == "dd.h5"
                items[k].dd = IMAS.hdf2imas(filename, h5path)
            elseif key == "dd.json"
                items[k].dd = IMAS.jstr2imas(H5_fid[h5path][])
            elseif key == "ini.h5"
                items[k].ini = SimulationParameters.hdf2par(H5_fid[h5path], ParametersInits())
            elseif key == "ini.json"
                items[k].ini = SimulationParameters.jstr2par(H5_fid[h5path][], ParametersInits())
            elseif key == "act.h5"
                items[k].act = SimulationParameters.hdf2par(H5_fid[h5path], ParametersActors())
            elseif key == "act.json"
                items[k].act = SimulationParameters.jstr2par(H5_fid[h5path][], ParametersActors())
            elseif key == "log.txt"
                items[k].log = H5_fid[h5path][]
            elseif key == "timer.txt"
                items[k].timer = H5_fid[h5path][]
            elseif key == "error.txt"
                items[k].error = H5_fid[h5path][]
            end
        end
    end

    close(H5_fid)

    return study_database(df, items)
end

"""
    sample_and_write_study_database(ori_DB_name::AbstractString, sampled_DB_name::AbstractString, conditions::Function)

Samples the database by filtering groups that satisfy the provided `conditions` function.

### Examples:
```julia
    sample_and_write_study_database("database.h5", "sampled.h5", x -> x.status == "fail")
    sample_and_write_study_database("database.h5", "sampled.h5", x -> x.R0>2 && x."<zeff>">1.5)
```
"""
function sample_and_write_study_database(ori_DB_name::AbstractString, sampled_DB_name::AbstractString, conditions::Function; kw...)
    @assert HDF5.ishdf5(ori_DB_name) "\"$ori_DB_name\" is not the HDF5 format"

    HDF5.h5open(ori_DB_name, "r") do H5_fid
        df = coalesce.(CSV.read(IOBuffer(H5_fid["/extract.csv"][]), DataFrame), NaN)
        parent_groups = filter(conditions, df)[!, :gparent]
        return sample_and_write_study_database(ori_DB_name, sampled_DB_name, parent_groups; kw...)
    end
end

"""
    sample_and_write_study_database(ori_DB_name::AbstractString, sampled_DB_name::AbstractString;
                              Nsamples::Union{Nothing,Int}=nothing, sampling_ratio::Union{Nothing,Float64}=nothing)

Samples the database by randomly selecting a subset of the original HDF5 file.
The sample size is determined by either a fixed number (`Nsamples`) or a fraction (`sampling_ratio`) of the
total number of rows.

### Examples:
```julia
    sample_and_write_study_database("database.h5", "sampled.h5"; sampling_ratio=0.2)
    sample_and_write_study_database("database.h5", "sampled.h5"; Nsamples=15)
```
"""
function sample_and_write_study_database(ori_DB_name::AbstractString, sampled_DB_name::AbstractString;
    Nsamples::Union{Nothing,Int}=nothing, sampling_ratio::Union{Nothing,Float64}=nothing, kw...)

    @assert HDF5.ishdf5(ori_DB_name) "\"$ori_DB_name\" is not the HDF5 format"
    if isnothing(Nsamples) && isnothing(sampling_ratio)
        error("Either Nsamples or sampling_ratio must be provided.")
    end

    HDF5.h5open(ori_DB_name, "r") do H5_fid
        ori_df = coalesce.(CSV.read(IOBuffer(H5_fid["/extract.csv"][]), DataFrame), NaN)
        if isnothing(Nsamples)
            Nsamples = clamp(ceil(Int, sampling_ratio * nrow(ori_df)), 1, nrow(ori_df))
        else
            Nsamples = clamp(Nsamples, 1, nrow(ori_df))
        end
        sampled_df = ori_df[Random.shuffle(1:nrow(ori_df))[1:Nsamples], :]
        return sample_and_write_study_database(ori_DB_name, sampled_DB_name, sampled_df[!, :gparent]; kw...)
    end
end

"""
    sample_and_write_study_database(ori_DB_name::AbstractString, sampled_DB_name::AbstractString, parent_group::AbstractString)

Sampling a single specific group from the original HDF5 file.

### Example:
```julia
    sample_and_write_study_database("database.h5", "sampled.h5", "/case01")
```
"""
function sample_and_write_study_database(ori_DB_name::AbstractString, sampled_DB_name::AbstractString, parent_group::AbstractString; kw...)
    return sample_and_write_study_database(ori_DB_name, sampled_DB_name, [parent_group]; kw...)
end

"""
    sample_and_write_study_database(ori_DB_name::AbstractString, sampled_DB_name::AbstractString, parent_groups::Vector{<:AbstractString})

Sampling the groups specified in `parent_groups` from the original HDF5 file.

### Example:
```julia
    df = sample_and_write_study_database("database.h5", "sampled.h5", ["/case01", "/case02"])
```
"""
function sample_and_write_study_database(ori_DB_name::AbstractString, sampled_DB_name::AbstractString, parent_groups::Vector{<:AbstractString}; sub_format::Union{Symbol,Nothing}=nothing)

    @assert HDF5.ishdf5(ori_DB_name) "\"$ori_DB_name\" is not the HDF5 format"
    @assert sub_format ∈ (:h5, :json, nothing) "sub_format must be either `:h5` or `:json" or `nothing`

    ori_fid = HDF5.h5open(ori_DB_name, "r")
    new_fid = HDF5.h5open(sampled_DB_name, "w")

    ori_df = coalesce.(CSV.read(IOBuffer(ori_fid["/extract.csv"][]), DataFrame), NaN)

    # normalize and unique the parent_groups
    parent_groups = IMAS.norm_hdf5_path.(parent_groups)
    unique!(parent_groups)

    # Check if any requested groups don't exist in the original database
    missing_groups = setdiff(parent_groups, ori_df.gparent)
    if !isempty(missing_groups)
        missing_list = join(["\n  [$i]: \"$group\"" for (i, group) in pairs(missing_groups)], "")
        @warn "Following $(length(missing_groups)) groups not found in original database:$missing_list"
    end

    sampled_df = filter(row -> string(row.gparent) in parent_groups, ori_df)
    sort!(sampled_df, "gparent")

    # prepare progressmeter
    ProgressMeter.ijulia_behavior(:clear)
    Ngparents = length(sampled_df.gparent)
    p = ProgressMeter.Progress(Ngparents; showspeed=true)

    if isnothing(sub_format)
        for (iter, gparent) in pairs(sampled_df.gparent)
            HDF5.copy_object(ori_fid, gparent, new_fid, gparent)
            ProgressMeter.next!(p; showvalues =[(:groups, "($iter/$Ngparents) \"$gparent\"")])
        end
    else
        for (iter, gparent) in pairs(sampled_df.gparent)
            for key in keys(ori_fid[gparent])

                ori_h5path = gparent * "/" * key

                if endswith(key, r"\.txt")
                    HDF5.copy_object(ori_fid, ori_h5path, new_fid, ori_h5path)
                elseif endswith(key, Regex(string(sub_format)))
                    HDF5.copy_object(ori_fid, ori_h5path, new_fid, ori_h5path)
                else
                    if startswith(key, r"dd")
                        if sub_format == :json
                            json_string = string(IMAS.hdf2imas(ori_DB_name, gparent * "/dd.h5"))
                            HDF5.write(new_fid, gparent * "/dd.json", json_string)
                        else
                            dd = IMAS.jstr2imas(ori_fid[ori_h5path][], IMAS.dd())
                            HDF5.create_group(new_fid, gparent * "/dd.h5")
                            IMAS.imas2hdf(dd, new_fid[gparent*"/dd.h5"])
                        end
                    elseif startswith(key, r"ini")
                        if sub_format == :json
                            ini = SimulationParameters.hdf2par(ori_fid[ori_h5path], ParametersInits())
                            json_string = SimulationParameters.par2jstr(ini)
                            HDF5.write(new_fid, gparent * "/ini.json", json_string)
                        else
                            ini = SimulationParameters.jstr2par(ori_fid[ori_h5path][], ParametersInits())
                            HDF5.create_group(new_fid, gparent * "/ini.h5")
                            SimulationParameters.par2hdf!(ini, new_fid[gparent*"/ini.h5"])
                        end
                    elseif startswith(key, r"act")
                        if sub_format == :json
                            act = SimulationParameters.hdf2par(ori_fid[ori_h5path], ParametersActors())
                            json_string = SimulationParameters.par2jstr(act)
                            HDF5.write(new_fid, gparent * "/act.json", json_string)
                        else
                            act = SimulationParameters.jstr2par(ori_fid[ori_h5path][], ParametersActors())
                            HDF5.create_group(new_fid, gparent * "/act.h5")
                            SimulationParameters.par2hdf!(act, new_fid[gparent*"/act.h5"])
                        end
                    end
                end
            end
            ProgressMeter.next!(p; showvalues =[(:groups, "($iter/$Ngparents) \"$gparent\"")])
        end
    end
    ProgressMeter.finish!(p)

    # write extract.csv into HDF5
    io_buffer = IOBuffer()
    CSV.write(io_buffer, sampled_df)
    csv_text = String(take!(io_buffer))
    HDF5.write(new_fid, "extract.csv", csv_text)
    attr = HDF5.attrs(new_fid["/extract.csv"])
    attr["date_time"] = Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS")
    attr["FUSE_version"] = string(pkgversion(FUSE))

    close(ori_fid)
    close(new_fid)

    return sampled_df
end


# ===================== #
# Equality comparisons  #
# ===================== #

"""
    isequal(x::study_database_item, y::study_database_item)

Content-wise comparison of two `study_database_item`s.

Rules:
- `dd`: uses `isequal` on IMAS.dd if both non-`nothing`; otherwise both must be `nothing`.
- `ini`/`act`: compares via `SimulationParameters.par2jstr` when both non-`nothing`; otherwise both must be `nothing`.
- `log`/`timer`/`error`: standard `isequal` on `Union{Nothing,String}`.
"""
function Base.isequal(x::study_database_item, y::study_database_item)
    
    name_eq = isequal(x.name, y.name)
    
    # dd: trust IMAS.dd equality when present; else both must be nothing
    dd_eq = (x.dd === nothing || y.dd === nothing) ? (x.dd === y.dd) : (x.dd == y.dd)

    # ini/act: use SimulationParameters.diff-based equality when both present
    ini_eq = if x.ini === nothing || y.ini === nothing
        x.ini === y.ini
    else
        !SimulationParameters.diff(x.ini, y.ini)
    end

    act_eq = if x.act === nothing || y.act === nothing
        x.act === y.act
    else
        !SimulationParameters.diff(x.act, y.act)
    end

    log_eq   = isequal(x.log,   y.log)
    timer_eq = isequal(x.timer, y.timer)
    err_eq   = isequal(x.error, y.error)

    return name_eq && dd_eq && ini_eq && act_eq && log_eq && timer_eq && err_eq
end

"""
    ==(x::study_database_item, y::study_database_item)

Defers to `isequal(x, y)` for convenience.
"""
Base.:(==)(x::study_database_item, y::study_database_item) = isequal(x, y)


"""
    isequal(x::study_database, y::study_database)

Content-wise comparison of two `study_database`s. Checks DataFrame equality and
pairwise `items` equality in order.
"""
function Base.isequal(x::study_database, y::study_database)
    if !isequal(x.df, y.df)
        return false
    end
    if length(x.items) != length(y.items)
        return false
    end
    @inbounds for i in eachindex(x.items)
        if !isequal(x.items[i], y.items[i])
            return false
        end
    end
    return true
end

"""
    ==(x::study_database, y::study_database)

Defers to `isequal(x, y)` for convenience.
"""
Base.:(==)(x::study_database, y::study_database) = isequal(x, y)
