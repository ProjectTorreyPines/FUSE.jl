module App
# == Packages ==
# set up Genie development environment. Use the Package Manager to install new packages
using GenieFramework
using FUSE
@genietools

struct NoDefault end
const NODEFAULT = NoDefault()
const EMPTY_SYMANY_DICT   = Dict{Symbol,Any}()
const EMPTY_STRING_DICT   = Dict{String,String}()

const case_methods = methods(FUSE.case_parameters)
const case_values = [Base.unwrap_unionall(m.sig).parameters[2] for m in case_methods]

# == Code import ==
# add your data analysis code here or in the lib folder. Code in lib/ will be
# automatically loaded
include("src/cases.jl")
include("src/session.jl")

# Format a value for display in the UI (Symbol -> :sym, String -> "str", etc.)
format_value(v::Symbol) = string(":", v)
format_value(v::AbstractString) = string("\"", v, "\"")
format_value(v) = string(v)

# Parse user input back to the appropriate type
function parse_value(::Type{Symbol}, s::AbstractString)
    s = strip(s)
    @assert !isempty(s) "Symbol string cannot be empty"
    @assert s[1] == ':' "Symbol string must start with ':'"
    @assert length(s) > 1 "Symbol string must have characters after ':'"
    return Symbol(s[2:end])
end

function parse_value(::Type{T}, s::AbstractString) where {T <: AbstractString}
    s = strip(s)
    @assert length(s) >= 2 "String must be quoted (e.g., \"value\")"
    @assert s[1] == '"' && s[end] == '"' "String must be quoted with \" (e.g., \"value\")"
    return String(s[2:end-1])
end

function parse_value(T::Type, s::AbstractString)
    s = strip(s)
    try
        return parse(T, s)
    catch e
        !(e isa MethodError) && rethrow(e)
        try
            return T(s)
        catch e2
            error("Failed to parse value \"$s\" as type $T")
        end
    end
end

"Given a kwarg spec like (T, default) or (T,), return (T, default_or_NODEFAULT)."
get_kw_type_and_default(spec::Tuple) =
    length(spec) == 1 ? (spec[1], NODEFAULT) : (spec[1], spec[2])


valid_cases = sort(unique([string(cv.parameters[1]) for cv in case_values if cv <: Val]), by=lowercase)

const default_usecase = :D3D
ms =  find_case_methods(default_usecase)
#mstrs = [tostring(m) for m in ms]

mstrs0, ms_dict0 = get_method_strdict(ms, default_usecase)

const default_usecase_method_strs = deepcopy(mstrs0)
const default_usecase_method_str = deepcopy(mstrs0[1])
const default_usecase_method = ms_dict0[default_usecase_method_str]
const default_usecase_arg_specs = get_args(default_usecase_method)
const default_usecase_kwarg_specs = get_kwargs(default_usecase_method)

# Convert arg specs to Vector{Dict} for UI iteration
function args_to_list(arg_specs::Dict{Symbol,<:Any})
    [Dict{String,Any}("name" => string(k), "type" => string(v)) for (k, v) in arg_specs]
end

# Convert kwarg specs to Vector{Dict} for UI iteration (includes default value)
function kwargs_to_list(kwarg_specs::Dict{Symbol,<:Tuple})
    [Dict{String,Any}(
        "name" => string(k),
        "type" => string(spec[1]),
        "default" => length(spec) >= 2 && spec[2] !== NODEFAULT ? format_value(spec[2]) : ""
    ) for (k, spec) in kwarg_specs]
end

const default_usecase_args_list = args_to_list(default_usecase_arg_specs)
const default_usecase_kwargs_list = kwargs_to_list(default_usecase_kwarg_specs)
const default_usecase_arg_values = Dict(string(k) => "" for (k, _) in default_usecase_arg_specs)
const default_usecase_kwarg_values = Dict(
    string(k) => begin
        T, d = get_kw_type_and_default(spec)
        d === NODEFAULT ? "" : format_value(d)
    end for (k, spec) in default_usecase_kwarg_specs
)

const default_ms_dict = ms_dict0

# == Reactive code ==
# add reactive code to make the UI interactive
@app begin
    #println(valid_cases)
    # == Reactive variables ==
    # reactive variables exist in both the Julia backend and the browser with two-way synchronization
    # @out variables can only be modified by the backend
    # @in variables can be modified by both the backend and the browser
    # variables must be initialized with constant values, or variables defined outside of the @app block
    @in usecase = default_usecase
    @out usecases = valid_cases
    @in usecase_method_str = default_usecase_method_str
    @out usecase_method_strs = default_usecase_method_strs
    @private usecase_method = default_usecase_method
    @private usecase_args = default_usecase_arg_specs
    @out usecase_args_list = default_usecase_args_list
    @private usecase_kwargs = default_usecase_kwarg_specs
    @out usecase_kwargs_list = default_usecase_kwargs_list
    @in usecase_arg_values   = default_usecase_arg_values      # Dict(String => String)
    @in usecase_kwarg_values = default_usecase_kwarg_values    # Dict(String => String)

    # Parsed, typed values (for actually calling FUSE.case_parameters)
    @out parsed_usecase_args   = EMPTY_SYMANY_DICT
    @out parsed_usecase_kwargs = EMPTY_SYMANY_DICT

    # Per-field error messages
    @out usecase_arg_errors   = EMPTY_STRING_DICT
    @out usecase_kwarg_errors = EMPTY_STRING_DICT

    @private ms_dict = default_ms_dict

    # Session state (holds initialized FUSE objects)
    @private session = WebGuiSession()

    # Initialize button
    @in initialize_clicked = false
    @out is_initializing = false
    @out init_status = ""
    @out init_error = ""

    # Log display
    @out log_text = ""
    @in clear_log_clicked = false


    # == Reactive handlers ==
    # reactive handlers watch a variable and execute a block of code when
    # its value changes
    @onchange usecase begin
        ms =  find_case_methods(usecase)
        mstrs, ms_dict = get_method_strdict(ms, usecase)
        usecase_method_strs = mstrs
        usecase_method_str = mstrs[1]
    end


    @onchange usecase_method_str begin
        usecase_method = ms_dict[usecase_method_str]
        usecase_args = get_args(usecase_method)
        usecase_args_list = args_to_list(usecase_args)
        usecase_kwargs = get_kwargs(usecase_method)
        usecase_kwargs_list = kwargs_to_list(usecase_kwargs)

        # reset UI-facing string values based on new specs
        usecase_arg_values = Dict(string(k) => "" for (k, _) in usecase_args)
        usecase_kwarg_values = Dict(
            string(k) => begin
                T, d = get_kw_type_and_default(spec)
                d === NODEFAULT ? "" : format_value(d)
            end for (k, spec) in usecase_kwargs
        )

        # clear parsed values and errors
        parsed_usecase_args   = EMPTY_SYMANY_DICT
        parsed_usecase_kwargs = EMPTY_SYMANY_DICT
        usecase_arg_errors    = EMPTY_STRING_DICT
        usecase_kwarg_errors  = EMPTY_STRING_DICT
    end

    @onchange usecase_arg_values begin
        new_parsed = Dict{Symbol,Any}()
        new_errors = Dict{String,String}()

        for (k_str, s) in usecase_arg_values
            k = Symbol(k_str)
            T = usecase_args[k]
            s_trim = strip(s)

            if isempty(s_trim)
                new_errors[k_str] = "Required"
                continue
            end

            try
                new_parsed[k] = parse_value(T, s_trim)
            catch e
                new_errors[k_str] = sprint(showerror, e)
            end
        end

        parsed_usecase_args = new_parsed
        usecase_arg_errors = new_errors
    end

    @onchange usecase_kwarg_values begin
        new_parsed = Dict{Symbol,Any}()
        new_errors = Dict{String,String}()

        for (k_str, s) in usecase_kwarg_values
            k    = Symbol(k_str)
            spec = usecase_kwargs[k]              # (Type, [default])
            T, default = get_kw_type_and_default(spec)
            s_trim = strip(s)

            if isempty(s_trim)
                if default === NODEFAULT
                    # no default and no user value; treat as omitted or required
                    # If required, uncomment:
                    # new_errors[k_str] = "Required"
                else
                    new_parsed[k] = default
                end
                continue
            end

            try
                new_parsed[k] = parse_value(T, s_trim)
            catch e
                new_errors[k_str] = sprint(showerror, e)
            end
        end

        parsed_usecase_kwargs = new_parsed
        usecase_kwarg_errors = new_errors
    end

    @onchange initialize_clicked begin
        if initialize_clicked
            # Check for errors before proceeding
            if !isempty(usecase_arg_errors) || !isempty(usecase_kwarg_errors)
                init_status = "error"
                init_error = "Please fix validation errors before initializing"
                initialize_clicked = false
                return
            end

            is_initializing = true
            init_status = ""
            init_error = ""

            # Build args tuple from parsed values (in order)
            arg_values = Tuple(parsed_usecase_args[k] for (k, _) in usecase_args)

            # Convert kwargs to Dict{Symbol,Any}
            kwargs = Dict{Symbol,Any}(parsed_usecase_kwargs)

            # Load the case using the session helper
            success = load_case!(session, usecase, arg_values, kwargs)

            if success
                init_status = "success"
            else
                init_status = "error"
                init_error = isempty(session.log_buffer) ? "Unknown error" : session.log_buffer[end]
            end

            # Update log display
            log_text = join(session.log_buffer, "\n")

            is_initializing = false
            initialize_clicked = false
        end
    end

    @onchange clear_log_clicked begin
        if clear_log_clicked
            clear_log!(session)
            log_text = ""
            clear_log_clicked = false
        end
    end

end

# == Pages ==
# register a new route and the page that will be loaded on access
@page("/", "app.jl.html")
end