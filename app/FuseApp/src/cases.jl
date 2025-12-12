#######################################################
# Utilities for getting case_parameter method options #
#######################################################
using CodeTracking

find_case_methods(val::Symbol) = case_methods[[Val{val} == cv for cv in case_values]]

tostring(m::Method) = string(CodeTracking.definition(m).args[1])

# Get argument names for a certain case_parameters method
function get_args(m::Method)
    names = Base.method_argnames(m)[3:end]
    types = Base.unwrap_unionall(m.sig).parameters[3:end]
    return Dict{Symbol,Type}(names .=> types)
end

function get_method_strdict(ms, usecase)
    margs = [get_args(m) for m in ms]
    fstr = "case_parameters(::Val{:$usecase}, "
    mstrs = [fstr * join(["$k::$v" for (k,v) in argdict], ", ") * ")" for argdict in margs]
    return mstrs, Dict(mstrs .=> ms)
end

eval_type(ex) = ex === :Any ? Any : Core.eval(@__MODULE__, ex)
eval_default(ex) = Core.eval(@__MODULE__, ex)

###################################
# Parsing argument spec from AST  #
###################################

function parse_arg_spec(node)
    if node isa Symbol
        # e.g. `x` or bare kw `alpha`
        return (node, :Any, false, nothing)

    elseif node isa Expr
        if node.head === :(::)
            # e.g. `x::Int` or `init_from::Symbol`
            var, ty = node.args
            var isa Symbol || error("Unsupported destructuring in arg $node")
            return (var, ty, false, nothing)

        elseif node.head === :(=)
            # e.g. `y = nothing` or `ne_setting::Symbol = :ne_ped`
            lhs, rhs = node.args
            name, ty, _, _ = parse_arg_spec(lhs)
            return (name, ty, true, rhs)

        elseif node.head === :kw
            # e.g. (kw beta 2.0)
            lhs, rhs = node.args
            if lhs isa Symbol
                return (lhs, :Any, true, rhs)
            else
                name, ty, _, _ = parse_arg_spec(lhs)
                return (name, ty, true, rhs)
            end
        else
            error("Unsupported argument node: $node")
        end
    else
        error("Unsupported argument node: $node")
    end
end

##########################################
# Split positional vs keyword arguments  #
##########################################

function split_positional_and_kw(call_expr::Expr)
    call_expr.head === :call || error("Signature must be a :call Expr")

    args = call_expr.args[2:end]  # drop function name
    idx = findfirst(a -> a isa Expr && (a::Expr).head === :parameters, args)

    if idx === nothing
        # no kwargs
        return Expr[]
    else
        params_expr = args[idx]::Expr  # Expr(:parameters, ...)
        kws = params_expr.args
        return kws
    end
end

function get_kwargs(m::Method)
    sig_expr = CodeTracking.definition(m).args[1]
    sig_expr.head === :call || error("Signature must be a :call Expr")

    kw_nodes = split_positional_and_kw(sig_expr)
    d = Dict{Symbol, Tuple}()

    for node in kw_nodes
        name, ty_ast, has_def, def_ast = parse_arg_spec(node)
        T = eval_type(ty_ast)
        if has_def
            d[name] = (T, eval_default(def_ast))
        else
            d[name] = (T,)
        end
    end

    return d
end