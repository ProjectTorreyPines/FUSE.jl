
get_types_list(k) = Dict(i => Symbol(s) for (i, s) in enumerate(join.(collect(Iterators.product(ntuple(_ -> 'A':'Z', k)...))[:])))

#generate a dict of 3 characters symbol 1 => :AAA, 2=> :BBB, =>:CCC
types_dic = get_types_list(3)

function _kwdef!(blk, params_args, call_args)
    for i in eachindex(blk.args)
        ei = blk.args[i]
        if ei isa Symbol
            #  var
            push!(params_args, ei)
            push!(call_args, ei)
        elseif ei isa Expr
            is_atomic = ei.head === :atomic
            ei = is_atomic ? first(ei.args) : ei # strip "@atomic" and add it back later
            is_const = ei.head === :const
            ei = is_const ? first(ei.args) : ei # strip "const" and add it back later
            # Note: `@atomic const ..` isn't valid, but reconstruct it anyway to serve a nice error
            if ei isa Symbol
                # const var
                push!(params_args, ei)
                push!(call_args, ei)
            elseif ei.head === :(=)
                lhs = ei.args[1]
                if lhs isa Symbol
                    #  var = defexpr
                    var = lhs
                elseif lhs isa Expr && lhs.head === :(::) && lhs.args[1] isa Symbol
                    #  var::T = defexpr
                    var = lhs.args[1]
                else
                    # something else, e.g. inline inner constructor
                    #   F(...) = ...
                    continue
                end
                defexpr = ei.args[2]  # defexpr
                push!(params_args, Expr(:kw, var, esc(defexpr)))
                push!(call_args, var)
                lhs = is_const ? Expr(:const, lhs) : lhs
                lhs = is_atomic ? Expr(:atomic, lhs) : lhs
                blk.args[i] = lhs # overrides arg
            elseif ei.head === :(::) && ei.args[1] isa Symbol
                # var::Typ
                var = ei.args[1]
                push!(params_args, var)
                push!(call_args, var)
            elseif ei.head === :block
                # can arise with use of @static inside type decl
                _kwdef!(ei, params_args, call_args)
            end
        end
    end
    blk
end

function _get_type!(blk)
    types = []

    for i in eachindex(blk.args)
        ei = blk.args[i]
        if ei isa Symbol

        elseif ei isa Expr
            is_atomic = ei.head === :atomic
            ei = is_atomic ? first(ei.args) : ei # strip "@atomic" and add it back later
            is_const = ei.head === :const
            ei = is_const ? first(ei.args) : ei # strip "const" and add it back later
            # Note: `@atomic const ..` isn't valid, but reconstruct it anyway to serve a nice error
            if ei isa Symbol

            elseif ei.head === :(=)
                lhs = ei.args[1]
                if lhs isa Symbol

                elseif lhs isa Expr && lhs.head === :(::) && lhs.args[1] isa Symbol

                    push!(types, lhs.args[2])
                else
                    # something else, e.g. inline inner constructor
                    #   F(...) = ...
                    continue
                end
            elseif ei.head === :(::) && ei.args[1] isa Symbol
                push!(types, ei.args[2])

            elseif ei.head === :block
                # can arise with use of @static inside type decl
                _get_type!(ei)
            end
        end
    end
    types
end


function get_new_type(types, existing_types)
    i = 1
    while (types_dic[i] ∈ types || types_dic[i] ∈ existing_types)
        i = i + 1
    end
    push!(types, types_dic[i])
    return types_dic[i]
end


function auto_type!(blk, types, existing_types)
    for i in eachindex(blk.args)
        ei = blk.args[i]
        if ei isa Symbol
            blk.args[i] = Expr(:(::), ei, get_new_type(types, existing_types))
        elseif ei isa Expr
            is_atomic = ei.head === :atomic
            ei = is_atomic ? first(ei.args) : ei # strip "@atomic" and add it back later
            is_const = ei.head === :const
            ei = is_const ? first(ei.args) : ei # strip "const" and add it back later
            # Note: `@atomic const ..` isn't valid, but reconstruct it anyway to serve a nice error
            if ei isa Symbol
                # const var
                blk.args[i] = Expr(:(::), ei, get_new_type(types, existing_types))
            elseif ei.head === :(=)
                lhs = ei.args[1]
                if lhs isa Symbol
                    #  var = defexpr
                    ei.args[1] = Expr(:(::), lhs, get_new_type(types, existing_types))
                elseif lhs isa Expr && lhs.head === :(::) && lhs.args[1] isa Symbol
                    #  var::T = defexpr
                else
                    # something else, e.g. inline inner constructor
                    #   F(...) = ...
                    continue
                end
            elseif ei.head === :(::) && ei.args[1] isa Symbol

            elseif ei.head === :block
                # can arise with use of @static inside type decl
                auto_type!(ei, types, existing_types)
            end
        end
    end
    blk
end

macro new_actor(expr)
    expr = macroexpand(__module__, expr) # to expand @static
    @assert !Base.isexpr(expr, :struct) "Invalid usage of @new_actor\n. Usage is @new_actor MyActorName"
    actor_name = expr
    actor_parameter_type = Symbol("FUSEparameters__" * string(actor_name))
    #make it mut
    expr1 = quote
        mutable struct $(actor_name){P<:$(actor_parameter_type)} <: ReactorAbstractActor
            dd::IMAS.dd
            par::P
        end
        
    end

    expr2 = quote    
        function $(actor_name)(dd::IMAS.dd, act::ParametersAllActors; kw...)
            par = act.$(actor_name)(kw...)
            actor = $(actor_name)(dd, par)
            step(actor)
            finalize(actor)
            return actor
        end
    end

    expr3 = quote
        function $(actor_name)(dd::IMAS.dd, par::$(actor_parameter_type); kw...)
            logging_actor_init($(actor_name))
            par = par(kw...)
            new{IMAS.dd,typeof(par)}(dd, par)
        end   
    end
    expr_final = quote
        $(esc(expr1))
        Base.@__doc__ $(esc(expr2))
        $(esc(expr3))
    end
    println("-----------\n", expr_final)
    return expr_final
end
macro test_kwdef(expr)
    expr = macroexpand(__module__, expr) # to expand @static
    Base.isexpr(expr, :struct) || error("Invalid usage of @kwdef")
    T = expr.args[2]
    if T isa Expr && T.head === :<:
        T = T.args[1]
    end

    params_ex = Expr(:parameters)
    call_args = Any[]

    _kwdef!(expr.args[3], params_ex.args, call_args)
    # Only define a constructor if the type has fields, otherwise we'll get a stack
    # overflow on construction
    if !isempty(params_ex.args)
        if T isa Symbol
            sig = :(($(esc(T)))($params_ex))
            call = :(($(esc(T)))($(call_args...)))
            body = Expr(:block, __source__, call)
            kwdefs = Expr(:function, sig, body)
        elseif Base.isexpr(T, :curly)
            # if T == S{A<:AA,B<:BB}, define two methods
            #   S(...) = ...
            #   S{A,B}(...) where {A<:AA,B<:BB} = ...
            S = T.args[1]
            P = T.args[2:end]
            Q = Any[Base.isexpr(U, :<:) ? U.args[1] : U for U in P]
            SQ = :($S{$(Q...)})
            body1 = Expr(:block, __source__, :(($(esc(S)))($(call_args...))))
            sig1 = :(($(esc(S)))($params_ex))
            def1 = Expr(:function, sig1, body1)
            body2 = Expr(:block, __source__, :(($(esc(SQ)))($(call_args...))))
            sig2 = :(($(esc(SQ)))($params_ex) where {$(esc.(P)...)})
            def2 = Expr(:function, sig2, body2)
            kwdefs = Expr(:block, def1, def2)
        else
            error("Invalid usage of @kwdef")
        end
    else
        kwdefs = nothing
    end
    expr_final = quote
        Base.@__doc__ $(esc(expr))
        $kwdefs
    end
    println("-----------\n", expr_final)
    return expr_final
end

macro actor_params(expr)
    expr = macroexpand(__module__, expr) # to expand @static
    Base.isexpr(expr, :struct) || error("Invalid usage of @kwdef")

    #make it mutable 
    expr.args[1] = true

    # add default parent type
    T0 = expr.args[2]
    if !(T0 isa Expr && T0.head === :<:)
        T0 = Expr(:<:, T0, :ParametersActor)
    else
        error("$(expr.args[2]) has already a supertype... The macro @actor_params automatically adds the supertype ParametersActor")
    end


    # Get struct name and   


    # add FUSEparameters__ prefix
    T = T0.args[1]
    if T isa Expr && T.head == :curly
        struct_name = string(T.args[1])
        full_name = "FUSEparameters__" * string(struct_name)
        T.args[1] = Symbol(full_name)
    else
        struct_name = string(T0.args[1])
        full_name = "FUSEparameters__" * string(T0.args[1])
        T0.args[1] = Symbol(full_name)
    end

    # replace struct header
    #println("struct name: ", struct_name)
    expr.args[2] = T0
    #println("expr: ", expr.args[2])

    #add _name and _parent fields
    insert!(expr.args[3].args, 1, :(_name::Symbol = :not_set))
    insert!(expr.args[3].args, 1, :(_parent::WeakRef = WeakRef(nothing)))
    #
    params_ex = Expr(:parameters)
    call_args = Any[]

    # collect current types 
    T = expr.args[2].args[1]
    existing_types = _get_type!(expr.args[3])
    if T isa Expr && T.head == :curly
        current_types = T.args[2:end]
    else
        current_types = []
    end
    # add types to fields without get_types_list
    auto_type!(expr.args[3], current_types, existing_types)

    # set the structure fields with added types  
    if length(current_types) > 0
        if T isa Symbol
            T = Expr(:curly, T)
        else
            while length(T.args) > 1
                deleteat!(T.args, length(T.args))
            end
        end
        for c in current_types
            push!(T.args, c)
        end
    end
    expr.args[2].args[1] = T

    # create constructors 
    Base._kwdef!(expr.args[3], params_ex.args, call_args)

    # Only define a constructor if the type has fields, otherwise we'll get a stack
    # overflow on construction
    if !isempty(params_ex.args)
        if T isa Symbol
            sig = :(($(esc(T)))($params_ex))
            call = :(($(esc(T)))($(call_args...)))
            body = Expr(:block, __source__, call)
            kwdefs = Expr(:function, sig, body)
        elseif Base.isexpr(T, :curly)
            # if T == S{A<:AA,B<:BB}, define two methods
            #   S(...) = ...
            #   S{A,B}(...) where {A<:AA,B<:BB} = ...
            S = T.args[1]
            P = T.args[2:end]
            Q = Any[Base.isexpr(U, :<:) ? U.args[1] : U for U in P]
            SQ = :($S{$(Q...)})
            body1 = Expr(:block, __source__, :(($(esc(S)))($(call_args...))))
            sig1 = :(($(esc(S)))($params_ex))
            def1 = Expr(:function, sig1, body1)
            body2 = Expr(:block, __source__, :(($(esc(SQ)))($(call_args...))))
            sig2 = :(($(esc(SQ)))($params_ex) where {$(esc.(P)...)})
            def2 = Expr(:function, sig2, body2)
            kwdefs = Expr(:block, def1, def2)
        else
            error("Invalid usage of @kwdef")
        end
    else
        kwdefs = nothing
    end

    #add constructor shortcurts MyStruct() = FUSEparameters__MyStruct() if not FuseParameter_ActorName to prevent overloading struct ActorName
    if !startswith(struct_name,"Actor")
        expr_f = quote 
        $(Symbol(struct_name))() = $(Symbol(full_name))()
        end
        push!(kwdefs.args, esc(expr_f))
    end

    #add doc string for struct
    expr_final = quote
        Base.@__doc__ $(esc(expr))
        $kwdefs
    end
    
    #println("-----------\n", expr_final)
    return expr_final
end

macro init_params(expr)
    expr = macroexpand(__module__, expr) # to expand @static
    Base.isexpr(expr, :struct) || error("Invalid usage of @kwdef")

    #make it mutable 
    expr.args[1] = true

    # add default parent type
    T0 = expr.args[2]
    if !(T0 isa Expr && T0.head === :<:)
        T0 = Expr(:<:, T0, ParametersInit)
    else
        error("$(expr.args[2]) has already a supertype... The macro @actor_params automatically adds the supertype ParametersActor")
    end


    # Get struct name and   


    # add FUSEparameters__ prefix
    T = T0.args[1]
    if T isa Expr && T.head == :curly
        struct_name = string(T.args[1])
        full_name = "FUSEparameters__" * string(struct_name)
        T.args[1] = Symbol(full_name)
    else
        struct_name = string(T0.args[1])
        full_name = "FUSEparameters__" * string(T0.args[1])
        T0.args[1] = Symbol(full_name)
    end

    # replace struct header
    #println("struct name: ", struct_name)
    expr.args[2] = T0
    #println("expr: ", expr.args[2])

    #add _name and _parent fields
    insert!(expr.args[3].args, 1, :(_name::Symbol = :not_set))
    insert!(expr.args[3].args, 1, :(_parent::WeakRef = WeakRef(nothing)))
    #
    params_ex = Expr(:parameters)
    call_args = Any[]

    # collect current types 
    T = expr.args[2].args[1]
    existing_types = _get_type!(expr.args[3])
    if T isa Expr && T.head == :curly
        current_types = T.args[2:end]
    else
        current_types = []
    end
    # add types to fields without get_types_list
    auto_type!(expr.args[3], current_types, existing_types)

    # set the structure fields with added types  
    if length(current_types) > 0
        if T isa Symbol
            T = Expr(:curly, T)
        else
            while length(T.args) > 1
                deleteat!(T.args, length(T.args))
            end
        end
        for c in current_types
            push!(T.args, c)
        end
    end
    expr.args[2].args[1] = T

    # create constructors 
    Base._kwdef!(expr.args[3], params_ex.args, call_args)

    # Only define a constructor if the type has fields, otherwise we'll get a stack
    # overflow on construction
    if !isempty(params_ex.args)
        if T isa Symbol
            sig = :(($(esc(T)))($params_ex))
            call = :(($(esc(T)))($(call_args...)))
            body = Expr(:block, __source__, call)
            kwdefs = Expr(:function, sig, body)
        elseif Base.isexpr(T, :curly)
            # if T == S{A<:AA,B<:BB}, define two methods
            #   S(...) = ...
            #   S{A,B}(...) where {A<:AA,B<:BB} = ...
            S = T.args[1]
            P = T.args[2:end]
            Q = Any[Base.isexpr(U, :<:) ? U.args[1] : U for U in P]
            SQ = :($S{$(Q...)})
            body1 = Expr(:block, __source__, :(($(esc(S)))($(call_args...))))
            sig1 = :(($(esc(S)))($params_ex))
            def1 = Expr(:function, sig1, body1)
            body2 = Expr(:block, __source__, :(($(esc(SQ)))($(call_args...))))
            sig2 = :(($(esc(SQ)))($params_ex) where {$(esc.(P)...)})
            def2 = Expr(:function, sig2, body2)
            kwdefs = Expr(:block, def1, def2)
        else
            error("Invalid usage of @kwdef")
        end
    else
        kwdefs = nothing
    end

    #add constructor shortcurts MyStruct() = FUSEparameters__MyStruct() if not FuseParameter_ActorName to prevent overloading struct ActorName
    if !startswith(struct_name,"Actor")
        expr_f = quote 
        $(Symbol(struct_name))() = $(Symbol(full_name))()
        end
        push!(kwdefs.args, esc(expr_f))
    end

    #add doc string for struct
    expr_final = quote
        Base.@__doc__ $(esc(expr))
        $kwdefs
    end
    
    println("-----------\n", expr_final)
    return expr_final
end