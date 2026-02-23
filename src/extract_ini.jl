import OrderedCollections

"""
    extract_ini(ini::ParametersAllInits)

Extract all leaf parameters from `ini` into an `OrderedDict{Symbol, IMAS.ExtractFunction}`,
compatible with `IMAS.extract(dd, :all)` output so results can be `merge!`'d before DataFrame conversion.

Skips leaves under `:build` (messy structure) and non-DataFrame-friendly value types
(`missing`, `Function`, `TimeData`, `Vector`, `AbstractRange`).

Column names follow the pattern `ini_equilibrium_R0`, `ini_nb_unit_1_power_launched`, etc.
"""
function extract_ini(ini::ParametersAllInits)
    xtract = OrderedCollections.OrderedDict{Symbol,IMAS.ExtractFunction}()

    for leaf in SimulationParameters.leaves(ini)
        p = SimulationParameters.path(leaf)

        # skip leaves under :build
        if :build in p
            continue
        end

        val = leaf.value

        # skip non-dataframe-friendly types
        if ismissing(val) || val isa Function || val isa TimeData || val isa AbstractVector || val isa AbstractRange
            continue
        end

        # keep Real (includes Bool, Int, Float64), String, Symbol
        if !(val isa Real || val isa AbstractString || val isa Symbol)
            continue
        end

        # convert Symbol values to String for DataFrame compatibility
        if val isa Symbol
            val = String(val)
        end

        col_name = Symbol("ini_" * join(string.(p[2:end]), "_"))
        xfun = IMAS.ExtractFunction(:ini, col_name, leaf.units, _ -> val)
        xfun.value = val
        xtract[col_name] = xfun
    end

    return xtract
end

"""
    extract_act(act::ParametersAllActors)

Extract `act.ActorTGLF.tglfnn_model` into an `OrderedDict{Symbol, IMAS.ExtractFunction}`,
compatible with `IMAS.extract(dd, :all)` output so results can be `merge!`'d before DataFrame conversion.
"""
function extract_act(act::ParametersAllActors)
    xtract = OrderedCollections.OrderedDict{Symbol,IMAS.ExtractFunction}()

    entry = getfield(act.ActorTGLF, :tglfnn_model)
    val = entry.value
    if !ismissing(val)
        if val isa Symbol
            val = String(val)
        end
        col_name = :act_ActorTGLF_tglfnn_model
        xfun = IMAS.ExtractFunction(:act, col_name, entry.units, _ -> val)
        xfun.value = val
        xtract[col_name] = xfun
    end

    return xtract
end
