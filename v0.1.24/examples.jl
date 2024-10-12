examples = [split(item[9:end], ".")[1] for item in readdir((@__DIR__)) if startswith(item, "example_") && endswith(item, ".md")]
dirs = unique([split(item, "__")[1] for item in examples])
popat!(dirs, findfirst(x -> x == "cases", dirs))
pushfirst!(dirs, "cases")

txt = ["""
# Worked examples

The following examples are available:
"""]

for dir in dirs
    push!(txt, "* **$(titlecase(dir))**")
    for example in examples
        if !startswith(example, dir * "__")
            continue
        end

        title = open("$(@__DIR__)/example_$example.md", "r") do io
            txt = read(io, String)
            try
                title = txt[findfirst(r"^# .*", txt)][3:end]
            catch
                error("example_$example.md must have header markdown cell (`# title`)")
            end
            return title
        end

        push!(txt, "  - [$title](example_$example.md)")
    end
end

open("$(@__DIR__)/examples.md", "w") do io
    return write(io, join(txt, "\n"))
end
