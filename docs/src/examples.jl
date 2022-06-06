examples = [split(item[9:end], ".")[1] for item in readdir(dirname(abspath(@__FILE__))) if startswith(item, "example_")]

txt = ["""
# examples

The following examples are available:
"""]

for example in examples

    title = open("src/example_$example.md", "r") do io
        txt = read(io, String)
        return txt[findfirst(r"^# .*", txt)][3:end]
    end

    push!(
        txt,
        """
* [$title](example_$example.md)
""",
    )
end

open("src/examples.md", "w") do io
    write(io, join(txt, "\n"))
end
