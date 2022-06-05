examples = [split(item[9:end], ".")[1] for item in readdir(dirname(abspath(@__FILE__))) if startswith(item,"example_")]
txt = ["""
# examples

The following examples are available:
"""]

for example in examples
    push!(
        txt,
        """
* [`$(replace(example,"_"=>" "))`](example_$example.md)
"""
    )
end

open("src/examples.md", "w") do io
    write(io, join(txt, "\n"))
end
