txt = open(joinpath(@__DIR__, "..", "..", "README.md"), "r") do io
    return read(io, String)
end

txt = replace(txt,
    "## Documentation\n\nFind the full documentation here: [https://fuse.help](https://fuse.help)" => "",
    "docs/src/assets/" => "assets/")

open(joinpath(@__DIR__, "index.md"), "w") do io
    write(io, txt)
    return write(io, "\nLast update on $(Dates.now())")
end
