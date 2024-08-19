txt = open(joinpath(@__DIR__, "..", "..", "README.md"), "r") do io
    return read(io)
end


open(joinpath(@__DIR__, "index.md"), "w") do io
    write(io, txt)
    return write(io, "\nLast update on $(Dates.now())")
end
