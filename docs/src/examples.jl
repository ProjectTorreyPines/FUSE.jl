assets = readdir(joinpath(dirname(abspath(@__FILE__)), "..", "build", "assets"))
examples = [split(item[9:end], ".")[1] for item in assets[findall(startswith("example_"), assets)]]
txt = ["""
  The following examples are available:
  """]
run(`pwd`) # make .html file 
for example in examples
    content = open("build/assets/example_$example.html", "r") do io
        read(io, String)
    end
    open("src/example_$example.md", "w") do io
        write(
            io,
            """
  # $example

  ```@raw html
  $content
  ```
  """
        )
    end

    push!(
        txt,
        """
* [$example](example_$example.md)
"""
    )
end

open("src/examples.md", "w") do io
    write(io, join(txt, "\n"))
end
