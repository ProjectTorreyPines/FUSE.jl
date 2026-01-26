using JSON

"""
    convert_notebook_to_litterate(ipynb_file::String, litterate_file::String)

Function to convert a Jupyter notebook (.ipynb) into a Litterate-compatible Julia script
"""
function convert_notebook_to_litterate(ipynb_file::String, litterate_file::String)
    # Read the Jupyter notebook file
    notebook_content = JSON.parsefile(ipynb_file)

    # Open the output Julia file for writing
    open(litterate_file, "w") do f
        # Loop through each cell in the notebook
        for cell in notebook_content["cells"]
            if cell["cell_type"] == "markdown"
                # Convert markdown cells to comments
                for line in cell["source"]
                    println(f, "# ", strip(line))
                end
            elseif cell["cell_type"] == "code"
                # Write code cells as-is
                for line in cell["source"]
                    println(f, strip(line))
                end
            end
            # Add a blank line between cells for readability
            println(f)
        end
    end
end

# ================================================================= #
# Command line execution
# ================================================================= #
if !isempty(ARGS)
    if length(ARGS) != 2
        println("Usage: julia script.jl <input.ipynb> <output.jl>")
        exit(1)
    end

    ipynb_file = ARGS[1]
    litterate_file = ARGS[2]

    # Check if input file exists
    if !isfile(ipynb_file)
        println("Error: Input file '$ipynb_file' not found")
        exit(1)
    end

    convert_notebook_to_litterate(ipynb_file, litterate_file)
    println("Converted '$ipynb_file' to '$litterate_file'")
end
