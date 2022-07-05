using Pkg
Pkg.activate("..")
using ProgressMeter

dirs = ["cases", "actors", "workflows", "tutorials"]

# Converts all notebooks in examples/cases to html and stores them in docs/build/assets
current_path = dirname(abspath(@__FILE__))

for dir in dirs
    example_folder = joinpath(current_path, "..", "examples", dir)
    files_to_convert = readdir(example_folder)[findall(x->endswith(x,".ipynb"), readdir(example_folder))]

    @showprogress for case in files_to_convert
        ipynb = joinpath(example_folder, case)
        casename = split(case, ".")[1]
        srcname = joinpath(example_folder, "$casename.md")
        srcfiles = joinpath(example_folder, casename * "_files")
        dstname = joinpath(current_path, "src", "example_$(dir)__$(casename).md")
        dstfiles = joinpath(current_path, "src", "assets", "$(casename)_files")
        
        if isfile(dstname)
            println("$dstname exists: skipping nbconvert")
        else
            run(`jupyter nbconvert --execute --to markdown $ipynb`)
            run(`rm -rf $dstfiles`)
            if isdir(srcfiles)
                run(`mv -f $srcfiles $dstfiles`)
            end
        end
        
        if isfile(srcname)
            run(`cp -f $srcname $dstname`)
            txt = open(dstname, "r") do io
                txt = read(io, String)
                txt = replace(txt, "$(casename)_files" => "assets/$(casename)_files")
                txt = replace(txt, "assets/assets/$(casename)_files" => "assets/$(casename)_files")
                txt = replace(txt, r"\[[0-9]+m" => "")
                return txt = replace(txt, r"```julia" => "```@julia")
            end
            open(dstname, "w") do io
                return write(io, txt)
            end
        end
    end

end
