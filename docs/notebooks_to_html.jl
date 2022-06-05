using Pkg
Pkg.activate("..")

using ProgressMeter
# Converts all notebooks in examples/cases to html and stores them in docs/build/assets
current_path = dirname(abspath(@__FILE__))
example_cases = joinpath(current_path, "..", "examples", "cases")
files_to_convert = readdir(example_cases)[findall(endswith(".ipynb"), readdir(example_cases))]

@showprogress for case in files_to_convert
    casename = split(case, ".")[1]
    srcname = joinpath(example_cases, "$casename.md")
    srcfiles = joinpath(example_cases, casename * "_files")
    dstname = joinpath(current_path, "src", "example_$casename.md")
    dstfiles = joinpath(current_path, "src", "assets", casename * "_files")
    if !isfile(dstname)
        run(`jupyter nbconvert --execute --to markdown $(joinpath(example_cases, case))`)
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
