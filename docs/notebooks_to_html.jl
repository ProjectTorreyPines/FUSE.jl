# Converts all notebooks in examples/cases to html and stores them in docs/build/assets
current_path= dirname(abspath(@__FILE__))
example_cases = joinpath(current_path,"..","examples","cases")
files_to_convert = readdir(example_cases)[findall(endswith(".ipynb"),readdir(example_cases))]

for case in files_to_convert
    run(`jupyter nbconvert --execute --to html $(joinpath(example_cases,case))`)
    run(`cp $(joinpath(example_cases,split(case,".")[1]*".html")) $(joinpath(current_path,"build","assets",string("example_",split(case,".")[1],".html")))`)
end