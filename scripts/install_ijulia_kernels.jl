#!/usr/bin/env julia
"""
Install IJulia kernels for FUSE without requiring `jupyter` on PATH.

Kernels are written under JUPYTER_DATA_DIR/kernels when set, otherwise
~/.local/share/jupyter/kernels. Pass specdir via IJulia.installkernel so
registration does not depend on the `jupyter` executable.
"""

using Pkg

const KERNEL_PACKAGES = ["Plots", "IJulia", "WebIO", "Interact"]

function kernels_parent_dir()
    if haskey(ENV, "JUPYTER_DATA_DIR") && !isempty(ENV["JUPYTER_DATA_DIR"])
        return joinpath(ENV["JUPYTER_DATA_DIR"], "kernels")
    end
    return joinpath(homedir(), ".local", "share", "jupyter", "kernels")
end

function kernel_slug(name::AbstractString)
    slug = lowercase(replace(strip(name), r"[^a-z0-9]+" => "-"))
    return isempty(slug) ? "julia" : slug
end

function install_one_kernel!(display_name::AbstractString, nthreads::AbstractString)
    import IJulia
    parent = kernels_parent_dir()
    specdir = joinpath(parent, kernel_slug(display_name))
    mkpath(specdir)
    IJulia.installkernel(
        display_name;
        specdir,
        env=Dict("JULIA_NUM_THREADS" => nthreads),
    )
    println("Installed kernel \"", display_name, "\" at ", specdir)
    return specdir
end

Pkg.add(KERNEL_PACKAGES)
Pkg.build("IJulia")
import IJulia

n = get(ENV, "JULIA_NUM_THREADS", string(length(Sys.cpu_info())))
install_one_kernel!("Julia (1 thread)", "1")
install_one_kernel!("Julia ($n threads)", n)

println()
println("Kernels directory: ", kernels_parent_dir())
println("List kernels with: python -m jupyter kernelspec list")
