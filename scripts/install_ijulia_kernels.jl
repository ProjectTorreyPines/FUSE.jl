#!/usr/bin/env julia
"""
Install IJulia kernels for FUSE without requiring `jupyter` on PATH.

Kernels are written to IJulia's kernel directory: JUPYTER_DATA_DIR/kernels
when set, otherwise the platform default (%APPDATA%\\jupyter\\kernels on
Windows, ~/Library/Jupyter/kernels on macOS, ~/.local/share/jupyter/kernels
elsewhere). Registration does not depend on the `jupyter` executable.
"""

using Pkg

const KERNEL_PACKAGES = ["Plots", "IJulia", "WebIO", "Interact"]

function kernel_slug(name::AbstractString)
    slug = lowercase(strip(name))
    slug = replace(slug, r"[^a-z0-9]+" => "-")
    slug = strip(slug, '-')
    return isempty(slug) ? "julia" : slug
end

function install_one_kernel!(display_name::AbstractString, nthreads::AbstractString)
    specdir = IJulia.installkernel(
        display_name;
        specname=kernel_slug(display_name),
        displayname=display_name,
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
println("Kernels directory: ", IJulia.kerneldir())
println("List kernels with: python -m jupyter kernelspec list")
