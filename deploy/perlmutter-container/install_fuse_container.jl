# Install FUSE + all ProjectTorreyPines packages into a container image and
# compile a PackageCompiler system image (sys_fuse.so).
#
# Adapted from deploy/perlmutter/install_fuse_environment.jl, with the Lmod
# module-file generation, IJulia kernelspec creation, and /global/common host
# paths removed. Those concerns are handled outside the container:
#   - the Jupyter kernel is a host-side kernel.json (see install_kernel.sh)
#   - module loading is replaced by `podman-hpc run`
#
# Driven entirely by environment variables set in the Containerfile:
#   FUSE_INSTALL_DIR   install/depot location          (default /opt/fuse)
#   JULIA_CPU_TARGET   sysimage CPU target              (required)
# Optional:
#   FUSE_MAKEFILE      path to FUSE Makefile to parse   (default /opt/build/Makefile)

@assert (Threads.nthreads() == 1) "Error: Installing FUSE sysimage requires running Julia with one thread"
@assert ("JULIA_CPU_TARGET" in keys(ENV)) "Error: Must define JULIA_CPU_TARGET environment variable"

install_dir = get(ENV, "FUSE_INSTALL_DIR", "/opt/fuse")
cpu_target = ENV["JULIA_CPU_TARGET"]
makefile = get(ENV, "FUSE_MAKEFILE", "/opt/build/Makefile")

import Pkg

println("### Setup main environment for installer")
Pkg.activate()
Pkg.Registry.add(Pkg.RegistrySpec(url="https://github.com/ProjectTorreyPines/FuseRegistry.jl.git"))
Pkg.Registry.add("General")
Pkg.add("PackageCompiler")
Pkg.update()  # also flips Pkg's "registry updated this session" flag so the
              # Pkg.add() below will not re-fetch (and revert) the registry edits
using PackageCompiler

println()
println("### Rewrite FuseRegistry SSH (git@github.com:) URLs to anonymous HTTPS")
# The FuseRegistry stores package repos as git@github.com: URLs. Julia's CLI git
# tries to clone those over SSH, which fails in a credential-free container
# (no ssh binary / no keys). Since all ProjectTorreyPines packages are public,
# rewrite the cloned registry's Package.toml repo fields to HTTPS so clones work
# anonymously. This is independent of git config / $HOME, so it always applies.
function rewrite_registry_ssh_to_https()
    n = 0
    for depot in DEPOT_PATH
        regroot = joinpath(depot, "registries")
        isdir(regroot) || continue
        for (root, _, files) in walkdir(regroot)
            for f in files
                f == "Package.toml" || continue
                path = joinpath(root, f)
                s = read(path, String)
                s2 = replace(s, "git@github.com:" => "https://github.com/")
                if s2 != s
                    chmod(path, 0o644)
                    write(path, s2)
                    n += 1
                end
            end
        end
    end
    return n
end
println("    rewrote ", rewrite_registry_ssh_to_https(), " Package.toml file(s)")

println()
println("### parse PTP packages from Makefile: ", makefile)
function get_packages_from_makefile(makefile)
    for line in eachline(makefile)
        if occursin(r"^FUSE_PACKAGES_MAKEFILE\s*:=", line)
            pkg_line = replace(line, r"^FUSE_PACKAGES_MAKEFILE\s*:=" => "")
            return split(strip(pkg_line))
        end
    end
    error("Could not find FUSE_PACKAGES_MAKEFILE in $makefile")
end
packages = get_packages_from_makefile(makefile)
pkgs_using = join(packages, ", ")
println("    ", packages)

println()
println("### Setup FUSE environment in ", install_dir)
Pkg.activate(install_dir)
Pkg.add([["FUSE", "Plots", "IJulia", "WebIO", "Interact", "EFIT", "ArgParse"]; packages])
Pkg.build("IJulia")

println()
println("### Eagerly install lazy artifacts (e.g. MKL) so the image runs offline")
# Pkg.add/instantiate skip artifacts marked `lazy = true`; those normally
# download on first @artifact_str access at runtime (observed: MKL). Walk every
# (Julia)Artifacts.toml in the depot and download all artifacts for this
# platform, including lazy ones, so they are baked into the image.
import Pkg.Artifacts: select_downloadable_artifacts, ensure_artifact_installed
function install_all_artifacts(depot)
    installed = 0
    pkgsdir = joinpath(depot, "packages")
    isdir(pkgsdir) || return installed
    for (root, _, files) in walkdir(pkgsdir)
        for f in files
            (f == "Artifacts.toml" || f == "JuliaArtifacts.toml") || continue
            toml = joinpath(root, f)
            local arts
            try
                arts = select_downloadable_artifacts(toml; include_lazy=true)
            catch err
                @warn "Could not parse $toml" exception = err
                continue
            end
            for (name, meta) in arts
                try
                    ensure_artifact_installed(name, meta, toml)
                    installed += 1
                catch err
                    @warn "Failed to install artifact $name from $toml" exception = err
                end
            end
        end
    end
    return installed
end
println("    installed/verified ", install_all_artifacts(first(DEPOT_PATH)), " downloadable artifact(s)")

println()
println("### Create precompile script")
precompile_execution_file = joinpath(install_dir, "precompile_script.jl")
precompile_cmds = """
using FUSE, EFIT, $pkgs_using
include(joinpath(pkgdir(FUSE), "docs", "src", "tutorial.jl"))
include(joinpath(pkgdir(FUSE), "test", "runtests.jl"))
"""
write(precompile_execution_file, precompile_cmds)

println()
println("### Precompile FUSE sys image")
sysimage_path = joinpath(install_dir, "sys_fuse.so")
create_sysimage(["FUSE"]; sysimage_path, precompile_execution_file, cpu_target)
chmod(sysimage_path, 0o555)

println()
println("### Record IJulia kernel.jl path for host-side kernelspec")
# The host kernel.json (see install_kernel.sh) needs the in-container path to
# IJulia's kernel.jl. Resolve and persist it next to the sysimage.
import IJulia
kernel_jl = joinpath(pkgdir(IJulia), "src", "kernel.jl")
@assert isfile(kernel_jl) "Could not locate IJulia kernel.jl at $kernel_jl"
write(joinpath(install_dir, "ijulia_kernel_path.txt"), kernel_jl)
println("    ", kernel_jl)

println()
println("### FUSE container install complete")
