#!/bin/bash
# Acceptance checks for the FUSE container on NERSC Perlmutter (podman-hpc).
#
#   FUSE_ENVIRONMENT=v1.2.0 ./deploy/perlmutter-container/acceptance.sh
#
# Runs against the migrated image in the per-user podman-hpc store, or a shared
# squash store when SQUASH_DIR is set. Checks mirror the omega suite
# (../omega-container/acceptance.sh) minus the SIF/launcher and /fusion-path
# checks that do not apply on NERSC.
#
# Optional:
#   SQUASH_DIR=/global/cfs/cdirs/m3739/shared_images   use the shared image store
#   SKIP_SLOW=1                                        skip the flux-matcher solve
#
# Each check prints PASS/FAIL independently; the summary exits nonzero if any
# failed. Checks assert the precondition they are actually testing, so a check
# cannot pass while silently exercising nothing.

set -uo pipefail

if [[ -n "${FUSE_ENVIRONMENT:-}" ]]; then
    version="$FUSE_ENVIRONMENT"
else
    version="$(curl -s https://api.github.com/repos/ProjectTorreyPines/FUSE.jl/releases/latest | jq -r .name)"
fi
[[ -n "$version" && "$version" != "null" ]] || { echo "ERROR: set FUSE_ENVIRONMENT" >&2; exit 1; }

command -v podman-hpc >/dev/null || { echo "ERROR: podman-hpc not found — run on Perlmutter" >&2; exit 1; }

image="localhost/fuse:$version"
out="${SCRATCH:?}/fuse_acceptance_${version}"
mkdir -p "$out"

podman=(podman-hpc)
[[ -n "${SQUASH_DIR:-}" ]] && podman=(podman-hpc --squash-dir "$SQUASH_DIR")

pass=0; fail=0; failed_names=()
check() {
    local name="$1"; shift
    echo "----- $name"
    if "$@"; then echo "PASS $name"; pass=$((pass+1))
    else echo "FAIL $name"; fail=$((fail+1)); failed_names+=("$name"); fi
}

echo "=== FUSE $version container acceptance — $(date -Is) on $(hostname)"
echo "=== image: $image${SQUASH_DIR:+ (squash-dir: $SQUASH_DIR)}"
echo

# --- identity ---------------------------------------------------------------
check "version-is-$version" bash -c "
    got=\$(${podman[*]} run --rm $image fuse -e 'import FUSE; print(pkgversion(FUSE))')
    echo \"  FUSE \$got\"; [ \"v\$got\" = \"$version\" ]"

check "sysimage-loaded" bash -c "
    p=\$(${podman[*]} run --rm $image fuse -e 'print(unsafe_string(Base.JLOptions().image_file))')
    echo \"  sysimage: \$p\"; [ \"\$p\" = /opt/fuse/sys_fuse.so ]"

# The `fuse` wrapper prepends $HOME/.julia_fuse_container when the baked depot
# is not writable (read-only squash rootfs). Either way, the effective first
# depot entry must be writable or runtime precompilation of new packages fails.
check "writable-depot" bash -c "
    ${podman[*]} run --rm $image fuse -e '
        d = first(DEPOT_PATH)
        println(\"  first depot: \", d)
        mkpath(d)
        p = joinpath(d, \".acceptance_write_test\")
        touch(p); rm(p)
        println(\"  depot is writable\")'"

check "core-package-versions" bash -c "
    ${podman[*]} run --rm $image fuse -e '
        import Pkg
        found = Dict(d.name => d.version for d in values(Pkg.dependencies()))
        for p in [\"FUSE\",\"IMAS\",\"IMASdd\",\"IMASggd\",\"TurbulentTransport\"]
            println(\"  \", p, \" => \", get(found, p, \"NOT INSTALLED\"))
        end'"

# --- image hygiene ----------------------------------------------------------
# Scan in Julia, not nested bash: through podman -> bash -c -> fuse -e the
# quoting is three deep, and a quoting slip reads as "the image still has LFS
# stubs". Collect into arrays -- `n += 1` in a top-level loop hits Julia soft
# scope and throws UndefVarError instead of counting.
check "no-LFS-pointer-stubs" bash -c "
    ${podman[*]} run --rm $image fuse -e '
        import TurbulentTransport
        root = joinpath(pkgdir(TurbulentTransport), \"models\")
        isdir(root) || (println(\"  no models dir at \", root); exit(1))
        files = String[]; stubs = String[]
        for (d, _, fs) in walkdir(root), f in fs
            p = joinpath(d, f); push!(files, p)
            startswith(String(open(io -> read(io, 40), p)), \"version https://git-lfs\") && push!(stubs, p)
        end
        isempty(files) && (println(\"  scanned nothing — models dir empty\"); exit(1))
        println(\"  scanned \", length(files), \" model file(s); \", length(stubs), \" still LFS stubs\")
        for s in first(stubs, 5); println(\"    STUB \", s); end
        exit(isempty(stubs) ? 0 : 1)'"

check "registries-world-readable" bash -c "
    bad=\$(${podman[*]} run --rm $image bash -c 'find /opt/fuse/.julia/registries ! -perm -o+r | wc -l')
    echo \"  non-world-readable entries: \$bad\"; [ \"\$bad\" = 0 ]"

# Revise is not a FUSE dependency: it is installed explicitly because
# FuseExamples notebooks open with `using Revise` and this image serves them
# through its Jupyter kernel.
check "revise-loads" bash -c "
    ${podman[*]} run --rm $image fuse -e 'using Revise; println(\"  Revise \", pkgversion(Revise))'"

check "host-tools-present" bash -c "
    ${podman[*]} run --rm $image bash -c 'for t in ps ssh rsync; do command -v \$t || { echo \"  missing \$t\"; exit 1; }; done'"

# --- physics smoke ----------------------------------------------------------
check "init-and-plot-png" bash -c "
    ${podman[*]} run --rm --volume $out:/out $image fuse -e '
        using FUSE, Plots
        ini, act = FUSE.case_parameters(:D3D, :L_mode)
        dd = FUSE.init(ini, act)
        plot(dd.equilibrium)
        savefig(\"/out/equilibrium.png\")' && ls -la $out/equilibrium.png"

if [[ "${SKIP_SLOW:-0}" != "1" ]]; then
    check "flux-matcher-offline" bash -c "
        ${podman[*]} run --rm --network none $image fuse -e '
            using FUSE
            ini, act = FUSE.case_parameters(:D3D, :L_mode)
            dd = FUSE.init(ini, act)
            FUSE.ActorFluxMatcher(dd, act)
            println(\"  ActorFluxMatcher completed with no network\")'"
fi

echo
echo "=== acceptance summary: $pass passed, $fail failed"
if [[ "$fail" -gt 0 ]]; then
    printf '    FAILED: %s\n' "${failed_names[@]}"
    exit 1
fi
echo "=== ALL CHECKS PASSED — $image is good on $(hostname)"
