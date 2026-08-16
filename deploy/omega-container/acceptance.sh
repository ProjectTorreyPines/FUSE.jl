#!/bin/bash
# Acceptance checks for a freshly built FUSE container, run before publishing.
#
#   FUSE_ENVIRONMENT=v1.2.0 ./deploy/omega-container/acceptance.sh
#
# Runs on the omega node that built the image (the podman store is node-local).
# Checks the podman image directly where possible, and the exported SIF through
# the fuse-container launcher for the parts that only matter read-only.
#
# Optional:
#   SIF=/path/to/fuse_<version>.sif   default: $SIF_DIR (or /local-scratch/$USER)
#   SKIP_SLOW=1                       skip the flux-matcher solve
#
# Each check prints PASS/FAIL independently; the summary exits nonzero if any
# failed. Checks assert the precondition they are actually testing, so a check
# cannot pass while silently exercising nothing.

set -uo pipefail

# non-interactive ssh may not set USER, and the default paths below are per-user
export USER="${USER:-$(id -un)}"

if [[ -r /etc/profile.d/00-modulepath.sh ]]; then
    source /etc/profile.d/00-modulepath.sh
    source /etc/profile.d/modules.sh
fi
module load singularity/3.11.3 2>/dev/null || true

if [[ -n "${FUSE_ENVIRONMENT:-}" ]]; then
    version="$FUSE_ENVIRONMENT"
else
    version="$(curl -s https://api.github.com/repos/ProjectTorreyPines/FUSE.jl/releases/latest | jq -r .name)"
fi
[[ -n "$version" && "$version" != "null" ]] || { echo "ERROR: set FUSE_ENVIRONMENT" >&2; exit 1; }

scriptdir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
scratch="/local-scratch/$USER"
shared="/fusion/projects/dt/fuse_containers"
sif="${SIF:-${SIF_DIR:-$scratch}/fuse_${version}.sif}"
image="localhost/fuse:$version"
launcher="$scriptdir/fuse-container"
out="$scratch/acceptance_${version}"
mkdir -p "$out"

podman=(podman)
[[ -d "$scratch" ]] && podman=(podman --root "$scratch/podman-storage" --runroot "$scratch/podman-run")

# The launcher needs squashfuse, which it looks for in bin/ next to the SIF. A
# freshly built SIF in /local-scratch has no bin/, so fall back to the copy
# published beside the shared images.
[[ -x "$(dirname "$sif")/bin/squashfuse" ]] || export PATH="$shared/bin:$PATH"

pass=0; fail=0; failed_names=()
check() {
    local name="$1"; shift
    echo "----- $name"
    if "$@"; then echo "PASS $name"; pass=$((pass+1))
    else echo "FAIL $name"; fail=$((fail+1)); failed_names+=("$name"); fi
}

echo "=== FUSE $version container acceptance — $(date -Is) on $(hostname)"
echo "=== image: $image"
echo "=== sif:   $sif"
echo

# --- identity ---------------------------------------------------------------
check "version-is-$version" bash -c "
    got=\$(${podman[*]} run --rm $image fuse -e 'import FUSE; print(pkgversion(FUSE))')
    echo \"  FUSE \$got\"; [ \"v\$got\" = \"$version\" ]"

check "sysimage-loaded" bash -c "
    p=\$(${podman[*]} run --rm $image fuse -e 'print(unsafe_string(Base.JLOptions().image_file))')
    echo \"  sysimage: \$p\"; [ \"\$p\" = /opt/fuse/sys_fuse.so ]"

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
# through its Jupyter kernel. v1.1.6 got it from a late Containerfile layer that
# shipped in that image but never merged to master, so v1.2.0 lost it — exactly
# the kind of released-vs-repo drift this script exists to catch.
check "revise-loads" bash -c "
    ${podman[*]} run --rm $image fuse -e 'using Revise; println(\"  Revise \", pkgversion(Revise))'"

check "host-tools-present" bash -c "
    ${podman[*]} run --rm $image bash -c 'for t in ps ssh rsync; do command -v \$t || { echo \"  missing \$t\"; exit 1; }; done'"

# --- physics smoke ----------------------------------------------------------
check "init-and-plot-png" bash -c "
    ${podman[*]} run --rm -v $out:/out:Z $image fuse -e '
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

# --- SIF + shared NN weight dirs --------------------------------------------
if [[ ! -f "$sif" ]]; then
    echo "----- SIF checks"
    echo "FAIL sif-exists ($sif not found)"; fail=$((fail+1)); failed_names+=("sif-exists")
else
    check "sif-smoke" bash -c "
        $launcher $sif -e '
            println(\"  sysimage: \", unsafe_string(Base.JLOptions().image_file))
            println(\"  cpu: \", Sys.CPU_NAME)
            using FUSE; println(\"  FUSE \", pkgversion(FUSE))'"

    # The shared pedestal weights are read ONLY on the ne_from=:nn_predictor
    # path. act.ActorPedestal.model=:EPED uses EPEDNN.jl's own baked weights and
    # never touches this directory -- driving the actor would pass while testing
    # nothing. Load the bundles directly, asserting the env var is honoured
    # first, which is exactly what the launcher's injection must make work.
    check "pedestal-nn-shared-dir" bash -c "
        FUSE_PEDESTAL_NN_DIR=$shared/pedestal_nn $launcher $sif -e '
            using FUSE
            dir = FUSE.resolve_pedestal_nn_dir()
            println(\"  resolved dir: \", dir)
            @assert dir == ENV[\"FUSE_PEDESTAL_NN_DIR\"] \"env var not honoured\"
            @assert FUSE.pedestal_nn_dir_complete(dir) \"shared pedestal-NN dir incomplete\"
            FUSE.load_pedestal_nn()
            println(\"  loaded PedestalNN bundles OK\")'"

    # SOLPS-NN runs through ActorSOL's model switch, not by calling ActorSOLPSNN
    # directly.
    check "solps-nn-shared-dir" bash -c "
        FUSE_SOLPS_NN_DIR=$shared/solps_nn_onnx $launcher $sif -e '
            using FUSE
            @assert ENV[\"FUSE_SOLPS_NN_DIR\"] == \"$shared/solps_nn_onnx\"
            ini, act = FUSE.case_parameters(:D3D, :L_mode)
            act.ActorSOL.model = :solps_nn
            dd = FUSE.init(ini, act)
            FUSE.ActorSOL(dd, act)
            println(\"  ActorSOL(:solps_nn) OK\")'"
fi

echo
echo "=== acceptance summary: $pass passed, $fail failed"
if [[ "$fail" -gt 0 ]]; then
    printf '    FAILED: %s\n' "${failed_names[@]}"
    exit 1
fi
echo "=== ALL CHECKS PASSED — safe to publish $version"
