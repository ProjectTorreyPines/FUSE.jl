#!/usr/bin/env bash
# Shared FUSE install logic for laptop and NERSC (executed, not sourced).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
INSTALL_DIR="${FUSE_WORK_DIR:-${PWD}}"
CONDA_ENV_NAME="${FUSE_CONDA_ENV:-fuse}"
MINICONDA_DIR="${MINICONDA_DIR:-${HOME}/.local/miniconda3}"
JULIA_MODULE="${FUSE_JULIA_MODULE:-julia/1.11.7}"

log() { echo "[install_fuse] $*"; }
die() { echo "[install_fuse] ERROR: $*" >&2; exit 1; }

need_cmd() {
    command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

ensure_julia() {
    if command -v julia >/dev/null 2>&1; then
        log "Julia: $(julia --version)"
        return 0
    fi
    die "julia is not on PATH after setup"
}

ensure_juliaup() {
    if command -v julia >/dev/null 2>&1; then
        return 0
    fi
    if [[ -x "${HOME}/.juliaup/bin/julia" ]]; then
        export PATH="${HOME}/.juliaup/bin:${PATH}"
        return 0
    fi
    log "julia not found — installing juliaup"
    need_cmd curl
    curl -fsSL https://install.julialang.org | sh -s -- -y
    export PATH="${HOME}/.juliaup/bin:${PATH}"
}

load_nersc_modules() {
    command -v module >/dev/null 2>&1 || die "Lmod is required on NERSC (module command not found)"
    log "Loading ${JULIA_MODULE}"
    module load "${JULIA_MODULE}"
    log "Loading conda"
    module load conda
    log "Julia: $(julia --version)"
}

miniconda_installer_url() {
    local os arch
    os="$(uname -s)"
    arch="$(uname -m)"
    case "${os}" in
        Linux)
            case "${arch}" in
                x86_64)  echo "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh" ;;
                aarch64|arm64) echo "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh" ;;
                *) die "Unsupported Linux architecture: ${arch}" ;;
            esac
            ;;
        Darwin)
            case "${arch}" in
                arm64)   echo "https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh" ;;
                x86_64)  echo "https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh" ;;
                *) die "Unsupported macOS architecture: ${arch}" ;;
            esac
            ;;
        *) die "Miniconda auto-install is only supported on Linux and macOS" ;;
    esac
}

install_miniconda() {
    local installer helper
    log "Installing Miniconda to ${MINICONDA_DIR}"
    need_cmd curl
    mkdir -p "$(dirname "${MINICONDA_DIR}")"
    installer="$(miniconda_installer_url)"
    helper="$(mktemp)"
    cat > "${helper}" <<EOF
#!/usr/bin/env bash
set -euo pipefail
# Constructor installers require the filename to end in ".sh" (they check \$0).
tmp="\$(mktemp "\${TMPDIR:-/tmp}/miniconda-installer-XXXXXX.sh")"
curl -fsSL "${installer}" -o "\${tmp}"
bash "\${tmp}" -b -p "${MINICONDA_DIR}"
rm -f "\${tmp}"
EOF
    chmod +x "${helper}"
    bash "${helper}"
    rm -f "${helper}"
    [[ -x "${MINICONDA_DIR}/bin/conda" ]] || die "Miniconda install failed (conda not found at ${MINICONDA_DIR}/bin/conda)"
}

# conda's shell functions (activation, `shell.bash hook`, sourced conda.sh)
# reference $PS1, which is unbound under `set -u` in a non-interactive shell.
# Run them with nounset temporarily disabled so the install does not abort with
# "PS1: unbound variable" (see conda/conda#3200).
conda_activate() {
    set +u
    conda activate "$@"
    set -u
}

activate_conda() {
  # shellcheck disable=SC1091
    set +u
    if [[ -f "${MINICONDA_DIR}/etc/profile.d/conda.sh" ]]; then
        source "${MINICONDA_DIR}/etc/profile.d/conda.sh"
    elif command -v conda >/dev/null 2>&1; then
        eval "$("$(command -v conda)" shell.bash hook 2>/dev/null)" || true
    else
        set -u
        die "conda is not available (install Miniconda or module load conda)"
    fi
    set -u
}

bootstrap_conda_channels() {
    conda config --set auto_activate_base false 2>/dev/null || true
    conda config --add channels conda-forge 2>/dev/null || true
    conda config --set channel_priority strict 2>/dev/null || true
    # Non-interactive installs (Miniconda 25+) require explicit TOS acceptance.
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main 2>/dev/null || true
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r 2>/dev/null || true
}

ensure_conda() {
    if command -v conda >/dev/null 2>&1; then
        log "conda: $(conda --version)"
        activate_conda
        return 0
    fi

    if [[ -x "${MINICONDA_DIR}/bin/conda" ]]; then
        log "Using existing Miniconda at ${MINICONDA_DIR}"
        activate_conda
        return 0
    fi

    if command -v module >/dev/null 2>&1; then
        log "Trying module load conda"
        module load conda 2>/dev/null || true
        if command -v conda >/dev/null 2>&1; then
            activate_conda
            return 0
        fi
    fi

    install_miniconda
    activate_conda
}

fuse_jupyter_yml() {
    julia -e '
        try
            using FUSE
        catch
            using Pkg
            Pkg.Registry.add(RegistrySpec(url="https://github.com/ProjectTorreyPines/FuseRegistry.jl.git"))
            Pkg.Registry.add("General")
            Pkg.add("FUSE")
            using FUSE
        end
        print(pkgdir(FUSE, "docs", "jupyter_environment.yml"))
    '
}

ensure_fuse_conda_env() {
    ensure_conda
    bootstrap_conda_channels
    local yml
    yml="$(fuse_jupyter_yml)"
    log "Jupyter environment file: ${yml}"
    if conda env list | awk '{print $1}' | grep -qx "${CONDA_ENV_NAME}"; then
        log "Conda env '${CONDA_ENV_NAME}' exists — updating"
        conda env update -n "${CONDA_ENV_NAME}" -f "${yml}" --prune
    else
        log "Creating conda env '${CONDA_ENV_NAME}'"
        conda env create -f "${yml}"
    fi
    conda_activate "${CONDA_ENV_NAME}"
    log "Active Python: $(python -c 'import sys; print(sys.executable)')"
}

fuse_pkg_dir() {
    julia -e 'using FUSE; print(pkgdir(FUSE))'
}

pick_fusebot_dir() {
    local candidates=()
    if command -v juliaup >/dev/null 2>&1; then
        candidates+=("$(dirname "$(command -v juliaup)")")
    fi
    if command -v julia >/dev/null 2>&1; then
        local julia_dir
        julia_dir="$(dirname "$(command -v julia)")"
        if [[ "${julia_dir}" == *juliaup* ]] || [[ -w "${julia_dir}" ]]; then
            candidates+=("${julia_dir}")
        fi
    fi
    candidates+=("${HOME}/.local/bin" "${HOME}/.local/shared/bin")

    local dir
    for dir in "${candidates[@]}"; do
        if mkdir -p "${dir}" 2>/dev/null && [[ -w "${dir}" ]]; then
            echo "${dir}"
            return 0
        fi
    done
    die "Could not find a writable directory for fusebot"
}

install_fusebot_cli() {
    local fusebot_dir
    fusebot_dir="$(pick_fusebot_dir)"
    log "fusebot install directory: ${fusebot_dir}"
    export FUSE_INSTALL_DIR="${fusebot_dir}"
    export PATH="${fusebot_dir}:${PATH}"
    if ! julia "${SCRIPT_DIR}/install_fuse_julia.jl" fusebot; then
        log "fusebot install via Julia failed — copying fusebot manually"
        local fuse_dir
        fuse_dir="$(fuse_pkg_dir)"
        cp "${fuse_dir}/fusebot" "${fusebot_dir}/fusebot"
        chmod +x "${fusebot_dir}/fusebot"
        if [[ "${FUSE_SETUP_SHELL:-false}" == "true" ]]; then
            julia -e "using FUSE; FUSE.setup_fusebot_shell!(\"${fusebot_dir}\")" || true
        fi
    fi
}

resolve_fuse_script() {
    local name="$1"
    if [[ -f "${SCRIPT_DIR}/${name}" ]]; then
        echo "${SCRIPT_DIR}/${name}"
        return 0
    fi
    local from_pkg=""
    if from_pkg="$(julia -e "using FUSE; print(joinpath(pkgdir(FUSE), \"scripts\", \"${name}\"))" 2>/dev/null)" \
        && [[ -n "${from_pkg}" && -f "${from_pkg}" ]]; then
        echo "${from_pkg}"
        return 0
    fi
    # Last resort for curl-bootstrap / older registry packages.
    local base_url="${FUSE_SCRIPT_BASE_URL:-https://raw.githubusercontent.com/ProjectTorreyPines/FUSE.jl/master/scripts}"
    local bundle_dir="${TMPDIR:-/tmp}/fuse-install-scripts-$$"
    mkdir -p "${bundle_dir}"
    curl -fsSL "${base_url}/${name}" -o "${bundle_dir}/${name}"
    if [[ "${name}" == *.sh ]]; then
        local companion="${name%.sh}.jl"
        curl -fsSL "${base_url}/${companion}" -o "${bundle_dir}/${companion}" 2>/dev/null || true
    fi
    [[ -f "${bundle_dir}/${name}" ]] || die "Could not locate scripts/${name}"
    echo "${bundle_dir}/${name}"
}

run_fusebot_or_make() {
    local target="$1"
    shift || true
    if command -v fusebot >/dev/null 2>&1; then
        log "fusebot ${target} $*"
        if fusebot "${target}" "$@"; then
            return 0
        fi
        log "fusebot ${target} failed — falling back to make"
    else
        log "fusebot not on PATH — falling back to make"
    fi

    local fuse_dir
    fuse_dir="$(fuse_pkg_dir)"
    log "Running make ${target} in ${fuse_dir}"
    (
        cd "${fuse_dir}"
        export PTP_ORIGINAL_DIR="${INSTALL_DIR}"
        make "${target}" "$@"
    )
}

install_ijulia_kernels() {
    # fusebot → make → direct install_ijulia.sh (no make required).
    # Keeps going when fusebot is missing/broken or make is not installed.
    if command -v fusebot >/dev/null 2>&1; then
        log "fusebot install_IJulia"
        if fusebot install_IJulia; then
            return 0
        fi
        log "fusebot install_IJulia failed — trying make / direct install"
    else
        log "fusebot not on PATH — trying make / direct install"
    fi

    local fuse_dir
    fuse_dir="$(fuse_pkg_dir)"

    if command -v make >/dev/null 2>&1; then
        log "make install_IJulia in ${fuse_dir}"
        if (
            cd "${fuse_dir}"
            export PTP_ORIGINAL_DIR="${INSTALL_DIR}"
            make install_IJulia
        ); then
            return 0
        fi
        log "make install_IJulia failed — running scripts/install_ijulia.sh directly"
    else
        log "make not found — running scripts/install_ijulia.sh directly"
    fi

    bash "${fuse_dir}/scripts/install_ijulia.sh"
}

clone_fuse_examples() {
    cd "${INSTALL_DIR}"
    if [[ -d FuseExamples/.git ]]; then
        log "Updating FuseExamples"
        git -C FuseExamples fetch origin
        git -C FuseExamples reset --hard origin/master
    else
        log "Cloning FuseExamples into ${INSTALL_DIR}"
        git clone https://github.com/ProjectTorreyPines/FuseExamples.git
    fi
}

verify_fluxmatcher_notebook() {
    local verify_sh
    verify_sh="$(resolve_fuse_script verify_fluxmatcher_notebook.sh)"
    log "Verifying fluxmatcher.ipynb cells 0–2 (often ~6 minutes on one thread the first time)"
    # Cell extraction needs Python from the fuse env.
    if ! command -v python >/dev/null 2>&1; then
        ensure_fuse_conda_env
    fi
    FUSE_WORK_DIR="${INSTALL_DIR}" bash "${verify_sh}"
}

run_julia_install() {
    local step="${1:-all}"
    julia "${SCRIPT_DIR}/install_fuse_julia.jl" "${step}"
}

install_fuse_stack() {
    cd "${INSTALL_DIR}"
    need_cmd git
    need_cmd curl
    ensure_julia

    log "Working directory: ${INSTALL_DIR}"
    run_julia_install packages
    run_julia_install revise
    install_fusebot_cli

    ensure_fuse_conda_env
    install_ijulia_kernels
    clone_fuse_examples
    run_julia_install smoke

    if [[ "${FUSE_SKIP_VERIFY:-0}" == "1" ]]; then
        export FUSE_VERIFY_FLUXMATCHER=false
    fi

    if [[ "${FUSE_VERIFY_FLUXMATCHER:-false}" == "true" ]]; then
        verify_fluxmatcher_notebook
        log "FUSE install complete (including fluxmatcher.ipynb cells 0–2)."
    else
        log "FUSE install complete."
        log "Step 2 — verify fluxmatcher.ipynb cells 0–2:"
        log "  bash $(resolve_fuse_script verify_fluxmatcher_notebook.sh)"
    fi
}

platform="${1:-}"
case "${platform}" in
    laptop)
        export FUSE_SETUP_SHELL="${FUSE_SETUP_SHELL:-false}"
        # Laptop one-command install finishes by running fluxmatcher cells 0–2.
        export FUSE_VERIFY_FLUXMATCHER="${FUSE_VERIFY_FLUXMATCHER:-true}"
        ensure_juliaup
        install_fuse_stack
        ;;
    nersc)
        export FUSE_SETUP_SHELL="${FUSE_SETUP_SHELL:-true}"
        # Same finish as laptop: fluxmatcher cells 0–2 (~6 minutes on 1 thread is fine on a login node).
        export FUSE_VERIFY_FLUXMATCHER="${FUSE_VERIFY_FLUXMATCHER:-true}"
        load_nersc_modules
        install_fuse_stack
        ;;
    *)
        die "Usage: $(basename "$0") laptop|nersc"
        ;;
esac
