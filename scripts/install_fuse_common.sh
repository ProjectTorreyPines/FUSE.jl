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
tmp="\$(mktemp)"
curl -fsSL "${installer}" -o "\${tmp}"
bash "\${tmp}" -b -p "${MINICONDA_DIR}"
rm -f "\${tmp}"
EOF
    chmod +x "${helper}"
    bash "${helper}"
    rm -f "${helper}"
    [[ -x "${MINICONDA_DIR}/bin/conda" ]] || die "Miniconda install failed (conda not found at ${MINICONDA_DIR}/bin/conda)"
}

activate_conda() {
  # shellcheck disable=SC1091
    if [[ -f "${MINICONDA_DIR}/etc/profile.d/conda.sh" ]]; then
        source "${MINICONDA_DIR}/etc/profile.d/conda.sh"
    elif command -v conda >/dev/null 2>&1; then
        eval "$("$(command -v conda)" shell.bash hook 2>/dev/null)" || true
    else
        die "conda is not available (install Miniconda or module load conda)"
    fi
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
    conda activate "${CONDA_ENV_NAME}"
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

run_fusebot_or_make() {
    local target="$1"
    shift || true
    if command -v fusebot >/dev/null 2>&1; then
        log "fusebot ${target} $*"
        fusebot "${target}" "$@"
        return 0
    fi

    local fuse_dir
    fuse_dir="$(fuse_pkg_dir)"
    log "fusebot not on PATH — running make ${target} in ${fuse_dir}"
    (
        cd "${fuse_dir}"
        export PTP_ORIGINAL_DIR="${INSTALL_DIR}"
        make "${target}" "$@"
    )
}

install_ijulia_kernels() {
    run_fusebot_or_make install_IJulia
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

    log "FUSE install complete."
    log "Step 2 — verify fluxmatcher.ipynb cells 0–2:"
    log "  bash ${SCRIPT_DIR}/verify_fluxmatcher_notebook.sh"
}

platform="${1:-}"
case "${platform}" in
    laptop)
        export FUSE_SETUP_SHELL="${FUSE_SETUP_SHELL:-false}"
        ensure_juliaup
        install_fuse_stack
        ;;
    nersc)
        export FUSE_SETUP_SHELL="${FUSE_SETUP_SHELL:-true}"
        load_nersc_modules
        install_fuse_stack
        ;;
    *)
        die "Usage: $(basename "$0") laptop|nersc"
        ;;
esac
