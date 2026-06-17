#!/bin/bash
set -euo pipefail

# Kairos-specific proxy configuration
export http_proxy="http://203.230.124.39:4128"
export https_proxy="http://203.230.124.39:4128"
export no_proxy="localhost,127.0.0.1,*.hpc.nfri.re.k"

# ===== Default values (can be overridden via environment variables) ===========
# To override, set environment variables before running this script:
#   export JULIA_VERSION=1.11.6
#   export JULIA_INSTALL_PREFIX=/custom/path
#   export JULIA_MODULE_DIR=/custom/modules
export JULIA_VERSION=${JULIA_VERSION:-1.11.6}
export JULIA_INSTALL_PREFIX=${JULIA_INSTALL_PREFIX:-/opt/nfri/glib/julia}
export JULIA_MODULE_DIR=${JULIA_MODULE_DIR:-/opt/nfri/glib/modulefiles}

need_var() { [[ -n "${!1:-}" ]] || { echo "ERROR: $1 must be set"; exit 1; }; }
need_bin() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found"; exit 1; }; }

need_var JULIA_VERSION
need_var JULIA_INSTALL_PREFIX
need_var JULIA_MODULE_DIR

need_bin curl
need_bin tar
need_bin awk

# ===== URL and paths (Generic Linux x86_64, glibc) ============================
MM="$(awk -F. '{print $1"."$2}' <<<"$JULIA_VERSION")"
TARBALL="julia-${JULIA_VERSION}-linux-x86_64.tar.gz"
URL="https://julialang-s3.julialang.org/bin/linux/x64/${MM}/${TARBALL}"

DEST_DIR="${JULIA_INSTALL_PREFIX}/julia-${JULIA_VERSION}"
MODULE_PATH="${JULIA_MODULE_DIR}/julia/${JULIA_VERSION}"

# Skip if already installed
if [[ -e "$DEST_DIR" ]]; then
  echo "[install-julia] Already present: $DEST_DIR"
  exit 0
fi

# ===== Download and extract ====================================================
TMPDIR="$(mktemp -d)"; trap 'rm -rf "$TMPDIR"' EXIT
echo "[install-julia] Downloading ${URL}"
curl -fL --retry 3 -o "${TMPDIR}/${TARBALL}" "${URL}"

mkdir -p "${JULIA_INSTALL_PREFIX}"
tar -xzf "${TMPDIR}/${TARBALL}" -C "${JULIA_INSTALL_PREFIX}"

# Normalize extracted directory name if it differs
if [[ ! -d "$DEST_DIR" ]]; then
  EXPANDED="$(find "${JULIA_INSTALL_PREFIX}" -maxdepth 1 -type d -name "julia-${JULIA_VERSION}*" | head -n1)"
  [[ -n "$EXPANDED" ]] || { echo "ERROR: extracted folder not found"; exit 1; }
  mv "${EXPANDED}" "${DEST_DIR}"
fi

# ===== Create a minimal Tcl modulefile ========================================
mkdir -p "$(dirname "${MODULE_PATH}")"
cat > "${MODULE_PATH}" <<EOF
#%Module1.0
proc ModulesHelp { } {
    puts stderr "Julia ${JULIA_VERSION} (tarball install)"
}
module-whatis "Julia programming language"

# Expose the real Julia binary (no juliaup launcher)
prepend-path PATH ${DEST_DIR}/bin

# Keep user's own depot first by default; do not hardcode foreign HOME.
# Optional site depot example (commented):
# setenv JULIA_DEPOT_PATH "\$::env(HOME)/.julia:/opt/site/julia_depot:"
EOF
chmod 0644 "${MODULE_PATH}"

echo "[install-julia] Installed: ${DEST_DIR}"
echo "[install-julia] Modulefile: ${MODULE_PATH}"
echo
echo "Usage:"
echo "  module use $(dirname "$(dirname "$MODULE_PATH")")"
echo "  module load julia/${JULIA_VERSION}"
echo "  which julia   # => ${DEST_DIR}/bin/julia"



unset http_proxy
unset https_proxy
unset no_proxy