# Julia Installation Guide for KAIROS Administrators

## Quick Setup

For KAIROS administrators, install Julia system-wide with:

```bash
export JULIA_VERSION=1.11.6
export JULIA_INSTALL_PREFIX=/opt/nfri/glib/julia
export JULIA_MODULE_DIR=/opt/nfri/glib/modulefiles

# Run with sudo if not logged in as root
sudo ./install_julia_binary.sh
```

## Installation Details

The script will:
1. Download Julia ${JULIA_VERSION} tarball from official repository
2. Extract to `/opt/nfri/glib/julia/julia-1.11.6/`
3. Create module file at `/opt/nfri/glib/modulefiles/julia/1.11.6`

## Usage

Users can load Julia with:
```bash
module use /opt/nfri/glib/modulefiles
module load julia/1.11.6
```

## Requirements
- curl, tar, awk commands available
- Root access or sudo privileges (for writing to `/opt/nfri/glib/`)
- Network access through KAIROS proxy (automatically configured)

## Notes
- Installation is skipped if Julia is already present at the target location
- No juliaup launcher is used - direct binary installation
- Proxy settings are automatically configured for KAIROS environment