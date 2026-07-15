# Step 2 after install: verify fluxmatcher.ipynb cells 0–2 from PowerShell.
#
#   .\scripts\verify_fluxmatcher_notebook.ps1

#Requires -Version 5.1
$ErrorActionPreference = "Stop"

$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$ScriptBaseUrl = if ($env:FUSE_SCRIPT_BASE_URL) { $env:FUSE_SCRIPT_BASE_URL } else {
    "https://raw.githubusercontent.com/ProjectTorreyPines/FUSE.jl/master/scripts"
}

if (-not (Test-Path (Join-Path $ScriptDir "verify_fluxmatcher_notebook.jl"))) {
    $bundleDir = Join-Path $env:TEMP "fuse-verify-$PID"
    New-Item -ItemType Directory -Force -Path $bundleDir | Out-Null
    Invoke-WebRequest -Uri "$ScriptBaseUrl/verify_fluxmatcher_notebook.jl" -OutFile (Join-Path $bundleDir "verify_fluxmatcher_notebook.jl")
    $ScriptDir = $bundleDir
}
$WorkDir = if ($env:FUSE_WORK_DIR) { $env:FUSE_WORK_DIR } else { (Get-Location).Path }
$Notebook = Join-Path $WorkDir "FuseExamples\fluxmatcher.ipynb"

function Write-VerifyLog {
    param([string]$Message)
    Write-Host "[verify_fluxmatcher] $Message"
}

$juliaupBin = Join-Path $env:USERPROFILE ".juliaup\bin"
if (Test-Path (Join-Path $juliaupBin "julia.exe")) {
    $env:Path = "$juliaupBin;$env:Path"
}

if (-not (Get-Command julia -ErrorAction SilentlyContinue)) {
    Write-Error "julia not on PATH"
}

if (-not (Test-Path $Notebook)) {
    Write-Error "$Notebook not found. Run the install script first, or set FUSE_WORK_DIR."
}

Set-Location $WorkDir

$pythonCandidates = @(
    (Join-Path $env:USERPROFILE ".local\miniconda3\envs\fuse\python.exe"),
    (Get-Command python -ErrorAction SilentlyContinue | Select-Object -ExpandProperty Source)
)
$Python = $pythonCandidates | Where-Object { $_ -and (Test-Path $_) } | Select-Object -First 1
if (-not $Python) {
    Write-Error "python not on PATH (activate the fuse conda env first)"
}

# Extract the code from the first three cells (cells 0 & 2; cell 1 is markdown)
# into a temp file that verify_fluxmatcher_notebook.jl runs, so the test stays
# faithful to whatever the notebook actually contains.
$env:FUSE_VERIFY_CELLS = [System.IO.Path]::GetTempFileName()

$pySource = @"
import json
import os
import sys
from pathlib import Path

out_path = Path(os.environ["FUSE_VERIFY_CELLS"])
path = Path(r"$Notebook")
nb = json.loads(path.read_text(encoding="utf-8"))
cells = nb["cells"]
if len(cells) < 3:
    sys.exit("fluxmatcher.ipynb must have at least 3 cells")
if cells[0]["cell_type"] != "code":
    sys.exit("Cell 0 must be a code cell")
if cells[1]["cell_type"] != "markdown":
    sys.exit("Cell 1 must be a markdown cell")
if cells[2]["cell_type"] != "code":
    sys.exit("Cell 2 must be a code cell")
text = "".join(cells[1].get("source", [])).lower()
if "flux-matcher" not in text and "flux matcher" not in text:
    sys.exit("Cell 1 markdown does not mention flux-matcher")

chunks = []
for i, cell in enumerate(cells[:3]):
    if cell["cell_type"] != "code":
        continue
    chunks.append("# --- fluxmatcher.ipynb cell %d ---" % i)
    chunks.append("".join(cell.get("source", [])))
out_path.write_text("\n".join(chunks) + "\n", encoding="utf-8")
print("[verify_fluxmatcher] Cells 0-2 present (cell 1 markdown, cells 0 & 2 code)")
"@

$pySource | & $Python -
if ($LASTEXITCODE -ne 0) {
    Remove-Item -Force $env:FUSE_VERIFY_CELLS -ErrorAction SilentlyContinue
    Write-Error "Notebook structure check failed"
}

$juliaVerify = Join-Path $ScriptDir "verify_fluxmatcher_notebook.jl"
& julia $juliaVerify
$juliaExit = $LASTEXITCODE
Remove-Item -Force $env:FUSE_VERIFY_CELLS -ErrorAction SilentlyContinue
if ($juliaExit -ne 0) {
    Write-Error "Julia verification failed"
}

Write-VerifyLog "fluxmatcher.ipynb cells 0–2 verification PASSED"
