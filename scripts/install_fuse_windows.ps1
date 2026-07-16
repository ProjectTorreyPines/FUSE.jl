# One-shot FUSE install for Windows laptops (juliaup + Miniconda if needed).
#
# Copy-paste (PowerShell, from any working directory):
#   winget install julia -s msstore --accept-source-agreements --accept-package-agreements --disable-interactivity; `
#   irm https://raw.githubusercontent.com/ProjectTorreyPines/FUSE.jl/master/scripts/install_fuse_windows.ps1 | iex
#
# Or download and run from a FUSE.jl clone:
#   .\scripts\install_fuse_windows.ps1
#
# Bootstrap-only (syntax / Julia + conda setup, no FUSE packages):
#   .\scripts\install_fuse_windows.ps1 -BootstrapOnly

#Requires -Version 5.1
[CmdletBinding()]
param(
    [switch]$BootstrapOnly
)

$ErrorActionPreference = "Stop"

$ScriptBaseUrl = if ($env:FUSE_SCRIPT_BASE_URL) { $env:FUSE_SCRIPT_BASE_URL } else {
    "https://raw.githubusercontent.com/ProjectTorreyPines/FUSE.jl/master/scripts"
}

function Resolve-ScriptDir {
    param([string]$SelfPath = $PSCommandPath)
    if ($selfPath) {
        $selfDir = Split-Path -Parent $selfPath
        if (Test-Path (Join-Path $selfDir "install_fuse_julia.jl")) {
            return $selfDir
        }
    }

    $bundleDir = Join-Path $env:TEMP "fuse-install-$PID"
    New-Item -ItemType Directory -Force -Path $bundleDir | Out-Null
    foreach ($file in @("install_fuse_julia.jl")) {
        $dest = Join-Path $bundleDir $file
        Invoke-WebRequest -Uri "$ScriptBaseUrl/$file" -OutFile $dest
    }
    return $bundleDir
}

$ScriptDir = Resolve-ScriptDir
$RepoRoot = Split-Path -Parent $ScriptDir
$InstallDir = if ($env:FUSE_WORK_DIR) { $env:FUSE_WORK_DIR } else { (Get-Location).Path }
$CondaEnvName = if ($env:FUSE_CONDA_ENV) { $env:FUSE_CONDA_ENV } else { "fuse" }
$MinicondaDir = if ($env:MINICONDA_DIR) { $env:MINICONDA_DIR } else { Join-Path $env:USERPROFILE ".local\miniconda3" }

function Write-InstallLog {
    param([string]$Message)
    Write-Host "[install_fuse] $Message"
}

function Write-InstallError {
    param([string]$Message)
    Write-Error "[install_fuse] ERROR: $Message"
}

function Test-Command {
    param([string]$Name)
    return [bool](Get-Command $Name -ErrorAction SilentlyContinue)
}

function Add-ToPath {
    param([string]$Directory)
    if (-not (Test-Path $Directory)) { return }
    $parts = $env:Path -split ';' | Where-Object { $_ -and $_ -ne $Directory }
    $env:Path = ($Directory, ($parts -join ';')) -join ';'
}

function Ensure-Julia {
    if (Test-Command julia) {
        Write-InstallLog "Julia: $(julia --version)"
        return
    }

    $juliaupBin = Join-Path $env:USERPROFILE ".juliaup\bin"
    $juliaupExe = Join-Path $juliaupBin "julia.exe"
    if (Test-Path $juliaupExe) {
        Add-ToPath $juliaupBin
        Write-InstallLog "Julia: $(julia --version)"
        return
    }

    Write-InstallLog "julia not found — installing via winget or Julia App Installer"
    if (Test-Command winget) {
        & winget install julia -s msstore --accept-source-agreements --accept-package-agreements --disable-interactivity
    }
    else {
        Add-AppxPackage -AppInstallerFile "https://install.julialang.org/Julia.appinstaller"
    }

    if (Test-Path $juliaupExe) {
        Add-ToPath $juliaupBin
    }

    # Refresh PATH from the user environment (winget / Store installs update the registry).
    $userPath = [Environment]::GetEnvironmentVariable("Path", "User")
    if ($userPath) {
        $env:Path = "$userPath;$env:Path"
    }

    if (-not (Test-Command julia)) {
        Write-InstallError "julia is not on PATH after setup. Open a new terminal and re-run this script."
    }
    Write-InstallLog "Julia: $(julia --version)"
}

function Install-Miniconda {
    Write-InstallLog "Installing Miniconda to $MinicondaDir"
    $installer = Join-Path $env:TEMP "Miniconda3-latest-Windows-x86_64.exe"
    Invoke-WebRequest -Uri "https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe" -OutFile $installer
    $args = @(
        "/InstallationType=JustMe",
        "/RegisterPython=0",
        "/S",
        "/D=$MinicondaDir"
    )
    $proc = Start-Process -FilePath $installer -ArgumentList $args -Wait -PassThru -NoNewWindow
    Remove-Item $installer -Force -ErrorAction SilentlyContinue
    if ($proc.ExitCode -ne 0) {
        Write-InstallError "Miniconda installer exited with code $($proc.ExitCode)"
    }
    $condaExe = Join-Path $MinicondaDir "Scripts\conda.exe"
    if (-not (Test-Path $condaExe)) {
        Write-InstallError "Miniconda install failed (conda not found at $condaExe)"
    }
}

function Initialize-Conda {
    $condaHook = Join-Path $MinicondaDir "shell\condabin\conda-hook.ps1"
    if (Test-Path $condaHook) {
        . $condaHook
        return
    }
    if (Test-Command conda) {
        (& (Get-Command conda).Source "shell.powershell" "hook") | Out-String | Invoke-Expression
        return
    }
    Write-InstallError "conda is not available (install Miniconda or add conda to PATH)"
}

function Set-CondaChannels {
    conda config --set auto_activate_base false 2>$null
    conda config --add channels conda-forge 2>$null
    conda config --set channel_priority strict 2>$null
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main 2>$null
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r 2>$null
}

function Ensure-Conda {
    if (Test-Command conda) {
        Write-InstallLog "conda: $(conda --version)"
        Initialize-Conda
        return
    }

    $condaExe = Join-Path $MinicondaDir "Scripts\conda.exe"
    if (Test-Path $condaExe) {
        Write-InstallLog "Using existing Miniconda at $MinicondaDir"
        Add-ToPath (Join-Path $MinicondaDir "Scripts")
        Add-ToPath (Join-Path $MinicondaDir "condabin")
        Initialize-Conda
        return
    }

    Install-Miniconda
    Add-ToPath (Join-Path $MinicondaDir "Scripts")
    Add-ToPath (Join-Path $MinicondaDir "condabin")
    Initialize-Conda
}

function Get-FuseJupyterYml {
    $juliaScript = @'
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
'@
    return julia -e $juliaScript
}

function Ensure-FuseCondaEnv {
    Ensure-Conda
    Set-CondaChannels
    $yml = Get-FuseJupyterYml
    Write-InstallLog "Jupyter environment file: $yml"
    $envNames = conda env list | ForEach-Object {
        if ($_ -match '^\s*(\S+)') { $Matches[1] }
    }
    if ($envNames -contains $CondaEnvName) {
        Write-InstallLog "Conda env '$CondaEnvName' exists — updating"
        conda env update -n $CondaEnvName -f $yml --prune
    }
    else {
        Write-InstallLog "Creating conda env '$CondaEnvName'"
        conda env create -f $yml
    }
    conda activate $CondaEnvName
    Write-InstallLog "Active Python: $(python -c 'import sys; print(sys.executable)')"
}

function Get-FusePkgDir {
    return julia -e 'using FUSE; print(pkgdir(FUSE))'
}

function Get-FusebotDir {
    $candidates = @()

    if (Test-Command juliaup) {
        $candidates += Split-Path (Get-Command juliaup).Source -Parent
    }
    if (Test-Command julia) {
        $juliaDir = Split-Path (Get-Command julia).Source -Parent
        if ($juliaDir -match 'juliaup' -or (Test-Path $juliaDir)) {
            $candidates += $juliaDir
        }
    }
    $candidates += @(
        (Join-Path $env:USERPROFILE ".local\bin"),
        (Join-Path $env:LOCALAPPDATA "Programs\Julia\bin")
    )

    foreach ($dir in $candidates) {
        if (-not $dir) { continue }
        try {
            New-Item -ItemType Directory -Force -Path $dir | Out-Null
            $testFile = Join-Path $dir ".write_test"
            New-Item -ItemType File -Force -Path $testFile | Out-Null
            Remove-Item $testFile -Force
            return $dir
        }
        catch {
            continue
        }
    }
    Write-InstallError "Could not find a writable directory for fusebot"
}

function Install-FusebotCli {
    $fusebotDir = Get-FusebotDir
    Write-InstallLog "fusebot install directory: $fusebotDir"
    $env:FUSE_INSTALL_DIR = $fusebotDir
    Add-ToPath $fusebotDir

    $juliaScript = Join-Path $ScriptDir "install_fuse_julia.jl"
    & julia $juliaScript fusebot
    if ($LASTEXITCODE -ne 0) {
        Write-InstallLog "fusebot install via Julia failed — copying fusebot manually"
        $fuseDir = Get-FusePkgDir
        Copy-Item (Join-Path $fuseDir "fusebot") (Join-Path $fusebotDir "fusebot") -Force
        if ($env:FUSE_SETUP_SHELL -eq "true") {
            julia -e "using FUSE; FUSE.setup_fusebot_shell!(\"$fusebotDir\")" 2>$null
        }
    }
}

function Invoke-FusebotOrMake {
    param(
        [Parameter(Mandatory = $true)][string]$Target,
        [Parameter(ValueFromRemainingArguments = $true)][string[]]$ExtraArgs
    )

    if (Test-Command fusebot) {
        Write-InstallLog "fusebot $Target $($ExtraArgs -join ' ')"
        & fusebot $Target @ExtraArgs
        return
    }

    $fuseDir = Get-FusePkgDir
    Write-InstallLog "fusebot not on PATH — running make $Target in $fuseDir"
    Push-Location $fuseDir
    try {
        $env:PTP_ORIGINAL_DIR = $InstallDir
        & make $Target @ExtraArgs
    }
    finally {
        Pop-Location
    }
}

function Install-IJuliaKernels {
    Invoke-FusebotOrMake install_IJulia
}

function Clone-FuseExamples {
    Set-Location $InstallDir
    $examplesDir = Join-Path $InstallDir "FuseExamples"
    if (Test-Path (Join-Path $examplesDir ".git")) {
        Write-InstallLog "Updating FuseExamples"
        git -C $examplesDir fetch origin
        git -C $examplesDir reset --hard origin/master
    }
    else {
        Write-InstallLog "Cloning FuseExamples into $InstallDir"
        git clone https://github.com/ProjectTorreyPines/FuseExamples.git
    }
}

function Invoke-JuliaInstall {
    param([string]$Step = "all")
    $juliaScript = Join-Path $ScriptDir "install_fuse_julia.jl"
    & julia $juliaScript $Step
    if ($LASTEXITCODE -ne 0) {
        Write-InstallError "Julia install step '$Step' failed"
    }
}

function Install-FuseStack {
    if (-not (Test-Command git)) {
        Write-InstallError "git is required but not found (install Git for Windows)"
    }

    Set-Location $InstallDir
    Ensure-Julia
    Write-InstallLog "Working directory: $InstallDir"

    Invoke-JuliaInstall packages
    Invoke-JuliaInstall revise
    Install-FusebotCli

    Ensure-FuseCondaEnv
    Install-IJuliaKernels
    Clone-FuseExamples

    if ($env:FUSE_SKIP_SMOKE -ne "1") {
        Invoke-JuliaInstall smoke
    }
    else {
        Write-InstallLog "Skipping smoke test (FUSE_SKIP_SMOKE=1)"
    }

    Write-InstallLog "FUSE install complete."
    Write-InstallLog "Step 2 — verify fluxmatcher.ipynb cells 0–2:"
    Write-InstallLog "  .\scripts\verify_fluxmatcher_notebook.ps1"
}

$env:FUSE_SETUP_SHELL = if ($env:FUSE_SETUP_SHELL) { $env:FUSE_SETUP_SHELL } else { "false" }

Ensure-Julia
Ensure-Conda

if ($BootstrapOnly) {
    Write-InstallLog "Bootstrap complete (Julia + conda available)."
    exit 0
}

Install-FuseStack
