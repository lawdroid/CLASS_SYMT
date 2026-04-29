# submit_phase3a.ps1 — launch Phase 3a baseline chain (beta = 0)
# Run from CLASS_SYMT/ on Windows PowerShell.
$ErrorActionPreference = 'Stop'

Set-Location $PSScriptRoot\..

# --- pre-flight: CLASS binary has afterglow patches? ---
if (-not (Test-Path .\class)) {
    throw "CLASS binary not found. Run 'make class' first."
}
$afterglow_in_binary = & strings.exe .\class 2>$null | Select-String "afterglow_on"
if (-not $afterglow_in_binary) {
    Write-Warning "Could not find 'afterglow_on' string in binary; verify with smoke test before launching long chain."
}

# --- pre-flight: data directories present? ---
$datasets = @(".\data\Planck_NPIPE", ".\data\DESI_DR2", ".\data\Pantheon+")
foreach ($d in $datasets) {
    if (-not (Test-Path $d)) {
        throw "Missing dataset folder: $d  (provision before launch)"
    }
}

# --- pre-flight: provenance tag exists? ---
$tag = & git tag --list "phase3-prereg-v1"
if (-not $tag) {
    Write-Warning "Provenance tag 'phase3-prereg-v1' not found. Per Phase3_MCMC_Windows_Execution.md §7 the tag should exist BEFORE launch."
}

# --- adjust nproc to your machine ---
$nproc = 8

Write-Host "Launching Phase 3a baseline (beta=0) with $nproc MPI ranks..."
mpiexec -n $nproc cobaya-run mcmc\cobaya_phase3a_baseline.yaml
