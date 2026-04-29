# submit_phase3b.ps1 — launch Phase 3b non-crossing chain (beta in (0, 2.4])
# Run from CLASS_SYMT/ on Windows PowerShell. Run only AFTER 3a converges.
$ErrorActionPreference = 'Stop'

Set-Location $PSScriptRoot\..

if (-not (Test-Path .\class)) {
    throw "CLASS binary not found. Run 'make class' first."
}

$datasets = @(".\data\Planck_NPIPE", ".\data\DESI_DR2", ".\data\Pantheon+")
foreach ($d in $datasets) {
    if (-not (Test-Path $d)) {
        throw "Missing dataset folder: $d"
    }
}

# Verify 3a chain has converged before starting 3b
$chain3a_converged = Test-Path .\chains\phase3a_baseline\afterglow_3a.converged
if (-not $chain3a_converged) {
    Write-Warning "3a convergence sentinel not found at chains\phase3a_baseline\afterglow_3a.converged. Per the master plan, 3b should run only after 3a converges."
}

$nproc = 8
Write-Host "Launching Phase 3b non-crossing (beta in (0, 2.4]) with $nproc MPI ranks..."
mpiexec -n $nproc cobaya-run mcmc\cobaya_phase3b_noncrossing.yaml
