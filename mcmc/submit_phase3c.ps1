# submit_phase3c.ps1 — launch Phase 3c crossing chain (beta in [2.65, 3.40])
# Run from CLASS_SYMT/ on Windows PowerShell. Run only AFTER 3b converges.
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

$chain3b_converged = Test-Path .\chains\phase3b_noncrossing\afterglow_3b.converged
if (-not $chain3b_converged) {
    Write-Warning "3b convergence sentinel not found. Per the master plan, 3c should run only after 3b converges."
}

$nproc = 8
Write-Host "Launching Phase 3c crossing (beta in [2.65, 3.40]) with $nproc MPI ranks..."
mpiexec -n $nproc cobaya-run mcmc\cobaya_phase3c_crossing.yaml
