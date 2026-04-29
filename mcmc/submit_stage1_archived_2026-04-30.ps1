# ==============================================================
# Stage-1 background MCMC launcher (Windows PowerShell)
# CLASS_SYMT afterglow — Martin & Koh (April 2026)
#
# Prereq: Patches 1-12 applied, CLASS built, Cobaya installed.
# Usage:  .\submit_stage1.ps1   [from mcmc/ directory]
# ==============================================================

$ErrorActionPreference = "Stop"

# --- sanity checks ---
if (-not (Test-Path "..\class")) {
    Write-Host "[ERROR] CLASS binary not found at ..\class." -ForegroundColor Red
    Write-Host "        Apply Patches 7-12, then run:  cd ..; make class" -ForegroundColor Yellow
    exit 1
}
if (-not (Get-Command "cobaya-run" -ErrorAction SilentlyContinue)) {
    Write-Host "[ERROR] Cobaya not installed." -ForegroundColor Red
    Write-Host "        Install:  pip install cobaya getdist" -ForegroundColor Yellow
    exit 1
}

# --- stage-1 config ---
$CHAINS     = 8
$YAML       = ".\cobaya_stage1.yaml"
$OUT        = ".\chains\stage1"
New-Item -ItemType Directory -Force -Path $OUT | Out-Null

Write-Host ""
Write-Host "================================================================" -ForegroundColor Cyan
Write-Host "  CLASS_SYMT Phase 3 — Stage 1 (background-only) MCMC" -ForegroundColor Cyan
Write-Host "  Datasets: Pantheon+SH0ES · DESI DR2 BAO · Planck 2018 cmp" -ForegroundColor Cyan
Write-Host "  Sampler:  Metropolis-Hastings × $CHAINS chains" -ForegroundColor Cyan
Write-Host "  Target:   R-1 < 0.02 (Gelman-Rubin)" -ForegroundColor Cyan
Write-Host "================================================================" -ForegroundColor Cyan
Write-Host ""

# --- launch MPI chains ---
# Windows MPI: use MS-MPI or WSL. Simplest portable path = separate python procs.
$jobs = @()
for ($i = 1; $i -le $CHAINS; $i++) {
    $logfile = "$OUT\chain_$i.log"
    $job = Start-Process -FilePath "cobaya-run" `
                         -ArgumentList $YAML, "--output-dir", "$OUT\c$i" `
                         -NoNewWindow -PassThru `
                         -RedirectStandardOutput $logfile
    $jobs += $job
    Write-Host "  Launched chain $i (PID $($job.Id)) -> $logfile"
    Start-Sleep -Seconds 2       # stagger to avoid IO collision at init
}

Write-Host ""
Write-Host "All $CHAINS chains running. Monitor with:" -ForegroundColor Green
Write-Host "  Get-Content $OUT\chain_1.log -Wait -Tail 20"
Write-Host ""
Write-Host "To check convergence (after ~1 hour):" -ForegroundColor Green
Write-Host "  python -c `"from getdist.mcsamples import loadMCSamples; s = loadMCSamples('$OUT\c1\afterglow_bg'); print(s.getConvergeTests())`""

# Wait for all
$jobs | ForEach-Object { $_ | Wait-Process }
Write-Host ""
Write-Host "Stage 1 complete. Analyze with:  python analyze_stage1.py" -ForegroundColor Green
