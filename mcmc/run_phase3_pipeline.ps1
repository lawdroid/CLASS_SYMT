# run_phase3_pipeline.ps1
# One-shot: fetch 3a/3b/3c chains to canonical paths and regenerate
# chains_summary.json with real data. Run from CLASS_SYMT/ root in PowerShell.
#
# Steps:
#   1. scp 3a + 3b from Ubuntu workstation -> chains/phase3a_baseline/, chains/phase3b_noncrossing/
#   2. scp 3c from Mac mini                 -> chains/phase3c_b1narrow/
#   3. python mcmc/export_chains_summary.py -> chains_summary.json with REAL data + pending_numbers
#
# Per CLAUDE.md §1, run from Windows PowerShell (not WSL) — WSL→LAN NAT
# blocks ssh to 192.168.200.*; Windows-side ssh works.
#
# Prerequisites:
#   - passwordless ssh to i@192.168.200.119 and ingko@192.168.200.173
#   - python on PATH with: numpy, pandas (only used if MCEvidence is installed)
#   - optional: pip install MCEvidence  (for logZ_3a, logZ_3b populated)
#   - optional: chains/lcdm_reference/ exists (for delta_logZ_3b_vs_lcdm)

$ErrorActionPreference = 'Continue'
$root = Split-Path -Parent $PSScriptRoot
Set-Location $root

Write-Host "==> Working in: $root"
Write-Host ""

# Canonical paths the export script reads
$paths = @{
    "phase3a_baseline"     = "i@192.168.200.119:/media/i/storage/CLASS_SYMT/chains/phase3a_baseline"
    "phase3b_noncrossing"  = "i@192.168.200.119:/media/i/storage/CLASS_SYMT/chains/phase3b_noncrossing"
    "phase3c_b1narrow"     = "ingko@192.168.200.173:/Volumes/AppsSSD/MCMC/CLASS_SYMT/chains/phase3c_b1narrow"
}

foreach ($name in $paths.Keys) {
    $remote = $paths[$name]
    $local  = "chains\$name"
    Write-Host "==> Fetching $name"
    Write-Host "    from: $remote"
    Write-Host "    to:   $local"
    New-Item -ItemType Directory -Path $local -Force | Out-Null
    scp -r -p "$remote/*" $local
    if ($LASTEXITCODE -ne 0) {
        Write-Host "    !! scp returned $LASTEXITCODE — chain $name may be missing or incomplete" -ForegroundColor Yellow
    } else {
        $count = (Get-ChildItem $local -Filter '*.txt' -ErrorAction SilentlyContinue).Count
        Write-Host "    OK ($count .txt files)"
    }
    Write-Host ""
}

Write-Host "==> Running export_chains_summary.py"
python mcmc/export_chains_summary.py
if ($LASTEXITCODE -ne 0) {
    Write-Host "    !! export failed. If 'python' is not on PATH, try 'py' or run from WSL:" -ForegroundColor Yellow
    Write-Host "       wsl python3 mcmc/export_chains_summary.py" -ForegroundColor Yellow
    exit $LASTEXITCODE
}

Write-Host ""
Write-Host "==> Done. Reload http://localhost:8000/phase3_summary.html to see real numbers."
Write-Host "    The banner will turn green and the three placeholder slots will populate."
Write-Host ""
Write-Host "    If logZ slots still show 'pending', install MCEvidence:"
Write-Host "        pip install MCEvidence"
Write-Host "    If delta_logZ still shows 'pending', a ΛCDM reference chain is needed at:"
Write-Host "        chains/lcdm_reference/"
