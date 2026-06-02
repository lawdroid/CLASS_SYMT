# view_phase3.ps1
# One command to view Phase 3 results with a confirmation step.
# Run from CLASS_SYMT/ root in PowerShell:
#   .\view_phase3.ps1

$ErrorActionPreference = 'Stop'
$root = $PSScriptRoot
Set-Location $root

$summary = Join-Path $root 'chains_summary.json'
if (-not (Test-Path $summary)) {
    Write-Host "ERROR: chains_summary.json not found at $summary" -ForegroundColor Red
    Write-Host "Run: .\mcmc\run_phase3_pipeline.ps1  (fetches chains, regenerates summary)"
    exit 1
}

# ---------- Confirmation ----------
$json = Get-Content $summary -Raw | ConvertFrom-Json
$src  = $json.data_source
$gen  = $json.generated_utc

Write-Host ""
Write-Host "==================== Phase 3 data ====================" -ForegroundColor Cyan
Write-Host "  Source     : $src"
Write-Host "  Generated  : $gen"
Write-Host "  Branches   :"
foreach ($p in $json.branches.PSObject.Properties) {
    $br  = $p.Name
    $b   = $p.Value
    $n   = $b.n_samples
    $R1  = if ($b.R1 -eq $null) { 'n/a' } else { '{0:F4}' -f $b.R1 }
    $lbl = $b.label
    Write-Host ("    {0}: n_samples={1}, R-1={2}" -f $br, $n, $R1)
    Write-Host ("         {0}" -f $lbl) -ForegroundColor DarkGray
}
Write-Host "======================================================" -ForegroundColor Cyan
Write-Host ""

if ($src -match 'MOCK') {
    Write-Host "WARNING: data_source still contains 'MOCK' — page will show fake numbers." -ForegroundColor Yellow
}
if ($src -match 'snapshot|missing|offline') {
    Write-Host "NOTE: provenance is mixed (partial real / snapshot). See data_source for details." -ForegroundColor Yellow
}

$reply = Read-Host "Open http://localhost:8000/phase3_summary.html in your browser? [Y/n]"
if ($reply -and $reply.ToLower().StartsWith('n')) {
    Write-Host "Cancelled. Nothing opened." -ForegroundColor Yellow
    exit 0
}

# ---------- Make sure something is serving the folder on :8000 ----------
$serverUp = $false
try {
    $probe = Test-NetConnection -ComputerName 'localhost' -Port 8000 -WarningAction SilentlyContinue
    $serverUp = $probe.TcpTestSucceeded
} catch { $serverUp = $false }

if (-not $serverUp) {
    Write-Host "Starting local HTTP server on :8000 (python -m http.server)" -ForegroundColor DarkGray
    Start-Process -FilePath python `
        -ArgumentList @('-m','http.server','8000','--bind','127.0.0.1') `
        -WorkingDirectory $root `
        -WindowStyle Minimized | Out-Null
    Start-Sleep -Seconds 2
}

Start-Process 'http://localhost:8000/phase3_summary.html'
Write-Host ""
Write-Host "Opened: http://localhost:8000/phase3_summary.html" -ForegroundColor Green
Write-Host "Banner color indicates provenance:"
Write-Host "  green = all real, all branches live"
Write-Host "  amber = partial / snapshot / missing branches"
Write-Host "  red   = mock"
