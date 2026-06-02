# upload_simple.ps1
# Simplest possible Drive upload: opens File Explorer at chains/ AND the
# Drive folder in your browser, you drag-and-drop. No rclone, no
# Drive desktop sync, no setup. Two confirmation steps.
#
# Run from CLASS_SYMT/ root:  .\upload_simple.ps1

$ErrorActionPreference = 'Stop'
$root = $PSScriptRoot
Set-Location $root

$driveUrl    = "https://drive.google.com/drive/folders/1NTaQAoOLeKXfCSWkLLYfu3LSLrkwppud"
$chainsPath  = Join-Path $root "chains"

if (-not (Test-Path $chainsPath)) {
    Write-Host "ERROR: $chainsPath not found." -ForegroundColor Red
    exit 1
}

# ---------- Confirmation 1: show what's being uploaded ----------
Write-Host ""
Write-Host "================ Phase 3 chains to upload ================" -ForegroundColor Cyan
$total = 0
foreach ($br in @('phase3a_baseline','phase3b_noncrossing','phase3c_b1narrow')) {
    $dir = Join-Path $chainsPath $br
    if (-not (Test-Path $dir)) {
        Write-Host ("  {0,-25} MISSING" -f $br) -ForegroundColor Yellow
        continue
    }
    $files = Get-ChildItem $dir -File
    $size  = [math]::Round(($files | Measure-Object Length -Sum).Sum / 1MB, 1)
    $total += $size
    Write-Host ("  {0,-25} {1,4} files   {2,7} MiB" -f $br, $files.Count, $size)
}
Write-Host ("  {0,-25}              {1,7} MiB" -f 'TOTAL', $total) -ForegroundColor Cyan
Write-Host "  Target: PhysicsPhD / Tom3abc (general-access)"
Write-Host "          $driveUrl"
Write-Host "==========================================================" -ForegroundColor Cyan
Write-Host ""

$reply = Read-Host "Open File Explorer + Drive folder in browser for drag-and-drop? [Y/n]"
if ($reply -and $reply.ToLower().StartsWith('n')) {
    Write-Host "Cancelled. Nothing opened." -ForegroundColor Yellow
    exit 0
}

# ---------- Open the two windows ----------
Start-Process explorer.exe $chainsPath
Start-Sleep -Milliseconds 800
Start-Process $driveUrl

Write-Host ""
Write-Host "Two windows just opened." -ForegroundColor Green
Write-Host ""
Write-Host "How to upload (~5-10 min on a 10 MB/s connection):"
Write-Host "  1. In File Explorer (chains\), select all three folders:"
Write-Host "     phase3a_baseline, phase3b_noncrossing, phase3c_b1narrow"
Write-Host "     (Ctrl+A, or click first + Shift-click last)"
Write-Host "  2. Drag them into the Drive browser tab"
Write-Host "  3. Drive will show an upload progress overlay at the bottom"
Write-Host ""
Write-Host "Tip: if you're signed into multiple Google accounts, the Drive tab"
Write-Host "needs to be on the same account that owns Tom3abc (ingyukoh2@gmail.com)."
Write-Host ""

# ---------- Confirmation 2: wait for completion ----------
Read-Host "Press Enter once the Drive UI shows 'uploads complete'"

Write-Host ""
Write-Host "Done. Verify in browser:" -ForegroundColor Green
Write-Host "  $driveUrl"
Write-Host ""
Write-Host "Optional sanity check (re-open in Drive):"
Write-Host "  - all 12 chain .txt files visible"
Write-Host "  - sizes match: 3a ~135 MB each, 3b ~210 MB each, 3c ~340 MB each"
Write-Host "  - checksums in MANIFEST_2026-06-02_complete.sha256 match"
