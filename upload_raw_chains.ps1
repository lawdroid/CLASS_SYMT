# upload_raw_chains.ps1
# Upload the 12 raw chain .txt files (~2.8 GB) to the Tom3abc Drive folder.
# Two paths: rclone (scripted), or Google Drive for desktop (manual sync).
#
# Run from CLASS_SYMT/ root in Windows PowerShell.

param(
    [ValidateSet('rclone','desktop-sync','dryrun')]
    [string]$Mode = 'dryrun',

    [string]$RcloneRemote = 'gdrive',          # rclone remote name
    [string]$RcloneFolderId = '1NTaQAoOLeKXfCSWkLLYfu3LSLrkwppud'  # Tom3abc
)

$ErrorActionPreference = 'Stop'
$root = $PSScriptRoot
Set-Location $root

$branches = @{
    '3a' = 'phase3a_baseline'
    '3b' = 'phase3b_noncrossing'
    '3c' = 'phase3c_b1narrow'
}

Write-Host ""
Write-Host "==> Phase 3 raw-chain upload, mode: $Mode" -ForegroundColor Cyan
Write-Host ""

# Inventory
$totalMiB = 0
foreach ($br in $branches.Keys | Sort-Object) {
    $dir = "chains\$($branches[$br])"
    $txt = Get-ChildItem -Path $dir -Filter "*.txt" -File -ErrorAction SilentlyContinue
    if (-not $txt) {
        Write-Host "  $br : NO .txt files in $dir — skipping" -ForegroundColor Yellow
        continue
    }
    $sizeMiB = [math]::Round(($txt | Measure-Object Length -Sum).Sum / 1MB, 1)
    $totalMiB += $sizeMiB
    Write-Host ("  {0}: {1,4} files, {2,7} MiB" -f $br, $txt.Count, $sizeMiB)
}
Write-Host ("  total: {0} MiB" -f $totalMiB)
Write-Host ""

if ($Mode -eq 'dryrun') {
    Write-Host "Dry run — no upload performed."
    Write-Host ""
    Write-Host "To actually upload, choose one of:"
    Write-Host ""
    Write-Host "  Option A (rclone, scripted):"
    Write-Host "    1. Install rclone:    winget install Rclone.Rclone"
    Write-Host "    2. One-time setup:    rclone config  (add a Google Drive remote named 'gdrive')"
    Write-Host "    3. Re-run this script: .\upload_raw_chains.ps1 -Mode rclone"
    Write-Host ""
    Write-Host "  Option B (Drive for desktop, manual):"
    Write-Host "    1. Install Google Drive for desktop"
    Write-Host "    2. Sign in with the same Google account that owns Tom3abc"
    Write-Host "    3. Re-run this script: .\upload_raw_chains.ps1 -Mode desktop-sync"
    Write-Host "       (it'll copy chains\* into the synced Drive folder for you)"
    return
}

if ($Mode -eq 'rclone') {
    # Verify rclone is installed and the remote exists
    $rclone = Get-Command rclone -ErrorAction SilentlyContinue
    if (-not $rclone) {
        Write-Host "ERROR: rclone not on PATH. Install: winget install Rclone.Rclone" -ForegroundColor Red
        exit 1
    }
    $remotes = & rclone listremotes
    if (-not ($remotes -match "^$RcloneRemote\:$")) {
        Write-Host "ERROR: rclone remote '$RcloneRemote' not configured. Run: rclone config" -ForegroundColor Red
        Write-Host "Configured remotes:"
        Write-Host $remotes
        exit 1
    }

    foreach ($br in $branches.Keys | Sort-Object) {
        $src = "chains\$($branches[$br])"
        $dst = "${RcloneRemote}:/${branches[$br]}"
        Write-Host ""
        Write-Host "==> Uploading $br: $src  ->  $dst" -ForegroundColor Green
        # Filter to chain artifacts only; skip locks
        rclone copy `
            --drive-root-folder-id $RcloneFolderId `
            --progress `
            --transfers 2 `
            --checkers 4 `
            --exclude "*.locked" `
            $src $dst
        if ($LASTEXITCODE -ne 0) {
            Write-Host "rclone exited with $LASTEXITCODE on branch $br" -ForegroundColor Red
            exit $LASTEXITCODE
        }
    }
    Write-Host ""
    Write-Host "==> All branches uploaded. Verify in Drive:" -ForegroundColor Green
    Write-Host "    https://drive.google.com/drive/folders/$RcloneFolderId"
    return
}

if ($Mode -eq 'desktop-sync') {
    Write-Host "Where is your Google Drive for desktop root? (e.g. G:\My Drive or %USERPROFILE%\Google Drive)"
    $driveRoot = Read-Host "Drive root path"
    if (-not (Test-Path $driveRoot)) {
        Write-Host "ERROR: $driveRoot doesn't exist" -ForegroundColor Red
        exit 1
    }
    $target = Join-Path $driveRoot "PhysicsPhD\Tom3abc"
    if (-not (Test-Path $target)) {
        Write-Host "ERROR: $target doesn't exist. Make sure PhysicsPhD\Tom3abc is synced." -ForegroundColor Red
        exit 1
    }
    foreach ($br in $branches.Keys | Sort-Object) {
        $src = "chains\$($branches[$br])"
        $dst = Join-Path $target $branches[$br]
        New-Item -ItemType Directory -Path $dst -Force | Out-Null
        Write-Host "==> Copying $br: $src  ->  $dst"
        Copy-Item -Path "$src\*" -Destination $dst -Recurse -Force -Exclude "*.locked"
    }
    Write-Host ""
    Write-Host "==> Copies done. Drive desktop sync will upload in background." -ForegroundColor Green
    return
}
