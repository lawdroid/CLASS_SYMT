# upload_chains_to_gdrive.ps1
# Upload the raw Phase 3 chains (too big for GitHub: 8 files > 100 MB) to the
# Google Drive folder that already holds MANIFEST.sha256.
#
# Target folder: https://drive.google.com/drive/folders/1NTaQAoOLeKXfCSWkLLYfu3LSLrkwppud
#
# WHY rclone and not the Claude/Drive connector: the connector carries content
# as base64 in a single call — fine for the ~2 KB manifest, impossible for
# 200 MB chain files. rclone streams natively and resumes on interruption.
#
# ONE-TIME SETUP (only needed once on this machine):
#   1. Install rclone:  winget install Rclone.Rclone     (or https://rclone.org/downloads/)
#   2. Create a Drive remote named 'gdrive':
#        rclone config
#        -> n (new remote) -> name: gdrive -> storage: drive
#        -> follow the browser OAuth; accept defaults; confirm.
#   That writes credentials to %APPDATA%\rclone\rclone.conf. Done once.
#
# Run from CLASS_SYMT\ :   .\mcmc\upload_chains_to_gdrive.ps1

param(
    [string]$Remote   = "gdrive",
    [string]$FolderId = "1NTaQAoOLeKXfCSWkLLYfu3LSLrkwppud"
)

$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $PSScriptRoot
Set-Location $root
Write-Host "==> Working in: $root"

# Fail early if rclone isn't installed or the remote isn't configured.
if (-not (Get-Command rclone -ErrorAction SilentlyContinue)) {
    Write-Host "!! rclone not found. Install it (see the SETUP block at the top of this script)." -ForegroundColor Red
    exit 1
}
$remotes = rclone listremotes
if ($remotes -notcontains "${Remote}:") {
    Write-Host "!! Remote '${Remote}:' not configured. Run 'rclone config' first (see SETUP block)." -ForegroundColor Red
    Write-Host "   Configured remotes: $remotes" -ForegroundColor Yellow
    exit 1
}

# Only the two real chains that today's export consumed — matches MANIFEST.sha256.
# (The old phase3a_FULL_BACKUP_* and the empty phase3c dir are deliberately skipped.)
$branches = @("phase3a_baseline", "phase3b_noncrossing")

foreach ($b in $branches) {
    $src = "chains\$b"
    if (-not (Test-Path $src)) {
        Write-Host "   skip $b (not present locally)" -ForegroundColor Yellow
        continue
    }
    Write-Host ""
    Write-Host "==> Uploading $src  ->  ${Remote}:/$b  (folder id $FolderId)"
    # --drive-root-folder-id makes the remote root = the shared folder, so the
    # destination path '/$b' becomes a subfolder inside it.
    rclone copy $src "${Remote}:/$b" `
        --drive-root-folder-id $FolderId `
        --progress --transfers 4 --checkers 8 `
        --drive-chunk-size 64M
}

# Refresh the manifest at the folder root too (rclone skips it if identical).
Write-Host ""
Write-Host "==> Uploading MANIFEST.sha256 to folder root"
rclone copy "chains\MANIFEST.sha256" "${Remote}:/" --drive-root-folder-id $FolderId

Write-Host ""
Write-Host "==> Done. Verify integrity after upload from any machine with the folder synced:"
Write-Host "      sha256sum -c MANIFEST.sha256        # Linux/Mac, run inside the chains dir"
Write-Host "      Get-FileHash <file> -Algorithm SHA256   # Windows, compare to the manifest"
