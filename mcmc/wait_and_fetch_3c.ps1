# wait_and_fetch_3c.ps1
# Persistent retry: poll the Mac mini until reachable, then fetch the 3c chain
# and re-export chains_summary.json with live 3c data.
#
# Run from Windows PowerShell (not WSL). Leave it running while you go wake
# the Mac mini; it will catch the host the moment it comes back online.
#
# Usage:
#   .\mcmc\wait_and_fetch_3c.ps1
#   .\mcmc\wait_and_fetch_3c.ps1 -PollSeconds 15 -MaxAttempts 480   # 2 hours

param(
    [string]$Host      = "192.168.200.173",
    [string]$User      = "ingko",
    [string]$RemotePath= "/Volumes/AppsSSD/MCMC/CLASS_SYMT/chains/phase3c_b1narrow",
    [int]$PollSeconds  = 15,
    [int]$MaxAttempts  = 240   # default ~1 hour at 15 s/attempt
)

$ErrorActionPreference = 'Continue'
$root = Split-Path -Parent $PSScriptRoot
Set-Location $root
$localDir = "chains\phase3c_b1narrow"
New-Item -ItemType Directory -Path $localDir -Force | Out-Null

Write-Host ""
Write-Host "==> wait_and_fetch_3c"
Write-Host "    target: $User@${Host}:$RemotePath"
Write-Host "    local : $root\$localDir"
Write-Host "    poll  : every $PollSeconds s, up to $MaxAttempts attempts"
Write-Host ""
Write-Host "    Go wake the Mac mini now. This script will fetch as soon as"
Write-Host "    port 22 is reachable. Ctrl-C to stop."
Write-Host ""

$attempt = 0
$reachable = $false
while ($attempt -lt $MaxAttempts -and -not $reachable) {
    $attempt++
    # WarningAction silences the noisy "TCP connect failed" yellow text on each miss
    $result = Test-NetConnection $Host -Port 22 -WarningAction SilentlyContinue
    if ($result.TcpTestSucceeded) {
        $reachable = $true
        Write-Host "[$attempt/$MaxAttempts] $Host:22 REACHABLE — fetching" -ForegroundColor Green
    } else {
        $stamp = Get-Date -Format 'HH:mm:ss'
        Write-Host "[$attempt/$MaxAttempts $stamp] $Host:22 unreachable ($($result.PingReplyDetails.Status)) — retry in ${PollSeconds}s"
        Start-Sleep -Seconds $PollSeconds
    }
}

if (-not $reachable) {
    Write-Host ""
    Write-Host "Gave up after $MaxAttempts attempts. Mac mini never came online." -ForegroundColor Red
    Write-Host "Possible causes:"
    Write-Host "  - Hardware off / unplugged"
    Write-Host "  - Different subnet / VPN / Wi-Fi network than this PC"
    Write-Host "  - sshd disabled (would show 'Connection refused', not Unreachable)"
    exit 1
}

Write-Host ""
Write-Host "==> scp -r $User@${Host}:$RemotePath/* $localDir\"
scp -r -p "${User}@${Host}:${RemotePath}/*" $localDir
if ($LASTEXITCODE -ne 0) {
    Write-Host "scp failed with exit code $LASTEXITCODE" -ForegroundColor Red
    exit $LASTEXITCODE
}

$count = (Get-ChildItem $localDir -Filter '*.txt' -ErrorAction SilentlyContinue).Count
$bytes = (Get-ChildItem $localDir -Recurse -File | Measure-Object -Property Length -Sum).Sum
$mib   = [math]::Round($bytes / 1MB, 1)
Write-Host ""
Write-Host "==> Fetched $count .txt files ($mib MiB) into $localDir" -ForegroundColor Green

Write-Host ""
Write-Host "==> Re-running export_chains_summary.py"
python mcmc/export_chains_summary.py
if ($LASTEXITCODE -ne 0) {
    Write-Host "export failed — try: wsl python3 mcmc/export_chains_summary.py" -ForegroundColor Yellow
    exit $LASTEXITCODE
}

Write-Host ""
Write-Host "==> Done. The 3c snapshot block has been replaced with live data."
Write-Host "    Reload http://localhost:8000/phase3_summary.html — banner flips green."
