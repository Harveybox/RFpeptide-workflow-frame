$ErrorActionPreference = "Stop"

$script = Join-Path $PSScriptRoot "workflow_gui.py"
$candidates = @(
    "py",
    "python3",
    "python",
    (Join-Path $env:USERPROFILE ".cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe")
)

foreach ($candidate in $candidates) {
    if ($candidate -match '[\\/]' -and -not (Test-Path $candidate)) {
        continue
    }

    $command = Get-Command $candidate -ErrorAction SilentlyContinue
    if ($command -or (Test-Path $candidate)) {
        & $candidate $script
        exit $LASTEXITCODE
    }
}

Write-Host "No Python interpreter found. Install Python 3 or run workflow_gui.py with an existing Python environment."
exit 1
