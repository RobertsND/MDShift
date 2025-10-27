# ==========================================
# Run_MDShift.ps1
# Copies MDShift.py into the specified data directory,
# executes it there, and deletes it afterward.
# Usage:
#     .\Run_MDShift.ps1 <DataDir>
# Examples:
#     .\Run_MDShift.ps1 "C:\Users\Henry Eyring\Documents\data"
#     .\Run_MDShift.ps1 cwd
# ==========================================

param(
    [string]$DataDir = "cwd"
)

# --- Define source script path ---
$ScriptPath = "D:\Grad school\Research\Python_scripts_for_data_analysis\MDShift\MDShift.py"

# --- Resolve DataDir ---
if ($DataDir -in @("cwd", "pwd", ".", "./")) {
    $DataDir = (Get-Location).Path
}

if (-not (Test-Path $DataDir)) {
    Write-Host "Error: Data directory not found at $DataDir" -ForegroundColor Red
    exit 1
}

$DestScript = Join-Path $DataDir "MDShift.py"

# --- Verify source exists ---
if (-Not (Test-Path $ScriptPath)) {
    Write-Host "Error: MDShift.py not found at $ScriptPath" -ForegroundColor Red
    exit 1
}

# --- Copy to data directory ---
Write-Host "Copying MDShift.py to $DataDir..."
Copy-Item -Path $ScriptPath -Destination $DestScript -Force

# --- Change to data directory ---
Push-Location $DataDir

# --- Run the Python script ---
Write-Host "Executing MDShift.py with arguments: label.txt 1.82..."
python "MDShift.py" "label.txt" "1.82"

# --- Check if Python completed successfully ---
if ($LASTEXITCODE -eq 0) {
    Write-Host "MDShift.py finished successfully."
} else {
    Write-Host "Warning: MDShift.py exited with code $LASTEXITCODE" -ForegroundColor Yellow
}

# --- Delete the copied Python script ---
Write-Host "Deleting MDShift.py from $DataDir..."
Remove-Item -Path $DestScript -Force -ErrorAction SilentlyContinue

# --- Return to original directory ---
Pop-Location

Write-Host "Run complete."
# ==========================================
