param([Parameter(Mandatory=$true)][string]$Jar)
$ErrorActionPreference = 'Stop'
$artifact = (Resolve-Path -LiteralPath $Jar).Path
Push-Location $PSScriptRoot
try {
    & java -jar $artifact study auto study.yaml --offline
    if ($LASTEXITCODE -ne 0) { throw "Study execution failed with exit code $LASTEXITCODE" }
    & java -jar $artifact study verify --project .
    if ($LASTEXITCODE -ne 0) { throw "Evidence verification failed with exit code $LASTEXITCODE" }
} finally { Pop-Location }
