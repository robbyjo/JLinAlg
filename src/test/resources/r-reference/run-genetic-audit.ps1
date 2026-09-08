# Run from repository root. Gradle resolves dependencies on a clean checkout.
[CmdletBinding()]
param(
    [switch]$RegenerateR,
    [string]$Rscript = 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe'
)

$ErrorActionPreference = 'Stop'
if (!(Test-Path -LiteralPath './gradlew.bat')) {
    throw 'Run this script from the JLinAlg repository root.'
}
if ($RegenerateR) {
    & $Rscript src/test/resources/r-reference/generate-genetic-audit-reference.R
    if ($LASTEXITCODE -ne 0) { throw 'R reference generation failed' }
    & $Rscript src/test/resources/r-reference/generate-genetic-final-fixes.R
    if ($LASTEXITCODE -ne 0) { throw 'Final independent-review R reference generation failed' }
}
& ./gradlew.bat check
if ($LASTEXITCODE -ne 0) { throw 'Gradle check failed' }
& ./gradlew.bat benchmarkGeneticAudit
if ($LASTEXITCODE -ne 0) { throw 'accuracy-gated benchmark failed' }
