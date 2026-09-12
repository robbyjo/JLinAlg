param(
    [string]$LinuxWorkspace = '/mnt/e/Projects/JLinAlg',
    [string]$Jar = 'build/cli/jlinalg-0.3.5.jar'
)
$ErrorActionPreference = 'Stop'
$benchmark = Join-Path (Get-Location) 'build/rare-meta-benchmark'
$run = [Guid]::NewGuid().ToString('N')
$results = @()
# Run outside the restricted tool filesystem sandbox. The same compressed
# score files are used; native RAREMETAL also writes its standard plots.
for ($i = 0; $i -lt 4; $i++) {
    $prefix = Join-Path $benchmark "java-$run-$i"
    $watch = [Diagnostics.Stopwatch]::StartNew()
    & java -jar $Jar rare-meta --cohorts (Join-Path $benchmark 'cohorts.tsv') --genome-build synthetic --test single --af-policy raremetal --out $prefix 2> "$prefix.stderr" | Out-Null
    if ($LASTEXITCODE -ne 0) { throw 'Java benchmark failed' }
    $watch.Stop()
    $results += [pscustomobject]@{program='JLinAlg'; run=$i; seconds=$watch.Elapsed.TotalSeconds}
}
for ($i = 0; $i -lt 4; $i++) {
    $prefix = "native-$run-$i"
    & wsl -d Ubuntu --cd "$LinuxWorkspace/build/rare-meta-benchmark" -- /usr/bin/time -f '%e' -o "$prefix.time" "$LinuxWorkspace/build/raremetal-reference/build/raremetal" --summaryFiles summaryfiles --covFiles covfiles --noPhoneHome --prefix $prefix > (Join-Path $benchmark "$prefix.console") 2>&1
    if ($LASTEXITCODE -ne 0) { throw 'RAREMETAL benchmark failed' }
    $seconds = [double]::Parse((Get-Content (Join-Path $benchmark "$prefix.time")),[Globalization.CultureInfo]::InvariantCulture)
    $results += [pscustomobject]@{program='RAREMETAL'; run=$i; seconds=$seconds}
}
$results | Export-Csv (Join-Path $benchmark 'timings.csv') -NoTypeInformation
$results | Format-Table
