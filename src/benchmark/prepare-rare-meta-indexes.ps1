param([string]$LinuxWorkspace = '/mnt/e/Projects/JLinAlg')
$ErrorActionPreference = 'Stop'
# Run after prepare_rare_meta_benchmark.R, with bgzip and tabix installed in Ubuntu.
# bgzip retains the plain input and refuses to replace existing compressed files.
foreach ($cohort in 1,2) {
    foreach ($kind in 'score','cov') {
        $file = "$LinuxWorkspace/build/rare-meta-benchmark/cohort$cohort.$kind"
        & wsl -d Ubuntu -- bgzip -k $file
        if ($LASTEXITCODE -ne 0) { throw "bgzip failed: $file" }
        & wsl -d Ubuntu -- tabix -c '#' -s 1 -b 2 -e 2 "$file.gz"
        if ($LASTEXITCODE -ne 0) { throw "tabix failed: $file.gz" }
    }
}
