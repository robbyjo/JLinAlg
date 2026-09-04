# JLinAlg LD reference format

Every downloadable linkage-disequilibrium reference is normalized to the same
directory contract:

    REFERENCE/
    |-- jlinalg-ld-reference.json
    \-- panels/
        |-- AFR/
        |   |-- genotypes.bed
        |   |-- genotypes.bim
        |   \-- genotypes.fam
        \-- EUR/
            |-- genotypes.bed
            |-- genotypes.bim
            \-- genotypes.fam

The set of panel identifiers depends on the database. File names and encoding
do not. LdReferenceLayout.panelPrefix(reference, panel) returns the prefix that
JLinAlg and external PLINK-compatible tools can consume.

## Version 1

- Genotypes are stored in variant-major PLINK 1 BED form: two bits per sample
  and direct byte-offset access to every variant.
- Variant metadata use the six-column PLINK BIM file associated with the BED.
  BIM allele 1 and allele 2 are preserved exactly from the source.
- Sample metadata use PLINK FAM.
- jlinalg-ld-reference.json records the format version, database ID, genome
  build, source URI and checksum, panel identities, canonical prefixes, and
  validated sample and variant counts.

During installation JLinAlg checks the archive checksum, normalizes source
names into this layout, verifies the BED magic bytes and variant-major flag,
and confirms that the BED length agrees exactly with the BIM and FAM counts.
The manifest is written last, so its presence marks a complete installation.
Failed installations remove their temporary and normalized data files.

## CLI

List supported sources:

    java -jar jlinalg-<version>.jar ld-db list

Install into the current directory:

    java -jar jlinalg-<version>.jar ld-db download --database 1000g-phase3

Install into a selected directory:

    java -jar jlinalg-<version>.jar ld-db download \
      --database 1000g-phase3 --location /data/ld/1000g-phase3

An omitted or empty --location means the current directory. Downloaded
archives are removed after successful normalization.

## LD clumping

Clump MR candidates against a panel from an installed database:

    java -jar jlinalg-<version>.jar clump \
      --database /data/ld/1000g-phase3 \
      --instrument exposure.mr.tsv \
      --population EUR \
      --ld-threshold 0.001 \
      --output exposure.clumped.tsv

The defaults match `TwoSampleMR::clump_data`: population `EUR`, a 10,000 kb
window, r-squared threshold 0.001, and index p-value threshold 1. Override the
last two non-r-squared settings with `--clump-kb` and `--p-threshold`.

The command recognizes `SNP`/`rsid`, uses `pval.exposure` before
`pval.outcome` and `pval`, and groups by `id.exposure` or `id`. Explicit column
overrides are available as `--snp-column`, `--pval-column`, and `--id-column`.
CSV, TSV, and gzip-compressed inputs are supported.

Within each multi-variant group, the lowest-p-value candidate is considered
first. Nearby candidates are removed when the founder-only haplotype maximum-
likelihood r-squared is strictly greater than the requested threshold, matching
PLINK 1.9 clumping semantics. Candidates absent from the reference panel are
removed. For compatibility with `ieugwasr::ld_clump`, a group containing only
one candidate is returned without consulting the panel. The output contains
the retained original rows and columns in input order.

Use `--overwrite` to replace an existing output. JLinAlg validates the database
manifest and BED/BIM/FAM sizes before clumping, and writes through a temporary
file so a failed run cannot leave a partial result.
