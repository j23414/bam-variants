# bam-variants

A modular workflow for generating a table of variants from alignment information.

## Usage

```
nextflow run j23414/bam-variants \
  --bam [path/*.bam] \
  --samplesheet [path/bam_samplesheet.csv] \
  --reference [path/reference.fasta] \
  --outdir "variant-results" \
  --depth 1000000 \
  --bcftools_mpileup_params "--per-sample-mF --count-orphans --no-BAQ -Ov --min-BQ 13 --min-MQ 30 --annotate AD,ADF,ADR,DP,SP " \
  --bcftools_call_params "-mv -Oz --prior-freqs AN,AC, -A --variants-only --keep-alts --keep-masked-ref" \
  --minor_variants_depth 20 \
  -profile stjude
```
