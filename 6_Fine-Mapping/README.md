# Meta-locus fine-mapping workflow

This repository contains fine-mapping workflow used for sex- and generation-stratified smoking-initiation analyses. 



```text
scripts/
  01_harmonize_metal_meta.py
  01_build_meta_locus_inputs.py
  02_make_joblists.py
  03_extract_locus_sumstats_metal.py
  03_susie_completion_qc.py
  04_combine_susie_outputs.py
  05_summarize_highpip_patterns.py
  06_meta_locus_evidence_tables.py
  07_exact_variant_overlap.py
  08_plot_moderated_maxpip.py
  run_susie_metal_weighted_ld_one_locus.R
  utils.py

bash/
  run_extract_locus_sumstats.sh
  submit_susie_arrays_template.sh
  run_weightedLD_susie_metal_array.slurm

config/
  example_config.yaml
  example_manifest.tsv
  example_ld_weights.tsv
```

## Required inputs

The workflow assumes these data types are prepared externally:

1. Stratified METAL GWAS summary statistics.
2. GRCh38 variant map files with at least `CHR`, `BP`, `variant_id`, `REF`, and `ALT`.
3. Sex- and generation-defined GWAS locus tables lifted to GRCh38.
4. 1000G GRCh38 PLINK2 PGEN reference panels by ancestry.
5. LD weight file with either `analysis_id`, `ancestry`, `weight` or `analysis_id`, `ancestry`, `N_ancestry`.

## Workflow

### 1. Harmonize METAL summary statistics

Create a manifest such as `config/example_manifest.tsv` with columns:

```text
analysis_id    gwas_file
sex_female_meta    /path/to/female_METAL.txt
sex_male_meta      /path/to/male_METAL.txt
```

Then run one analysis at a time:

```bash
python scripts/01_harmonize_metal_meta.py \
  --analysis-id sex_female_meta \
  --manifest config/example_manifest.tsv \
  --map-dir /path/to/variant_maps \
  --out-dir outputs/sex_female_meta/harmonized_META_N80 \
  --chunksize 100000
```

The script writes chromosome-level files named:

```text
smoking_META_chr<CHR>_harmonized.tsv
```

### 2. Build meta-locus coordinates and buffered union file

```bash
python scripts/01_build_meta_locus_inputs.py \
  --sex-loci-hg38 path/to/sex_loci_main_and_moderated.hg38.tsv \
  --generation-loci-hg38 path/to/generation_loci_main_and_moderated.hg38.tsv \
  --meta-locus-mapping path/to/sex_generation_original_loci_to_meta_loci_mapping.tsv \
  --meta-source-summary path/to/meta_source_summary.csv \
  --out-dir outputs/meta_locus \
  --buffer-bp 400000
```

### 3. Extract locus-level summary statistics

```bash
python scripts/03_extract_locus_sumstats_metal.py \
  --analysis-id sex_female_meta \
  --harm-dir outputs/sex_female_meta/harmonized_META_N80 \
  --union-file outputs/meta_locus/meta_loci_union_hg38_buffered_for_sumstats.tsv \
  --out-dir outputs/meta_locus/sex_female_meta/locus_sumstats \
  --min-snps 10
```

A convenience template is included:

```bash
bash bash/run_extract_locus_sumstats.sh
```

### 4. Create A/B/C SuSiE joblists

```bash
python scripts/02_make_joblists.py \
  --summary outputs/meta_locus/sex_female_meta/locus_sumstats/locus_sumstats_extraction_summary.tsv \
  --out-dir outputs/sex_female_meta/meta_locus_400kb/joblists
```

### 5. Run weighted-LD SuSiE on a SLURM job

Submit the array template using environment variables:

```bash
sbatch \
  --array=1-500%50 \
  --export=ALL,PROJECT=/path/to/project,PGEN_DIR=/path/to/1000G_GRCh38_pgen,ANALYSIS_ID=sex_female_meta,METAL_MODE=meta_locus_400kb,JOB_GROUP=A,RSCRIPT=scripts/run_susie_metal_weighted_ld_one_locus.R,WEIGHTS_FILE=config/example_ld_weights.tsv \
  bash/run_weightedLD_susie_metal_array.slurm
```

Or use the chunked submit template:

```bash
ANALYSIS_ID=sex_female_meta \
JOB_GROUP=A \
WEIGHTS_FILE=config/example_ld_weights.tsv \
JOBLIST=outputs/sex_female_meta/meta_locus_400kb/joblists/jobs_A_nsnps_le5000_minimal.tsv \
bash bash/submit_susie_arrays_template.sh
```

### 6. Check completion

```bash
python scripts/03_susie_completion_qc.py \
  --expected-loci outputs/meta_locus/meta_loci_union_hg38_buffered_for_sumstats.tsv \
  --result-root outputs \
  --metal-mode meta_locus_400kb \
  --out-dir outputs/meta_locus/completion_qc
```

### 7. Combine SuSiE outputs

```bash
python scripts/04_combine_susie_outputs.py \
  --result-root outputs \
  --metal-mode meta_locus_400kb \
  --union-file outputs/meta_locus/meta_loci_union_hg38_buffered_for_sumstats.tsv \
  --out-dir outputs/meta_locus/combined_susie_results
```

### 8. Summarize high-PIP patterns and core/buffer support

```bash
python scripts/05_summarize_highpip_patterns.py \
  --pip-wide outputs/meta_locus/combined_susie_results/all_analyses_pip_wide.tsv \
  --pip-long outputs/meta_locus/combined_susie_results/all_analyses_pip_long_with_core_buffer.tsv \
  --out-dir outputs/meta_locus/combined_susie_results
```

### 9. Create evidence, exact-overlap, and plotting tables

```bash
python scripts/06_meta_locus_evidence_tables.py --help
python scripts/07_exact_variant_overlap.py --help
python scripts/08_plot_moderated_maxpip.py --help
```

## Definitions

- High-confidence variant: `PIP >= 0.8`.
- Primary support: high-PIP variant inside the core collapsed meta-locus interval.
- Buffer-only variants are excluded from primary summaries.
- A meta-locus can contain main-effect and moderated evidence; scripts can use moderated-priority classification where a mutually exclusive class is needed.
