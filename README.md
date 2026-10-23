# Automated rare disease variant pipeline

`auto_pipeline.R` takes a patient's HPO phenotype terms and a whole-genome VCF (GRCh38). It ranks the variants most likely to explain the condition and writes up to 10 rows to a submission file. It runs the same way every time, needs no manual curation, and makes no network requests.

The output is a ranked list of research candidates. The scores are hand-set heuristics, not calibrated probabilities, and a high score does not confirm a diagnosis.

## Requirements

- Linux, with R 4.4 or later. The pipeline was written and run on Linux. Windows users can run it in WSL2 (Windows Subsystem for Linux) with a Linux distribution such as Ubuntu, and install R, the Bioconductor packages, and HTSlib inside it. It has not been tested on native Windows or macOS.
- These Bioconductor packages installed:
  - `VariantAnnotation`
  - `GenomicScores`
  - `TxDb.Hsapiens.UCSC.hg38.knownGene`
  - `BSgenome.Hsapiens.UCSC.hg38`
  - `MafH5.gnomAD.v4.0.GRCh38`
  - `org.Hs.eg.db`
- A bgzipped VCF with a tabix index (`.vcf.gz` and `.vcf.gz.tbi`) in the same folder. It must contain one sample.
- The HPO file `phenotype_to_genes.txt` (from the HPO project downloads).

## Quick start

1. Edit `SEED_GENES` in the CONFIG section of `auto_pipeline.R` if you suspect a particular disorder.
2. Run it with your file paths and the patient's HPO terms:

```bash
Rscript auto_pipeline.R \
  --vcf /path/to/input/proband_sample.vcf.gz \
  --hpo /path/to/hpo/phenotype_to_genes.txt \
  --terms HP:0002859,HP:0000121,HP:0004322 \
  --out /path/to/output/Martins_auto_pipeline_submission.csv \
  --audit /path/to/output/audit_variants.tsv
```

The script prints progress for stages 1 to 6, then a summary table of the findings. Run `Rscript auto_pipeline.R --help` to see the options.

### Command line options

| Option | Meaning |
| --- | --- |
| `--vcf` | Input VCF. Required. Must have a `.tbi` index next to it. |
| `--hpo` | HPO file `phenotype_to_genes.txt`. Required. |
| `--terms` | The patient's HPO term IDs, separated by commas, for example `HP:0002859,HP:0000121`. Defaults to the proband's eight terms. The script stops if an ID is not in the form `HP:` plus seven digits, and warns if the HPO file has no genes for a term. |
| `--out` | Where to write the submission file. Defaults to `Martins_auto_pipeline_submission.csv` in the current folder. |
| `--audit` | Where to write the audit file. Defaults to `audit_variants.tsv` in the current folder. |

The script stops with a message if a required option or input file is missing.

## Configuration

File paths are command line options (see above). The other settings are at the top of the script.

| Setting | What it does |
| --- | --- |
| `PROBAND_ID` | Value written to the `proband_id` column in every row (`PROBAND01`). The challenge scorer accepts only this identifier. |
| `SEED_GENES` | Genes that are always searched, for a suspected disorder. Use `character(0)` for a purely phenotype driven run. |
| `N_PANEL` | How many HPO ranked genes to search (500). |
| `PAD` | Bases of flanking sequence read around each gene (5000). |
| `MAX_FINDINGS` | Number of rows in the submission (10). |
| `N_PRIMARY` | How many top rows are labelled `primary`. The rest are `secondary`. |
| `MIN_DP`, `MIN_GQ`, `MIN_AB` | Minimum read depth, genotype quality, and allele balance for heterozygous calls. |
| `MAX_AF_DOM`, `MAX_AF_REC` | Highest gnomAD allele frequency allowed for dominant and recessive models. |

## What the script does

1. **Gene panel.** Weights each HPO term by how specific it is, scores genes by the terms they match, and keeps the top `N_PANEL` genes plus the seed genes.
2. **Read the VCF.** Reads only the panel genes and their flanking sequence, one row per alternate allele.
3. **Quality filters.** Keeps non-reference calls that pass the filter, depth, genotype quality, and allele balance thresholds.
4. **Annotation.** Adds genomic region, coding effect, amino acid change, and gnomAD v4 frequency. A variant missing from gnomAD is treated as frequency zero. Writes the audit file.
5. **Findings.** Builds three kinds of finding from protein altering variants:
   - `homozygous`: one homozygous variant, frequency at most `MAX_AF_REC`.
   - `compound_het`: two heterozygous variants in one gene. Pairs proven to be on the same copy (in cis) by GATK phasing are dropped. Pairs proven to be on opposite copies (in trans) score higher. Pairs without phasing score slightly lower.
   - `heterozygous_lof`: one loss of function variant, frequency at most `MAX_AF_DOM`.
6. **Score and write.** Ranks findings by score, converts the score to a 0 to 1 value, and writes the top `MAX_FINDINGS` rows.

## Output

### Submission file

Written to the `--out` path (default `Martins_auto_pipeline_submission.csv`). It is comma separated, with one row per finding. Text fields are in double quotes, so commas inside `notes` do not split the row. The script sets the file so that only your user can read it.

| Column | Meaning |
| --- | --- |
| `proband_id` | Proband identifier from the dataset, the same in every row (`PROBAND01`). It is set by `PROBAND_ID` in the script. The sample name in the VCF is not used. |
| `chrom_1`, `pos_1`, `ref_1`, `alt_1` | The first variant. |
| `chrom_2`, `pos_2`, `ref_2`, `alt_2` | The second variant of a compound heterozygous pair. Empty for single variant findings. |
| `epcr` | Ranking score from 0.02 to 0.95. Use it to compare findings against each other. It is not a probability. |
| `finding_type` | `primary` or `secondary`. |
| `notes` | Plain text summary that starts with the gene symbol, then the consequence, amino acid change, frequency, quality values, and phasing status. |

### Audit file

Written to the `--audit` path (default `audit_variants.tsv`). It is tab separated, and like the submission file it is readable only by your user. It contains patient genotypes, so do not share it. It lists every annotated variant in the panel genes that passed quality filtering, with gene, position, genotype, consequence, frequency, depth, genotype quality, allele balance, and phasing fields. Use it to check why a variant was or was not reported.

## Limitations

- Only small variants in the VCF are used. Structural variants, copy number changes, and deep intronic variants are not detected.
- Only the panel genes are searched. A causal variant in a gene outside the panel is never seen.
- Score weights are chosen by hand and are not fitted to data.
- Without phasing information, a compound heterozygous pair may actually be on the same copy of the gene. Parental testing is needed to confirm it.

## License

Copyright 2026 Kiley Martin.

This work is licensed under the Creative Commons Attribution 4.0 International License (CC BY 4.0). See `LICENSE`.
