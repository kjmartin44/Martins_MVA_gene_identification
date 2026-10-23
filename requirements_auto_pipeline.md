# Requirements for auto_pipeline.R

## System

- Linux. The pipeline was written and run on Linux. Windows users can run it in WSL2 (Windows Subsystem for Linux) with a Linux distribution such as Ubuntu, and install R, the Bioconductor packages, and HTSlib inside it. It has not been tested on native Windows or macOS.
- R 4.4 or later. The script adds `/home/martiki/R/rocker-rstudio/4.4` to the library path on its first line. Change or remove that line if your packages are installed elsewhere.
- Bioconductor, installed through `BiocManager`. Use the Bioconductor release that matches your R version.
- `tabix` and `bgzip` (part of HTSlib) to compress and index the VCF if it is not already indexed.
- Enough disk space and memory for the genome and annotation packages below. They are large, and the gnomAD v4 data is the largest.

## R packages

The script loads these six packages directly.

| Package | Used for |
| --- | --- |
| `VariantAnnotation` | Reading the VCF, locating variants, and predicting coding effects. |
| `GenomicScores` | Looking up gnomAD allele frequencies (`gscores`). |
| `TxDb.Hsapiens.UCSC.hg38.knownGene` | Gene and transcript coordinates for hg38. |
| `BSgenome.Hsapiens.UCSC.hg38` | The hg38 reference sequence, used to predict amino acid changes. |
| `MafH5.gnomAD.v4.0.GRCh38` | gnomAD v4.0 allele frequencies. |
| `org.Hs.eg.db` | Converting gene symbols to Entrez IDs. |

These come in automatically as dependencies of the packages above, and the script uses them without loading them by name: `GenomicRanges`, `GenomicFeatures`, `Biostrings`, `Rsamtools`, `AnnotationDbi`, `IRanges`, and `S4Vectors`.

## Install

Run this once in R. It needs internet access.

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "VariantAnnotation",
  "GenomicScores",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "BSgenome.Hsapiens.UCSC.hg38",
  "MafH5.gnomAD.v4.0.GRCh38",
  "org.Hs.eg.db"
))
```

To install into the library the script expects, set the path first:

```r
.libPaths(c("/home/martiki/R/rocker-rstudio/4.4", .libPaths()))
```

## Check the install

This should print without errors:

```r
suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(GenomicScores)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(MafH5.gnomAD.v4.0.GRCh38)
  library(org.Hs.eg.db)
})
```

## Input data

| File | Notes |
| --- | --- |
| VCF | Bgzipped, with a `.tbi` index next to it. One sample, GRCh38. Fields used: GT, DP, GQ, AD, and optionally PID and PGT for phasing. |
| `phenotype_to_genes.txt` | From the HPO project downloads. Needs a term ID column named `hpo_id` and a gene column named `gene_symbol` or `entrez_gene_symbol`. |

## Notes

- Internet access is only needed to install packages. The first call to the gnomAD package may also need to download its data. After that, the script itself makes no network requests.
- Package versions are not pinned here. To record the versions you used, run `sessionInfo()` after a successful run and save the output with your results.
- If the `MafH5` package fails to install, check that the `rhdf5` dependency installed. It sometimes needs system HDF5 libraries.
