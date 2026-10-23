#!/usr/bin/env Rscript
# =====================================================================
# AUTOMATED RARE-DISEASE PIPELINE -> challenge submission
#
#   HPO terms + VCF  ->  ranked findings  ->  submission.tsv
#
# No manual curation. Deterministic. Runs offline start to finish.
#
#   Rscript auto_pipeline.R --vcf sample.vcf.gz --hpo phenotype_to_genes.txt
#
# Stages
#   1  build a phenotype-weighted gene panel from HPO annotations
#   2  read only those genes from the VCF (fast, avoids a genome pass)
#   3  genotype and quality filters
#   4  annotate: region, coding consequence, gnomAD v4 frequency
#   5  assemble recessive and dominant findings, phased where possible
#   6  score each finding, convert to a 0-1 epcr, rank, write the table
#
# All data stays on this machine.
# =====================================================================

# Optional: uncomment and edit if your R packages live in a personal library.
# .libPaths(c("/path/to/your/R/package/library", .libPaths()))

suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(GenomicScores)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(MafH5.gnomAD.v4.0.GRCh38)
  library(org.Hs.eg.db)
})

# ---------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------
PROBAND_ID  <- "PROBAND01"
# The example rows in the spec show a constant sample ID across loci, but if
# the challenge actually wants the gene symbol in this column, set TRUE.
ID_IS_GENE  <- TRUE

# File paths come from the command line:
#   --vcf    bgzipped, tabix-indexed, single-sample GRCh38 VCF (the .tbi
#            index must sit next to it)                         [required]
#   --hpo    HPO file mapping phenotype terms to genes
#            (phenotype_to_genes.txt)                           [required]
#   --out    ranked findings table to submit  [default: submission.tsv]
#   --audit  every annotated variant that passed QC; contains patient
#            genotypes, do not share          [default: audit_variants.tsv]
USAGE <- paste(
  "Usage: Rscript auto_pipeline.R --vcf <sample.vcf.gz> --hpo <phenotype_to_genes.txt>",
  "                               [--out <submission.tsv>] [--audit <audit_variants.tsv>]",
  sep = "\n")

parse_args <- function(argv) {
  opts <- list(vcf = NULL, hpo = NULL,
               out = "submission.tsv", audit = "audit_variants.tsv")
  if (any(argv %in% c("-h", "--help"))) { cat(USAGE, "\n"); quit(status = 0) }
  i <- 1L
  while (i <= length(argv)) {
    a <- argv[i]
    if (!grepl("^--", a)) stop("unexpected argument: ", a, "\n", USAGE, call. = FALSE)
    key <- sub("^--", "", sub("=.*$", "", a))
    if (!key %in% names(opts)) stop("unknown option: --", key, "\n", USAGE, call. = FALSE)
    if (grepl("=", a, fixed = TRUE)) {
      val <- sub("^[^=]*=", "", a)
    } else {
      i <- i + 1L
      if (i > length(argv)) stop("missing value for --", key, "\n", USAGE, call. = FALSE)
      val <- argv[i]
    }
    opts[[key]] <- val
    i <- i + 1L
  }
  opts
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(args$vcf) || is.null(args$hpo))
  stop("--vcf and --hpo are required\n", USAGE, call. = FALSE)
if (!file.exists(args$vcf)) stop("VCF not found: ", args$vcf, call. = FALSE)
if (!file.exists(paste0(args$vcf, ".tbi")))
  stop("tabix index not found: ", args$vcf, ".tbi", call. = FALSE)
if (!file.exists(args$hpo)) stop("HPO file not found: ", args$hpo, call. = FALSE)

VCF_PATH  <- args$vcf
HPO_ANNOT <- args$hpo
OUT_TSV   <- args$out
AUDIT_TSV <- args$audit

HPO_TERMS <- c("HP:0002859", "HP:0000121", "HP:0004322", "HP:0001508",
               "HP:0003202", "HP:0001622", "HP:0001518", "HP:0200067")  # the proband's terms

# Genes always searched regardless of HPO ranking, for a suspected disorder.
# Leave as character(0) for a purely phenotype-driven run.
SEED_GENES <- c("BUB1B", "CEP57", "TRIP13", "BUB1", "BUB3", "CEP192")

N_PANEL      <- 500L      # HPO-ranked genes to search
PAD          <- 5000L     # bp around each gene (promoter, UTR, splice)
MAX_FINDINGS <- 10L       # rows in the submission
N_PRIMARY    <- 1L        # top N labelled "primary", rest "secondary"

# Scoring rewards ranking only, and secondary findings carry no penalty, so
# these are deliberately permissive: fill all 10 slots and let the ranking do
# the work rather than filtering down to a handful.
MIN_DP <- 10; MIN_GQ <- 20; MIN_AB <- 0.20
MAX_AF_DOM <- 0.005       # dominant / de novo models
MAX_AF_REC <- 0.05        # recessive: carrier alleles can be commoner

txdb  <- TxDb.Hsapiens.UCSC.hg38.knownGene
bsg   <- BSgenome.Hsapiens.UCSC.hg38
mafdb <- MafH5.gnomAD.v4.0.GRCh38

LOF <- c("nonsense", "frameshift")
MISSENSE <- "nonsynonymous"

# ---------------------------------------------------------------------
# Seqlevel handling. Done explicitly: seqlevelsStyle() fails silently
# when alt contigs are present, which yields zero overlaps and no error.
# ---------------------------------------------------------------------
to_ucsc <- function(gr) {
  sl <- seqlevels(gr)
  new <- ifelse(grepl("^chr", sl), sl, paste0("chr", sub("^MT$", "M", sl)))
  if (anyDuplicated(new) == 0) seqlevels(gr) <- new
  gr
}
to_vcf <- function(gr) {
  sl <- seqlevels(gr)
  new <- ifelse(grepl("^chr", sl), sub("^chr", "", sl), sl)
  if (anyDuplicated(new) == 0) seqlevels(gr) <- new
  gr
}

# ---------------------------------------------------------------------
# STAGE 1 -- phenotype-weighted gene panel
# ---------------------------------------------------------------------
message("[1/6] building HPO gene panel")

hpo_weight <- setNames(rep(0, 0), character(0))
panel <- toupper(SEED_GENES)

if (file.exists(HPO_ANNOT) && length(HPO_TERMS)) {
  ann <- read.delim(HPO_ANNOT, stringsAsFactors = FALSE, comment.char = "")
  names(ann) <- tolower(sub("^#", "", names(ann)))
  gcol <- if ("gene_symbol" %in% names(ann)) "gene_symbol" else "entrez_gene_symbol"
  ann$sym <- toupper(ann[[gcol]])

  total <- length(unique(ann$sym))
  # Specificity weight: a term annotated to 20 genes carries far more signal
  # than one annotated to 3000. Log-scaled so the rarest do not dominate.
  w <- sapply(HPO_TERMS, function(t) {
    n <- length(unique(ann$sym[ann$hpo_id == t]))
    if (n > 0) log(total / n) else 0
  })
  sub_ann <- ann[ann$hpo_id %in% HPO_TERMS, ]
  if (nrow(sub_ann)) {
    sc <- tapply(sub_ann$hpo_id, sub_ann$sym,
                 function(x) sum(w[unique(x)], na.rm = TRUE))
    sc <- sort(sc, decreasing = TRUE)
    hpo_weight <- sc / max(sc)                  # normalised 0-1
    panel <- unique(c(panel, names(sc)[seq_len(min(N_PANEL, length(sc)))]))
  }
} else {
  warning("HPO annotations not found -- running on SEED_GENES only")
}
message("  panel: ", length(panel), " genes")
if (!length(panel)) stop("empty gene panel")

# ---------------------------------------------------------------------
# STAGE 2 -- read those genes from the VCF
# ---------------------------------------------------------------------
message("[2/6] reading VCF over panel regions")

eg <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db, keys = panel, keytype = "SYMBOL", columns = "ENTREZID"))
eg <- eg[!is.na(eg$ENTREZID), ]

gg <- suppressWarnings(genes(txdb, single.strand.genes.only = TRUE))
gr <- gg[names(gg) %in% eg$ENTREZID]
gr$symbol <- eg$SYMBOL[match(names(gr), eg$ENTREZID)]
gr <- keepSeqlevels(gr, seqlevelsInUse(gr), pruning.mode = "coarse")
gr <- gr + PAD
message("  ", length(gr), " genes with coordinates")

v <- readVcf(TabixFile(VCF_PATH), "hg38", param = ScanVcfParam(which = to_vcf(gr)))
message("  ", nrow(v), " records")
if (nrow(v) == 0) stop("zero records -- check contig naming against the VCF")

v  <- expand(v)                 # one row per ALT allele
rr <- rowRanges(v)

# ---------------------------------------------------------------------
# STAGE 3 -- genotype and quality filters
# ---------------------------------------------------------------------
message("[3/6] quality filtering")

gt <- as.character(geno(v)$GT[, 1])
carried <- !(gt %in% c("0/0", "0|0", "./.", ".|."))

dp <- if (!is.null(geno(v)$DP)) suppressWarnings(as.numeric(geno(v)$DP[, 1])) else rep(Inf, nrow(v))
gq <- if (!is.null(geno(v)$GQ)) suppressWarnings(as.numeric(geno(v)$GQ[, 1])) else rep(Inf, nrow(v))
dp[is.na(dp)] <- 0; gq[is.na(gq)] <- 0

ab <- rep(NA_real_, nrow(v))
ad <- geno(v)$AD
if (!is.null(ad) && length(dim(ad)) == 3L) {
  tot <- ad[, 1, 1] + ad[, 1, 2]
  ab <- ifelse(!is.na(tot) & tot > 0, ad[, 1, 2] / tot, NA_real_)
}

zyg <- ifelse(gt %in% c("1/1", "1|1"), "hom",
       ifelse(gt %in% c("0/1", "1/0", "0|1", "1|0"), "het", "other"))

filt_ok <- lengths(filt(v)) == 0 | filt(v) %in% c("PASS", ".", "")
qual_ok <- carried & filt_ok & dp >= MIN_DP & gq >= MIN_GQ &
           !(zyg == "het" & !is.na(ab) & ab < MIN_AB)
message("  ", sum(qual_ok), " carried variants pass QC")
if (!any(qual_ok)) stop("no variants survived QC")

# ---------------------------------------------------------------------
# STAGE 4 -- annotation
# ---------------------------------------------------------------------
message("[4/6] annotating")

idx <- which(qual_ok)
rr_k <- rr[idx]

hitidx <- findOverlaps(to_ucsc(rr_k), gr, select = "first")
symbol <- gr$symbol[hitidx]

rr_u <- to_ucsc(rr_k)
suppressWarnings(seqlevels(rr_u, pruning.mode = "coarse") <-
                   intersect(seqlevels(rr_u), seqlevels(txdb)))

location <- rep(NA_character_, length(rr_k))
consequence <- rep(NA_character_, length(rr_k))
aa <- rep(NA_character_, length(rr_k))

if (length(rr_u)) {
  loc <- suppressWarnings(locateVariants(rr_u, txdb, AllVariants()))
  lb <- loc[!duplicated(loc$QUERYID), ]
  location[lb$QUERYID] <- as.character(lb$LOCATION)

  cod <- suppressWarnings(suppressMessages(tryCatch(
    predictCoding(rr_u, txdb, seqSource = bsg, varAllele = DNAStringSet(rr_u$ALT)),
    error = function(e) { message("  predictCoding: ", conditionMessage(e)); NULL })))
  if (!is.null(cod) && length(cod)) {
    cb <- cod[!duplicated(cod$QUERYID), ]
    consequence[cb$QUERYID] <- as.character(cb$CONSEQUENCE)
    # PROTEINLOC is an IntegerList: an indel spanning a codon boundary has
    # more than one position, so take the first rather than coercing.
    ploc <- vapply(cb$PROTEINLOC,
                   function(z) if (length(z)) as.character(z[1]) else "",
                   character(1))
    aa[cb$QUERYID] <- paste0(as.character(cb$REFAA), ploc,
                             as.character(cb$VARAA))
  }
}

# gnomAD v4. Query must be clean width-1 UCSC ranges; guard the length so a
# short or empty return becomes NA rather than breaking the data frame.
q <- GRanges(paste0("chr", sub("^chr", "", as.character(seqnames(rr_k)))),
             IRanges(start(rr_k), width = 1))
af <- tryCatch({
  gs <- gscores(mafdb, q, pop = "AF")
  # score() looks for a column named "score"/"default". With pop= set, the
  # column is named after the population ("AF"), so score() returns nothing
  # and every frequency silently becomes NA. Read mcols directly instead.
  col <- if ("AF" %in% names(mcols(gs))) "AF" else names(mcols(gs))[1]
  s <- as.numeric(mcols(gs)[[col]])
  if (length(s) == length(rr_k)) s else rep(NA_real_, length(rr_k))
}, error = function(e) {
  message("  gscores failed (", conditionMessage(e), ") -- frequencies unavailable")
  rep(NA_real_, length(rr_k))
})
message("  gnomAD frequencies for ", sum(!is.na(af)), "/", length(af), " variants")

getg <- function(nm) {
  if (!is.null(geno(v)[[nm]])) as.character(geno(v)[[nm]][idx, 1])
  else rep(NA_character_, length(idx))
}

d <- data.frame(
  gene = symbol,
  chrom = as.character(seqnames(rr_k)), pos = start(rr_k),
  ref = as.character(rr_k$REF), alt = as.character(rr_k$ALT),
  zyg = zyg[idx], location = location, consequence = consequence, aa = aa,
  af = af, dp = dp[idx], gq = gq[idx], ab = round(ab[idx], 3),
  pid = getg("PID"), pgt = getg("PGT"),
  stringsAsFactors = FALSE)
d <- d[!is.na(d$gene), ]

d$is_lof <- (!is.na(d$consequence) & d$consequence %in% LOF) |
            (!is.na(d$location) & d$location == "spliceSite")
d$is_mis <- !is.na(d$consequence) & d$consequence == MISSENSE
# predictCoding is brittle and can return nothing, which would leave both
# flags FALSE and produce an empty findings list. locateVariants is far more
# robust, so its region call also counts. Less precise -- "coding" does not
# separate missense from synonymous -- but with only 10 submission slots and
# no penalty for weak guesses, over-generating and ranking beats filtering
# down to zero.
d$impactful <- d$is_lof | d$is_mis |
               (!is.na(d$location) & d$location %in% c("coding", "spliceSite"))

# NA frequency means absent from gnomAD, which is the interesting case.
d$af_eff <- ifelse(is.na(d$af), 0, d$af)

write.table(d, AUDIT_TSV, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
Sys.chmod(AUDIT_TSV, "600")

cand <- d[d$impactful, ]
message("  ", nrow(cand), " rare/protein-altering candidates")

# ---------------------------------------------------------------------
# STAGE 5 -- assemble findings
# ---------------------------------------------------------------------
message("[5/6] assembling findings")

gw <- function(g) if (g %in% names(hpo_weight)) as.numeric(hpo_weight[[g]]) else 0
seeded <- function(g) g %in% toupper(SEED_GENES)

# Gene-level evidence. Deliberately kept OUT of the phase adjustment: how
# well a gene matches the phenotype does not become less true because two
# variants in it could not be phased.
gene_term <- function(g) 2.2 * gw(g) + 0.8 * seeded(g)

# Severity falls back to the locateVariants region call when predictCoding
# returns nothing, so the term still discriminates instead of flattening to
# a single constant for every variant.
sev <- function(r) {
  if (isTRUE(r$is_lof)) 1.0
  else if (isTRUE(r$is_mis)) 0.6
  else if (!is.na(r$location) && r$location == "spliceSite") 0.9
  else if (!is.na(r$location) && r$location == "coding") 0.55
  else 0.3
}
rarity <- function(a) if (a == 0) 1.0 else if (a < 1e-4) 0.8 else
                      if (a < 1e-3) 0.55 else 0.3

# Allele evidence, scaled so that the two recessive models are directly
# comparable: a homozygote and a compound heterozygote with the same
# severity and rarity on both alleles score identically. Previously the
# homozygous coefficient (3.0 on one variant) competed against 1.6+1.6 on
# two, and the whole compound-het sum was then multiplied by 0.75.
allele_term <- function(rs) {
  sevs <- vapply(rs, sev, numeric(1))
  rars <- vapply(rs, function(r) rarity(r$af_eff), numeric(1))
  3.0 * mean(sevs) + 1.6 * mean(rars)
}

# Phase is an additive adjustment, not a multiplier. Two rare variants in
# one gene are more often in trans than in cis, since inheriting both from
# one parent requires them to already share a haplotype -- so the penalty
# for "unknown" is modest.
PHASE_TRANS    <-  0.35
PHASE_UNPHASED <- -0.25

findings <- list()
addf <- function(v1, v2, score, type, notes) {
  findings[[length(findings) + 1L]] <<- list(v1 = v1, v2 = v2,
                                             score = score, model = type,
                                             notes = notes)
}

for (g in unique(cand$gene)) {
  sub <- cand[cand$gene == g, ]
  sub <- sub[order(sub$pos), ]
  wg <- gw(g)

  # --- homozygous: recessive model satisfied by one variant ----------
  hom <- sub[sub$zyg == "hom" & sub$af_eff <= MAX_AF_REC, ]
  for (i in seq_len(nrow(hom))) {
    r <- hom[i, ]
    s <- allele_term(list(r)) + gene_term(g)
    addf(r, NULL, s, "homozygous",
         sprintf("%s homozygous %s%s; gnomAD AF %s; DP=%g GQ=%g AB=%s",
                 g, ifelse(is.na(r$consequence), r$location, r$consequence),
                 ifelse(is.na(r$aa) || !nzchar(r$aa), "", paste0(" ", r$aa)),
                 ifelse(is.na(r$af), "absent", format(r$af, digits = 3)),
                 r$dp, r$gq, ifelse(is.na(r$ab), "NA", r$ab)))
  }

  # --- compound heterozygote ----------------------------------------
  het <- sub[sub$zyg == "het" & sub$af_eff <= MAX_AF_REC, ]
  if (nrow(het) >= 2) {
    for (a in 1:(nrow(het) - 1)) for (b in (a + 1):nrow(het)) {
      r1 <- het[a, ]; r2 <- het[b, ]
      phase <- "phase unresolved (no shared GATK phase group); parental testing required"
      padj <- PHASE_UNPHASED
      if (!is.na(r1$pid) && !is.na(r2$pid) && r1$pid == r2$pid) {
        if (r1$pgt != r2$pgt) {
          phase <- "in trans, confirmed by GATK physical phasing"
          padj <- PHASE_TRANS
        } else next   # in cis: other allele intact, not a compound het
      }
      s <- allele_term(list(r1, r2)) + gene_term(g) + padj
      addf(r1, r2, s, "compound_het",
           sprintf("%s compound heterozygote; %s / %s; gnomAD %s / %s; %s",
                   g,
                   ifelse(is.na(r1$consequence), r1$location, r1$consequence),
                   ifelse(is.na(r2$consequence), r2$location, r2$consequence),
                   ifelse(is.na(r1$af), "absent", format(r1$af, digits = 3)),
                   ifelse(is.na(r2$af), "absent", format(r2$af, digits = 3)),
                   phase))
    }
  }

  # --- dominant / de novo candidate: single very rare LoF -------------
  dom <- sub[sub$zyg == "het" & sub$is_lof & sub$af_eff <= MAX_AF_DOM, ]
  for (i in seq_len(nrow(dom))) {
    r <- dom[i, ]
    s <- 0.55 * allele_term(list(r)) + gene_term(g) - 0.5
    addf(r, NULL, s, "heterozygous_lof",
         sprintf(paste("%s heterozygous %s; gnomAD %s;",
                       "dominant or unresolved recessive - second allele may be",
                       "deep-intronic, splice-region or a CNV, none of which are",
                       "visible in a small-variant VCF"),
                 g, ifelse(is.na(r$consequence), r$location, r$consequence),
                 ifelse(is.na(r$af), "absent", format(r$af, digits = 3))))
  }
}

if (!length(findings)) {
  message("  no findings assembled")
}

# ---------------------------------------------------------------------
# STAGE 6 -- score, rank, write
# ---------------------------------------------------------------------
message("[6/6] scoring and writing submission")

COLS <- c("proband_id", "chrom_1", "pos_1", "ref_1", "alt_1",
          "chrom_2", "pos_2", "ref_2", "alt_2",
          "epcr", "finding_type", "notes")

chrfmt <- function(x) paste0("chr", sub("^chr", "", as.character(x)))

out <- data.frame(matrix(ncol = length(COLS), nrow = 0), stringsAsFactors = FALSE)
names(out) <- COLS

if (length(findings)) {
  sc <- vapply(findings, function(f) f$score, numeric(1))
  # Logistic squash to a bounded probability. Centre and slope chosen so a
  # homozygous LoF in a top-ranked phenotype gene lands near 0.9 and a lone
  # heterozygote near 0.25. Calibrated by construction, not from data.
  epcr <- 1 / (1 + exp(-(sc - 4.2) / 1.3))
  epcr <- pmin(0.95, pmax(0.02, epcr))

  ord <- order(sc, decreasing = TRUE)
  ord <- ord[seq_len(min(MAX_FINDINGS, length(ord)))]

  rows <- lapply(seq_along(ord), function(k) {
    f <- findings[[ord[k]]]
    v1 <- f$v1; v2 <- f$v2
    data.frame(
      proband_id = if (ID_IS_GENE) as.character(v1$gene) else PROBAND_ID,
      chrom_1 = chrfmt(v1$chrom), pos_1 = as.character(v1$pos),
      ref_1 = v1$ref, alt_1 = v1$alt,
      chrom_2 = if (is.null(v2)) "" else chrfmt(v2$chrom),
      pos_2   = if (is.null(v2)) "" else as.character(v2$pos),
      ref_2   = if (is.null(v2)) "" else v2$ref,
      alt_2   = if (is.null(v2)) "" else v2$alt,
      epcr = round(epcr[ord[k]], 3),
      finding_type = if (k <= N_PRIMARY) "primary" else "secondary",
      notes = f$notes, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)[, COLS]
}

write.table(out, OUT_TSV, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
Sys.chmod(OUT_TSV, "600")

cat("\n=== submission ===\n")
if (nrow(out)) {
  print(out[, c("chrom_1", "pos_1", "chrom_2", "pos_2", "epcr", "finding_type")],
        row.names = FALSE)
} else cat("no findings\n")
cat("\n", nrow(out), " rows -> ", OUT_TSV, "\n", sep = "")
cat("audit trail (all annotated variants) -> ", AUDIT_TSV, "\n", sep = "")

invisible(out)
