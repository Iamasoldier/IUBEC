############################################################
## example_usage.R
##
## Purpose:
##   A minimal runnable example for IUBEC_github.R.
##   This script creates a small toy isoform count matrix,
##   tx2gene mapping, and sampleData, then runs fitDTU_IUBEC().
##
## How to run:
##   Put this file in the same directory as IUBEC_github.R or IUBEC.R,
##   then run:
##     Rscript example_usage.R
############################################################

## ---------------- 1. Load required packages ----------------
required_pkgs <- c("SummarizedExperiment", "S4Vectors", "satuRn")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0L) {
  stop(
    "Missing required R packages: ", paste(missing_pkgs, collapse = ", "),
    "\nPlease install them before running this example. For example:\n",
    "  if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')\n",
    "  BiocManager::install(c('SummarizedExperiment', 'S4Vectors', 'satuRn'))\n",
    call. = FALSE
  )
}

library(SummarizedExperiment)
library(S4Vectors)
library(satuRn)

## ---------------- 2. Source IUBEC code ----------------
## Use IUBEC.R if you renamed the GitHub file.
## Otherwise use IUBEC_github.R.
if (file.exists("IUBEC.R")) {
  source("IUBEC.R")
} else if (file.exists("IUBEC_github.R")) {
  source("IUBEC_github.R")
} else {
  stop(
    "Cannot find IUBEC.R or IUBEC_github.R in the current working directory.\n",
    "Please put example_usage.R and IUBEC_github.R in the same folder, or rename IUBEC_github.R to IUBEC.R.",
    call. = FALSE
  )
}

## ---------------- 3. Create toy input data ----------------
## In real use, replace this section with read.csv() commands.
## The three required objects are:
##   countData  : transcripts x samples count matrix
##   tx2gene    : transcript-to-gene mapping table
##   sampleData : sample annotation table with group and optional batch

set.seed(123)

## Sample annotation.
## Important: group must exist.
## batch is optional, but if provided for IUBEC correction it should have at least two levels.
sampleData <- data.frame(
  sample = paste0("sample_", 1:8),
  group  = rep(c("control", "case"), each = 4),
  batch  = rep(c("batch_1", "batch_2"), times = 4),
  stringsAsFactors = FALSE
)
rownames(sampleData) <- sampleData$sample
sampleData$group <- factor(sampleData$group, levels = c("control", "case"))
sampleData$batch <- factor(sampleData$batch)

## Transcript-to-gene mapping.
n_genes <- 12
isoforms_per_gene <- 3
n_tx <- n_genes * isoforms_per_gene

tx2gene <- data.frame(
  TXNAME  = paste0("tx_", seq_len(n_tx)),
  gene_id = rep(paste0("gene_", seq_len(n_genes)), each = isoforms_per_gene),
  stringsAsFactors = FALSE
)

## Simulated isoform count matrix.
countData <- matrix(
  0,
  nrow = n_tx,
  ncol = nrow(sampleData),
  dimnames = list(tx2gene$TXNAME, rownames(sampleData))
)

for (g in seq_len(n_genes)) {
  rows <- which(tx2gene$gene_id == paste0("gene_", g))

  for (j in seq_len(nrow(sampleData))) {
    gene_total <- rnbinom(1, mu = 200 + 10 * g, size = 20)

    ## Baseline within-gene isoform usage proportions.
    prob <- c(0.55, 0.30, 0.15)

    ## Add a simple batch-related isoform usage shift.
    if (sampleData$batch[j] == "batch_2") {
      prob <- prob + c(-0.10, 0.10, 0.00)
    }

    ## Add a small condition-related isoform usage shift for the first 3 genes.
    if (sampleData$group[j] == "case" && g <= 3) {
      prob <- prob + c(-0.12, 0.12, 0.00)
    }

    prob <- pmax(prob, 0.001)
    prob <- prob / sum(prob)

    countData[rows, j] <- as.vector(rmultinom(1, size = gene_total, prob = prob))
  }
}

## Optional checks before running IUBEC.
stopifnot(!is.null(rownames(countData)))
stopifnot(!is.null(colnames(countData)))
stopifnot(all(colnames(countData) %in% rownames(sampleData)))
stopifnot(all(rownames(countData) %in% tx2gene$TXNAME))
stopifnot("group" %in% colnames(sampleData))

cat("Input check passed.\n")
cat("countData dimension:", paste(dim(countData), collapse = " x "), "\n")
cat("Number of genes:", length(unique(tx2gene$gene_id)), "\n")
cat("Group x batch table:\n")
print(table(sampleData$group, sampleData$batch))

## ---------------- 4. Run IUBEC + satuRn fitting ----------------
sm <- fitDTU_IUBEC(
  countData  = countData,
  tx2gene    = tx2gene,
  sampleData = sampleData,
  quiet      = FALSE
)

cat("IUBEC + satuRn fitting finished.\n")

## ---------------- 5. Export corrected counts ----------------
## fitDTU_IUBEC() returns a SummarizedExperiment-like object from satuRn.
## The corrected counts are stored in the counts assay.
corrected_counts <- SummarizedExperiment::assay(sm, "counts")

write.csv(corrected_counts, "IUBEC_example_corrected_counts.csv")
saveRDS(sm, "IUBEC_example_fitDTU_result.rds")

cat("Output files written:\n")
cat("  IUBEC_example_corrected_counts.csv\n")
cat("  IUBEC_example_fitDTU_result.rds\n")

cat("First rows of corrected counts:\n")
print(round(corrected_counts[1:6, 1:4], 3))

## ---------------- 6. Real data template ----------------
## For real data, comment out the toy-data section above and use this template:
##
## countData <- read.csv("counts.csv", row.names = 1, check.names = FALSE)
## tx2gene <- read.csv("tx2gene.csv", stringsAsFactors = FALSE)
## sampleData <- read.csv("sampleData.csv", stringsAsFactors = FALSE)
## rownames(sampleData) <- sampleData$sample
## sampleData$group <- factor(sampleData$group, levels = c("control", "case"))
## sampleData$batch <- factor(sampleData$batch)
##
## sm <- fitDTU_IUBEC(
##   countData = countData,
##   tx2gene = tx2gene,
##   sampleData = sampleData,
##   quiet = FALSE
## )
