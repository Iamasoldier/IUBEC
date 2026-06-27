# IUBEC Usage

## 1. Dependencies

The core IUBEC correction logic mainly uses base R. The current ready-to-run script connects the corrected counts to `satuRn::fitDTU()`, so the following packages are required:

```r
library(SummarizedExperiment)
library(S4Vectors)
library(satuRn)
```

The `compositions` package is optional. If it is not installed, the code uses an internal ILR fallback.

## 2. Required input data

IUBEC requires three objects: `countData`, `tx2gene`, and `sampleData`.

### `countData`

`countData` is an isoform count matrix.

- Rows are isoforms.
- Columns are samples.
- `rownames(countData)` must be isoform IDs.
- `colnames(countData)` must be sample IDs.

Example:

```r
countData[1:3, 1:3]
#            sample_1 sample_2 sample_3
# isoform_1       100      120       95
# isoform_2        30       25       40
# isoform_3        10       12        8
```

### `tx2gene`

`tx2gene` is the isoform-to-gene mapping table.

- It must contain one isoform ID column.
- It must contain one gene ID column.
- Isoform IDs must cover `rownames(countData)`.

Example:

```r
tx2gene <- data.frame(
  TXNAME = c("isoform_1", "isoform_2", "isoform_3"),
  gene_id = c("gene_1", "gene_1", "gene_2")
)
```

### `sampleData`

`sampleData` is the sample information table.

- It must contain `group`. This is the biological group tested by the downstream DIU analysis.
- It can contain the covariate to be removed by IUBEC, such as `batch`.
- Sample IDs in `sampleData` must be alignable to `colnames(countData)`.

Example:

```r
sampleData <- data.frame(
  sample = c("sample_1", "sample_2", "sample_3", "sample_4"),
  group = c("control", "control", "case", "case"),
  batch = c("batch_1", "batch_2", "batch_1", "batch_2")
)
rownames(sampleData) <- sampleData$sample
```

It is recommended to set the order of `group` levels before running IUBEC:

```r
sampleData$group <- factor(sampleData$group, levels = c("control", "case"))
```

## 3. How to run `example_usage.R`

Put the following two files in the same folder:

```text
IUBEC.R
example_usage.R
```

Then run:

```bash
Rscript example_usage.R
```

No extra input files are required for testing `example_usage.R`. The script generates example `countData`, `tx2gene`, and `sampleData` internally.

After the script finishes, it writes example output files, such as corrected counts and the fitted `satuRn` object.

## 4. What data does the example generate?

`example_usage.R` automatically generates:

```text
countData   # simulated isoform count matrix
tx2gene     # simulated isoform-to-gene mapping table
sampleData  # simulated sample information with group and batch
```

It then calls:

```r
sm <- fitDTU_IUBEC(
  countData  = countData,
  tx2gene    = tx2gene,
  sampleData = sampleData,
  quiet      = FALSE
)
```

This example is only used to check whether the code can run successfully. It does not represent a real biological analysis result.

## 7. How IUBEC connects to the downstream `satuRn::fitDTU()` model

`IUBEC.R` first applies IUBEC correction to isoform usage proportions and reconstructs corrected counts. It then builds a `SummarizedExperiment` object and calls:

```r
sm <- satuRn::fitDTU(
  object  = sumExp,
  formula = ~ 0 + group
)
```

The `group` column is required.

If IUBEC has removed a covariate, such as `batch`, do not include the same covariate again in the downstream DTU model. That is, do not use:

```r
formula = ~ 0 + group + batch
```

Use:

```r
formula = ~ 0 + group
```

## 8. Notes

1. `countData` must have isoform IDs as row names and sample IDs as column names.
2. `sampleData` must be alignable to the columns of `countData`.
3. `sampleData` must contain `group`.
4. If IUBEC correction is used, `sampleData` should contain the covariate to remove, such as `batch`.
5. The covariate removed by IUBEC must have at least two levels. Two levels are preferred.
6. If `batch` and `group` are fully confounded, for example all case samples are in batch 1 and all control samples are in batch 2, IUBEC may remove true biological signals.
7. Before running IUBEC, it is recommended to check the relationship between `group` and the covariate:

```r
table(sampleData$group, sampleData$batch)
```

8. Corrected counts from IUBEC can be non-integers. This is expected because IUBEC first corrects proportions and then reconstructs counts using the original gene total counts. Rounding corrected counts is not recommended by default.
9. If a gene has only one isoform, or if the total count of a gene-sample pair is 0, the code keeps the original counts and does not force ILR correction.
