#' fitDTU_IUBEC: IUBEC_satuRn fitting entry point
#'
#' This function applies IUBEC correction to isoform-level count data and then
#' fits satuRn on the corrected counts. It keeps the original validation and
#' alignment checks from the working script, but removes code that is only used
#' for internal evaluation, logging, or ablation experiments.
#'
#' IUBEC correction follows four steps:
#'   1) within-gene isoform usage proportion construction;
#'   2) ILR transformation;
#'   3) batch-related direction estimation and projection removal in ILR space;
#'   4) inverse ILR transformation and corrected count reconstruction.
#'
#' @param countData matrix/data.frame. Rows are isoforms or transcripts, and columns are samples.
#' @param tx2gene data.frame. It must contain one transcript ID column and one gene ID column.
#' @param sampleData data.frame. It must contain group or condition. It may contain batch.
#' @param quiet logical. If TRUE, messages are suppressed where possible.
#' @param tx_col optional character. Column name for transcript IDs in tx2gene.
#' @param gene_col optional character. Column name for gene IDs in tx2gene.
#' @return A satuRn StatModel object.
#' @export
fitDTU_IUBEC <- function(countData,
                         tx2gene,
                         sampleData,
                         quiet = TRUE,
                         tx_col = NULL,
                         gene_col = NULL
) {
  ## ----------------- 基本检查 -----------------
  countData <- as.matrix(countData)
  if (is.null(rownames(countData))) {
    stop("fitDTU_IUBEC: countData 必须有转录本 ID 作为行名(rownames)。")
  }
  if (is.null(colnames(countData))) {
    stop("fitDTU_IUBEC: countData 必须有样本 ID 作为列名(colnames)。")
  }
  if (!is.data.frame(sampleData)) {
    stop("fitDTU_IUBEC: sampleData 应该是 data.frame，当前类型为 ",
         paste(class(sampleData), collapse = "/"))
  }

  ## ----------------- sampleData 与 countData 显式对齐（修复：不再按位置） -----------------
  align_sampleData <- function(sampleData, countData) {
    cn <- colnames(countData)

    ## 优先：rownames(sampleData) 直接覆盖 colnames(countData)
    if (!is.null(rownames(sampleData)) && all(cn %in% rownames(sampleData))) {
      sampleData2 <- sampleData[cn, , drop = FALSE]
      rownames(sampleData2) <- cn
      return(sampleData2)
    }

    ## 次优：sample_id / sample / sampleID 等列
    id_candidates <- c("sample_id", "sample", "sampleID", "SampleID", "Sample", "id", "ID")
    id_col <- intersect(id_candidates, colnames(sampleData))[1]
    if (!is.na(id_col)) {
      ids <- as.character(sampleData[[id_col]])
      m <- match(cn, ids)
      if (any(is.na(m))) {
        missing <- cn[is.na(m)]
        stop("fitDTU_IUBEC: sampleData 无法与 countData 对齐：以下样本在 sampleData 的列 ", id_col,
             " 中找不到：", paste(missing, collapse = ", "),
             "\n请确保 sampleData 的 rownames 或 sample_id 列能覆盖 countData 的所有列名。")
      }
      sampleData2 <- sampleData[m, , drop = FALSE]
      rownames(sampleData2) <- cn
      return(sampleData2)
    }

    stop("fitDTU_IUBEC: 无法对齐 sampleData 与 countData。\n",
         "需要满足其一：\n",
         "  (a) rownames(sampleData) 覆盖并包含 countData 的 colnames；或\n",
         "  (b) sampleData 存在 sample_id/sample/sampleID 等列用于匹配。\n")
  }

  sampleData <- align_sampleData(sampleData, countData)

  ## ----------------- group / case 标签 -----------------
  if (!("group" %in% colnames(sampleData))) {
    if ("condition" %in% colnames(sampleData)) {
      sampleData$group <- as.character(sampleData$condition)
    } else {
      stop("fitDTU_IUBEC: sampleData 需要包含 group 或 condition 列。")
    }
  }
  grp <- factor(sampleData$group)

  ## ----------------- batch 方向估计使用所有样本 -----------------
  batch_fit_idx <- seq_len(ncol(countData))
  ## ----------------- batch（多水平不二值化） -----------------
  have_batch <- ("batch" %in% colnames(sampleData))
  batch <- NULL
  y_batch01 <- NULL
  batch_levels <- NULL
  n_batch_levels <- 0L

  if (have_batch) {
    batch <- factor(sampleData$batch)
    n_batch_levels <- nlevels(batch)
    batch_levels <- levels(batch)
    if (n_batch_levels < 2L) {
      have_batch <- FALSE
    } else if (n_batch_levels == 2L) {
      y_batch01 <- as.integer(batch != batch_levels[1])
    }
  }

  ## ----------------- tx2gene：严格识别 tx/gene 列（修复：不再前两列回退） -----------------
  tx_cols_all <- colnames(tx2gene)
  tx_col_candidates   <- c("TXNAME", "tx_id", "txname", "transcript_id", "isoform_id", "feature_id", "feature")
  gene_col_candidates <- c("gene_id", "GENEID", "gene", "GENE", "geneID")

  if (is.null(tx_col)) {
    tx_col <- intersect(tx_col_candidates, tx_cols_all)[1]
  }
  if (is.null(gene_col)) {
    gene_col <- intersect(gene_col_candidates, tx_cols_all)[1]
  }

  
if (is.na(tx_col) || is.na(gene_col) || is.null(tx_col) || is.null(gene_col)) {
  ## 与旧代码保持“可回退”的行为，但增加两个硬约束：
  ## 1) 必须能识别出一对 tx/gene 列；2) 该 tx 列必须覆盖全部 rownames(countData)
  if (length(tx_cols_all) >= 2L) {
    tx_col_try <- tx_cols_all[1]
    gene_col_try <- tx_cols_all[2]
    tx_ids_try <- as.character(tx2gene[[tx_col_try]])
    if (all(rownames(countData) %in% tx_ids_try)) {
      tx_col <- tx_col_try
      gene_col <- gene_col_try
      if (!quiet) {
        message("fitDTU_IUBEC: tx2gene 列名未匹配候选列表，安全回退使用前两列: ",
                tx_col, " / ", gene_col)
      }
    } else {
      stop("fitDTU_IUBEC: 无法在 tx2gene 中自动识别 tx/gene 列，且前两列也无法覆盖全部转录本。",
     "\n当前 tx2gene 列名为：", paste(tx_cols_all, collapse = ", "),
     "\n请通过 tx_col / gene_col 显式指定正确列名。")
    }
  } else {
    stop("fitDTU_IUBEC: tx2gene 列数不足 2，无法识别 tx/gene 映射列。请提供规范 tx2gene 或显式指定列名。")
  }
}
  if (!(tx_col %in% tx_cols_all) || !(gene_col %in% tx_cols_all)) {
    stop("fitDTU_IUBEC: tx_col 或 gene_col 不存在于 tx2gene 列名中。")
  }

  tx_ids_all   <- as.character(tx2gene[[tx_col]])
  gene_ids_all <- as.character(tx2gene[[gene_col]])

  idx <- match(rownames(countData), tx_ids_all)
  if (any(is.na(idx))) {
    missing <- rownames(countData)[is.na(idx)]
    stop("fitDTU_IUBEC: 有 rownames(countData) 在 tx2gene 的转录本列中找不到。示例：",
         paste(head(missing, 10), collapse = ", "),
         if (length(missing) > 10) " ..." else "")
  }

  TXNAME  <- tx_ids_all[idx]
  gene_id <- gene_ids_all[idx]

  tx2gene_model <- data.frame(
    TXNAME  = as.character(TXNAME),
    gene_id = as.character(gene_id),
    stringsAsFactors = FALSE
  )
  rownames(tx2gene_model) <- tx2gene_model$TXNAME

  ## 基因内索引
  tx_index <- split(seq_len(nrow(tx2gene_model)), tx2gene_model$gene_id)

  ## ----------------- ILR / 逆 ILR 工具（修复 fallback） -----------------
  have_comp <- requireNamespace("compositions", quietly = TRUE)

  ## 标准 Helmert ILR basis：K x (K-1)，列正交，且与 1 向量正交
  helmert_ilr_basis <- function(K) {
    if (K < 2L) stop("ILR basis requires K>=2")
    V <- matrix(0, nrow = K, ncol = K - 1)
    for (j in 1:(K - 1)) {
      denom <- sqrt(j * (j + 1))
      V[1:j, j] <- 1 / denom
      V[j + 1, j] <- -j / denom
    }
    ## 安全检查（极小开销）：确保 V^T 1 ≈ 0
    if (max(abs(colSums(V))) > 1e-8) {
      stop("fitDTU_IUBEC: ILR fallback basis check failed: max|1^T V| > 1e-8")
    }
    V
  }

  ilr_vec <- function(p) {
    K <- length(p)
    if (have_comp) {
      as.numeric(compositions::ilr(compositions::acomp(p)))
    } else {
      V <- helmert_ilr_basis(K)          # K x (K-1)
      as.numeric(log(p) %*% V)           # 1 x (K-1)
    }
  }

  ilr_inv_vec <- function(z) {
    K1 <- length(z); K <- K1 + 1
    if (have_comp) {
      as.numeric(compositions::ilrInv(z))
    } else {
      V <- helmert_ilr_basis(K)          # K x (K-1)
      logp <- as.numeric(V %*% z)        # K x 1
      p <- exp(logp)
      p / sum(p)
    }
  }

  ## ----------------- 每个基因做 IUBEC + 逆变换 -----------------
  corrected_counts <- matrix(0, nrow = nrow(countData), ncol = ncol(countData),
                             dimnames = dimnames(countData))

  ## [v3-patch] eps 仅保留给数值稳定（如投影范数/求解时的安全下界），
  ## 不再直接作为所有 gene/sample 共用的固定伪计数。
  eps <- 1e-8

  ## 预先准备多水平 batch 的 one-hot（仅在需要时）
  ## batch 方向估计默认使用所有样本。
  if (have_batch && n_batch_levels > 2L) {
    Yb_full_all <- stats::model.matrix(~ 0 + batch)  # n x L
    colnames(Yb_full_all) <- batch_levels
  } else {
    Yb_full_all <- NULL
  }
  for (gid in names(tx_index)) {
    rows <- tx_index[[gid]]
    K <- length(rows)

    if (K < 2L) {
      corrected_counts[rows, ] <- countData[rows, , drop = FALSE]
      next
    }

    Xg <- countData[rows, , drop = FALSE]          # K x n
    Xg_sum0 <- colSums(Xg)                         # n
    zero_cols <- (Xg_sum0 <= 0)

    ## [v3-patch-4] 对极稀疏的 gene-sample 列跳过 IUBEC：
    ## 若该样本在该基因下 <2 个 transcript 非零，则不做 ILR/去 batch，后面直接保留原始 counts。
    nonzero_tx_per_sample <- colSums(Xg > 0)
    skip_cols_local <- (!zero_cols) & (nonzero_tx_per_sample < 2L)

    ## [v3-patch-2] 记录原始 transcript-level 0 的位置；
    ## 仅对“gene total > 0 且未被 skip”的列，后续才尝试做回零。
    orig_zero_mask <- (Xg == 0) &
      matrix(rep((!zero_cols) & (!skip_cols_local), each = K), nrow = K)

    ## 比例矩阵（严格闭合；全 0 列/skip 列先给一个占位比例，后续会将 Z 置 0，避免污染评分）
    Pg <- matrix(0, nrow = K, ncol = ncol(countData))
    inactive_cols <- zero_cols | skip_cols_local
    Pg[, inactive_cols] <- 1 / K

    ## [v3-patch-1] 自适应伪计数：
    ## 不再对所有 gene/sample 一律加固定 eps，而是随原始 gene total 调整伪计数大小。
    active_cols <- which(!inactive_cols)
    if (length(active_cols) > 0L) {
      eps_pc_vec <- pmin(1e-4, 0.5 / pmax(Xg_sum0[active_cols], 1) / K)
      eps_pc_mat <- matrix(rep(eps_pc_vec, each = K), nrow = K)
      denom_mat <- matrix(rep(Xg_sum0[active_cols] + K * eps_pc_vec, each = K), nrow = K)
      Pg[, active_cols] <- (Xg[, active_cols, drop = FALSE] + eps_pc_mat) / denom_mat
    }

    ## ILR：对每个样本列向量做 ilr，得到 n x (K-1)
    ## 注意：当 K=2 (K-1=1) 时，apply() 会把结果简化成向量（无 dim），
    ## 进而导致 t() 变成 1 x n（行/列颠倒），后续按样本行索引会报
    ## “(subscript) logical subscript too long”。这里强制恢复矩阵形状。
    Z_tmp <- apply(Pg, 2, ilr_vec)                 # (K-1) x n  (K=2 时可能退化为 length-n 向量)
    if (is.null(dim(Z_tmp))) {
      Z_tmp <- matrix(Z_tmp, nrow = K - 1L, ncol = ncol(Pg))
    }
    Z <- t(Z_tmp)                                  # n x (K-1)

    ## [v3-patch-4] 对全 0 列和极稀疏 skip 列，强制 Z=0，
    ## 避免这些列参与 batch/disease 方向学习与评分。
    if (any(inactive_cols)) {
      Z[inactive_cols, ] <- 0
    }

    ## -------- batch 对抗投影 --------
    if (have_batch) {
      if (n_batch_levels == 2L) {
        ## 二分类 batch：估计“最区分 batch 的单一方向”，并从所有样本中投影删除
        ## 方向估计使用所有可用样本，并只在当前 gene 的 active 列内学习 batch 方向
        idx_fit <- batch_fit_idx

        ## [v3-patch-4] 仅在当前 gene 的 active 列内学习 batch 方向；
        ## [v3-patch-3] 同时按 gene total 加权，降低低深度/高稀疏样本的影响。
        idx_fit <- intersect(idx_fit, which(!inactive_cols))

        ## 若子集中 batch 只有单一水平，或可用 active 样本过少，则无法估计方向：直接跳过
        if (length(idx_fit) < 2L || length(unique(y_batch01[idx_fit])) < 2L) {
          b <- rep(0, ncol(Z))
          Z_hat <- rep(0, nrow(Z))
          Z_adj <- Z
        } else {
          Z_fit <- Z[idx_fit, , drop = FALSE]
          y_fit <- y_batch01[idx_fit]
          w_fit <- pmax(Xg_sum0[idx_fit], 1)
          sw <- sqrt(w_fit)
          Z_fit_w <- Z_fit * sw
          y_fit_w <- y_fit * sw
          ZZ <- crossprod(Z_fit_w)
          b <- tryCatch(
            solve(ZZ + diag(ncol(Z)), crossprod(Z_fit_w, y_fit_w)),
            error = function(e) rep(0, ncol(Z))
          )
          nb2 <- sum(b * b) + eps
          Z_hat <- as.vector(Z %*% b)
          Z_adj <- Z - tcrossprod(Z_hat / nb2, b)
        }

      } else {
        ## 多水平 batch：one-hot 回归得到多方向，并删除其张成的子空间
        ## 方向估计使用所有可用样本，并只在当前 gene 的 active 列内学习 batch 方向
        idx_fit <- batch_fit_idx
        Yb_fit_full <- Yb_full_all[idx_fit, , drop = FALSE]

        ## [v3-patch-4] 仅在当前 gene 的 active 列内学习 batch 方向；
        ## [v3-patch-3] 同时按 gene total 加权，降低低深度/高稀疏样本的影响。
        idx_fit <- intersect(idx_fit, which(!inactive_cols))
        Yb_fit_full <- Yb_full_all[idx_fit, , drop = FALSE]

        ## 子集中若仅包含 1 个 batch 水平，或可用 active 样本过少，则跳过对抗投影
        n_present <- sum(colSums(Yb_fit_full) > 0)
        if (length(idx_fit) < 2L || n_present < 2L) {
          B <- matrix(0, nrow = ncol(Z), ncol = n_batch_levels,
                      dimnames = list(NULL, batch_levels))
          Z_adj <- Z
        } else {
          ## 去中心化（在“用于估计方向的子集”内去中心化）
          Yb_fit <- scale(Yb_fit_full, center = TRUE, scale = FALSE)
          colnames(Yb_fit) <- batch_levels
          Z_fit <- Z[idx_fit, , drop = FALSE]

          w_fit <- pmax(Xg_sum0[idx_fit], 1)
          sw <- sqrt(w_fit)
          Z_fit_w <- Z_fit * sw
          Yb_fit_w <- Yb_fit * sw

          ## B: (K-1) x L
          B <- tryCatch(
            solve(crossprod(Z_fit_w) + diag(ncol(Z)), crossprod(Z_fit_w, Yb_fit_w)),
            error = function(e) matrix(0, nrow = ncol(Z), ncol = ncol(Yb_fit))
          )
          colnames(B) <- batch_levels

          ## 删除 span(B)（用正交基 U）
          qrB <- qr(B)
          r <- qrB$rank
          if (r > 0) {
            U <- qr.Q(qrB, complete = FALSE)[, seq_len(r), drop = FALSE]   # (K-1) x r
            Z_adj <- Z - (Z %*% U) %*% t(U)
          } else {
            Z_adj <- Z
          }
        }
      }

    } else {
      Z_adj <- Z
    }


    ## IUBEC 只删除 batch 相关方向；不估计疾病方向，也不做额外增强。
    Z_used <- Z_adj
    ## 逆 ILR -> 比例 -> counts（按原基因总量 Xg_sum0 缩放，保证守恒）
    P_corr <- t(apply(Z_used, 1, ilr_inv_vec))      # n x K
    ## 数值安全：再次闭合
    P_corr <- t(apply(P_corr, 1, function(p) p / sum(p)))

    ## [v3-patch-2] 对“原来为 0 且逆变换后仍极小”的位置做回零，然后再归一化。
    ## 目的：尽量保留 transcript-level 结构性 0，避免它们被伪计数 + ILR + 线性投影永久抬成非零。
    active_cols2 <- which((!zero_cols) & (!skip_cols_local))
    if (length(active_cols2) > 0L) {
      for (j in active_cols2) {
        p <- P_corr[j, ]
        tau <- min(1e-6, 0.5 / max(Xg_sum0[j], 1))
        zmask <- orig_zero_mask[, j]
        p[zmask & (p < tau)] <- 0
        s <- sum(p)
        if (s > 0) {
          p <- p / s
        }
        P_corr[j, ] <- p
      }
    }

    X_corr <- t(P_corr) * matrix(Xg_sum0,
                                 nrow = K,
                                 ncol = length(Xg_sum0),
                                 byrow = TRUE)

    ## [v3-patch-4] 对极稀疏 skip 列直接保留原始 counts，不做 IUBEC 校正。
    if (any(skip_cols_local)) {
      X_corr[, skip_cols_local] <- Xg[, skip_cols_local, drop = FALSE]
    }

    corrected_counts[rows, ] <- X_corr
  }
  ## counts 下界裁剪（保持原逻辑）
  corrected_counts[corrected_counts < 0] <- 0

  ## ----------------- 拟合 satuRn GLM（不含 batch） -----------------
  if (!quiet) message("[fitDTU_IUBEC] 调用原生 fitDTU() （只用 group，不含 batch）...")

  sampleData_fit <- as.data.frame(sampleData)
  if ("batch" %in% colnames(sampleData_fit)) {
    sampleData_fit$batch <- NULL
  }
  ## 只保留 group（和可选的 sample 列），供 ~0 + group 使用
  keep_cols <- intersect(c("sample", "group"), colnames(sampleData_fit))
  sampleData_fit <- sampleData_fit[, keep_cols, drop = FALSE]

  ## 为 SummarizedExperiment 补 satuRn 常用列名（保持原逻辑）
  tx2gene_se <- tx2gene_model
  tx2gene_se$isoform_id <- tx2gene_se$TXNAME
  tx2gene_se$GENEID     <- tx2gene_se$gene_id

  sumExp <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(counts = corrected_counts),
    colData = S4Vectors::DataFrame(sampleData_fit),
    rowData = S4Vectors::DataFrame(tx2gene_se)
  )

  sm <- satuRn::fitDTU(
    object  = sumExp,
    formula = ~ 0 + group
  )

  if (!quiet) message("[fitDTU_IUBEC] 完成。")

  sm
}
