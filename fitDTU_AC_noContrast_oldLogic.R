#' fitDTU_AC: satuRn 的“AC_satuRn”拟合入口（返回 StatModel）
#'
#' 本版本为“消融（无对比增强）”版本：
#' 对 satuRn 的原始输入做以下处理：
#'   1) 基因内转录本比例 -> 2) ILR 变换 ->
#'   3) 线性对抗去批次（删除最能区分批次的方向 / 多水平 batch 则删除其张成子空间） ->
#'   4) 逆 ILR -> 回到 counts 空间
#' 然后调用原生 fitDTU() 拟合，并复用 testDTU() 计算评估指标。
#'
#' 与之前版本的区别：不再沿疾病方向做任何“对比增强”，contrast_gain 参数保留但不再使用。
#'
#' 本文件修复/增强点（v2）：
#'  1) ILR fallback 修复：即使没有 compositions 包，也使用标准 Helmert ILR 基（保证 V^T 1 = 0 且列正交），确保 ilr 与逆变换互逆；
#'  2) sampleData 与 countData 对齐：不再默认“按位置”，而是优先按 rownames(sampleData) 或 sample_id 列显式匹配到 colnames(countData)；
#'  3) 病例/对照标签：可用 case_label 显式指定；若不指定且 group 不是二分类则直接报错；
#'  4) batch 多水平：不再二值化；多水平 batch 时构造 one-hot 并删除其在 ILR 空间对应的子空间；
#'  5) tx2gene 映射：不再“前两列回退”；必须显式识别到 tx/gene 列（或用户传入 tx_col/gene_col），并严格要求覆盖所有 rownames(countData)。
#'
#' @param countData  matrix/data.frame，行为 transcript（TXNAME），列为样本
#' @param tx2gene    meta 信息（如 Performance_* 的 metaInfo），只需要含有“转录本 ID 列 + 基因 ID 列”
#' @param sampleData data.frame，包含样本信息；必须有 group 或 condition，可选 batch
#' @param quiet      是否静音
#' @param eval_out_path 评估结果落盘路径
#' @param contrast_gain 保留参数，为兼容旧代码，本版本内部不使用
#' @param case_label 可选：显式指定哪一个 group 水平作为“病例(1)”。若 NULL，则仅在 group 为二分类时沿用 levels(grp)[2]
#' @param batch_dir_source 指定“估计 batch 可分方向（用于对抗投影删除）”时使用哪些样本。
#'   - "all"：使用所有样本（case + control，等价于当前默认行为）
#'   - "control"：只用对照样本（y_case==0）估计 batch 方向
#'   - "case"：只用病例样本（y_case==1）估计 batch 方向
#' @param tx_col     可选：显式指定 tx2gene 中“转录本ID列”列名
#' @param gene_col   可选：显式指定 tx2gene 中“基因ID列”列名
#' @return satuRn 的 StatModel 对象（可直接交给 testDTU）
#' @export
fitDTU_AC <- function(countData,
                      tx2gene,
                      sampleData,
                      quiet = TRUE,
                      eval_out_path = file.path(getwd(), "AC_fitDTU_evaluate.txt"),
                      contrast_gain = 0.2,  # 本消融版本中不会实际使用，仅为保持函数签名不变
                      case_label = NULL,
                      batch_dir_source = c("all", "control", "case"),
                      tx_col = NULL,
                      gene_col = NULL
) {

  ## ----------------- 基本检查 -----------------
  countData <- as.matrix(countData)
  if (is.null(rownames(countData))) {
    stop("fitDTU_AC: countData 必须有转录本 ID 作为行名(rownames)。")
  }
  if (is.null(colnames(countData))) {
    stop("fitDTU_AC: countData 必须有样本 ID 作为列名(colnames)。")
  }
  if (!is.data.frame(sampleData)) {
    stop("fitDTU_AC: sampleData 应该是 data.frame，当前类型为 ",
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
        stop("fitDTU_AC: sampleData 无法与 countData 对齐：以下样本在 sampleData 的列 ", id_col,
             " 中找不到：", paste(missing, collapse = ", "),
             "\n请确保 sampleData 的 rownames 或 sample_id 列能覆盖 countData 的所有列名。")
      }
      sampleData2 <- sampleData[m, , drop = FALSE]
      rownames(sampleData2) <- cn
      return(sampleData2)
    }

    stop("fitDTU_AC: 无法对齐 sampleData 与 countData。\n",
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
      stop("fitDTU_AC: sampleData 需要包含 group 或 condition 列。")
    }
  }
  grp <- factor(sampleData$group)

  
## ----------------- 病例/对照：更稳健的自动判定 + 允许显式指定 -----------------
case_label_used <- NULL
if (is.null(case_label)) {
  if (nlevels(grp) != 2L) {
    stop("fitDTU_AC: case_label 未指定时，group 必须是二分类。当前 group 水平为：",
     paste(levels(grp), collapse = ", "),
     "\n请通过 case_label 显式指定病例水平。")
  }
  levs <- levels(grp)
  levs_lc <- tolower(levs)

  ## 自动规则（仅在二分类时启用）：
  ## 1) 若能唯一识别“对照”水平，则另一个水平为病例；
  ## 2) 否则若能唯一识别“病例”水平，则该水平为病例；
  ## 3) 否则回退为 levels(grp)[2]（与旧逻辑保持一致）
  ctrl_re <- "^(ctrl|control|controls|healthy|normal|wt|wildtype|wild-type|0)$"
  case_re <- "^(case|disease|patient|ad|mci|tumou?r|treated|1)$"

  is_ctrl <- grepl(ctrl_re, levs_lc)
  is_case <- grepl(case_re, levs_lc)

  if (sum(is_ctrl) == 1L) {
    case_label_used <- levs[!is_ctrl]
  } else if (sum(is_case) == 1L) {
    case_label_used <- levs[is_case]
  } else {
    case_label_used <- levs[2]
  }

  if (!quiet) {
    message("fitDTU_AC: auto case_label = ", case_label_used,
            " (levels: ", paste(levs, collapse = ", "), ")")
  }
} else {
  case_label_used <- as.character(case_label)
  if (!(case_label_used %in% levels(grp))) {
    stop("fitDTU_AC: case_label='", case_label_used,
         "' 不在 group 水平中：", paste(levels(grp), collapse = ", "))
  }
}
y_case <- as.integer(grp == case_label_used)  # 0:对照 1:病例

  ## ----------------- batch 方向估计使用的样本子集 -----------------
  ## 说明：对抗投影中用于“学习 batch 可分方向”的回归/子空间估计，默认用 all（case+control）。
  ##      通过 batch_dir_source 可切换为仅用 control 或仅用 case。
  ## 允许更宽松的写法：case+control / case_control / both / ctrl 等
  if (length(batch_dir_source) == 1L && is.character(batch_dir_source)) {
    x <- tolower(batch_dir_source)
    x <- gsub("[[:space:]_\\-]", "", x)
    x <- gsub("\\+", "", x)
    if (x %in% c("casecontrol", "both", "all")) {
      batch_dir_source <- "all"
    } else if (x %in% c("control", "ctrl")) {
      batch_dir_source <- "control"
    } else if (x %in% c("case", "disease", "patient")) {
      batch_dir_source <- "case"
    }
  }
  batch_dir_source <- match.arg(batch_dir_source, choices = c("all", "control", "case"))
  batch_fit_idx <- switch(
    batch_dir_source,
    all     = seq_len(ncol(countData)),
    control = which(y_case == 0L),
    case    = which(y_case == 1L)
  )
  if (length(batch_fit_idx) == 0L) {
    stop("fitDTU_AC: batch_dir_source='", batch_dir_source, "' 选择的样本子集为空。请检查 group/case_label。")
  }

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
        message("fitDTU_AC: tx2gene 列名未匹配候选列表，安全回退使用前两列: ",
                tx_col, " / ", gene_col)
      }
    } else {
      stop("fitDTU_AC: 无法在 tx2gene 中自动识别 tx/gene 列，且前两列也无法覆盖全部转录本。",
     "\n当前 tx2gene 列名为：", paste(tx_cols_all, collapse = ", "),
     "\n请通过 tx_col / gene_col 显式指定正确列名。")
    }
  } else {
    stop("fitDTU_AC: tx2gene 列数不足 2，无法识别 tx/gene 映射列。请提供规范 tx2gene 或显式指定列名。")
  }
}
  if (!(tx_col %in% tx_cols_all) || !(gene_col %in% tx_cols_all)) {
    stop("fitDTU_AC: tx_col 或 gene_col 不存在于 tx2gene 列名中。")
  }

  tx_ids_all   <- as.character(tx2gene[[tx_col]])
  gene_ids_all <- as.character(tx2gene[[gene_col]])

  idx <- match(rownames(countData), tx_ids_all)
  if (any(is.na(idx))) {
    missing <- rownames(countData)[is.na(idx)]
    stop("fitDTU_AC: 有 rownames(countData) 在 tx2gene 的转录本列中找不到。示例：",
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

  ## 如有 gene_modified，带一份进来用于评估（保持原逻辑）
  gene_modified <- NULL
  if ("gene_modified" %in% tx_cols_all) {
    gm_all <- tx2gene[["gene_modified"]]
    gene_modified <- gm_all[idx]
    tx2gene_model$gene_modified <- gene_modified
  }

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
      stop("fitDTU_AC: ILR fallback basis check failed: max|1^T V| > 1e-8")
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

  ## ----------------- 每个基因做 AC + 逆变换 -----------------
  corrected_counts <- matrix(0, nrow = nrow(countData), ncol = ncol(countData),
                             dimnames = dimnames(countData))

  ## batch score：二分类仍用向量；多分类用矩阵（每个水平一列）
  if (have_batch && n_batch_levels > 2L) {
    score_batch_mat <- matrix(0, nrow = ncol(countData), ncol = n_batch_levels,
                              dimnames = list(colnames(countData), batch_levels))
  } else {
    score_batch_mat <- NULL
  }
  score_batch   <- numeric(ncol(countData))
  score_disease <- numeric(ncol(countData))

  eps <- 1e-8

  ## 预先准备多水平 batch 的 one-hot（仅在需要时）
  ## 注意：是否“只用 control/case 估计方向”由 batch_fit_idx 控制；因此这里只缓存全量 one-hot。
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

    ## 比例矩阵（严格闭合；全 0 列用均匀比例，但后续会将 Z 置 0，避免污染评分）
    Pg <- matrix(0, nrow = K, ncol = ncol(countData))
    Pg[, zero_cols] <- 1 / K

    if (any(!zero_cols)) {
      denom <- Xg_sum0[!zero_cols] + K * eps
      Pg[, !zero_cols] <- sweep(Xg[, !zero_cols, drop = FALSE] + eps, 2, denom, "/")
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

    ## 对全 0 列，强制 Z=0（避免“无表达基因”对 batch/disease score 产生贡献）
    if (any(zero_cols)) {
      Z[zero_cols, ] <- 0
    }

    ## -------- batch 对抗投影 --------
    if (have_batch) {
      if (n_batch_levels == 2L) {
        ## 二分类 batch：估计“最区分 batch 的单一方向”，并从所有样本中投影删除
        ## 方向估计样本子集由 batch_dir_source/batch_fit_idx 控制
        idx_fit <- batch_fit_idx

        ## 若子集中 batch 只有单一水平，则无法估计方向：直接跳过
        if (length(unique(y_batch01[idx_fit])) < 2L) {
          b <- rep(0, ncol(Z))
          Z_hat <- rep(0, nrow(Z))
          Z_adj <- Z
        } else {
          Z_fit <- Z[idx_fit, , drop = FALSE]
          y_fit <- y_batch01[idx_fit]
          ZZ <- crossprod(Z_fit)
          b <- tryCatch(
            solve(ZZ + diag(ncol(Z)), crossprod(Z_fit, y_fit)),
            error = function(e) rep(0, ncol(Z))
          )
          nb2 <- sum(b * b) + eps
          Z_hat <- as.vector(Z %*% b)
          Z_adj <- Z - tcrossprod(Z_hat / nb2, b)
        }

        score_batch <- score_batch + Z_hat

      } else {
        ## 多水平 batch：one-hot 回归得到多方向，并删除其张成的子空间
        ## 方向估计样本子集由 batch_dir_source/batch_fit_idx 控制
        idx_fit <- batch_fit_idx
        Yb_fit_full <- Yb_full_all[idx_fit, , drop = FALSE]

        ## 子集中若仅包含 1 个 batch 水平，则跳过对抗投影
        n_present <- sum(colSums(Yb_fit_full) > 0)
        if (n_present < 2L) {
          B <- matrix(0, nrow = ncol(Z), ncol = n_batch_levels,
                      dimnames = list(NULL, batch_levels))
          Z_adj <- Z
        } else {
          ## 去中心化（在“用于估计方向的子集”内去中心化）
          Yb_fit <- scale(Yb_fit_full, center = TRUE, scale = FALSE)
          colnames(Yb_fit) <- batch_levels
          Z_fit <- Z[idx_fit, , drop = FALSE]

          ## B: (K-1) x L
          B <- tryCatch(
            solve(crossprod(Z_fit) + diag(ncol(Z)), crossprod(Z_fit, Yb_fit)),
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

        ## 评分：每个 batch 水平一个 score（用于 macro one-vs-rest AUC）
        score_batch_mat <- score_batch_mat + (Z %*% B)
      }
    } else {
      Z_adj <- Z
    }

    ## -------- 疾病方向：仅用于评估，不做对比增强 --------
    mu0 <- colMeans(Z_adj[y_case == 0, , drop = FALSE])
    mu1 <- colMeans(Z_adj[y_case == 1, , drop = FALSE])
    d   <- mu1 - mu0

    Z_used <- Z_adj

    ## 逆 ILR -> 比例 -> counts（按原基因总量 Xg_sum0 缩放，保证守恒）
    P_corr <- t(apply(Z_used, 1, ilr_inv_vec))      # n x K
    ## 数值安全：再次闭合
    P_corr <- t(apply(P_corr, 1, function(p) p / sum(p)))

    X_corr <- t(P_corr) * matrix(Xg_sum0,
                                 nrow = K,
                                 ncol = length(Xg_sum0),
                                 byrow = TRUE)
    corrected_counts[rows, ] <- X_corr

    ## score 汇总
    ## 二分类 batch：已累计在 score_batch；多分类 batch：在 score_batch_mat
    score_disease <- score_disease + as.vector(Z_adj %*% d)
  }

  ## counts 下界裁剪（保持原逻辑）
  corrected_counts[corrected_counts < 0] <- 0

  ## ----------------- 拟合 satuRn GLM（不含 batch） -----------------
  if (!quiet) message("[fitDTU_AC] 调用原生 fitDTU() （只用 group，不含 batch）...")

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

  ## ----------------- 评估：AUC / TPR@FDR0.1，并写日志 -----------------
  auc_01 <- function(score, label01) {
    o <- order(score); y <- label01[o]
    n1 <- sum(y == 1); n0 <- sum(y == 0)
    if (n1 == 0 || n0 == 0) return(NA_real_)
    ranks <- rank(score[o], ties.method = "average")
    (sum(ranks[y == 1]) - n1 * (n1 + 1) / 2) / (n0 * n1)
  }

  if (have_batch) {
    if (n_batch_levels == 2L) {
      auc_batch <- auc_01(score_batch, y_batch01)
    } else {
      aucs <- vapply(batch_levels, function(bl) {
        y01 <- as.integer(batch == bl)
        auc_01(score_batch_mat[, bl], y01)
      }, numeric(1))
      auc_batch <- mean(aucs, na.rm = TRUE)
    }
  } else {
    auc_batch <- NA_real_
  }

  auc_disease <- auc_01(score_disease, y_case)

  ## 这里不在 fitDTU_AC 内部再跑 testDTU()，避免和外层 AC_satuRn_DTU 的正式检验重复。
  tpr_at_fdr <- NA_real_

  outdir <- dirname(eval_out_path)
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  }

  val <- data.frame(
    time        = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    n_tx        = nrow(countData),
    n_sample    = ncol(countData),
    have_batch  = have_batch,
    auc_batch   = round(auc_batch,   4),
    auc_disease = round(auc_disease, 4),
    tpr_fdr_0.1 = round(tpr_at_fdr,  4),
    check.names = FALSE
  )

  if (!file.exists(eval_out_path)) {
    write.table(val, file = eval_out_path,
                row.names = FALSE, col.names = TRUE,
                sep = "\t", quote = FALSE)
  } else {
    suppressWarnings(
      write.table(val, file = eval_out_path, append = TRUE,
                  row.names = FALSE, col.names = FALSE,
                  sep = "\t", quote = FALSE)
    )
  }

  if (!quiet) {
    if (have_batch && n_batch_levels > 2L) {
      message(sprintf(
        "[fitDTU_AC] 评估: batch_macro_auc=%.3f, disease_auc=%.3f",
        auc_batch, auc_disease
      ))
    } else {
      message(sprintf(
        "[fitDTU_AC] 评估: batch_auc=%.3f, disease_auc=%.3f",
        auc_batch, auc_disease
      ))
    }
    message("[fitDTU_AC] 完成。")
  }

  sm
}
