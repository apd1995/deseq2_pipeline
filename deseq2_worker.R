suppressPackageStartupMessages({
  library(DESeq2)
  library(glmGamPoi)
  library(data.table)
})

args       <- commandArgs(trailingOnly = TRUE)
input_path <- args[1]
out_path   <- args[2]
worker_id  <- args[3]

df <- data.table::fread(input_path)

meta       <- df[, .(perturbation, replicate, condition, n_pert_cells, n_ctrl_cells)]
counts_dt  <- df[, !c("perturbation", "replicate", "condition",
                       "n_pert_cells", "n_ctrl_cells")]
gene_names <- colnames(counts_dt)
pert       <- unique(meta$perturbation)[1]

cat(sprintf("[worker %s] starting: %s\n", worker_id, pert))
t0 <- proc.time()["elapsed"]

idx       <- which(meta$perturbation == pert)
count_mat <- t(as.matrix(counts_dt[idx, ]))
colnames(count_mat) <- paste0("s", seq_len(ncol(count_mat)))

# DESeq2 requires integer counts — round and cast
# pseudobulk sums from sparse float32 matrices can have tiny floating point residuals (e.g. 4.9999999) so rounding before casting is essential
if (!is.integer(count_mat)) {
  count_mat <- round(count_mat)
  storage.mode(count_mat) <- "integer"
}

col_data <- data.frame(
  condition = factor(meta$condition[idx], levels = c("ctrl", "pert"))
)

result <- tryCatch({
  dds <- DESeqDataSetFromMatrix(count_mat, col_data, ~ condition)
  dds <- DESeq(dds, fitType = "glmGamPoi", quiet = TRUE)
  res <- results(dds, contrast = c("condition", "pert", "ctrl"))
  data.frame(
    gene           = gene_names,
    perturbation   = pert,
    baseMean       = res$baseMean,
    log2FoldChange = res$log2FoldChange,
    fold_change    = 2^res$log2FoldChange,   # linear fold change
    lfcSE          = res$lfcSE,
    stat           = res$stat,
    pvalue         = res$pvalue,
    fdr            = res$padj               # BH-adjusted p-value
  )
}, error = function(e) {
  cat(sprintf("[worker %s] ERROR for %s: %s\n", worker_id, pert, e$message))
  NULL
})

elapsed <- proc.time()["elapsed"] - t0
cat(sprintf("[worker %s] done: %s  (%.2fs)\n", worker_id, pert, elapsed))

if (!is.null(result)) {
  data.table::fwrite(result, out_path)
}
