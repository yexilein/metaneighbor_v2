
bplot <- function(nv_mat) {
    Celltype = rep(colnames(nv_mat),each=dim(nv_mat)[1])
    ROCValues = unlist(lapply(seq_len(dim(nv_mat)[2]), function(i) nv_mat[,i]))
    beanplot::beanplot(ROCValues ~ Celltype,
                       border="NA",
                       col="gray",
                       ylab="AUROC",
                       what=c(0,1,1,1),
                       frame.plot = FALSE,
		       las = 2)
}

plot_NV_heatmap <- function(
  dat, reorder_entries = TRUE, breaks = seq(0, 1, length = 21),
  label_size = 0.2, norm = "", row_colors, col_colors
) {
  cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(20))
  if (reorder_entries) {
    reorder_entries <- as.dendrogram(hclust(as.dist(1-dat), method="average"))
  }
  if (norm == "rank") {
     dat <- rank_normalize(dat)
  } else if (norm == "log") {
     dat <- log_normalize(dat)
  }
  gplots::heatmap.2(
    dat, margins = c(11,11),
    key = TRUE, keysize = 1, key.xlab="AUROC", key.title="NULL",
    labRow = NA, offsetRow=0.1, cexRow = label_size,
    labCol = NA, offsetCol=0.1, cexCol = label_size,
    trace = "none", density.info = "none", col = cols, breaks = breaks,
    Rowv = reorder_entries, Colv = reorder_entries, dendrogram = "row",
    RowSideColors = row_colors$colors, ColSideColors = col_colors$colors
  )
  par(lend = 1)
  legend("topright", inset = c(0, .2),
         legend = col_colors$legend,
         col = col_colors$color_scale, pt.cex = 1, cex = 0.5, lwd = 10, bty="n")
  legend("topright", inset = c(.05, .3),
         legend = row_colors$legend,
         col = row_colors$color_scale, cex = 0.5, lwd = 10, bty="n")
}

rank_normalize <- function(matrix_) {
  row_labels <- rownames(result)
  col_labels <- colnames(result)
  result <- matrix(rank(matrix_), nrow = nrow(matrix_))
  rownames(result) <- row_labels
  colnames(result) <- col_labels
  return(result / max(result))
}

log_normalize <- function(matrix_) {
  result <- matrix_
  high <- result[result <= 0.5]
  low <- result[result > 0.5]
  result[result <= 0.5] <- log(high * (1 - high))
  result[result > 0.5] <- -log(low * (1 - low))
  result <- result - min(result)
  return(result / max(result))
}
