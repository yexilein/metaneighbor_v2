#' Runs unsupervised version of MetaNeighbor
#'
#' When it is difficult to know how cell type labels compare across datasets this
#' function helps users to make an educated guess about the overlaps without
#' requiring in-depth knowledge of marker genes
#'
#' @param var_genes vector of high variance genes.
#' @param dat SummarizedExperiment object containing gene-by-sample
#' expression matrix.
#' @param i default value 1; non-zero index value of assay containing the matrix
#' data
#' @param study_id a vector that lists the Study (dataset) ID for each sample
#' @param cell_type a vector that lists the cell type of each sample
#'
#' @return The output is a cell type-by-cell type mean AUROC matrix, which is
#' built by treating each pair of cell types as testing and training data for
#' MetaNeighbor, then taking the average AUROC for each pair (NB scores will not
#' be identical because each test cell type is scored out of its own dataset,
#' and the differential heterogeneity of datasets will influence scores).
#'
#' @examples
#' data(mn_data)
#' var_genes = variableGenes(dat = mn_data, exp_labels = mn_data$study_id)
#' celltype_NV = MetaNeighborUS(var_genes = var_genes,
#'                              dat = mn_data,
#'                              study_id = mn_data$study_id,
#'                              cell_type = mn_data$cell_type)
#' celltype_NV
#'
#' @export
#'

MetaNeighborUS <- function(dat, study_id, cell_type, ranked = TRUE, n_centroids = 0) {
  check_input(dat, study_id, cell_type)
  dat <- normalize_cols(dat)
  colnames(dat) <- paste(study_id, cell_type, sep = "|")
  if (n_centroids > 0) {
    centroids <- compute_centroids(dat, n_centroids)
  }

  result <- create_result_matrix(colnames(dat))
  studies <- unique(study_id)
  for (study_A_index in seq_along(studies)) {
    study_A <- studies[study_A_index]
    for (study_B in studies[study_A_index:length(studies)]) {
      # study B votes for study A
      candidates <- get_subset(dat, study_A)
      if (n_centroids > 0) {
        voters <- get_subset(centroids, study_B)
      } else {
        voters <- get_subset(dat, study_B)
      }
      network <- build_network(candidates, voters, ranked = ranked)
      aurocs <- compute_aurocs(network)
      result[rownames(aurocs), colnames(aurocs)] <- aurocs

      # study A votes for study B
      if (n_centroids > 0) {
        network <- build_network(get_subset(dat, study_B),
                                 get_subset(centroids, study_A),
                                 ranked = ranked)
      } else {
        network <- t(network)
      }
      aurocs <- compute_aurocs(network)
      result[rownames(aurocs), colnames(aurocs)] <- aurocs
    }
  }
  return(result)
}

check_input <- function(dat, study_id, cell_type) {
  if(length(study_id) != ncol(dat)){
    stop('study_id length does not match number of samples')
  }
  if(length(cell_type) != ncol(dat)){
    stop('cell_type length does not match number of samples')
  }
  if(length(unique(study_id)) < 2){
      stop('Found only 1 unique experiment_label. Please use data from more than 1 study!')
  }
}

normalize_cols <- function(M, ranked = TRUE) {
  if (ranked) { M <- apply(M, 2, rank) }
  M <- t(M) - apply(M, 2, mean)
  return(t(M / apply(M, 1, norm, type = "2")))
}

compute_centroids <- function(dat, n_centroids = 1) {
  result <- lapply(unique(colnames(dat)),
                   function(c) k_mean_centroids(dat[, colnames(dat) == c], c, n_centroids))
  return(do.call(cbind, result))
}

k_mean_centroids <- function(dat, label, n_centroids) {
  # careful: dat can be reduced to a single vector
  dat <- as.matrix(dat)
  if (ncol(dat) > n_centroids) {
    result <- t(kmeans(t(dat), n_centroids)$centers)
  } else {
    result <- dat
  }
  colnames(result) <- rep(label, ncol(result))
  return(result)
}

get_subset <- function(dat, study_name) {
  study <- get_study_id(colnames(dat))
  return(dat[, study == study_name])
}

build_network <- function(set_A, set_B, ranked = TRUE) {
  result <- crossprod(set_A, set_B)
  if (!ranked) { return(result) }
  A_labels <- rownames(result)
  B_labels <- colnames(result)
  result <- matrix(pseudo_rank(result), nrow = nrow(result))
  rownames(result) <- A_labels
  colnames(result) <- B_labels
  return(result)
}

pseudo_rank <- function(x, breaks = 1000) {
  m <- min(x)
  M <- max(x)
  bins <- floor((x-m) / ((1+1e-10)*(M-m)) * breaks) + 1
  num_per_bin <- rep(0, breaks)
  for (row in seq_len(nrow(x))) {
    for (val in x[row,]) { num_per_bin[val] <- num_per_bin[val] }
  }
  rank_per_bin <- (c(0, cumsum(num_per_bin)[-length(num_per_bin)]) +
                   (num_per_bin+1)/2) / length(x)
  return(rank_per_bin[as.numeric(bins)])
}

compute_aurocs <- function(network) {
  voter_labels <- colnames(network)
  candidate_labels <- rownames(network)
  votes <- network %*% design_matrix(voter_labels) / rowSums(network)
  positives <- design_matrix(candidate_labels)
  n_positives <- colSums(positives)
  n_negatives <- nrow(positives) - n_positives
  sum_of_positive_ranks <- t(positives) %*% apply(abs(votes), MARGIN = 2, FUN = rank)
  result <- (sum_of_positive_ranks / n_positives - (n_positives+1)/2) / n_negatives
  return(result)
}

design_matrix <- function(cell_type) {
  result <- model.matrix(~cell_type-1)
  colnames(result) <- levels(as.factor(cell_type))
  return(result)
}

create_result_matrix <- function(cell_type) {
  unique_cell_type <- unique(cell_type)
  result <- matrix(0, nrow = length(unique_cell_type), ncol = length(unique_cell_type))
  rownames(result) <- unique_cell_type
  colnames(result) <- unique_cell_type
  return(result)
}

get_study_id <- function(cluster_name) {
  return(sapply(strsplit(cluster_name, "\\|"), head, 1))
}

get_cell_type <- function(cluster_name) {
  return(sapply(strsplit(cluster_name, "\\|"), tail, 1))
}

plot_NV_heatmap <- function(
  dat, reorder_entries = TRUE, breaks = seq(0, 1, length = 21),
  label_size = 0.3, norm = "", row_colors, col_colors
) {
  cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(20))
  if (reorder_entries) {
    reorder_entries <- as.dendrogram(hclust(as.dist(1-dat)))
  }
  if (norm == "rank") {
     dat <- rank_normalize(dat)
  } else if (norm == "log") {
     dat <- log_normalize(dat)
  }
  gplots::heatmap.2(
    dat, margins = c(1,11),
    key = TRUE, keysize = 1, key.xlab="AUROC", key.title="NULL",
    labRow = NA, labCol = NA,
    trace = "none", density.info = "none", col = cols, breaks = breaks,
#    offsetRow=0.1, offsetCol=0.1, cexRow = label_size, cexCol = label_size,
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
