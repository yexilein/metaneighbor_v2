#' Runs MetaNeighbor
#'
#' For each gene set of interest, the function builds a network of rank
#' correlations between all cells. Next,It builds a network of rank correlations
#' between all cells for a gene set. Next, the neighbor voting predictor
#' produces a weighted matrix of predicted labels by performing matrix
#' multiplication between the network and the binary vector indicating cell type
#' membership, then dividing each element by the null predictor (i.e., node
#' degree). That is, each cell is given a score equal to the fraction of its
#' neighbors (including itself), which are part of a given cell type. For
#' cross-validation, we permute through all possible combinations of
#' leave-one-dataset-out cross-validation, and we report how well we can recover
#' cells of the same type as area under the receiver operator characteristic
#' curve (AUROC). This is repeated for all folds of cross-validation, and the
#' mean AUROC across folds is reported. Calls
#' \code{\link{neighborVoting}}.
#'
#' @param dat A SummarizedExperiment object containing gene-by-sample
#' expression matrix.
#' @param i default value 1; non-zero index value of assay containing the matrix
#' data
#' @param experiment_labels A numerical vector that indicates the source of each
#' sample.
#' @param celltype_labels A matrix that indicates the cell type of each sample.
#' @param genesets Gene sets of interest provided as a list of vectors.
#' @param bplot default true, beanplot is generated
#' @return A matrix of AUROC scores representing the mean for each gene set
#' tested for each celltype is returned directly (see \code{\link{neighborVoting}}).
#'
#' @seealso \code{\link{neighborVoting}}
#' @examples
#' data("mn_data")
#' data("GOmouse")
#' library(SummarizedExperiment)
#' AUROC_scores = MetaNeighbor(dat = mn_data,
#'                             experiment_labels = as.numeric(factor(mn_data$study_id)),
#'                             celltype_labels = metadata(colData(mn_data))[["cell_labels"]],
#'                             genesets = GOmouse,
#'                             bplot = TRUE)
#' @export
#'
MetaNeighbor <- function(dat, study_id, cell_type, genesets) {
  check_input(dat, study_id, cell_type)
  check_genesets(dat, genesets)
  colnames(dat) <- cell_type

  nv_mat <- matrix(0, ncol = length(unique(cell_type)), nrow = length(genesets))
  rownames(nv_mat) <- names(genesets)
  colnames(nv_mat) <- levels(as.factor(cell_type))
  for (l in seq_along(genesets)) {
      print(names(genesets)[l])
      geneset <- genesets[[l]]
      geneset_dat <- dat[!is.na(match(rownames(dat), geneset)), ]
      geneset_dat <- normalize_cols(geneset_dat)
      aurocs <- c()
      for (study in unique(study_id)) {
        network <- build_network(geneset_dat[, study_id == study],
                                 geneset_dat[, study_id != study])
        aurocs <- rbind(aurocs, diag(compute_aurocs(network)))
      }
      nv_mat[l,] <- colMeans(aurocs)
  }
  return(nv_mat)
}

check_genesets <- function(dat, genesets) {
  #check genesets matches more than 1 genes in gene_matrix
  genes_in_geneset <- as.character(unlist(genesets))
  genes_in_matrix <- rownames(dat)
  if(length(intersect(genes_in_geneset,genes_in_matrix)) < 1)
      stop('No matching genes between genesets and gene_matrix')
}

bplot <- function(nv_mat) {
    Celltype = rep(colnames(nv_mat),each=dim(nv_mat)[1])
    ROCValues = unlist(lapply(seq_len(dim(nv_mat)[2]), function(i) nv_mat[,i]))
    beanplot::beanplot(ROCValues ~ Celltype,
                       border="NA",
                       col="gray",
                       ylab="AUROC",
                       what=c(0,1,1,1),
                       frame.plot = FALSE)
}
