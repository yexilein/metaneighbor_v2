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
MetaNeighborUS <- function(dat, study_id, cell_type, ranked = TRUE, n_centroids = 0, fast = TRUE) {
  check_input(dat, study_id, cell_type)
  dat <- normalize_cols(dat)
  colnames(dat) <- paste(study_id, cell_type, sep = "|")
  studies <- unique(study_id)
  data_subsets <- find_subsets(study_id, studies)
  if (n_centroids > 0) {
    centroids <- compute_centroids(dat, n_centroids)
    centroid_subsets <- find_subsets(get_study_id(colnames(centroids)), studies)
  }

  result <- create_result_matrix(colnames(dat))
  for (study_A_index in seq_along(studies)) {
    study_A <- dat[, data_subsets[, study_A_index]]
    if (n_centroids > 0) {
      centroids_A <- centroids[, centroid_subsets[, study_A_index]]
    }
    for (study_B_index in study_A_index:length(studies)) {
      study_B <- dat[, data_subsets[, study_B_index]]
      if (n_centroids > 0) {
        centroids_B <- centroids[, centroid_subsets[, study_B_index]]
      }

      # study B votes for study A
      if (fast) {
        if (n_centroids > 0) {
          votes <- compute_votes_without_network(study_A, centroids_B)
        } else {
          votes <- compute_votes_without_network(study_A, study_B)
        }
      } else {
        if (n_centroids > 0) {
          network <- build_network(study_A, centroids_B, ranked = ranked)
        } else {
          network <- build_network(study_A, study_B, ranked = ranked)
        }
        votes <- compute_votes_from_network(network)
      }
      aurocs <- compute_aurocs(votes)
      result[rownames(aurocs), colnames(aurocs)] <- aurocs

      # study A votes for study B
      if (fast) {
        if (n_centroids > 0) {
          votes <- compute_votes_without_network(study_B, centroids_A)
        } else {
          votes <- compute_votes_without_network(study_B, study_A)
        }
      } else {
        if (n_centroids > 0) {
          network <- build_network(study_B, centroids_A, ranked = ranked)
        } else {
          network <- t(network)
        }
        votes <- compute_votes_from_network(network)
      }
      aurocs <- compute_aurocs(votes)
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
  M <- scale(M, scale=FALSE)
  return(sweep(M, 2, apply(M, 2, function(c) sqrt(crossprod(c,c))), FUN="/"))
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

find_subsets <- function(full_list, list_names) {
  return(sapply(list_names, function(name) full_list == name))
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

pseudo_rank <- function(x, breaks = 1000, depth = 1000) {
  m <- min(x)
  M <- max(x)
  bins <- floor((x-m) / ((1+1e-10)*(M-m)) * breaks) + 1
  num_per_bin <- rep(0, breaks)
  if (is.null(depth)) {
    for (row in seq_len(nrow(bins))) {
      for (val in bins[row,]) { num_per_bin[val] <- num_per_bin[val] + 1 }
    }
  } else {
    indices <- sample.int(length(bins), breaks*depth)
    for (i in indices) {
      val <- bins[i]
      num_per_bin[val] <- num_per_bin[val] + 1
    }
  }
  rank_per_bin <- (c(0, cumsum(num_per_bin)[-length(num_per_bin)]) +
                   (num_per_bin+1)/2) / sum(num_per_bin)
  return(rank_per_bin[bins])
}

compute_aurocs <- function(votes) {
  candidate_labels <- rownames(votes)
  positives <- design_matrix(candidate_labels)
  n_positives <- colSums(positives)
  n_negatives <- nrow(positives) - n_positives
  sum_of_positive_ranks <- t(positives) %*% apply(abs(votes), MARGIN = 2, FUN = rank)
  result <- (sum_of_positive_ranks / n_positives - (n_positives+1)/2) / n_negatives
  return(result)
}

compute_votes_from_network <- function(network) {
  voter_labels <- colnames(network)
  return(network %*% design_matrix(voter_labels) / rowSums(network))
}

compute_votes_without_network <- function(candidates, voters) {
  voter_identity <- design_matrix(colnames(voters))
  raw_votes <- crossprod(candidates, voters %*% voter_identity)
  return(sweep(raw_votes, 2, colSums(voter_identity), FUN = "+")) /
         (c(crossprod(candidates, rowSums(voters))) + rep(ncol(voters), ncol(candidates)))
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
