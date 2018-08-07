#' Identify a highly variable gene set
#'
#' Identifies genes with high variance compared to their median expression
#' (top quartile) within each experimentCertain function
#'
#' @param dat SummarizedExperiment object containing gene-by-sample
#' expression matrix.
#' @param i default value 1; non-zero index value of assay containing the matrix data
#' @param exp_labels character vector that denotes the source (Study ID) of
#' each sample.
#'
#' @return The output is a vector of gene names that are highly variable in
#' every experiment (intersect)
#'
#' @examples
#' data(mn_data)
#' var_genes = variableGenes(dat = mn_data, exp_labels = mn_data$study_id)
#' var_genes
#'
#' @export
#'

variableGenes <- function(dat, i = 1, exp_labels) {
  dat <- SummarizedExperiment::assay(dat, i = i)
  check_variable_gene_input(dat, exp_labels)
  result <- vector("list")
  for (exp in unique(exp_labels)) {
    gene_list <- variable_genes_from_experiment(dat[ , exp_labels == exp])
    # we remove the last bin, corresponding to highly expressed genes
    result[[length(result) + 1]] <- unlist(gene_list[-length(gene_list)])
  }
  return(result)
}

check_variable_gene_input <- function(dat, exp_labels) {
  if(length(exp_labels) != length(colnames(dat))){
      stop('experiment_labels length does not match number of samples')
  }
  if(length(unique(exp_labels)) < 2){
      stop('Found only 1 unique exp_labels. Please use data from more than 1 study!')
  }
}

variable_genes_from_experiment <- function(data_subset) {
  quantile_data <- get_deciles(data_subset[Matrix::rowSums(data_subset) > 0, ])
  variance_data <- apply(data_subset, MARGIN = 1, FUN = stats::var)
  result <- vector("list", length = length(levels(quantile_data)))
  for(i in seq_along(levels(quantile_data))) {
    var_temp <- variance_data[quantile_data == levels(quantile_data)[i]]
    quant_var <- stats::quantile(var_temp, na.rm = TRUE)
    result[[i]] <- names(var_temp)[var_temp > quant_var[4]]
  }
  return(result)
}

get_deciles <- function(data_subset) {
  geom_mean <- apply(data_subset, MARGIN = 1, FUN = function(x) sum(log1p(x)))
  quants <- stats::quantile(geom_mean,
                            probs = seq(from = 0, to = 1, length = 11),
                            type = 5)
  return(cut(geom_mean, quants, include.lowest = TRUE))
}

extract_variable_genes <- function(dat, var_genes) {
  matching_vargenes <- match(rownames(dat), var_genes)
  matching_vargenes_count <- sum(!is.na(matching_vargenes))
  if(matching_vargenes_count < 2){
    stop("matching_vargenes should have more than 1 matching genes!",
         call. = TRUE)
  } else if(matching_vargenes_count < 5) {
    warning("matching_vargenes should have more matching genes!",
            immediate. = TRUE)
  }
  return(dat[!is.na(matching_vargenes),])
}
