
iterative_metaneighbor <- function(dataset, depth=3, output_dir, suffix = "") {
  aurocs <- MetaNeighborUS(dataset$data, dataset$study_id, dataset$cell_type)
  test_plots(aurocs, file.path(output_dir, paste0("aurocs", suffix, ".pdf")))
  if (depth > 1) {
    new_sets <- divide_dataset(dataset, aurocs)
    iterative_metaneighbor(new_sets[[1]], depth-1, output_dir, paste0(suffix, 1))
    iterative_metaneighbor(new_sets[[2]], depth-1, output_dir, paste0(suffix, 2))
  }
}

divide_dataset <- function(dataset_, aurocs) {
  group_labels <- divide_labels(aurocs)
  group_1 <- colnames(dataset_$data) %in% group_labels[[1]]
  group_2 <- colnames(dataset_$data) %in% group_labels[[2]]
  return(list(dataset_subset(dataset_, group_1),
              dataset_subset(dataset_, group_2)))
}

divide_labels <- function(auroc_table) {
  groups <- cutree(hclust(auroc_as_distance(auroc_table)), 2)
  return(list(names(groups)[groups == 1],
              names(groups)[groups == 2]))
}

auroc_as_distance <- function(auroc_table) {
  symmetric_auroc <- pmin(auroc_table, t(auroc_table))
  return(as.dist(1-symmetric_auroc))
}

dataset_subset <- function(dataset_, column_selector) {
  return(list(data = dataset_$data[, column_selector],
              study_id = dataset_$study_id[column_selector],
              cell_type = dataset_$cell_type[column_selector]))
}
