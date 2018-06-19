
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
