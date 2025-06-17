


library(CSCORE)
library(Seurat)
library(abind)


csCORE_pval <- function(ourdf, genecols, metacols, groupcol, logcountcol, B){

  grouplevels <- as.character(levels(ourdf[, groupcol]))
  
  if (sum(ourdf[, groupcol] == grouplevels[1]) <= sum(ourdf[, groupcol] == grouplevels[2])){
    smallgroup <- grouplevels[1]
    biggroup <- grouplevels[2]    
  } else {
    biggroup <- grouplevels[1]
    smallgroup <- grouplevels[2]
  }
  
  nsmall <- sum(ourdf[, groupcol] == smallgroup)
  nbig <- sum(ourdf[, groupcol] == biggroup)
  
  if (nsmall < nbig){
    drop_idcs <- sample(which(ourdf[, groupcol] == biggroup), size=nbig - nsmall)
    ourdf <- ourdf[-drop_idcs, ]
  }
  
  ourdf_pre <- ourdf[ourdf[, groupcol] == grouplevels[1], ]
  ourdf_post <- ourdf[ourdf[, groupcol] == grouplevels[2], ]
  
  seurat_pre <- CreateSeuratObject(counts = as(t(as.matrix(ourdf_pre[, genecols])), "dgCMatrix"), meta.data = ourdf_pre[, metacols])
  
  cscore_pre <- CSCORE(object = seurat_pre, genes = rownames(seurat_pre), seq_depth = exp(seurat_pre@meta.data[[logcountcol]]))
  
  seurat_post <- CreateSeuratObject(counts = as(t(as.matrix(ourdf_post[, genecols])), "dgCMatrix"), meta.data = ourdf_post[, metacols])
  
  cscore_post <- CSCORE(object = seurat_post, genes = rownames(seurat_post), seq_depth = exp(seurat_post@meta.data[[logcountcol]]))
  
  corrdiff <- cscore_post$est - cscore_pre$est
  
  
  mat_list <- list()
  for (i in 1:B) {
    # now let's shuffle pre cells and post cells
    ourdf_temp <- ourdf
    ourdf_temp[, groupcol] <- ourdf_temp[sample(1:nrow(ourdf_temp), size=nrow(ourdf_temp)), groupcol]
    
    ourdf_pre_temp <- ourdf_temp[ourdf_temp[, groupcol] == grouplevels[1], ]
    ourdf_post_temp <- ourdf_temp[ourdf_temp[, groupcol] == grouplevels[2], ]
    
    seurat_pre_temp <- CreateSeuratObject(counts = as(t(as.matrix(ourdf_pre_temp[, genecols])), "dgCMatrix"), meta.data = ourdf_pre_temp[, metacols])
    
    cscore_pre_temp <- CSCORE(object = seurat_pre_temp, genes = rownames(seurat_pre_temp), seq_depth = exp(seurat_pre_temp@meta.data[[logcountcol]]))
    
    seurat_post_temp <- CreateSeuratObject(counts = as(t(as.matrix(ourdf_post_temp[, genecols])), "dgCMatrix"), meta.data = ourdf_post_temp[, metacols])
    
    cscore_post_temp <- CSCORE(object = seurat_post_temp, genes = rownames(seurat_post_temp), seq_depth = exp(seurat_post_temp@meta.data[[logcountcol]]))
    
    corrdiff_temp <- cscore_post_temp$est - cscore_pre_temp$est
    
    mat_list[[i]] <- corrdiff_temp
  }
  
  stacked_array <- do.call(abind, c(mat_list, along = 3))
  
  dim(stacked_array)
  
  n_row <- dim(stacked_array)[1]
  n_col <- dim(stacked_array)[2]
  
  quantile_map <- matrix(NA, nrow = n_row, ncol = n_col)
  rownames(quantile_map) <- rownames(corrdiff)
  colnames(quantile_map) <- colnames(corrdiff)
  
  for (i in 1:n_row) {
    for (j in 1:n_col) {
      quantile_map[i, j] <- mean(abs(stacked_array[i, j, ]) >= abs(corrdiff[i, j]))
    }
  }
  
  lt_idx <- which(lower.tri(quantile_map), arr.ind = TRUE)
  
  gene_pairs <- apply(lt_idx, 1, function(idx) {
    paste(rownames(quantile_map)[idx[1]], "_", colnames(quantile_map)[idx[2]], sep = "")
  })
  
  pvals <- quantile_map[lower.tri(quantile_map)]
  names(pvals) <- gene_pairs

  return(pvals)
  
}









