
library(Seurat)
library(sctransform)
library(abind)

scrho_pval <- function(ourdf, genecols, metacols, groupcol, logcountcol, B){
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

  umidf_pre <- data.frame(log_umi = ourdf_pre$logn_counts)
  rownames(umidf_pre) <- rownames(ourdf_pre)
  
  umidf_post <- data.frame(log_umi = ourdf_post$logn_counts)
  rownames(umidf_post) <- rownames(ourdf_post)
  
  sct_pre <- vst(
    umi = as(t(as.matrix(ourdf_pre[, genecols])), "dgCMatrix"),
    cell_attr = umidf_pre,
    vst.flavor = "v2",
    verbosity=0
  )
  
  cor_pre <- cor(t(sct_pre$y))
  
  sct_post <- vst(
    umi = as(t(as.matrix(ourdf_post[, genecols])), "dgCMatrix"),
    cell_attr = umidf_post,
    vst.flavor = "v2",
    verbosity=0
  )
  
  cor_post <- cor(t(sct_post$y))
  
  
  corrdiff <- cor_post - cor_pre
  
  
  mat_list <- list()
  for (i in 1:B) {
    # now let's shuffle pre cells and post cells
    ourdf_temp <- ourdf
    ourdf_temp[, groupcol] <- ourdf_temp[sample(1:nrow(ourdf_temp), size=nrow(ourdf_temp)), groupcol]
    
    ourdf_pre_temp <- ourdf_temp[ourdf_temp[, groupcol] == grouplevels[1], ]
    ourdf_post_temp <- ourdf_temp[ourdf_temp[, groupcol] == grouplevels[2], ]
    
    umidf_pre_temp <- data.frame(log_umi = ourdf_pre_temp$logn_counts)
    rownames(umidf_pre_temp) <- rownames(ourdf_pre_temp)
    
    umidf_post_temp <- data.frame(log_umi = ourdf_post_temp$logn_counts)
    rownames(umidf_post_temp) <- rownames(ourdf_post_temp)
    
    sct_pre_temp <- vst(
      umi = as(t(as.matrix(ourdf_pre_temp[, genecols])), "dgCMatrix"),
      cell_attr = umidf_pre_temp,
      vst.flavor = "v2",
      verbosity=0
    )

    cor_pre_temp <- cor(t(sct_pre_temp$y))
    
    sct_post_temp <- vst(
      umi = as(t(as.matrix(ourdf_post_temp[, genecols])), "dgCMatrix"),
      cell_attr = umidf_post_temp,
      vst.flavor = "v2",
      verbosity=0
    )
    
    cor_post_temp <- cor(t(sct_post_temp$y))
  
    corrdiff_temp <- cor_post_temp - cor_pre_temp
    
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










