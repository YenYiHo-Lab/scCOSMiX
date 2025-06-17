
library(dplyr)

source("C:/Users/abussing/Downloads/scCOSMiX_functions.R")


ourdf <- readRDS("C:/Users/abussing/Downloads/gse108989_filtered_data.rds")


eq1 <- y1 ~ Group + s(Patient, bs="re", by=Group)
eq2 <- y2 ~ Group + s(Patient, bs="re", by=Group)
eq3 <- ~ Group
eq4 <- ~ Group
eq5 <- ~ Group + s(Patient, bs="re", by=Group)
eq6 <- ~ Group + s(Patient, bs="re")
eq7 <- ~ Group + s(Patient, bs="re")
eqlist <- list(eq1, eq2, eq3, eq4, eq5, eq6, eq7)



gene_pairs <- t(combn(colnames(ourdf)[1:276],2))


gene_pair_string <- apply(gene_pairs, 1, function(a) paste(a[1], "_", a[2], sep=""))

retmat <- data.frame(rho_CD4 = rep(0, nrow(gene_pairs)), rho_CD8 = 0, pval = 0)
rownames(retmat) <- gene_pair_string

for (k in 1:nrow(gene_pairs)){
  
  gene1 <- gene_pairs[k, 1]
  gene2 <- gene_pairs[k, 2]
  
  tempdf <- ourdf[, c(gene1, gene2, "Patient_ID", "celltype_group", "nCount_pre")]
  colnames(tempdf) <- c("y1", "y2", "Patient", "Group", "n_counts")
  
  # check which patients we need to drop
  numeric_cols <- names(tempdf)[1:2]
  factor_cols <- c("Patient", "Group")
  
  
  zinf_prop <- tempdf %>%
    group_by(across(all_of(factor_cols))) %>%
    summarise(across(all_of(numeric_cols), ~ mean(.x > 0)))
  
  zinf_prop_CD4 <- zinf_prop[(1:nrow(zinf_prop) %% 2)==1, ]
  zinf_prop_CD8 <- zinf_prop[(1:nrow(zinf_prop) %% 2)==0, ]
  
  zinf_count <- tempdf %>%
    group_by(across(all_of(factor_cols))) %>%
    summarise(across(all_of(numeric_cols), ~ sum(.x > 0)))
  
  zinf_count_CD4 <- zinf_count[(1:nrow(zinf_count) %% 2)==1, ]
  zinf_count_CD8 <- zinf_count[(1:nrow(zinf_count) %% 2)==0, ]
  
  zero_count <- tempdf %>%
    group_by(across(all_of(factor_cols))) %>%
    summarise(across(all_of(numeric_cols), ~ sum(.x == 0)))
  
  zero_count_CD4 <- zero_count[(1:nrow(zero_count) %% 2)==1, ]
  zero_count_CD8 <- zero_count[(1:nrow(zero_count) %% 2)==0, ]
  
  combo_satisfies <- (zinf_prop_CD4[, 3:ncol(zinf_prop_CD4)] > 0.3) &
    (zinf_prop_CD8[, 3:ncol(zinf_prop_CD8)] > 0.3) &
    (zinf_count_CD4[, 3:ncol(zinf_count_CD4)] > 30) &
    (zinf_count_CD8[, 3:ncol(zinf_count_CD8)] > 30) &
    (zero_count_CD4[, 3:ncol(zero_count_CD4)] > 3) &
    (zero_count_CD8[, 3:ncol(zero_count_CD8)] > 3)
  
  patients_keep <- as.character(zinf_prop_CD4$Patient[rowSums(combo_satisfies) == 2])
  
  tempdf <- tempdf[tempdf$Patient %in% patients_keep, ]
  tempdf$Patient <- droplevels(tempdf$Patient)
  
  
  res_out <- sccosmix(formula = eqlist, data=tempdf, offset1 = log(tempdf$n_counts), offset2 = log(tempdf$n_counts))
  
  retmat[gene_pair_string[k], "rho_CD4"] <- tanh(summary(res_out)$tableP5[1,1])
  retmat[gene_pair_string[k], "rho_CD8"] <- tanh(sum(summary(res_out)$tableP5[1:2,1]))
  retmat[gene_pair_string[k], "pval"] <- summary(res_out)$tableP5[2,4]
  
}



retmat$pval_BH <- p.adjust(retmat$pval, method = "BH")


retmat[retmat$pval_BH < 0.05, ]




