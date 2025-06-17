
library(dplyr)

source("C:/Users/abussing/Downloads/scCOSMiX_functions.R")


ourdf <- readRDS("C:/Users/abussing/Downloads/gse266919_filtered_data.rds")


eq1 <- y1 ~ Group + s(Patient, bs="re", by=Group)
eq2 <- y2 ~ Group + s(Patient, bs="re", by=Group)
eq3 <- ~ Group
eq4 <- ~ Group
eq5 <- ~ Group + s(Patient, bs="re", by=Group)
eq6 <- ~ Group + s(Patient, bs="re")
eq7 <- ~ Group + s(Patient, bs="re")
eqlist <- list(eq1, eq2, eq3, eq4, eq5, eq6, eq7)



gene_pairs <- t(combn(colnames(ourdf)[1:251],2))


gene_pair_string <- apply(gene_pairs, 1, function(a) paste(a[1], "_", a[2], sep=""))

retmat <- data.frame(rho_pre = rep(0, nrow(gene_pairs)), rho_post = 0, pval = 0)
rownames(retmat) <- gene_pair_string

for (k in 1:nrow(gene_pairs)){
  print(k)
  
  gene1 <- gene_pairs[k, 1]
  gene2 <- gene_pairs[k, 2]
  
  tempdf <- ourdf[, c(gene1, gene2, "Patient", "Group", "n_counts")]
  colnames(tempdf) <- c("y1", "y2", "Patient", "Group", "n_counts")
  
  numeric_cols <- names(tempdf)[1:2]
  factor_cols <- c("Patient", "Group")
  
  
  zinf_prop <- tempdf %>%
    group_by(across(all_of(factor_cols))) %>%
    summarise(across(all_of(numeric_cols), ~ mean(.x > 0)))
  
  zinf_prop_pre <- zinf_prop[(1:nrow(zinf_prop) %% 2)==1, ]
  zinf_prop_post <- zinf_prop[(1:nrow(zinf_prop) %% 2)==0, ]
  
  zinf_count <- tempdf %>%
    group_by(across(all_of(factor_cols))) %>%
    summarise(across(all_of(numeric_cols), ~ sum(.x > 0)))
  
  zinf_count_pre <- zinf_count[(1:nrow(zinf_count) %% 2)==1, ]
  zinf_count_post <- zinf_count[(1:nrow(zinf_count) %% 2)==0, ]
  
  zero_count <- tempdf %>%
    group_by(across(all_of(factor_cols))) %>%
    summarise(across(all_of(numeric_cols), ~ sum(.x == 0)))
  
  zero_count_pre <- zero_count[(1:nrow(zero_count) %% 2)==1, ]
  zero_count_post <- zero_count[(1:nrow(zero_count) %% 2)==0, ]
  
  combo_satisfies <- (zinf_prop_pre[, 3:ncol(zinf_prop_pre)] > 0.3) &
    (zinf_prop_post[, 3:ncol(zinf_prop_post)] > 0.3) &
    (zinf_count_pre[, 3:ncol(zinf_count_pre)] > 25) &
    (zinf_count_post[, 3:ncol(zinf_count_post)] > 25) &
    (zero_count_pre[, 3:ncol(zero_count_pre)] > 2) &
    (zero_count_post[, 3:ncol(zero_count_post)] > 2)
  
  patients_keep <- as.character(zinf_prop_pre$Patient[rowSums(combo_satisfies) == 2])
  
  
  tempdf <- tempdf[tempdf$Patient %in% patients_keep, ]
  tempdf$Patient <- droplevels(tempdf$Patient)
  
  
  
  res_out <- sccosmix(formula = eqlist, data=tempdf, offset1 = log(tempdf$n_counts), offset2 = log(tempdf$n_counts))

  retmat[gene_pair_string[k], "rho_pre"] <- tanh(summary(res_out)$tableP5[1,1])
  retmat[gene_pair_string[k], "rho_post"] <- tanh(sum(summary(res_out)$tableP5[1:2,1]))
  retmat[gene_pair_string[k], "pval"] <- summary(res_out)$tableP5[2,4]
  
}



retmat$pval_BH <- p.adjust(retmat$pval, method = "BH")


retmat[retmat$pval_BH < 0.05, ]




