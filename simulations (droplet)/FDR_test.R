
source(".../functions/scDesign3_functions.R")
source(".../functions/CSCORE_functions.R")
source(".../functions/sctransformrho_functions.R")
source(".../functions/CoCoA_functions.R")
source(".../functions/CNM_functions.R")
library(GJRM)
library(scDECO)

library(SingleCellExperiment)
library(dplyr)
library(Matrix)


methods_use <-  c("CSCORE","sctransform-rho","CoCoA","CNM.full","GJRM","scDECO","scCOSMiX")



bigdf <- readRDS(".../data/gse266919_filtered_data.rds")

bigdf$Group <- factor(as.character(bigdf$Group), levels = c("Pre", "Post"))

table(bigdf$Group, bigdf$Patient)

bigdf$logn_counts <- log(bigdf$n_counts)

bigdf$group_x_patient <- interaction(bigdf$Group, bigdf$Patient, sep = "_")
bigdf$group_x_patient <- as.factor(bigdf$group_x_patient)

colnames(bigdf)

metacols <- colnames(bigdf)[(ncol(bigdf)-14):ncol(bigdf)]


genes_keep <- c("PNRC1", "TGFBI", "MS4A7", "DBI", "TXN", "HLA-DMB", "TSPO", "ARHGDIB", "GPNMB", "FCGR2A", "SGK1", "FCGRT", "MS4A6A", "LITAF", "PKM",     
             "CTSL", "DDX5", "ANXA5", "TMEM176A", "BTG2", "IFITM3", "ITGB2", "LILRB4", "ANXA2", "HSPB1", "IGSF6", "PLEK", "ZFP36", "RASSF4", "SERPING1",
             "SYNGR2", "GLUL", "ASAH1", "ATP5MC2", "LGALS3", "CTSH", "TMEM176B", "SERPINA1", "MAFB", "FCGR3A", "STAT1", "TAGLN2", "FYB1", "TREM2", "SPI1",    
             "CYBB", "LGMN", "JUNB", "GSTP1", "CD68")


gene_pairs <- t(combn(genes_keep, 2))
gene_pairs <- apply(gene_pairs, 1, function(a) paste(a[2], "_", a[1], sep=""))


pats_keep <- c("P038", "P041", "P042", "P049", "P051", "P056", "P057", "P058")



subdf <- bigdf[bigdf$Patient %in% pats_keep, ]
subdf$Patient <- droplevels(subdf$Patient)

subdf$group_x_patient <- droplevels(subdf$group_x_patient)

genesdf <- subdf[, genes_keep]
metadf <- subdf[, metacols]


sce <- SingleCellExperiment(
  assays = list(counts = t(as.matrix(genesdf))),
  colData = metadf
)


BATCH_data <- construct_data(
  sce = sce,
  assay_use = "counts",
  celltype = "Group",
  pseudotime = NULL,
  spatial = NULL,
  other_covariates = c("logn_counts", "group_x_patient", "Patient"),
  corr_by = "group_x_patient"
)


BATCH_marginal_zinf <- fit_marginal(
  data = BATCH_data,
  predictor = "gene",
  mu_formula = "s(Patient, bs='re')", # "Group + Group:Patient + offset(logn_counts)",
  sigma_formula = "1", 
  family_use = "zinb",
  n_cores = 1,
  usebam = FALSE
)


BATCH_marginal <- fit_marginal(
  data = BATCH_data,
  predictor = "gene",
  mu_formula = "Group + Group:Patient + offset(logn_counts)",
  sigma_formula = "Group", 
  family_use = "nb",
  n_cores = 1,
  usebam = FALSE
)


BATCH_data_ignoring_group <- construct_data(
  sce = sce,
  assay_use = "counts",
  celltype = "Group",
  pseudotime = NULL,
  spatial = NULL,
  other_covariates = c("logn_counts", "group_x_patient", "Patient"),
  corr_by = "Patient"
)

for (i in 1:length(BATCH_marginal)){

  BATCH_marginal[[i]]$fit$mu.coefficients <- BATCH_marginal_zinf[[i]]$fit$mu.coefficients
  BATCH_marginal[[i]]$fit$sigma.coefficients <- BATCH_marginal_zinf[[i]]$fit$sigma.coefficients
  BATCH_marginal[[i]]$fit$residuals <- BATCH_marginal_zinf[[i]]$fit$residuals
  
}


BATCH_copula <- fit_copula(
  sce = sce,
  assay_use = "counts",
  marginal_list = BATCH_marginal,
  family_use = "nb",
  copula = "gaussian",
  DT = FALSE,
  n_cores = 1,
  input_data = BATCH_data$dat
)

BATCH_copula_ignoring_group <- fit_copula(
  sce = sce,
  assay_use = "counts",
  marginal_list = BATCH_marginal,
  family_use = "nb",
  copula = "gaussian",
  DT = FALSE,
  n_cores = 1,
  input_data = BATCH_data_ignoring_group$dat
)


example_para <- extract_para(
  sce = sce,
  marginal_list = BATCH_marginal,
  n_cores = 1,
  family_use = "nb",
  new_covariate = BATCH_data$newCovariate,
  data = BATCH_data$dat
)

example_para_zinf <- extract_para(
  sce = sce,
  marginal_list = BATCH_marginal_zinf,
  n_cores = 1,
  family_use = "zinb",
  new_covariate = BATCH_data$newCovariate,
  data = BATCH_data$dat
)


B <- 1 # 100


retmat <- expand.grid(method = methods_use, b = 1:B)
for (k in 1:length(gene_pairs)){
  retmat[[gene_pairs[k]]] <- 0
}

dataset_list <- list()


eq1 <- y1 ~ Group + s(Patient, bs="re", by=Group)
eq2 <- y2 ~ Group + s(Patient, bs="re", by=Group)
eq3 <- ~ Group
eq4 <- ~ Group
eq5 <- ~ Group + s(Patient, bs="re", by=Group)
eq6 <- ~ Group + s(Patient, bs="re")
eq7 <- ~ Group + s(Patient, bs="re")
eqlist <- list(eq1, eq2, eq3, eq4, eq5)


for (b in 1:B){
  
  for (i in 1:length(BATCH_copula$copula_list)){
    
    BATCH_copula$copula_list[[i]] <- BATCH_copula_ignoring_group$copula_list[[ceiling(i/2)]] + matrix(rnorm(50*50, 0, 0.05), nrow=50, ncol=50)
    
  }
  
  for (i in 1:length(BATCH_copula$copula_list)){
    BATCH_copula$copula_list[[i]] <- as.matrix(nearPD(BATCH_copula$copula_list[[i]], corr=TRUE)$mat)
  }
  
  
  new_dat <- simu_new(
    sce=sce,
    assay_use = "counts",
    mean_mat=example_para$mean_mat,
    sigma_mat=example_para$sigma_mat,
    zero_mat=example_para$zero_mat,
    quantile_mat = NULL,
    copula_list=BATCH_copula$copula_list,
    n_cores = 1,
    family_use = "nb",
    input_data = BATCH_data$dat,
    new_covariate = BATCH_data$newCovariate,
    filtered_gene = BATCH_data$filtered_gene
  )
  
  new_dat <- as.data.frame(t(new_dat))
  
  new_dat$Group <- BATCH_data$newCovariate$Group
  new_dat$Patient <- BATCH_data$newCovariate$Patient
  new_dat$logn_counts <- BATCH_data$newCovariate$logn_counts
  
  new_dat[, 1:50] <- new_dat[, 1:50]*rbinom(n=length(example_para_zinf$zero_mat), size=1, prob=1-as.vector(example_para_zinf$zero_mat))

  dataset_list[[paste("rep_", b, sep="")]] <- new_dat

  for (method_idx in 1:(length(methods_use)-1)){

    if (methods_use[method_idx] == "CSCORE"){

      what <- csCORE_pval(
        ourdf = new_dat,
        genecols = colnames(new_dat)[1:50],
        metacols = c("Group", "Patient", "logn_counts"),
        groupcol = "Group",
        logcountcol = "logn_counts",
        B = 100
      )

      retmat[(retmat$method == "CSCORE") & (retmat$b == b), names(what)] <- what

    }

    if (methods_use[method_idx] == "sctransform-rho"){

      what <- scrho_pval(
        ourdf = new_dat,
        genecols = colnames(new_dat)[1:50],
        metacols = c("Group", "Patient", "logn_counts"),
        groupcol = "Group",
        logcountcol = "logn_counts",
        B = 100
      )

      retmat[(retmat$method == "sctransform-rho") & (retmat$b == b), names(what)] <- what

    }

    if (methods_use[method_idx] == "CoCoA"){

      tempvec <- rep(0, length(gene_pairs))
      names(tempvec) <- gene_pairs

      for (k in 1:length(gene_pairs)){

        gene1 <- strsplit(gene_pairs[k], "_", 2)[[1]][1]
        gene2 <- strsplit(gene_pairs[k], "_", 2)[[1]][2]

        tempdf <- new_dat[, c(gene1, gene2, "Group", "Patient", "logn_counts")]
        colnames(tempdf) <- c("y1", "y2", "Group", "Patient", "logn_counts")

        tempdf$Group_bern <- 1*(tempdf$Group == "Post")

        tempdf$y1.norm <- gamlss::gamlss(formula = y1 ~ Group + offset(logn_counts),
                                         sigma.formula = ~ Group,
                                         data = tempdf, family = gamlss.dist::NBI, trace=FALSE)$residuals

        tempdf$y2.norm <- gamlss::gamlss(formula = y2 ~ Group + offset(logn_counts),
                                         sigma.formula = ~ Group,
                                         data = tempdf, family = gamlss.dist::NBI, trace=FALSE)$residuals

        tempdf$Tpar <- NA
        tempdf$ZTpar <- NA

        cocoa_df = data.frame(X = tempdf$y1.norm, Y = tempdf$y2.norm, Z=tempdf$Group_bern)


        what <- get_params_mle(dat_xyz = cocoa_df)

        tempvec[gene_pairs[k]] <- what$out2$p[nrow(what$out2)]

      }

      retmat[(retmat$method == "CoCoA") & (retmat$b == b), names(tempvec)] <- tempvec

    }

    if (methods_use[method_idx] == "CNM.full"){

      tempvec <- rep(0, length(gene_pairs))
      names(tempvec) <- gene_pairs

      for (k in 1:length(gene_pairs)){

        gene1 <- strsplit(gene_pairs[k], "_", 2)[[1]][1]
        gene2 <- strsplit(gene_pairs[k], "_", 2)[[1]][2]

        tempdf <- new_dat[, c(gene1, gene2, "Group", "Patient", "logn_counts")]
        colnames(tempdf) <- c("y1", "y2", "Group", "Patient", "logn_counts")

        tempdf$Group_bern <- 1*(tempdf$Group == "Post")

        tempdf$y1.norm <- gamlss::gamlss(formula = y1 ~ Group + offset(logn_counts),
                                         sigma.formula = ~ Group,
                                         data = tempdf, family = gamlss.dist::NBI, trace=FALSE)$residuals

        tempdf$y2.norm <- gamlss::gamlss(formula = y2 ~ Group + offset(logn_counts),
                                         sigma.formula = ~ Group,
                                         data = tempdf, family = gamlss.dist::NBI, trace=FALSE)$residuals

        tempdf$rownum <- 1:nrow(tempdf)

        ourdat_long <- tempdf %>%
          pivot_longer(cols = c("y1.norm", "y2.norm"),
                       names_to = "Margin",
                       values_to = "Value") %>%
          mutate(Margin = ifelse(Margin == "y1.norm", "Margin 1", "Margin 2"))

        zcor_matrix <- model.matrix(~ Group_bern, data = tempdf)

        cnm_res <- summary(geese(Value ~ Margin,
                                 id = rownum,
                                 data = ourdat_long,
                                 family = gaussian,
                                 corstr = "exchangeable",
                                 sformula = ~ Margin,
                                 zcor = zcor_matrix,
                                 cor.link = "fisherz"))$correlation

        tempvec[gene_pairs[k]] <- cnm_res$p[2]

      }

      retmat[(retmat$method == "CNM.full") & (retmat$b == b), names(tempvec)] <- tempvec

    }


    if (methods_use[method_idx] == "GJRM"){

      tempvec <- rep(0, length(gene_pairs))
      names(tempvec) <- gene_pairs

      for (k in 1:length(gene_pairs)){

        print(k)

        gene1 <- strsplit(gene_pairs[k], "_", 2)[[1]][1]
        gene2 <- strsplit(gene_pairs[k], "_", 2)[[1]][2]

        tempdf <- new_dat[, c(gene1, gene2, "Group", "Patient", "logn_counts")]
        colnames(tempdf) <- c("y1", "y2", "Group", "Patient", "logn_counts")

        what <- gjrm(formula = eqlist, data=tempdf, margins = c("NBI", "NBI"), model = "B", offset1 = tempdf$logn_counts, offset2 = tempdf$logn_counts)

        tempvec[gene_pairs[k]] <- summary(what)$tableP5[2,4]

      }

      retmat[(retmat$method == "GJRM") & (retmat$b == b), names(tempvec)] <- tempvec

    }

    if (methods_use[method_idx] == "scDECO"){

      tempvec <- rep(0, length(gene_pairs))
      names(tempvec) <- gene_pairs

      for (k in 1:length(gene_pairs)){
        print(k)
        gene1 <- strsplit(gene_pairs[k], "_", 2)[[1]][1]
        gene2 <- strsplit(gene_pairs[k], "_", 2)[[1]][2]

        tempdf <- new_dat[, c(gene1, gene2, "Group", "Patient", "logn_counts")]
        colnames(tempdf) <- c("y1", "y2", "Group", "Patient", "logn_counts")

        mcmc.out <- scdeco.cop(y=tempdf[, c("y1", "y2")], x=1*(tempdf$Group == "Post"), marginals=c("ZINB", "ZINB"),
                               n.mcmc=100, burn=0, thin=1, offset1 = tempdf$logn_counts, offset2 = tempdf$logn_counts) # update offset to 6000, burn to 1000 for real test

        pval <- 2*(min(mean(mcmc.out[,"tau1"] < 0), mean(mcmc.out[, "tau1"] > 0)))

        tempvec[gene_pairs[k]] <- pval

      }

      retmat[(retmat$method == "scDECO") & (retmat$b == b), names(tempvec)] <- tempvec

    }
    
  }
    
}
    
    
    



library(stringr)

source(".../functions/scCOSMiX_functions.R")

eqlist <- list(eq1, eq2, eq3, eq4, eq5, eq6, eq7)

for (j in 1:length(dataset_list)){
  
  new_dat <- dataset_list[[j]]
  
  dataset_name <- names(dataset_list)[j]
  
  print(dataset_name)
  
  b <- as.integer(str_extract_all(dataset_name, "\\d+")[[1]])
  
  tempvec <- rep(0, length(gene_pairs))
  names(tempvec) <- gene_pairs
  
  for (k in 1:length(gene_pairs)){
    
    print(k)
    
    gene1 <- strsplit(gene_pairs[k], "_", 2)[[1]][1]
    gene2 <- strsplit(gene_pairs[k], "_", 2)[[1]][2]
    
    tempdf <- new_dat[, c(gene1, gene2, "Group", "Patient", "logn_counts")]
    colnames(tempdf) <- c("y1", "y2", "Group", "Patient", "logn_counts")
    
    what <- sccosmix(formula = eqlist, data=tempdf, offset1 = tempdf$logn_counts, offset2 = tempdf$logn_counts)
    
    tempvec[gene_pairs[k]] <- summary(what)$tableP5[2,4]
    
  }
  
  retmat[(retmat$method == "scCOSMiX") & (retmat$b == b), names(tempvec)] <- tempvec

}


retmat



