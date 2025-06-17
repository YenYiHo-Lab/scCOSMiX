
library(SingleCellExperiment)
library(dplyr)
library(biomaRt)


# Download GSE266919_Myeloid from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266919

sce_obj <- readRDS(".../GSE266919_Myeloid.rds")


sce_sub <- sce_obj[, (sce_obj$Treatment %in% c("Nab-PTX+Anti-PD-L1", "Nab-PTX")) &
                     (sce_obj$Response == "R") & 
                     (sce_obj$Tissue == "breast") &
                     (sce_obj$Subset %in% c("Macro-CXCL9", "Macro-CCL3L1", "Macro-FOLR2", "Macro-IL1B", "Mono-1L1B", "Macro-ISG15", "Macro-CCL2", "Macro-FTH1"))]

ourtab <- table(sce_sub$Patient, sce_sub$Group)
ourtab

patients_keep <- rownames(ourtab)[(ourtab[,1] > 30) & (ourtab[,2] > 30)]

sce_sub <- sce_sub[, sce_sub$Patient %in% patients_keep]


count.mat <- t(as.matrix(assay(sce_sub, "counts")))

meta.dat <- as.data.frame(colData(sce_sub))


subdf <- cbind(count.mat[rownames(meta.dat), ], meta.dat)
subdf$Patient <- as.factor(subdf$Patient)
subdf$Group <- as.factor(subdf$Group)



colnames(subdf)[(ncol(subdf) - 13):ncol(subdf)]
levels(subdf$Group)
levels(subdf$Patient)


m <- ncol(subdf) - 13
numeric_cols <- names(subdf)[1:m]
factor_cols <- c("Patient", "Group")

zinf_prop <- subdf %>%
  group_by(across(all_of(factor_cols))) %>%
  summarise(across(all_of(numeric_cols), ~ mean(.x > 0)))

zinf_prop_Teff <- zinf_prop[(1:nrow(zinf_prop) %% 2)==1, ]
zinf_prop_Tsem <- zinf_prop[(1:nrow(zinf_prop) %% 2)==0, ]

zinf_count <- subdf %>%
  group_by(across(all_of(factor_cols))) %>%
  summarise(across(all_of(numeric_cols), ~ sum(.x > 0)))

zinf_count_Teff <- zinf_count[(1:nrow(zinf_count) %% 2)==1, ]
zinf_count_Tsem <- zinf_count[(1:nrow(zinf_count) %% 2)==0, ]

zero_count <- subdf %>%
  group_by(across(all_of(factor_cols))) %>%
  summarise(across(all_of(numeric_cols), ~ sum(.x == 0)))

zero_count_Teff <- zero_count[(1:nrow(zero_count) %% 2)==1, ]
zero_count_Tsem <- zero_count[(1:nrow(zero_count) %% 2)==0, ]


combo_satisfies <- (zinf_prop_Teff[, 3:ncol(zinf_prop_Teff)] > 0.3) &
  (zinf_prop_Tsem[, 3:ncol(zinf_prop_Tsem)] > 0.3) &
  (zinf_count_Teff[, 3:ncol(zinf_count_Teff)] > 25) &
  (zinf_count_Tsem[, 3:ncol(zinf_count_Tsem)] > 25) &
  (zero_count_Teff[, 3:ncol(zero_count_Teff)] > 2) &
  (zero_count_Tsem[, 3:ncol(zero_count_Tsem)] > 2)


numpatients_satisfied <- colSums(combo_satisfies)


geneskeep <- names(numpatients_satisfied)[numpatients_satisfied > 5]


subdf <- subdf[, c(geneskeep, colnames(subdf)[(ncol(subdf)-12):ncol(subdf)])]


listMarts()
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")  # Adjust if needed

annot <- getBM(
  attributes = c("hgnc_symbol", "gene_biotype", "ensembl_gene_id"),
  filters    = "hgnc_symbol",
  values     = geneskeep,
  mart       = mart
)

geneskeep <- unique(annot[annot$gene_biotype == "protein_coding", "hgnc_symbol"])

# remove housekeeping genes
hk_genes <- read.table("https://www.tau.ac.il/~elieis/HKG/HK_genes.txt", header = FALSE, stringsAsFactors = FALSE)
hk_genes <- hk_genes$V1

geneskeep <- setdiff(geneskeep, hk_genes)

# second set of housekeeping genes comes from https://housekeeping.unicamp.br/?download
load(".../Housekeeping_GenesHuman.RData")

geneskeep <- setdiff(geneskeep, c(Housekeeping_Genes$Gene.name))


# remove ribosomal genes
geneskeep <- geneskeep[!grepl("^RPL|^RPS", geneskeep)]

# remove mitochondrial genes
geneskeep <- geneskeep[-grep("^MT-", geneskeep)]



subdf <- subdf[, c(geneskeep, colnames(subdf)[(ncol(subdf)-12):ncol(subdf)])]


m <- ncol(subdf) - 13
numeric_cols <- names(subdf)[1:m]
factor_cols <- c("Patient", "Group")

quant_sub <- subdf %>%
  group_by(across(all_of(factor_cols))) %>%
  summarise(across(all_of(numeric_cols), ~ quantile(.x , 0.975)))


quant_sub_Teff <- quant_sub[(1:nrow(quant_sub) %% 2)==1, ]
quant_sub_Tsem <- quant_sub[(1:nrow(quant_sub) %% 2)==0, ]

combo_satisfies <- (quant_sub_Teff[, 3:ncol(quant_sub_Teff)] > 4) &
  (quant_sub_Tsem[, 3:ncol(quant_sub_Tsem)] > 4)



numpatients_satisfied <- colSums(combo_satisfies)


geneskeep <- names(numpatients_satisfied)[numpatients_satisfied > 5]

subdf <- subdf[, c(geneskeep, colnames(subdf)[(ncol(subdf)-12):ncol(subdf)])]

subdf$Group <- ifelse(subdf$Group == "Pre-treatment", "Pre", "Post")

subdf$Group <- as.factor(subdf$Group)
subdf$Group <- relevel(subdf$Group, ref = "Pre")


genepairs <- t(combn(geneskeep, 2))

saveRDS(subdf, ".../data/gse266919_filtered_data.rds")






