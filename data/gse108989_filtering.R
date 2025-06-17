
library(dplyr)
library(SingleCellExperiment)
library(biomaRt)
library(Seurat)

url <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc=GSE108989&id=166188&db=GeoDb_blob163"
download.file(url, destfile = "C:/Users/abussing/Downloads/GSE108989_metadata.txt", mode = "wb")

# then we clean it up in excel and save as csv

meta.data <- read.csv(".../GSE108989_metadata.csv", header = TRUE)

# there is a UniqueCell_ID called "" which we should remove
meta.data <- meta.data[!(meta.data$UniqueCell_ID %in% unique(meta.data$UniqueCell_ID[duplicated(meta.data$UniqueCell_ID)])), ]

rownames(meta.data) <- meta.data$UniqueCell_ID

meta.data <- meta.data[, -1]

# download from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108989
count.data <- read.csv(".../GSE108989_CRC.TCell.S11138.count.csv")

# remove duplicated gene names (there are only 3 and they are meaningless garble)
count.data <- count.data[!(count.data$symbol %in% unique(count.data$symbol[duplicated(count.data$symbol)])), ]

rownames(count.data) <- count.data$symbol

count.data <- count.data[, -c(1,2)]


# meta.data uses - but count.data uses .    let's change them both to _
rownames(meta.data) <- gsub("-", "_", rownames(meta.data))
colnames(count.data) <- gsub("\\.", "_", colnames(count.data))


setdiff(rownames(meta.data), colnames(count.data))


count.data <- t(count.data)

bigdf <- cbind(count.data, meta.data[rownames(count.data), ])


# let's cut out the cell and gene names which the author cut in their curated (normalized) data
# download from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108989
author_corrected_data <- read.table(".../GSE108989_CRC.TCell.S10805.norm.centered.txt", 
                                    header = TRUE, sep = "\t", stringsAsFactors = FALSE)


colnames(author_corrected_data) <- gsub("\\.", "_", colnames(author_corrected_data))

# they have NA as a gene name. Let's take that out
author_corrected_data <- author_corrected_data[-which(author_corrected_data$geneSymbol %in% unique(author_corrected_data$geneSymbol[duplicated(author_corrected_data$geneSymbol)])), ]


rownames(author_corrected_data) <- author_corrected_data$geneSymbol


bigdf$nCount_pre <- rowSums(bigdf[, 1:(ncol(bigdf)-3)])

# # remove those huge nCounts
bigdf <- bigdf[bigdf$nCount < 2.35e06, ]

bigdf <- bigdf[intersect(rownames(bigdf), colnames(author_corrected_data)), ]

bigdf <- bigdf[, c(intersect(colnames(bigdf), rownames(author_corrected_data)), "Patient_ID", "majorCluster", "sampleType", "nCount_pre")]




# subset down to just Tumor cells
subdf <- bigdf[(bigdf$sampleType %in% c("TP7", "TTC", "TTH", "TTR", "TTY")), ]

table(subdf[, c("Patient_ID", "majorCluster")])


# we must remove patient P1207 because they have no CD4 cells
subdf <- subdf[subdf$Patient_ID != "P1207", ]
subdf$Patient_ID <- as.factor(subdf$Patient_ID)
subdf$Patient_ID <- droplevels(subdf$Patient_ID)

# break it down into CD4 and CD8 Tcells
subdf$celltype_group <- "ok"
subdf[grepl("^CD4", subdf$majorCluster), "celltype_group"] <- "CD4"
subdf[grepl("^CD8", subdf$majorCluster), "celltype_group"] <- "CD8"

table(subdf$celltype_group)

subdf$celltype_group <- as.factor(subdf$celltype_group)

subdf <- subdf[subdf$celltype_group != "ok", ]

subdf$celltype_group <- droplevels(subdf$celltype_group)

table(subdf[, c("celltype_group", "Patient_ID")])

colnames(subdf)[(ncol(subdf) - 5):ncol(subdf)]
levels(subdf$celltype_group)
levels(subdf$Patient_ID)




# Extract gene expression matrix (all columns except the last 5)
gene_counts <- t(subdf[, 1:(ncol(subdf) - 5)])
rownames(gene_counts) <- colnames(subdf)[1:(ncol(subdf) - 5)]  # Set gene names
colnames(gene_counts) <- rownames(subdf)  # Set cell barcodes as column names

# Extract metadata (last 5 columns)
metadata <- subdf[, (ncol(subdf) - 4):ncol(subdf)]
rownames(metadata) <- colnames(gene_counts)  # Ensure rownames match cell barcodes

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = gene_counts, meta.data = metadata)

# Overwrite nCount_RNA with pre-existing values from metadata
seurat_obj$nCount_RNA <- metadata$nCount_pre


levels(seurat_obj@meta.data$Patient_ID)
levels(seurat_obj@meta.data$celltype_group)


seurat_obj <- SCTransform(seurat_obj, return.only.var.genes = FALSE, verbose = TRUE)

seurat_df <- t(as.data.frame(seurat_obj@assays$SCT$scale.data))
seurat_df <- cbind(seurat_df, seurat_obj@meta.data)








m <- ncol(subdf) - 5
numeric_cols <- names(subdf)[1:m]
factor_cols <- c("Patient_ID", "celltype_group")


# we will require 9/11 of the patients to have less than 70% zinf for both cell types
# furthermore, we will require the overall zinf of the gene be under 70%
# furthermore, we will require 9/11 of the patients to have more than 30 obs for both cell types

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
  (zinf_count_Teff[, 3:ncol(zinf_count_Teff)] > 30) &
  (zinf_count_Tsem[, 3:ncol(zinf_count_Tsem)] > 30) &
  (zero_count_Teff[, 3:ncol(zero_count_Teff)] > 3) &
  (zero_count_Tsem[, 3:ncol(zero_count_Tsem)] > 3)

numpatients_satisfied <- colSums(combo_satisfies)


geneskeep <- names(numpatients_satisfied)[numpatients_satisfied > 8]






# now take out housekeeping etc genes

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")



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



subdf <- subdf[, c(geneskeep, colnames(subdf)[(ncol(subdf)-4):ncol(subdf)])]






seurat_df <- seurat_df[, c(geneskeep, c("Patient_ID", "majorCluster", "sampleType", "nCount_pre", "celltype_group"))]
colnames(seurat_df)[colnames(seurat_df)=="nCount_pre"] <- "nCount"


# now check the variance across patients of the different genes

m <- ncol(seurat_df) - 5
numeric_cols <- names(seurat_df)[1:m]
factor_cols <- c("Patient_ID", "celltype_group")

patvar <- seurat_df %>%
  group_by(across(all_of(factor_cols))) %>%
  summarise(across(all_of(numeric_cols), ~ sd(.x)))

patvar_CD4 <- patvar[(1:nrow(patvar) %% 2)==1, ]
patvar_CD8 <- patvar[(1:nrow(patvar) %% 2)==0, ]

var_of_var_CD4 <- apply(patvar_CD4[, 3:m], 2, FUN = function(x) sd(x))
var_of_var_CD8 <- apply(patvar_CD8[, 3:m], 2, FUN = function(x) sd(x))

mean_of_var_CD4 <- apply(patvar_CD4[, 3:m], 2, FUN = function(x) mean(x))
mean_of_var_CD8 <- apply(patvar_CD8[, 3:m], 2, FUN = function(x) mean(x))

mean_of_var_diff <- abs(mean_of_var_CD4 - mean_of_var_CD8)


hist(var_of_var_CD4, breaks=50)
hist(var_of_var_CD8, breaks=50)
hist(mean_of_var_diff, breaks=50)


# filter genes who have var_of_var < 0.4 and have mean_of_var_diff > 0.1
whichgenes_seurat <- names(var_of_var_CD4)[(var_of_var_CD4 < 0.4) & (var_of_var_CD8 < 0.4) & (mean_of_var_diff > 0.1)]


subdf <- subdf[, c(whichgenes_seurat, colnames(subdf)[(ncol(subdf)-4):ncol(subdf)])]


saveRDS(subdf, ".../data/gse108989_filtered_data.rds")


