# Merging SoupX corrected, final gene model datasets

# Incorporating with all other datasets.

# I can do both uncorrected and corrected datasets
# Uncorrected

library(scater)
library(Seurat)
library(monocle) # This is monocle 2.14.0
library(monocle3)
library(expss)
library(dplyr)

setwd("~/Dropbox (VU Basic Sciences)/Miller lab/10X Genomics")

pan.1.c <- readRDS("./final_gene_model/1806-ST-1/SoupX/020620_pan_1_sce_new_SoupX_corrected.rds")
pan.2.c <- readRDS("./final_gene_model/1806-ST-2/SoupX/020620_pan_2_sce_new_SoupX_corrected.rds")
unc86.c <- readRDS("./final_gene_model/3697-ST-1/SoupX/020620_unc86_sce_new_SoupX_corrected.rds")
unc47.1.c <- readRDS("./final_gene_model/2658-ST-1/SoupX/020620_unc_47_1_sce_new_SoupX_corrected.rds")
unc47.2.c <- readRDS("./final_gene_model/3131-ST-1/SoupX/020620_unc_47_2_sce_new_SoupX_corrected.rds")
nmr1.c <- readRDS("./final_gene_model/2966-ST-1/SoupX/020620_nmr_1_sce_new_SoupX_corrected.rds")
eat4.c <- readRDS("./final_gene_model/3070-ST-1/SoupX/020620_eat_4_sce_new_SoupX_corrected.rds")
unc3.c <- readRDS("./final_gene_model/3183-ST-1/SoupX/020620_unc_3_sce_new_SoupX_corrected.rds")
ift20.c <- readRDS("./final_gene_model/3239-ST-1/SoupX/020620_ift_20_sce_new_SoupX_corrected.rds")
cho1.1.c <- readRDS("./final_gene_model/3441-ST-1/SoupX/020620_cho_1_1_sce_new_SoupX_corrected.rds")
cho1.2.c <- readRDS("./final_gene_model/3441-ST-2/SoupX/020620_cho_1_2_sce_new_SoupX_corrected.rds")
tph1.c <- readRDS("./final_gene_model/3465-ST-1/SoupX/020620_tph_1_sce_new_SoupX_corrected.rds")
acr2.c <- readRDS("./final_gene_model/3495-ST-1/SoupX/020620_acr_2_sce_new_SoupX_corrected.rds")
ceh34.c <- readRDS("./final_gene_model/3503-ST-1/SoupX/020620_ceh_34_sce_new_SoupX_corrected.rds")
unc53.c <- readRDS("./final_gene_model/3831-ST-1/SoupX/020620_unc_53_sce_new_SoupX_corrected.rds")
nlp13.c <- readRDS("./final_gene_model/4138-ST-1/SoupX/020620_nlp_13_sce_new_SoupX_corrected.rds")
ceh28.c <- readRDS("./final_gene_model/4170-ST-1/SoupX2/020620_ceh_28_sce_new_SoupX_corrected.rds")

pan.1.sc <- as.Seurat(pan.1.c)
pan.2.sc <- as.Seurat(pan.2.c)
unc86.sc <- as.Seurat(unc86.c)
unc47.1.sc <- as.Seurat(unc47.1.c)
unc47.2.sc <- as.Seurat(unc47.2.c)
nmr1.sc <- as.Seurat(nmr1.c)
eat4.sc <- as.Seurat(eat4.c)
unc3.sc <- as.Seurat(unc3.c)
ift20.sc <- as.Seurat(ift20.c)
cho1.1.sc <- as.Seurat(cho1.1.c)
cho1.2.sc <- as.Seurat(cho1.2.c)
tph1.sc <- as.Seurat(tph1.c)
acr2.sc <- as.Seurat(acr2.c)
ceh34.sc <- as.Seurat(ceh34.c)
unc53.sc <- as.Seurat(unc53.c)
nlp13.sc <- as.Seurat(nlp13.c)
ceh28.sc <- as.Seurat(ceh28.c)

pan.1.sc$Experiment <- "Pan-1"
pan.2.sc$Experiment <- "Pan-2"
unc86.sc$Experiment <- "unc-86"
unc53.sc$Experiment <- "unc-53"
unc47.2.sc$Experiment <- "unc-47_2"
unc47.1.sc$Experiment <- "unc-47_1"
unc3.sc$Experiment <- "unc-3"
tph1.sc$Experiment <- "tph-1_ceh-10"
nmr1.sc$Experiment <- "nmr-1"
ift20.sc$Experiment <- "ift-20"
eat4.sc$Experiment <- "eat-4"
cho1.2.sc$Experiment <- "cho-1_2"
cho1.1.sc$Experiment <- "cho-1_1"
ceh34.sc$Experiment <- "ceh-34"
acr2.sc$Experiment <- "acr-2"
nlp13.sc$Experiment <- "nlp-13_ceh-2"
ceh28.sc$Experiment <- "ceh-28_dat-1"

Idents(pan.1.sc) <- pan.2.sc@meta.data$Experiment
Idents(pan.2.sc) <- pan.1.sc@meta.data$Experiment
Idents(unc86.sc) <- unc86.sc@meta.data$Experiment
Idents(unc47.2.sc) <- unc47.2.sc@meta.data$Experiment
Idents(unc47.1.sc) <- unc47.1.sc@meta.data$Experiment
Idents(unc53.sc) <- unc53.sc@meta.data$Experiment
Idents(unc3.sc) <- unc3.sc@meta.data$Experiment
Idents(tph1.sc) <- tph1.sc@meta.data$Experiment
Idents(nmr1.sc) <- nmr1.sc@meta.data$Experiment
Idents(ift20.sc) <- ift20.sc@meta.data$Experiment
Idents(eat4.sc) <- eat4.sc@meta.data$Experiment
Idents(cho1.2.sc) <- cho1.2.sc@meta.data$Experiment
Idents(cho1.1.sc) <- cho1.1.sc@meta.data$Experiment
Idents(ceh34.sc) <- ceh34.sc@meta.data$Experiment
Idents(acr2.sc) <- acr2.sc@meta.data$Experiment
Idents(nlp13.sc) <- nlp13.sc@meta.data$Experiment
Idents(ceh28.sc) <- ceh28.sc@meta.data$Experiment

# Making cell barcode names consistent with previous versions
unc47.1.sc@meta.data$Barcode <- paste(unc47.1.sc@meta.data$Barcode, "GABA", sep = "_")
pan.1.sc@meta.data$Barcode <- paste("Pan", pan.1.sc@meta.data$Barcode, "1", sep = "_")
pan.2.sc@meta.data$Barcode <- paste("Pan", pan.2.sc@meta.data$Barcode, "2", sep = "_")
nmr1.sc@meta.data$Barcode <- paste(nmr1.sc@meta.data$Barcode, "nmr1", sep = "_")

# Merging all experiments into one object
whole_L4.sc <- merge(pan.1.sc, y = c(pan.2.sc, unc47.1.sc, nmr1.sc, unc47.2.sc, eat4.sc, ift20.sc, unc3.sc, tph1.sc, cho1.1.sc, cho1.2.sc, acr2.sc, ceh34.sc, unc86.sc, unc53.sc, nlp13.sc, ceh28.sc), add.cell.ids = c('Pan-1', "Pan-2", 'GABA', "nmr1", "G2", "eat4", "ci", "u3", "tph1_ceh10", 'cho1-1', "cho1-2", "acr2", 'ceh34', "unc86", "unc53", "nc", "c28"))

table(whole_L4.sc$pct_counts_Mito <= 20)
# FALSE  TRUE 
# 5191 105613 

rownames(whole_L4.sc@meta.data) <- ifelse(
  whole_L4.sc@meta.data$Experiment == "Pan-1" | whole_L4.sc@meta.data$Experiment == "Pan-2",
  as.character(whole_L4.sc@meta.data$Barcode),
  as.character(rownames(whole_L4.sc@meta.data))
)

rownames(whole_L4.sc@meta.data) <- ifelse(
  whole_L4.sc@meta.data$Experiment == "nmr-1",
  paste(rownames(whole_L4.sc@meta.data), 'nmr1', sep = "_"),
  as.character(rownames(whole_L4.sc@meta.data))
)

rownames(whole_L4.sc@meta.data) <- ifelse(
  whole_L4.sc@meta.data$Experiment == "unc-47_1",
  paste(rownames(whole_L4.sc@meta.data), 'GABA', sep = "_"),
  as.character(rownames(whole_L4.sc@meta.data))
)

# Loading in annotations from previous versions of the data
whole_L4.sc@meta.data$Neuron.type <- vlookup(whole_L4.sc@meta.data$Barcode, all.LUT, result_column = 2, lookup_column = 1)

DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size=1000e6)

# Generating monocle3 alpha object for finding variable genes

data.sc <- GetAssayData(whole_L4.sc, slot = "counts")

data.sc[1:6,1:6]

pD <- new("AnnotatedDataFrame", data = whole_L4.sc@meta.data)
gene_list <- read.table("./final_gene_model/genes.txt")
head(gene_list)
fData.sc <- data.frame(row.names = rownames(data.sc), id = rownames(data.sc), gene_short_name = vlookup(rownames(data.sc), gene_list, result_column = 2, lookup_column = 1))
head(fData.sc)
fD <- new("AnnotatedDataFrame", data = fData.sc)

head(rownames(pD))
head(colnames(data.sc))
colnames(data.sc) <- rownames(pD)
data.sc[1:6,1:6]
L4.sc <- newCellDataSet(data.sc, phenoData = pD, featureData = fD)
library(monocle)
L4.sc
# 46911 features, 110804 samples
head(Biobase::fData(L4.sc))
L4.sc <- detectGenes(L4.sc)
table(Biobase::pData(L4.sc)$pct_counts_Mito <= 20)
# FALSE TRUE
#  5191 105613

table(Biobase::fData(L4.sc)$num_cells_expressed > 0)
# FALSE TRUE
# 24372 22539

L4.20.sc <- L4.sc[, Biobase::pData(L4.sc)$pct_counts_Mito <= 20]
L4.20.sc <- detectGenes(L4.20.sc)
L4.20.sc <- L4.20.sc[Biobase::fData(L4.20.sc)$num_cells_expressed > 0,]
L4.20.sc
# 22525 features, 105612 samples

L4.20.sc <- estimateSizeFactors(L4.20.sc)
L4.20.sc <- estimateDispersions(L4.20.sc)

dT <- dispersionTable(L4.20.sc)
table(dT$mean_expression > 0.001 & dT$dispersion_empirical > 50)
# FALSE TRUE
# 13822  7691 

clust_genes <- subset(dT, mean_expression > 0.001 & dispersion_empirical > 50)
L4.20.sc <- setOrderingFilter(L4.20.sc, clust_genes$gene_id)
ambient <- readRDS("./SoupX/112519_ambient_genes_nlp-13_ceh-28.rds")
ambient.L4 <- intersect(ambient, Biobase::fData(L4.20.sc)$id)
fData(L4.20.sc)[ambient.L4, 4] <- FALSE
clust_genes_ids <- as.character(Biobase::fData(L4.20.sc)[which(Biobase::fData(L4.20.sc)$use_for_ordering == TRUE),]$id)

plot_ordering_genes(L4.20.sc)

saveRDS(L4.20.sc, "./final_gene_model/021320_L4_all_cells_g20_Mito_SoupX2_corrected_monocle3alpha_cds.rds")
saveRDS(L4.sc, "./final_gene_model/021320_L4_all_EmptyDrops_cells_ceh28_SoupX2_monocle3alpha_cds.rds")

# Monocle3 object construction
L4.20.m3 <- new_cell_data_set(L4.20.sc@assayData$exprs, cell_metadata = Biobase::pData(L4.20.sc), gene_metadata = Biobase::fData(L4.20.sc))

L4.20.m3 <- preprocess_cds(L4.20.m3, method = "PCA", num_dim = 135, use_genes = clust_genes_ids)
monocle3::plot_pc_variance_explained(L4.20.m3)

L4.20.m3 <- reduce_dimension(L4.20.m3, reduction_method = "UMAP", preprocess_method = "PCA", umap.min_dist = 0.3, umap.n_neighbors = 75)

plot_cells(L4.20.m3, color_cells_by = "Neuron.type", cell_size = 0.7, group_label_size = 3)

# Here I realized that I still have doublets included.  I will use the old cds that had doublets annotated to create a new LUT

table(Biobase::pData(L4.20.sc)$Neuron.type)
table(Biobase::pData(L4.sc)$Neuron.type)

old.LUT <- data.frame(Barcode = Biobase::pData(L4.sc.old)$Barcode, Neuron.type = Biobase::pData(L4.sc.old)$Neuron.type)
library(dplyr)
new.LUT <- full_join(all.LUT, old.LUT, by = "Barcode")

new.LUT$Neuron.type <- ifelse(
  is.na(new.LUT$Neuron.type.x), 
  as.character(new.LUT$Neuron.type.y),
  as.character(new.LUT$Neuron.type.x)
)

table(new.LUT$Neuron.type, exclude = NULL)

new.LUT$Neuron.type.x <- NULL
new.LUT$Neuron.type.y <- NULL

colData(L4.20.m3)$Neuron.type <- as.character(vlookup(colData(L4.20.m3)$Barcode, new.LUT))
table(colData(L4.20.m3)$Neuron.type, exclude = NULL)

Non.neurons <- as.character(unique(colData(L4.20.m3.olds)[which(colData(L4.20.m3.olds)$Tissue == "Non-neuronal"),]$Neuron.type))
Doublets <- c("Coelomocyte_doublets", "Doublets", "Epidermal_doublets", "GABA_doublet", "Muscle_doublets", "RIA_AUA_doublets", "RIA_doublets?", "RIA_M1_doublets", "RIC_doublets", "RIC_RIA doublets")
table(colData(L4.20.m3)$Tissue)
colData(L4.20.m3)$Tissue <- ifelse(
  colData(L4.20.m3)$Neuron.type %in% Non.neurons,
  "Non-neuronal",
  "Neuron"
)

colData(L4.20.m3)[which(colData(L4.20.m3)$Neuron.type %in% Doublets),]$Tissue <- "Doublet"

table(colData(L4.20.m3)$Tissue, exclude = NULL)

saveRDS(L4.20.m3, "./final_gene_model/021320_L4_all_cells_g20_Mito_with_doublets_monocle3_cds.rds")

# Below is sample code for how neuron type annotations were made for the last added dataset
# based on clustering with previously annotated cells.

# Because the UMAP dimensionality reduction is not identical each time it is run, the cluster and partition
# numbers as shown below will not necessarily be reproducible. Code is shown for illustrative purposes.

# Removing annotated doublets.
Doublets <- c("Coelomocyte_doublets", "Doublets", "Epidermal_doublets", "GABA_doublet", "Muscle_doublets", "RIA_AUA_doublets", "RIA_doublets?", "RIA_M1_doublets", "RIC_doublets", "RIC_RIA doublets")
table(colData(L4.20.m3)$Tissue)
colData(L4.20.m3)$Tissue <- ifelse(
  colData(L4.20.m3)$Neuron.type %in% Non.neurons,
  "Non-neuronal",
  "Neuron"
)

colData(L4.20.m3)[which(colData(L4.20.m3)$Neuron.type %in% Doublets),]$Tissue <- "Doublet"

L4.20.m3 <- L4.20.m3[,colData(L4.20.m3)$Tissue %in% c("Neuron", "Non-neuronal")]
L4.20.m3 <- detect_genes(L4.20.m3)
L4.20.m3 <- L4.20.m3[rowData(L4.20.m3)$num_cells_expressed > 0,]
L4.20.m3
# 22469 features, 101147 cells

L4.20.m3 <- preprocess_cds(L4.20.m3, method = "PCA", num_dim = 135, use_genes = clust_genes_ids)
monocle3::plot_pc_variance_explained(L4.20.m3)
L4.20.m3 <- reduce_dimension(L4.20.m3, reduction_method = "UMAP", preprocess_method = "PCA", umap.min_dist = 0.3, umap.n_neighbors = 75)
L4.20.m3 <- cluster_cells(L4.20.m3, resolution = 3e-3, k = 18)

# Storing the UMAP coordinates, partitions and clusters in the cell metadata
colData(L4.20.m3)$UMAP_1 <- L4.20.m3@reducedDims$UMAP[,1]
colData(L4.20.m3)$UMAP_2 <- L4.20.m3@reducedDims$UMAP[,2]

colData(L4.20.m3)$partition <- as.character(monocle3::partitions(L4.20.m3))
colData(L4.20.m3)$UMAP.clusters <- as.character(monocle3::clusters(L4.20.m3))

plot_cells(L4.20.m3, color_cells_by = "Neuron.type", cell_size = 0.7, group_label_size = 3, label_groups_by_cluster = FALSE)
plot_cells(L4.20.m3, color_cells_by = "partition", cell_size = 0.7, group_label_size = 3, label_groups_by_cluster = FALSE)
plot_cells(L4.20.m3, color_cells_by = "cluster", cell_size = 0.7, group_label_size = 3)

# Visualizing the location of cells without annotations
ggplot(as.data.frame(colData(L4.20.m3)), aes(UMAP_1, UMAP_2)) + 
  geom_point(aes(color = ifelse(is.na(Neuron.type), "red", "blue")), size = 0.7) + 
  theme(legend.position = "none")

# Setting annotations based on the already annotated cells in a partition and/or cluster
colData(L4.20.m3)$Neuron.type <- ifelse(
  is.na(colData(L4.20.m3)$Neuron.type) & colData(L4.20.m3)$partition == "25",
  "Pharyngeal_muscle",
  as.character(colData(L4.20.m3)$Neuron.type)
)

colData(L4.20.m3)$Neuron.type <- ifelse(
  is.na(colData(L4.20.m3)$Neuron.type) & colData(L4.20.m3)$partition == "14",
  "ASJ",
  as.character(colData(L4.20.m3)$Neuron.type)
)

colData(L4.20.m3)$Neuron.type <- ifelse(
  is.na(colData(L4.20.m3)$Neuron.type) & colData(L4.20.m3)$partition == "18",
  "ADF",
  as.character(colData(L4.20.m3)$Neuron.type)
)

colData(L4.20.m3)$Neuron.type <- ifelse(
  is.na(colData(L4.20.m3)$Neuron.type) & colData(L4.20.m3)$partition == "45",
  "Excretory_cell",
  as.character(colData(L4.20.m3)$Neuron.type)
)

# Separating the non-neuronal and neuronal sets

nn.cds <- L4.20.m3[, colData(L4.20.m3)$Tissue == "Non-neuronal"]
neuron.cds <- L4.20.m3[, colData(L4.20.m3)$Tissue == "Neuron"]

neuron.cds <- detect_genes(neuron.cds)
neuron.cds <- neuron.cds[rowData(neuron.cds)$num_cells_expressed > 0,]
neuron.cds
# 21576 73577 

# After more rounds of annotation and finding some additional doublets, we reached a dataset of 73528 cells,
# 70296 are identified neurons, 3232 cells are Unannotated. 

# Full cell type and tissue level annotations for cells with < 20% of UMIs from mitochondrial genes can be 
# found in the full.data.LUT.rds file.

colData(neuron.cds)$Neuron.type <- as.character(vlookup(colData(neuron.cds)$Barcode, full.data.LUT, 
                                                        result_column = "Cell.type", lookup_column = "Barcode"))

colData(neuron.cds)$Tissue <- as.character(vlookup(colData(neuron.cds)$Barcode, full.data.LUT, 
                                                        result_column = "Tissue", lookup_column = "Barcode"))

# Keeping only neurons and the unannotated cells
neuron.cds <- neuron.cds[, colData(neuron.cds)$Tissue %in% c("Neuron", "Unannotated")]

# Removing genes expressed in fewer than 5 single cells
neuron.cds <- detect_genes(neuron.cds)

neuron.cds <- neuron.cds[rowData(neuron.cds)$num_cells_expressed > 4,]

# Generating aggregate expression (transcripts per million, TPM), proportion of cells expressing, 
# nUMIs, and number of cells expressing for each neuron type

colData(neuron.cds)$n.umi <- Matrix::colSums(exprs(neuron.cds))

# This function was written by the monocle lab developers

get.expr.info.by.facet.3 <-   function(cds, colData.column) {
  cell.assignments = colData(cds)[, colData.column]
  
  if (!(class(cell.assignments) %in% c("factor", "integer", "character"))) {
    message("class(colData(cds)[, colData.column]) must be factor, integer, or character")
  }
  
  if (class(cell.assignments) == "factor") {
    unique.assignments = levels(cell.assignments)
  } else {
    unique.assignments = sort(setdiff(unique(cell.assignments), NA))
  }
  
  norm.expr = exprs(cds)
  norm.expr@x = norm.expr@x / rep.int(colData(cds)$Size_Factor, diff(norm.expr@p))
  
  facet.norm.expr = list()
  for (facet in unique.assignments) {
    facet.norm.expr[[facet]] = norm.expr[, !is.na(cell.assignments) & cell.assignments == facet]
  }
  
  message("Computing gene expression statistics")
  
  norm.means = sapply(unique.assignments, function(facet) {
    Matrix::rowMeans(facet.norm.expr[[facet]])
  })
  
  tpm = sweep(norm.means, 2, colSums(norm.means), "/") * 1000000
  
  prop.cells.expr = sapply(unique.assignments, function(facet) {
    apply(facet.norm.expr[[facet]], 1, function(x) sum(x > 0) / length(x))
  })
  
  n.cells = sapply(unique.assignments, function(facet) {
    apply(facet.norm.expr[[facet]], 1, function(x) sum(x > 0))
  })
  
  n.umi = sapply(unique.assignments, function(facet) {
    Matrix::rowSums(exprs(cds)[, !is.na(cell.assignments) & cell.assignments == facet])
  })
  
  total.n.umi.for.facet = sapply(unique.assignments, function(facet) {
    sum(colData(cds)$n.umi[!is.na(cell.assignments) & cell.assignments == facet])
  })
  
  return(list(
    tpm = tpm,
    prop.cells.expr = prop.cells.expr,
    n.umi = n.umi,
    n.cells = n.cells,
    total.n.umi.for.facet = total.n.umi.for.facet,
    norm.expr = norm.expr,
    gene.annotations = rowData(cds),
    facet = colData.column
  ))
}

L4.neuron.expr <- get.expr.info.by.facet.3(neuron.cds, "Neuron.type")

L4.TPM <- L4.neuron.expr[[1]]
L4.prop <- L4.neuron.expr[[2]]
L4.UMI <- L4.neuron.expr[[3]]
L4.ncells <- L4.neuron.expr[[4]]

saveRDS(L4.TPM, "L4_neuron_aggregate_TPM_expression.rds")
saveRDS(L4.prop, "L4_neuron_aggregate_proportion_expressing.rds")
saveRDS(L4.UMI, "L4_neuron_aggregate_UMIs.rds")
saveRDS(L4.ncells, "L4_neuron_aggregate_ncells_expressing.rds")

saveRDS(L4.neuron.expr, "L4_neuron_aggregate_expression_info.rds")
saveRDS(neuron.cds, "L4_neuron_only_cds.rds")
saveRDS(nn.cds, "L4_non-neuronal_cds.rds")
saveRDS(L4.20.m3, "L4_all_cells_passQC.rds")

### New R Session for generating final UMAPs

# Generating images for figures

# Generating final UMAPS 

library(monocle3)
library(expss)
library(dplyr)
library(ggplot2)

source("Monocle3Functions.txt")
source("MonocleFunctions.txt")

L4.cds <- readRDS("L4_neuron_only_cds.rds")
L4.test.cds <- L4.cds
L4.test.cds <- preprocess_cds(L4.test.cds, num_dim = 125, method = "PCA", norm_method = "log")

# Using the align_cds function to correct for batch effects (correcting by Experiment)
L4.test.cds <- align_cds(L4.test.cds, preprocess_method = "PCA", alignment_group = "Experiment", alignment_k = 5)
L4.test.cds <- reduce_dimension(L4.test.cds, reduction_method = "UMAP", preprocess_method = "Aligned", umap.min_dist = 0.3, umap.n_neighbors = 75)
L4.test.cds <- cluster_cells(L4.test.cds, res = 3e-3)

# Storing cluster, partition information and UMAP coordinates in the metadata for custom functions
colData(L4.test.cds)$UMAP.clusters <- monocle3::clusters(L4.test.cds)
colData(L4.test.cds)$partition <- monocle3::partitions(L4.test.cds)

colData(L4.test.cds)$UMAP_1 <- reducedDims(L4.test.cds)[["UMAP"]][,1]
colData(L4.test.cds)$UMAP_2 <- reducedDims(L4.test.cds)[["UMAP"]][,2]

# Subsetting the large central cluster of cells

L4.test.large.sub <- L4.test.cds[, colData(L4.test.cds)$partition %in% c(2)]
L4.test.large.sub <- detect_genes(L4.test.large.sub)
table(rowData(L4.test.large.sub)$num_cells_expressed > 0)

L4.test.large.sub <- L4.test.large.sub[rowData(L4.test.large.sub)$num_cells_expressed > 0,]
L4.test.large.sub
L4.test.large.sub <- preprocess_cds(L4.test.large.sub, num_dim = 110, method = "PCA", norm_method = "log")
plot_pc_variance_explained(L4.test.large.sub)
# Looks ok with 100 PCs

L4.test.large.sub <- align_cds(L4.test.large.sub, alignment_group = "Experiment", alignment_k = 50)
L4.test.large.sub <- reduce_dimension(L4.test.large.sub, reduction_method = "UMAP", preprocess_method = "Aligned", umap.min_dist = 0.2, umap.n_neighbors = 50)
L4.test.large.sub <- cluster_cells(L4.test.large.sub, res = 3e-4, k = 18)

colData(L4.test.large.sub)$UMAP.clusters <- monocle3::clusters(L4.test.large.sub)
colData(L4.test.large.sub)$partition <- monocle3::partitions(L4.test.large.sub)

colData(L4.test.large.sub)$UMAP_1 <- reducedDims(L4.test.large.sub)[["UMAP"]][,1]
colData(L4.test.large.sub)$UMAP_2 <- reducedDims(L4.test.large.sub)[["UMAP"]][,2]

plot_cells(L4.test.large.sub, color_cells_by = "Cell.type", group_label_size = 3, label_groups_by_cluster = FALSE, cell_size = 0.5)

saveRDS(L4.test.cds, "./final_gene_model/frozen_data_sets/050520_L4_neuron_m3_aligned_cds.rds")
saveRDS(L4.test.large.sub, "./final_gene_model/frozen_data_sets/050520_L4_neuron_motor_neuron_subset_aligned_cds.rds")

# all cells

L4.all.cells <- readRDS("L4_all_cells_passQC.rds")
L4.all.cells <- align_cds(L4.all.cells, preprocess_method = "PCA", alignment_group = "Experiment")
L4.all.cells <- reduce_dimension(L4.all.cells, preprocess_method = "PCA", reduction_method = "UMAP", umap.min_dist = 0.3, umap.n_neighbors = 75)
L4.all.cells <- cluster_cells(L4.all.cells, res = 3e-4)

colData(L4.all.cells)$UMAP.clusters <- monocle3::clusters(L4.all.cells)
colData(L4.all.cells)$partition <- monocle3::partitions(L4.all.cells)
colData(L4.all.cells)$UMAP_1 <- reducedDims(L4.all.cells)[["UMAP"]][,1]
colData(L4.all.cells)$UMAP_2 <- reducedDims(L4.all.cells)[["UMAP"]][,2]

plot_cells(L4.all.cells, color_cells_by = "Cell.type", group_label_size = 3.5, label_groups_by_cluster = FALSE, cell_size = 0.5) + theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12), axis.line.x.bottom = element_line(size = 1), axis.line.y.left = element_line(size = 1), axis.ticks = element_line(size = 1))

saveRDS(L4.all.cells, "L4_all_cells_aligned.rds")




