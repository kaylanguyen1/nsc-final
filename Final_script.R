
# Dataset: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE305314
# Paper used for referencing methods: https://pmc.ncbi.nlm.nih.gov/articles/PMC12623948/ 

# Install only once
#install.packages("pak")
#BiocManager::install("dittoSeq") #cell proportions across conditions
#install.packages("harmony") #integration
#install.packages('devtools')
#devtools::install_github('immunogenomics/presto') #speeds up differential expression testing
#install.packages('Seurat') #single-cell clustering
#devtools::install_github("constantAmateur/SoupX",ref='devel')
#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', force = TRUE)
#BiocManager::install("scDblFinder")
#BiocManager::install('glmGamPoi')
#install.packages("hdf5r")
#install.packages('Rfast2')
#pak::pak("jinworks/CellChat")
#pak::pak("rpolicastro/scProportionTest")

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(harmony)
library(future)
library(patchwork)
#library(DropletUtils)
#library(SoupX)
#library(DoubletFinder)
library(SingleCellExperiment)
library(scDblFinder)
library(dittoSeq)
library(scProportionTest)
library(CellChat)
library(jsonlite)

#-------------------------SINGLE CELL------------------------------------------------------

# Set working directory
setwd("D:/path_to_folder")

# Set path to parent directory containing single-cell sample folders
data_dir <- "D:/path_to_folder/dataset/snRNA-seq"

tau10_samples <- c("Tau_10mo_1", "Tau_10mo_2", "Tau_10mo_3")
tau20_samples <- c("Tau_20mo_1", "Tau_20mo_2", "Tau_20mo_3")
wt10_samples <- c("WT_10mo_1", "WT_10mo_2", "WT_10mo_3")
wt20_samples <- c("WT_20mo_1", "WT_20mo_2", "WT_20mo_3")

tau10_list <- list()
tau20_list <- list()
wt10_list <- list()
wt20_list <- list()

# ------------LOAD ALL SINGLE-CELL SAMPLES AND PREPROCESS----------------------------

# Load samples for tau condition at 10 months and preprocess/remove doublets
for (sample in tau10_samples) {
  sample_path <- file.path(data_dir, sample)
  sc_data <- Read10X(data.dir = sample_path)
  sc_data <- sc_data[["Gene Expression"]]
  
  # Create Seurat object from raw count matrix
  seurat_obj <- CreateSeuratObject(counts = sc_data, project = sample, min.cells = 3, min.features = 200)
  
  # Add metadata: sample ID, experimental condition, and age
  seurat_obj$samples <- sample
  seurat_obj$condition <- "tau"
  seurat_obj$age <- "10mo"
  
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  
  sce <- as.SingleCellExperiment(seurat_obj)
  sce <- scDblFinder(sce)
  seurat_obj$doublet_score <- sce$scDblFinder.score
  seurat_obj$doublet_class <- sce$scDblFinder.class
  table(seurat_obj$doublet_class)
  seurat_singlets <- subset(seurat_obj, subset = doublet_class == "singlet")
  
  seurat_obj <- NormalizeData(seurat_singlets)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  
  tau10_list[[sample]] <- seurat_obj
  rm(seurat_obj, sc_data, seurat_singlets, sce)
}
gc()

# Load samples for tau condition at 20 months and preprocess/remove doublets
for (sample in tau20_samples) {
  sample_path <- file.path(data_dir, sample)
  sc_data <- Read10X(data.dir = sample_path)
  sc_data <- sc_data[["Gene Expression"]]
  
  # Create Seurat object from raw count matrix
  seurat_obj <- CreateSeuratObject(counts = sc_data, project = sample, min.cells = 3, min.features = 200)
  
  # Add metadata: sample ID, experimental condition, and age
  seurat_obj$samples <- sample
  seurat_obj$condition <- "tau"
  seurat_obj$age <- "20mo"
  
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  
  sce <- as.SingleCellExperiment(seurat_obj)
  sce <- scDblFinder(sce)
  seurat_obj$doublet_score <- sce$scDblFinder.score
  seurat_obj$doublet_class <- sce$scDblFinder.class
  table(seurat_obj$doublet_class)
  seurat_singlets <- subset(seurat_obj, subset = doublet_class == "singlet")
  
  seurat_obj <- NormalizeData(seurat_singlets)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  
  tau20_list[[sample]] <- seurat_obj
  rm(seurat_obj, sc_data, seurat_singlets, sce)
}
gc()

# Load samples for wt condition at 10 months and preprocess/remove doublets
for (sample in wt10_samples) {
  sample_path <- file.path(data_dir, sample)
  sc_data <- Read10X(data.dir = sample_path)
  sc_data <- sc_data[["Gene Expression"]]
  
  # Create Seurat object from raw count matrix
  seurat_obj <- CreateSeuratObject(counts = sc_data, project = sample, min.cells = 3, min.features = 200)
  
  # Add metadata: sample ID, experimental condition, and age
  seurat_obj$samples <- sample
  seurat_obj$condition <- "wt"
  seurat_obj$age <- "10mo"
  
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  
  sce <- as.SingleCellExperiment(seurat_obj)
  sce <- scDblFinder(sce)
  seurat_obj$doublet_score <- sce$scDblFinder.score
  seurat_obj$doublet_class <- sce$scDblFinder.class
  table(seurat_obj$doublet_class)
  seurat_singlets <- subset(seurat_obj, subset = doublet_class == "singlet")
  
  seurat_obj <- NormalizeData(seurat_singlets)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  
  wt10_list[[sample]] <- seurat_obj
  rm(seurat_obj, sc_data, seurat_singlets, sce)
}
gc()

# Load samples for wt condition at 20 months and preprocess/remove doublets
for (sample in wt20_samples) {
  sample_path <- file.path(data_dir, sample)
  sc_data <- Read10X(data.dir = sample_path)
  sc_data <- sc_data[["Gene Expression"]]
  
  # Create Seurat object from raw count matrix
  seurat_obj <- CreateSeuratObject(counts = sc_data, project = sample, min.cells = 3, min.features = 200)
  
  # Add metadata: sample ID, experimental condition, and age
  seurat_obj$samples <- sample
  seurat_obj$condition <- "wt"
  seurat_obj$age <- "20mo"
  
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  
  sce <- as.SingleCellExperiment(seurat_obj)
  sce <- scDblFinder(sce)
  seurat_obj$doublet_score <- sce$scDblFinder.score
  seurat_obj$doublet_class <- sce$scDblFinder.class
  table(seurat_obj$doublet_class)
  seurat_singlets <- subset(seurat_obj, subset = doublet_class == "singlet")
  
  seurat_obj <- NormalizeData(seurat_singlets)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  
  wt20_list[[sample]] <- seurat_obj
  rm(seurat_obj, sc_data, seurat_singlets, sce)
}
gc()

#--------------MERGE SINGLE-CELL SAMPLES----------------------------

# Merge all wt samples into one Seurat object
seurat_wt10 <- merge(wt10_list[[1]], y = wt10_list[-1], add.cell.ids = wt10_samples, project = "WT10")
seurat_wt20 <- merge(wt20_list[[1]], y = wt20_list[-1], add.cell.ids = wt20_samples, project = "WT20")
merged_wt <- merge(seurat_wt10, y = seurat_wt20, add.cell.ids = c("WT10", "WT20"))
rm(wt10_list, wt20_list, seurat_wt10, seurat_wt20)
gc()

# Merge all tau samples into one Seurat object
seurat_tau10 <- merge(tau10_list[[1]], y = tau10_list[-1], add.cell.ids = tau10_samples, project = "Tau10")
seurat_tau20 <- merge(tau20_list[[1]], y = tau20_list[-1], add.cell.ids = tau20_samples, project = "Tau20")
merged_tau <- merge(seurat_tau10, y = seurat_tau20, add.cell.ids = c("Tau10", "Tau20"))
rm(tau10_list, tau20_list, seurat_tau10, seurat_tau20)
gc()

# Merge wt and tau samples into one Seurat object
merged_obj <- merge(merged_wt, y = merged_tau, add.cell.ids = c("WT", "Tau"))
merged_obj[["percent.mt"]] <- PercentageFeatureSet(merged_obj, pattern = "^Mt-")
rm(merged_wt, merged_tau)
gc()

# Visualize QC metrics as a violin plot individually so more legible
VlnPlot(merged_obj, features = "nFeature_RNA", group.by = "samples")
VlnPlot(merged_obj, features = "nCount_RNA", group.by = "samples")
VlnPlot(merged_obj, features = "percent.mt", group.by = "samples")

merged_obj <- subset(merged_obj, subset = nFeature_RNA < 5000 & percent.mt < 10)
saveRDS(merged_obj, "merged_obj.rds")
merged_obj <- readRDS("merged_obj.rds")

# Preprocessing 
# Normalize, find variable features, scale, and run PCA globally
merged_obj <- NormalizeData(merged_obj)
merged_obj <- FindVariableFeatures(merged_obj)
merged_obj <- ScaleData(merged_obj)
merged_obj <- RunPCA(merged_obj)

ElbowPlot(merged_obj, ndims = 50)

merged_obj <- FindNeighbors(merged_obj, dims = 1:19)
merged_obj <- FindClusters(merged_obj, resolution = 0.3)
merged_obj <- RunUMAP(merged_obj, dims = 1:19)

#--------------EXPORT PRE HARMONY PLOT--------------------------------
#png("sc_pre_harmony_clusters.png")
plt1 <- DimPlot(merged_obj, group.by = "condition")
plt2 <- DimPlot(merged_obj, label = TRUE)
plt1+plt2
gc()
#plot(plt1+plt2)
#dev.off()

#---------------- HARMONY INTEGRATION ------------------------------------
merged_obj_har <- IntegrateLayers(object = merged_obj, method = HarmonyIntegration, 
                                  orig.reduction = "pca", new.reduction = "harmony",
                                  verbose = FALSE)
merged_obj_har[["RNA"]] <- JoinLayers(merged_obj_har[["RNA"]])
rm(merged_obj)
gc()

merged_obj_har <- FindNeighbors(merged_obj_har, reduction = "harmony", dims = 1:19)
merged_obj_har <- FindClusters(merged_obj_har, resolution = 0.25)
merged_obj_har <- RunUMAP(merged_obj_har, reduction = "harmony", dims = 1:19)

#----------EXPORT POST HARMONY PLOT
# Visualize UMAP grouped by condition after Harmony batch correction
#png("sc_post_harmony_clusters.png")
plt1<- DimPlot(merged_obj_har, reduction = "umap", label = TRUE)
plt2<- DimPlot(merged_obj_har, reduction = "umap", label = TRUE, group.by = c("condition"))
plt1+plt2
#plot(plt1+plt2)
#dev.off()

#----------------------- CELL-TYPE ANNOTATION ---------------------------------
#--------EXPORT CSV------------------
# Find gene markers in each cluster
markers2 <- FindAllMarkers(merged_obj_har, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) #much faster after installing presto library
write.csv(markers2, "markers_single_cell_clusters_harmony_r0.1.csv")

top20 = matrix(ncol=7)
colnames(top20) = colnames(markers2)
for (i in 0:22) {
  top20 = rbind(top20, head(markers2[which(markers2$cluster==i),], 20))
}
top20 = top20[-1,]
write.csv(top20, "sc_top20.csv")

# Plot markers for each cluster's cell type
FeaturePlot(merged_obj_har, c("Mog", "Mag", "Cldn11"), label = TRUE)                     # oligodendrocyte markers high in cluster 0
FeaturePlot(merged_obj_har, c("Prox1", "C1ql3", "Cdh9"), label = TRUE)                   # DG/CA3 Glut markers in cluster 1
FeaturePlot(merged_obj_har, c("Aqp4", "Aldh1l1", "Slc1a2", "Rorb"), label = TRUE)        # astrocyte markers high in cluster 2
FeaturePlot(merged_obj_har, c("Pdgfra", "Cspg4"), label = TRUE)                          # OPC markers high in cluster 3
FeaturePlot(merged_obj_har, c("Wfs1", "Mpped1", "Fibcd1"), label = TRUE)                 # CA1 Glut markers high in cluster 4
FeaturePlot(merged_obj_har, c("Foxp2", "Grik3", "Vwc2l"), label = TRUE)                  # CT/NP SUB Glut markers high in cluster 5
FeaturePlot(merged_obj_har, c("Ntng1", "Tle4", "Bcl11b"), label = TRUE)                  # SUB-ProS Glut markers in cluster 6
FeaturePlot(merged_obj_har, c("Gad1"), label = TRUE)                                     # General GABA markers in cluster 7
FeaturePlot(merged_obj_har, c("C1qa", "Cfh", "Cped1"), label = TRUE)                     # Microglia markers high in cluster 8
FeaturePlot(merged_obj_har, c("Satb2", "Tshz2", "Tshz3", "Tox", "Vwc2l"), label = TRUE)  # L2/3 IT PPP markers high in cluster 9
FeaturePlot(merged_obj_har, c("Gad1", "Sst", "Npy", "Kcnc2"), label = TRUE)              # SST/Pvalb GABA markers high in cluster 10
FeaturePlot(merged_obj_har, c("Grik4", "Itgbl1"), label = TRUE)                          # CA3/DG Glut markers high in cluster 11
FeaturePlot(merged_obj_har, c("Vip", "Cck", "Ndnf"), label = TRUE)                       # VIP/RHP-COA Ndnf GABA high in cluster 12
FeaturePlot(merged_obj_har, c("Hcn1", "Ndst4", "Ntng1", "Nts"), label = TRUE)            # CA1-ProS Glut markers high in cluster 13
FeaturePlot(merged_obj_har, c("Col1a2", "Cped1", "Itgbl1"), label = TRUE)                # VLMC markers high in cluster 14
FeaturePlot(merged_obj_har, c("Lamp5", "Lhx6", "Gad1"), label = TRUE)                    # Lamp5 Lhx6 markers high in cluster 15
FeaturePlot(merged_obj_har, c("Slc4a5", "Ccdc170", "Htr2c"), label = TRUE)               # Chor markers high in cluster 16
FeaturePlot(merged_obj_har, c("Ebf1", "Slc2a1"), label = TRUE)                           # Endothelial markers high in cluster 17
FeaturePlot(merged_obj_har, c("Ebf1", "Rbms3"), label = TRUE)                            # Pericytes markers high in cluster 18
FeaturePlot(merged_obj_har, c("Ccdc170", "Cdhr4"), label = TRUE)                         # Ependymal markers high in cluster 19
FeaturePlot(merged_obj_har, c("Cdh9", "Cfh"), label = TRUE)                              # Other Immune cell markers in cluster 20
FeaturePlot(merged_obj_har, c("Cd163", "Pf4", "Cd79a"), label = TRUE)                    # BAM/Lymph markers in cluster 21

# New labels for clusters based off of each cluster's top markers
new.cluster.ids <- c('0' = "Oligo",
                     '1' = "DG/CA3 Glut",
                     '2' = "Astrocytes",
                     '3' = "OPC",
                     '4' = "CA1 Glut",
                     '5' = "CT/NP SUB Glut",
                     '6' = "SUB-ProS Glut",
                     '7' = "GABA", 
                     '8' = "Micro",
                     '9' = "L2/3 IT PPP", 
                     '10' = "SST/Pvalb GABA", 
                     '11' = "CA3/DG Glut", 
                     '12' = "VIP/RHP-COA Ndnf GABA",
                     '13' = "CA1-ProS Glut",
                     '14' = "VLMC",
                     '15' = "Lamp5 Lhx6 GABA",
                     '16' = "Chor",
                     '17' = "Endo",
                     '18' = "Peri", 
                     '19' = "Ependymal",
                     '20' = "Immune",
                     '21' = "BAM/Lymph")

# Relabel each cluster with its cell type
merged_obj_har <- RenameIdents(merged_obj_har, new.cluster.ids)

#----------------EXPORT ANNOTATED CLUSTERS-------------------
png("sc_annotated_clusters.png")
DimPlot(merged_obj_har, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
plot(DimPlot(merged_obj_har, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend())
dev.off()

new.cluster.order <- c("Oligo",
                       "DG/CA3 Glut",
                       "Astrocytes",
                       "OPC",
                       "CA1 Glut",
                       "CT/NP SUB Glut",
                       "SUB-ProS Glut",
                       "GABA",
                       "Micro",
                       "L2/3 IT PPP",
                       "SST/Pvalb GABA",
                       "CA3/DG Glut",
                       "VIP/RHP-COA Ndnf GABA",
                       "CA1-ProS Glut",
                       "VLMC",
                       "Lamp5 Lhx6 GABA",
                       "Chor",
                       "Endo",
                       "Peri",
                       "Ependymal",
                       "Immune",
                       "BAM/Lymph")

new.feature.order <- c("Mog",
                       "Prox1",
                       "Aqp4",
                       "Pdgfra",
                       "Wfs1",
                       "Grik3",
                       "Ntng1",
                       "Gad1",
                       "C1qa",
                       "Satb2",
                       "Sst",
                       "Itgbl1",
                       "Vip", 
                       "Ndst4",
                       "Col1a2",
                       "Lamp5",
                       "Ccdc170",
                       "Ebf1", "Slc2a1",
                       "Rbms3", "Cdhr4"
)

# Plot top features of clusters 
# ---------MAKE PLOT PRETTIER AND EXPORT-------------------------
Idents(merged_obj_har) <- factor(Idents(merged_obj_har), levels = new.cluster.order)
DotPlot(merged_obj_har, features = new.feature.order) + RotatedAxis() + 
  scale_color_gradientn(
    colors = c("khaki1", "yellow", "lightsalmon", "salmon", "palevioletred3", "magenta4", "purple4"))

markers2 %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 2) %>%
  ungroup() -> top2
DoHeatmap(merged_obj_har, features = top2$gene) + NoLegend()

#---------------------COMPARE PROPORTION OF CELL TYPES -------------------------------------
# Visualize proportion of cell types between condition and age groups
merged_obj_har$CellType <- Idents(merged_obj_har)
dittoBarPlot(merged_obj_har, "CellType", group.by = "condition")
dittoBarPlot(merged_obj_har, "CellType", group.by = "age")

# Statistical analysis of cell type proportions
prop_test <- sc_utils(merged_obj_har)
prop_test <- permutation_test(
  prop_test,
  cluster_identity = "CellType",
  sample_identity = "condition",
  sample_1 = "tau",
  sample_2 = "wt",         
  n_permutations = 1000
)

permutation_plot(prop_test)

prop_test <- sc_utils(merged_obj_har)
prop_test <- permutation_test(
  prop_test,
  cluster_identity = "CellType",
  sample_identity = "age",
  sample_1 = "10mo",
  sample_2 = "20mo",         
  n_permutations = 1000
)

permutation_plot(prop_test)

#--------------- DIFFERENTIAL EXPRESSION BETWEEN CONDITION --------------------------------------
# Add to metadata cell type and condition
merged_obj_har$celltype.condition <- paste(merged_obj_har$CellType, merged_obj_har$condition, sep = "_")
Idents(merged_obj_har) <- "celltype.condition"

# Oligodendrocyte cluster 
oligo.response <- FindMarkers(merged_obj_har, ident.1 = "Oligo_wt", ident.2 = "Oligo_tau", verbose = FALSE) 
head(oligo.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Grip2"), split.by = "condition", cols = c("grey","red"), reduction = "umap", label = TRUE)

# DG/CA3 Glut cluster
dg_ca3_glut.response <- FindMarkers(merged_obj_har, ident.1 = "DG/CA3 Glut_wt", ident.2 = "DG/CA3 Glut_tau", verbose = FALSE) 
head(dg_ca3_glut.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Chrm2"), split.by = "condition", cols = c("grey","red"), reduction = "umap", label = TRUE)

# Astrocyte cluster
astro.response <- FindMarkers(merged_obj_har, ident.1 = "Astrocytes_wt", ident.2 = "Astrocytes_tau", verbose = FALSE) 
head(astro.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Ptgds"), split.by = "condition", cols = c("grey","red"), reduction = "umap", label = TRUE)

# OPC cluster
opc.response <- FindMarkers(merged_obj_har, ident.1 = "OPC_wt", ident.2 = "OPC_tau", verbose = FALSE) 
head(opc.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Ttr"), split.by = "condition", cols = c("grey","red"), reduction = "umap", label = TRUE)

# CA1 Glut cluster
ca1_glut.response <- FindMarkers(merged_obj_har, ident.1 = "CA1 Glut_wt", ident.2 = "CA1 Glut_tau", verbose = FALSE) 
head(ca1_glut.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Ttr"), split.by = "condition", cols = c("grey","red"), reduction = "umap", label = TRUE)

# CT/NP SUB Glut cluster
ctnp_glut.response <- FindMarkers(merged_obj_har, ident.1 = "CT/NP SUB Glut_wt", ident.2 = "CT/NP SUB Glut_tau", verbose = FALSE) 
head(ctnp_glut.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Ttr"), split.by = "condition", cols = c("grey","red"), reduction = "umap", label = TRUE)

# SUB-ProS Glut cluster
subpro_glut.response <- FindMarkers(merged_obj_har, ident.1 = "SUB-ProS Glut_wt", ident.2 = "SUB-ProS Glut_tau", verbose = FALSE) 
head(subpro_glut.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Ttr"), split.by = "condition", cols = c("grey","red"), reduction = "umap", label = TRUE)

# GABA cluster
gaba.response <- FindMarkers(merged_obj_har, ident.1 = "GABA_wt", ident.2 = "GABA_tau", verbose = FALSE) 
head(gaba.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Ttr"), split.by = "condition", cols = c("grey","red"), reduction = "umap", label = TRUE)

# Microglia cluster
micro.response <- FindMarkers(merged_obj_har, ident.1 = "Micro_wt", ident.2 = "Micro_tau", verbose = FALSE) 
head(micro.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Ttr"), split.by = "condition", cols = c("grey","red"), reduction = "umap", label = TRUE)

# L2/3 IT PPP cluster
l23_it.response <- FindMarkers(merged_obj_har, ident.1 = "L2/3 IT PPP_wt", ident.2 = "L2/3 IT PPP_tau", verbose = FALSE) 
head(l23_it.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Ttr"), split.by = "condition", cols = c("grey","red"), reduction = "umap")

# SST/Pvalb cluster
sst_gaba.response <- FindMarkers(merged_obj_har, ident.1 = "SST/Pvalb GABA_wt", ident.2 = "SST/Pvalb GABA_tau", verbose = FALSE) 
head(sst_gaba.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Ttr"), split.by = "condition", cols = c("grey","red"), reduction = "umap")

# CA3/DG Glut cluster
dg_glut.response <- FindMarkers(merged_obj_har, ident.1 = "CA3/DG Glut_wt", ident.2 = "CA3/DG Glut_tau", verbose = FALSE) 
head(dg_glut.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Ttr"), split.by = "condition", cols = c("grey","red"), reduction = "umap")

# VIP/Ndnf GABA cluster
vip_gaba.response <- FindMarkers(merged_obj_har, ident.1 = "VIP/RHP-COA Ndnf GABA_wt", ident.2 = "VIP/RHP-COA Ndnf GABA_tau", verbose = FALSE) 
head(vip_gaba.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Ttr"), split.by = "condition", cols = c("grey","red"), reduction = "umap")

# CA1-ProS Glut cluster
ca1pro_glut.response <- FindMarkers(merged_obj_har, ident.1 = "CA1-ProS Glut_wt", ident.2 = "CA1-ProS Glut_tau", verbose = FALSE) 
head(ca1pro_glut.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Ttr"), split.by = "condition", cols = c("grey","red"), reduction = "umap")

# VLMC cluster
vlmc.response <- FindMarkers(merged_obj_har, ident.1 = "VLMC_wt", ident.2 = "VLMC_tau", verbose = FALSE) 
head(vlmc.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Sncg"), split.by = "condition", cols = c("grey","red"), reduction = "umap")

# Lamp5 Lhx6 cluster
lamp5_gaba.response <- FindMarkers(merged_obj_har, ident.1 = "Lamp5 Lhx6 GABA_wt", ident.2 = "Lamp5 Lhx6 GABA_tau", verbose = FALSE) 
head(lamp5_gaba.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Ttr"), split.by = "condition", cols = c("grey","red"), reduction = "umap")

# Chor cluster
chor.response <- FindMarkers(merged_obj_har, ident.1 = "Chor_wt", ident.2 = "Chor_tau", verbose = FALSE) 
head(chor.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Cst3"), split.by = "condition", cols = c("grey","red"), reduction = "umap")

# Endo cluster
endo.response <- FindMarkers(merged_obj_har, ident.1 = "Endo_wt", ident.2 = "Endo_tau", verbose = FALSE) 
head(endo.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Chrm2"), split.by = "condition", cols = c("grey","red"), reduction = "umap")

# Peri cluster
peri.response <- FindMarkers(merged_obj_har, ident.1 = "Peri_wt", ident.2 = "Peri_tau", verbose = FALSE) 
head(peri.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Chrm2"), split.by = "condition", cols = c("grey","red"), reduction = "umap")

# Ependymal cluster
epend.response <- FindMarkers(merged_obj_har, ident.1 = "Ependymal_wt", ident.2 = "Ependymal_tau", verbose = FALSE) 
head(epend.response, n = 15)
FeaturePlot(merged_obj_har, features = c("S1pr3"), split.by = "condition", cols = c("grey","red"), reduction = "umap")

# Immune cluster
imm.response <- FindMarkers(merged_obj_har, ident.1 = "Immune_wt", ident.2 = "Immune_tau", verbose = FALSE) 
head(imm.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Chrm2"), split.by = "condition", cols = c("grey","red"), reduction = "umap")

# BAM/Lymph cluster
bam.response <- FindMarkers(merged_obj_har, ident.1 = "BAM/Lymph_wt", ident.2 = "BAM/Lymph_tau", verbose = FALSE) 
head(bam.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Magi2"), split.by = "condition", cols = c("grey","red"), reduction = "umap", label = TRUE)

rm(astro.response, bam.response, ca1_glut.response, ca1pro_glut.response, chor.response, ctnp_glut.response,
   dg_ca3_glut.response, dg_glut.response, endo.response, epend.response, gaba.response, imm.response,
   l23_it.response, lamp5_gaba.response, micro.response, oligo.response, peri.response, opc.response,
   sst_gaba.response, subpro_glut.response, vip_gaba.response, vlmc.response
   )
gc()

#------------------ DIFFERENTIAL EXPRESSION BETWEEN AGES --------------------------------------
# Add to metadata cell type and age
merged_obj_har$celltype.age <-  paste(merged_obj_har$CellType, merged_obj_har$age, sep = "_")
Idents(merged_obj_har) <- "celltype.age"

# Oligodendrocyte cluster 
oligo.response <- FindMarkers(merged_obj_har, ident.1 = "Oligo_10mo", ident.2 = "Oligo_20mo", verbose = FALSE) #specify stat test
head(oligo.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Ttr"), split.by = "age", cols = c("grey","red"), reduction = "umap", label = TRUE)

# DG/CA3 Glut cluster
dg_ca3_glut.response <- FindMarkers(merged_obj_har, ident.1 = "DG/CA3 Glut_10mo", ident.2 = "DG/CA3 Glut_20mo", verbose = FALSE) 
head(dg_ca3_glut.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Ttr"), split.by = "age", cols = c("grey","red"), reduction = "umap")

# Astrocyte cluster
astro.response <- FindMarkers(merged_obj_har, ident.1 = "Astrocytes_10mo", ident.2 = "Astrocytes_20mo", verbose = FALSE) 
head(astro.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Ttr"), split.by = "age", cols = c("grey","red"), reduction = "umap")

# OPC cluster
opc.response <- FindMarkers(merged_obj_har, ident.1 = "OPC_10mo", ident.2 = "OPC_20mo", verbose = FALSE) 
head(opc.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Ttr"), split.by = "age", cols = c("grey","red"), reduction = "umap")

# CA1 Glut cluster
ca1_glut.response <- FindMarkers(merged_obj_har, ident.1 = "CA1 Glut_10mo", ident.2 = "CA1 Glut_20mo", verbose = FALSE) 
head(ca1_glut.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Ttr"), split.by = "age", cols = c("grey","red"), reduction = "umap")

# CT/NP SUB Glut cluster
ctnp_glut.response <- FindMarkers(merged_obj_har, ident.1 = "CT/NP SUB Glut_10mo", ident.2 = "CT/NP SUB Glut_20mo", verbose = FALSE) 
head(ctnp_glut.response, n = 15)
FeaturePlot(merged_obj_har, features = c("RGD1565158"), split.by = "age", cols = c("grey","red"), reduction = "umap")

# SUB-ProS Glut cluster
subpro_glut.response <- FindMarkers(merged_obj_har, ident.1 = "SUB-ProS Glut_10mo", ident.2 = "SUB-ProS Glut_20mo", verbose = FALSE) 
head(subpro_glut.response, n = 15)
FeaturePlot(merged_obj_har, features = c("RGD1565158"), split.by = "age", cols = c("grey","red"), reduction = "umap")

# GABA cluster
gaba.response <- FindMarkers(merged_obj_har, ident.1 = "GABA_10mo", ident.2 = "GABA_20mo", verbose = FALSE) 
head(gaba.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Ttr"), split.by = "age", cols = c("grey","red"), reduction = "umap")

# Microglia cluster
micro.response <- FindMarkers(merged_obj_har, ident.1 = "Micro_10mo", ident.2 = "Micro_20mo", verbose = FALSE) 
head(micro.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Grxcr1"), split.by = "age", cols = c("grey","red"), reduction = "umap")

# L2/3 IT PPP cluster
l23_it.response <- FindMarkers(merged_obj_har, ident.1 = "L2/3 IT PPP_10mo", ident.2 = "L2/3 IT PPP_20mo", verbose = FALSE) 
head(l23_it.response, n = 15)
FeaturePlot(merged_obj_har, features = c("RGD1565158"), split.by = "age", cols = c("grey","red"), reduction = "umap")

# SST/Pvalb cluster
sst_gaba.response <- FindMarkers(merged_obj_har, ident.1 = "SST/Pvalb GABA_10mo", ident.2 = "SST/Pvalb GABA_20mo", verbose = FALSE) 
head(sst_gaba.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Ttr"), split.by = "age", cols = c("grey","red"), reduction = "umap")

# CA3/DG Glut cluster
dg_glut.response <- FindMarkers(merged_obj_har, ident.1 = "CA3/DG Glut_10mo", ident.2 = "CA3/DG Glut_20mo", verbose = FALSE) 
head(dg_glut.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Ttr"), split.by = "age", cols = c("grey","red"), reduction = "umap")

# VIP/Ndnf GABA cluster
vip_gaba.response <- FindMarkers(merged_obj_har, ident.1 = "VIP/RHP-COA Ndnf GABA_10mo", ident.2 = "VIP/RHP-COA Ndnf GABA_20mo", verbose = FALSE) 
head(vip_gaba.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Ttr"), split.by = "age", cols = c("grey","red"), reduction = "umap")

# CA1-ProS Glut cluster
ca1pro_glut.response <- FindMarkers(merged_obj_har, ident.1 = "CA1-ProS Glut_10mo", ident.2 = "CA1-ProS Glut_20mo", verbose = FALSE) 
head(ca1pro_glut.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Ttr"), split.by = "age", cols = c("grey","red"), reduction = "umap")

# VLMC cluster
vlmc.response <- FindMarkers(merged_obj_har, ident.1 = "VLMC_10mo", ident.2 = "VLMC_20mo", verbose = FALSE) 
head(vlmc.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Ttr"), split.by = "age", cols = c("grey","red"), reduction = "umap")

# Lamp5 Lhx6 cluster
lamp5_gaba.response <- FindMarkers(merged_obj_har, ident.1 = "Lamp5 Lhx6 GABA_10mo", ident.2 = "Lamp5 Lhx6 GABA_20mo", verbose = FALSE) 
head(lamp5_gaba.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Ttr"), split.by = "age", cols = c("grey","red"), reduction = "umap")

# Chor cluster
chor.response <- FindMarkers(merged_obj_har, ident.1 = "Chor_10mo", ident.2 = "Chor_20mo", verbose = FALSE) 
head(chor.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Cst3"), split.by = "age", cols = c("grey","red"), reduction = "umap")

# Endo cluster
endo.response <- FindMarkers(merged_obj_har, ident.1 = "Endo_10mo", ident.2 = "Endo_20mo", verbose = FALSE) 
head(endo.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Adgrl3"), split.by = "age", cols = c("grey","red"), reduction = "umap")

# Peri cluster
peri.response <- FindMarkers(merged_obj_har, ident.1 = "Peri_10mo", ident.2 = "Peri_20mo", verbose = FALSE) 
head(peri.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Ttr"), split.by = "age", cols = c("grey","red"), reduction = "umap")

# Ependymal cluster
epend.response <- FindMarkers(merged_obj_har, ident.1 = "Ependymal_10mo", ident.2 = "Ependymal_20mo", verbose = FALSE) 
head(epend.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Lsamp"), split.by = "age", cols = c("grey","red"), reduction = "umap")

# Immune cluster
imm.response <- FindMarkers(merged_obj_har, ident.1 = "Immune_10mo", ident.2 = "Immune_20mo", verbose = FALSE) 
head(imm.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Ttr"), split.by = "age", cols = c("grey","red"), reduction = "umap")

# BAM/Lymph cluster
bam.response <- FindMarkers(merged_obj_har, ident.1 = "BAM/Lymph_10mo", ident.2 = "BAM/Lymph_20mo", verbose = FALSE) 
head(bam.response, n = 15)
FeaturePlot(merged_obj_har, features = c("Ttr"), split.by = "age", cols = c("grey","red"), reduction = "umap", label = TRUE)


rm(astro.response, bam.response, ca1_glut.response, ca1pro_glut.response, chor.response, ctnp_glut.response,
   dg_ca3_glut.response, dg_glut.response, endo.response, epend.response, gaba.response, imm.response,
   l23_it.response, lamp5_gaba.response, micro.response, oligo.response, peri.response, opc.response,
   sst_gaba.response, subpro_glut.response, vip_gaba.response, vlmc.response
)
rm(merged_obj_har)
gc()

#---------------SPATIAL----------------------------------------------------
# Set file path
setwd("D:/path_to_folder")

data_dir <- "D:/path_to_folder/dataset/spatial"

samples <- c("visium_tau_10mo", "visium_tau_20mo", "visium_wt_10mo", "visium_wt_20mo")
spatial_list <- list()

# If loading in samples from online dataset, uncomment below
# Convert individual .tar.gz files into h5 file to load into Load10X_Spatial
#for (sample in samples){
#  fp <- file.path(data_dir, sample, "filtered_feature_bc_matrix")
#  write_path <- file.path(data_dir, sample, "filtered_feature_bc_matrix.h5")

#  filter_matrix <- Read10X(fp)
#  write10xCounts(write_path, filter_matrix, type = "HDF5", version = "3", overwrite = TRUE)
#  rm(filter_matrix)
#}

# Load Visium data for each sample and preprocess
for (sample in samples) {
  sample_path <- file.path(data_dir, sample)
  list.files(file.path(sample_path, "spatial"))
  obj <- Load10X_Spatial(data.dir = sample_path, filename = "filtered_feature_bc_matrix.h5", assay = "Spatial")
  
  names(obj@images) <- sample
  
  obj$sample <- sample
  obj$condition <- ifelse(grepl("tau", sample), "tau", "wt")
  obj$age <- ifelse(grepl("10mo", sample), "10mo", "20mo")
  
  plt1 <- VlnPlot(obj, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
  plt2 <- SpatialFeaturePlot(obj, features = "nCount_Spatial") + theme(legend.position = "right")
  print(wrap_plots(plt1, plt2))
  print(VlnPlot(obj, features = c("nFeature_Spatial", "nCount_Spatial"), ncol = 2, pt.size = 0.1))
  
  DefaultAssay(obj) <- "Spatial"
  
  obj <- SCTransform(obj, assay = "Spatial", new.assay.name = "SCT", return.only.var.genes = FALSE, verbose = FALSE)
  spatial_list[[sample]] <- obj
  rm(obj, plt1, plt2)
}

# Merge samples into one seurat object
spatial_merged <- merge(spatial_list[[1]], y = spatial_list[-1], add.cell.ids = samples, project = "Spatial")
SpatialFeaturePlot(spatial_merged, features = c("Hpca", "Gfap"))

p1 <- SpatialFeaturePlot(spatial_merged, features = "Hpca", pt.size.factor = 1)
p2 <- SpatialFeaturePlot(spatial_merged, features = "Hpca", alpha = c(0.1, 1))
p1+p2

DefaultAssay(spatial_merged) <- "SCT"
features <- SelectIntegrationFeatures(object.list = spatial_list, nfeatures = 3000)
VariableFeatures(spatial_merged) <- features

# Dimensionality reduction, clustering, and visualization
spatial_merged <- RunPCA(spatial_merged, assay = "SCT", features = features, verbose = FALSE)
ElbowPlot(spatial_merged, ndims = 50)

spatial_merged <- FindNeighbors(spatial_merged, reduction = "pca", dims = 1:40)
spatial_merged <- FindClusters(spatial_merged, resolution = 0.5, verbose = FALSE)
spatial_merged <- RunUMAP(spatial_merged, reduction = "pca", dims = 1:40)

p1 <- DimPlot(spatial_merged, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(spatial_merged, label = TRUE, label.size = 3)
p1 + p2

for (sample in samples) {
    print(
      SpatialDimPlot(spatial_merged, images = sample, group.by = "seurat_clusters", label = TRUE, label.size = 3
      ) + ggtitle(sample)
    )
}

# Identify genes with spatially patterned expression
spatial_merged <- FindSpatiallyVariableFeatures(spatial_merged, assay = "SCT", features = VariableFeatures(spatial_merged)[1:1000],
                                                selection.method = "moransi")
top.features <- head(SpatiallyVariableFeatures(spatial_merged, method = "moransi"), 6)
SpatialFeaturePlot(spatial_merged, features = top.features, ncol = 3, alpha = c(0.1,1))

DefaultAssay(spatial_merged) <- "SCT"
Idents(spatial_merged) <- "seurat_clusters"
spatial_merged <- PrepSCTFindMarkers(spatial_merged)

table(Idents(spatial_merged))
spatial_markers <- FindAllMarkers(spatial_merged, assay = "SCT", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top20_spatial = matrix(ncol=7)
colnames(top20_spatial) = colnames(spatial_markers)
for (i in 0:20) {
  top20_spatial = rbind(top20_spatial, head(spatial_markers[which(spatial_markers$cluster==i),], 20))
}
top20_spatial = top20_spatial[-1,]
write.csv(top20_spatial, "spatial_top20.csv")

#---------------------------- CLUSTER CELL/REGION ANNOTATION ------------------------------------------------------------------
# Marker genes from spatial_markers paired with region/cell of interest
FeaturePlot(spatial_merged, features = "Mag", label = TRUE)                           # marker gene for cluster 0,17 = Fiber Tracts
FeaturePlot(spatial_merged, features = "Stac2", label = TRUE)                         # marker gene for cluster 1,15 = L5/6 Cortex
FeaturePlot(spatial_merged, features = "Prkcd", label = TRUE)                         # marker gene for cluster 2,9,12 = Thalamus
FeaturePlot(spatial_merged, features = "Arc", label = TRUE)                           # marker gene for cluster 3 = L2/3 Cortex
FeaturePlot(spatial_merged, features = "Ecel1", label = TRUE)                         # marker gene for cluster 4 = Hypothalamus
FeaturePlot(spatial_merged, features = "Olfm1", label = TRUE)                         # marker gene for cluster 5,6 = Olfactory Areas
FeaturePlot(spatial_merged, features = "Six3", label = TRUE)                          # marker gene for cluster 7,10 = Striatum
FeaturePlot(spatial_merged, features = "Cabp7", label = TRUE)                         # marker gene for cluster 8 = Dentate Gyrus
FeaturePlot(spatial_merged, features = c("Doc2b", "Cabp7", "Gabra5"), label = TRUE)   # marker gene for cluster 11 = CA1
FeaturePlot(spatial_merged, features = c("Cabp7", "Gabra5"), label = TRUE)            # marker gene for cluster 13 = CA2/3
FeaturePlot(spatial_merged, features = "Gfra2", label = TRUE)                         # marker gene for cluster 18 = L4 Cortex

# Add merged region/cell to cluster identity
spatial_merged$original_clusters <- spatial_merged$seurat_clusters
spatial_merged$merged_region <- as.character(spatial_merged$seurat_clusters)
spatial_merged$merged_region[spatial_merged$seurat_clusters %in% c("0", "17")] <- "Fiber Tracts"
spatial_merged$merged_region[spatial_merged$seurat_clusters %in% c("1", "15")] <- "L5/6 Cortex"
spatial_merged$merged_region[spatial_merged$seurat_clusters %in% c("2", "9", "12")] <- "Thalamus"
spatial_merged$merged_region[spatial_merged$seurat_clusters %in% c("3")] <- "L2/3 Cortex"
spatial_merged$merged_region[spatial_merged$seurat_clusters %in% c("4")] <- "Hypothalamus"
spatial_merged$merged_region[spatial_merged$seurat_clusters %in% c("5", "6")] <- "Olfactory Areas"
spatial_merged$merged_region[spatial_merged$seurat_clusters %in% c("7", "10")] <- "Striatum"
spatial_merged$merged_region[spatial_merged$seurat_clusters %in% c("8")] <- "Dentate Gyrus"
spatial_merged$merged_region[spatial_merged$seurat_clusters %in% c("11")] <- "CA1"
spatial_merged$merged_region[spatial_merged$seurat_clusters %in% c("13")] <- "CA2/3"
spatial_merged$merged_region[spatial_merged$seurat_clusters %in% c("18")] <- "L4 Cortex"
spatial_merged$merged_region[spatial_merged$seurat_clusters %in% c("14", "16")] <- "Other"

new.cluster.order <- c("Fiber Tracts", "Thalamus", "L5/6 Cortex", "Hypothalamus", "L2/3 Cortex", "Olfactory Areas", "Striatum", "L4 Cortex", "CA1", "Dentate Gyrus", "CA2/3", "Other")
new.feature.order <-  c("Mag", "Prkcd", "Stac2", "Ecel1", "Arc", "Olfm1", "Six3", "Gfra2", "Doc2b", "Cabp7", "Gabra5")

# Organize regions and their features in fixed position
spatial_merged$merged_region <- factor(spatial_merged$merged_region, levels = new.cluster.order)
Idents(spatial_merged) <- "merged_region"

features.use <- intersect(new.feature.order, rownames(spatial_merged))

# Dot plot of features with merged regions against their marker gene
DotPlot(spatial_merged, features = features.use, group.by = "merged_region") + 
  RotatedAxis() + 
  scale_color_gradientn(
    colors = c("khaki1", "yellow", "lightsalmon", "salmon", "palevioletred3", "magenta4", "purple4")
  )

# Plot of annotated clusters
DimPlot(spatial_merged, reduction = "umap", group.by = "merged_region", label = TRUE, repel = TRUE) + NoLegend()

# Plot annotated regions on each sample
for (sample in samples) {
  print(
    SpatialDimPlot(spatial_merged, images = sample, group.by = "merged_region", label = TRUE, label.size = 4) +
    ggtitle(sample)
  )
}

# Heatmap of top 2 genes for each region
spatial_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 2) %>%
  ungroup() -> top2
DoHeatmap(spatial_merged, features = top2$gene) + NoLegend()

#--------------------------- DIFFERENTIAL EXPRESSION -------------------------------------------------
spatial_merged$region.condition.age <- paste(spatial_merged$merged_region, spatial_merged$condition, spatial_merged$age, sep = "_")
Idents(spatial_merged) <- "region.condition.age"

# DIFFERENTIAL EXPRESSION IN FIBER TRACTS
# Wt 10mo vs tau 10mo
fiber.response <- FindMarkers(spatial_merged, ident.1 = "Fiber Tracts_wt_10mo", ident.2 = "Fiber Tracts_tau_10mo", verbose = FALSE) 
head(fiber.response, n = 15)
FeaturePlot(spatial_merged, features = c("Cd81"), split.by = "condition", cols = c("grey","red"), reduction = "umap")

# Wt 20mo vs tau 20mo
fiber.response <- FindMarkers(spatial_merged, ident.1 = "Fiber Tracts_wt_20mo", ident.2 = "Fiber Tracts_tau_20mo", verbose = FALSE) 
head(fiber.response, n = 15)
FeaturePlot(spatial_merged, features = c("Gfap"), split.by = "condition", cols = c("grey","red"), reduction = "umap")

# Tau 10mo vs tau 20mo
fiber.response <- FindMarkers(spatial_merged, ident.1 = "Fiber Tracts_tau_10mo", ident.2 = "Fiber Tracts_tau_20mo", verbose = FALSE) 
head(fiber.response, n = 15)
FeaturePlot(spatial_merged, features = c("Aplp1"), split.by = "age", cols = c("grey","red"), reduction = "umap")

# DIFFERENTIAL EXPRESSION IN THALAMUS
# WT 10mo vs 10mo
thalamus.response <- FindMarkers(spatial_merged, ident.1 = "Thalamus_wt_10mo", ident.2 = "Thalamus_tau_10mo", verbose = FALSE) 
head(thalamus.response, n = 15)
FeaturePlot(spatial_merged, features = c("Cd81"), split.by = "condition", cols = c("grey","red"), reduction = "umap")

# Wt 20mo vs tau 20mo
thalamus.response <- FindMarkers(spatial_merged, ident.1 = "Thalamus_wt_20mo", ident.2 = "Thalamus_tau_20mo", verbose = FALSE) 
head(thalamus.response, n = 15)
FeaturePlot(spatial_merged, features = c("Slc6a11"), split.by = "condition", cols = c("grey","red"), reduction = "umap")

# Tau 10mo vs tau 20mo
thalamus.response <- FindMarkers(spatial_merged, ident.1 = "Thalamus_tau_10mo", ident.2 = "Thalamus_tau_20mo", verbose = FALSE) 
head(thalamus.response, n = 15)
FeaturePlot(spatial_merged, features = c("Arf1"), split.by = "age", cols = c("grey","red"), reduction = "umap")

# DIFFERENTIAL EXPRESSION IN L5/6 CORTEX
# WT 10mo vs 10mo
l56_cortex.response <- FindMarkers(spatial_merged, ident.1 = "L5/6 Cortex_wt_10mo", ident.2 = "L5/6 Cortex_tau_10mo", verbose = FALSE) 
head(l56_cortex.response, n = 15)
FeaturePlot(spatial_merged, features = c("Cd81"), split.by = "condition", cols = c("grey","red"), reduction = "umap")

# Wt 20mo vs tau 20mo
l56_cortex.response <- FindMarkers(spatial_merged, ident.1 = "L5/6 Cortex_wt_20mo", ident.2 = "L5/6 Cortex_tau_20mo", verbose = FALSE) 
head(l56_cortex.response, n = 15)
FeaturePlot(spatial_merged, features = c("Rtn1"), split.by = "condition", cols = c("grey","red"), reduction = "umap")

# Tau 10mo vs tau 20mo
l56_cortex.response <- FindMarkers(spatial_merged, ident.1 = "L5/6 Cortex_tau_10mo", ident.2 = "L5/6 Cortex_tau_20mo", verbose = FALSE) 
head(l56_cortex.response, n = 15)
FeaturePlot(spatial_merged, features = c("Aplp1"), split.by = "age", cols = c("grey","red"), reduction = "umap")

# DIFFERENTIAL EXPRESSION IN HYPOTHALAMUS
# WT 10mo vs 10mo
hypo.response <- FindMarkers(spatial_merged, ident.1 = "Hypothalamus_wt_10mo", ident.2 = "Hypothalamus_tau_10mo", verbose = FALSE) 
head(hypo.response, n = 15)
FeaturePlot(spatial_merged, features = c("Cd81"), split.by = "condition", cols = c("grey","red"), reduction = "umap")

# Wt 20mo vs tau 20mo
hypo.response <- FindMarkers(spatial_merged, ident.1 = "Hypothalamus_wt_20mo", ident.2 = "Hypothalamus_tau_20mo", verbose = FALSE) 
head(hypo.response, n = 15)
FeaturePlot(spatial_merged, features = c("Slc6a11"), split.by = "condition", cols = c("grey","red"), reduction = "umap")

# Tau 10mo vs tau 20mo
hypo.response <- FindMarkers(spatial_merged, ident.1 = "Hypothalamus_tau_10mo", ident.2 = "Hypothalamus_tau_20mo", verbose = FALSE) 
head(hypo.response, n = 15)
FeaturePlot(spatial_merged, features = c("Aplp1"), split.by = "age", cols = c("grey","red"), reduction = "umap")

# DIFFERENTIAL EXPRESSION IN L2/3 CORTEX
# WT 10mo vs 10mo
l23_cortex.response <- FindMarkers(spatial_merged, ident.1 = "L2/3 Cortex_wt_10mo", ident.2 = "L2/3 Cortex_tau_10mo", verbose = FALSE) 
head(l23_cortex.response, n = 15)
FeaturePlot(spatial_merged, features = c("Gfap"), split.by = "condition", cols = c("grey","red"), reduction = "umap")

# Wt 20mo vs tau 20mo
l23_cortex.response <- FindMarkers(spatial_merged, ident.1 = "L2/3 Cortex_wt_20mo", ident.2 = "L2/3 Cortex_tau_20mo", verbose = FALSE) 
head(l23_cortex.response, n = 15)
FeaturePlot(spatial_merged, features = c("Lingo1"), split.by = "condition", cols = c("grey","red"), reduction = "umap")

# Tau 10mo vs tau 20mo
l23_cortex.response <- FindMarkers(spatial_merged, ident.1 = "L2/3 Cortex_tau_10mo", ident.2 = "L2/3 Cortex_tau_20mo", verbose = FALSE) 
head(l23_cortex.response, n = 15)
FeaturePlot(spatial_merged, features = c("Arf3"), split.by = "age", cols = c("grey","red"), reduction = "umap")

# No differential expression for L4 Cortex because not enough cells

# DIFFERENTIAL EXPRESSION IN OLFACTORY AREAS
# WT 10mo vs 10mo
olf.response <- FindMarkers(spatial_merged, ident.1 = "Olfactory Areas_wt_10mo", ident.2 = "Olfactory Areas_tau_10mo", verbose = FALSE) 
head(olf.response, n = 15)
FeaturePlot(spatial_merged, features = c("Cd81"), split.by = "condition", cols = c("grey","red"), reduction = "umap")

# Wt 20mo vs tau 20mo
olf.response <- FindMarkers(spatial_merged, ident.1 = "Olfactory Areas_wt_20mo", ident.2 = "Olfactory Areas_tau_20mo", verbose = FALSE) 
head(olf.response, n = 15)
FeaturePlot(spatial_merged, features = c("Gfap"), split.by = "condition", cols = c("grey","red"), reduction = "umap")

# Tau 10mo vs tau 20mo
olf.response <- FindMarkers(spatial_merged, ident.1 = "Olfactory Areas_tau_10mo", ident.2 = "Olfactory Areas_tau_20mo", verbose = FALSE) 
head(olf.response, n = 15)
FeaturePlot(spatial_merged, features = c("Arf3"), split.by = "age", cols = c("grey","red"), reduction = "umap")

# DIFFERENTIAL EXPRESSION IN Striatum
# WT 10mo vs 10mo
striatum.response <- FindMarkers(spatial_merged, ident.1 = "Striatum_wt_10mo", ident.2 = "Striatum_tau_10mo", verbose = FALSE) 
head(striatum.response, n = 15)
FeaturePlot(spatial_merged, features = c("Eif2s3y"), split.by = "condition", cols = c("grey","red"), reduction = "umap")

# Wt 20mo vs tau 20mo
striatum.response <- FindMarkers(spatial_merged, ident.1 = "Striatum_wt_20mo", ident.2 = "Striatum_tau_20mo", verbose = FALSE) 
head(striatum.response, n = 15)
FeaturePlot(spatial_merged, features = c("Snap25"), split.by = "condition", cols = c("grey","red"), reduction = "umap")

# Tau 10mo vs tau 20mo
striatum.response <- FindMarkers(spatial_merged, ident.1 = "Striatum_tau_10mo", ident.2 = "Striatum_tau_20mo", verbose = FALSE) 
head(striatum.response, n = 15)
FeaturePlot(spatial_merged, features = c("Arf3"), split.by = "age", cols = c("grey","red"), reduction = "umap")

# DIFFERENTIAL EXPRESSION IN CA1
# WT 10mo vs 10mo
ca1.response <- FindMarkers(spatial_merged, ident.1 = "CA1_wt_10mo", ident.2 = "CA1_tau_10mo", verbose = FALSE) 
head(ca1.response, n = 15)
FeaturePlot(spatial_merged, features = c("Eif2s3y"), split.by = "condition", cols = c("grey","red"), reduction = "umap")

# Wt 20mo vs tau 20mo
ca1.response <- FindMarkers(spatial_merged, ident.1 = "CA1_wt_20mo", ident.2 = "CA1_tau_20mo", verbose = FALSE) 
head(ca1.response, n = 15)
FeaturePlot(spatial_merged, features = c("Jph4"), split.by = "condition", cols = c("grey","red"), reduction = "umap")

# Tau 10mo vs tau 20mo
ca1.response <- FindMarkers(spatial_merged, ident.1 = "CA1_tau_10mo", ident.2 = "CA1_tau_20mo", verbose = FALSE) 
head(ca1.response, n = 15)
FeaturePlot(spatial_merged, features = c("Arf3"), split.by = "age", cols = c("grey","red"), reduction = "umap")

# DIFFERENTIAL EXPRESSION IN DG
# WT 10mo vs 10mo
dg.response <- FindMarkers(spatial_merged, ident.1 = "Dentate Gyrus_wt_10mo", ident.2 = "Dentate Gyrus_tau_10mo", verbose = FALSE) 
head(dg.response, n = 15)
FeaturePlot(spatial_merged, features = c("Cd81"), split.by = "condition", cols = c("grey","red"), reduction = "umap")

# Wt 20mo vs tau 20mo
dg.response <- FindMarkers(spatial_merged, ident.1 = "Dentate Gyrus_wt_20mo", ident.2 = "Dentate Gyrus_tau_20mo", verbose = FALSE) 
head(dg.response, n = 15)
FeaturePlot(spatial_merged, features = c("Prnp"), split.by = "condition", cols = c("grey","red"), reduction = "umap")

# Tau 10mo vs tau 20mo
dg.response <- FindMarkers(spatial_merged, ident.1 = "Dentate Gyrus_tau_10mo", ident.2 = "Dentate Gyrus_tau_20mo", verbose = FALSE) 
head(dg.response, n = 15)
FeaturePlot(spatial_merged, features = c("Arf3"), split.by = "age", cols = c("grey","red"), reduction = "umap")

# DIFFERENTIAL EXPRESSION IN CA2/3
# WT 10mo vs 10mo
ca23.response <- FindMarkers(spatial_merged, ident.1 = "CA2/3_wt_10mo", ident.2 = "CA2/3_tau_10mo", verbose = FALSE) 
head(ca23.response, n = 15)
FeaturePlot(spatial_merged, features = c("Eif2s3y"), split.by = "condition", cols = c("grey","red"), reduction = "umap")

# Wt 20mo vs tau 20mo
ca23.response <- FindMarkers(spatial_merged, ident.1 = "CA2/3_wt_20mo", ident.2 = "CA2/3_tau_20mo", verbose = FALSE) 
head(ca23.response, n = 15)
FeaturePlot(spatial_merged, features = c("Ghitm"), split.by = "condition", cols = c("grey","red"), reduction = "umap")

# Tau 10mo vs tau 20mo
ca23.response <- FindMarkers(spatial_merged, ident.1 = "CA2/3_tau_10mo", ident.2 = "CA2/3_tau_20mo", verbose = FALSE) 
head(ca23.response, n = 15)
FeaturePlot(spatial_merged, features = c("Aplp1"), split.by = "age", cols = c("grey","red"), reduction = "umap")

rm(ca1.response, ca23.response, dg.response, fiber.response, hypo.response, l23_cortex.response, l56_cortex.response, olf.response, striatum.response, thalamus.response)
gc()

#--------------------------CELLCHAT-----------------------------------------------------------------
spatial_merged$condition_age <- paste(spatial_merged$condition, spatial_merged$age, sep = "_")
spatial_merged$cellchat_region <- make.names(as.character(spatial_merged$merged_region))

# Use subset of CellChatDB for cell-cell communciation analysis
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")

cellchat_list <- list()

for (sample in samples) {
  cells.use <- colnames(spatial_merged)[spatial_merged$sample == sample]
  obj <- subset(spatial_merged, cells = cells.use)
  
  cellchat <- createCellChat(obj)
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  cellchat_list[[sample]] <- cellchat
  rm(cellchat, obj)
}

#---------Visualize interactions and pathways for tau 10mo----------------
cellchat_tau_10 <- cellchat_list$visium_tau_10mo
gg1 <- netVisual_heatmap(cellchat_tau_10)
gg2 <- netVisual_heatmap(cellchat_tau_10, measure = "weight")
gg1+gg2

cellchat_tau_10@netP$pathways
par(mfrow=c(1,1))

# Visualize individual signaling pathways in each sample
netVisual_heatmap(cellchat_tau_10, signaling = "PSAP", color.heatmap = "Reds")
netAnalysis_contribution(cellchat_tau_10, signaling = "PSAP")

netVisual_heatmap(cellchat_tau_10, signaling = "PTN", color.heatmap = "Reds")
netAnalysis_contribution(cellchat_tau_10, signaling = "PTN")

netVisual_heatmap(cellchat_tau_10, signaling = "SLITRK", color.heatmap = "Reds")
netAnalysis_contribution(cellchat_tau_10, signaling = "SLITRK")

netVisual_heatmap(cellchat_tau_10, signaling = "NT", color.heatmap = "Reds")
netAnalysis_contribution(cellchat_tau_10, signaling = "NT")

netVisual_heatmap(cellchat_tau_10, signaling = "NPY", color.heatmap = "Reds")
netAnalysis_contribution(cellchat_tau_10, signaling = "NPY")

netVisual_heatmap(cellchat_tau_10, signaling = "WNT", color.heatmap = "Reds")
netAnalysis_contribution(cellchat_tau_10, signaling = "WNT")

netVisual_heatmap(cellchat_tau_10, signaling = "SOMATOSTATIN", color.heatmap = "Reds")
netAnalysis_contribution(cellchat_tau_10, signaling = "SOMATOSTATIN")

#---------Visualize interactions and pathways for tau 20mo----------------
cellchat_tau_20 <- cellchat_list$visium_tau_20mo
gg1 <- netVisual_heatmap(cellchat_tau_20)
gg2 <- netVisual_heatmap(cellchat_tau_20, measure = "weight")
gg1+gg2

cellchat_tau_20@netP$pathways
par(mfrow=c(1,1))

# Visualize individual signaling pathways in each sample
netVisual_heatmap(cellchat_tau_20, signaling = "PSAP", color.heatmap = "Reds")
netAnalysis_contribution(cellchat_tau_20, signaling = "PSAP")

netVisual_heatmap(cellchat_tau_20, signaling = "PTN", color.heatmap = "Reds")
netAnalysis_contribution(cellchat_tau_20, signaling = "PTN")

netVisual_heatmap(cellchat_tau_20, signaling = "SLITRK", color.heatmap = "Reds")
netAnalysis_contribution(cellchat_tau_20, signaling = "SLITRK")

netVisual_heatmap(cellchat_tau_20, signaling = "NT", color.heatmap = "Reds")
netAnalysis_contribution(cellchat_tau_20, signaling = "NT")

netVisual_heatmap(cellchat_tau_20, signaling = "FGF", color.heatmap = "Reds")
netAnalysis_contribution(cellchat_tau_20, signaling = "FGF")

netVisual_heatmap(cellchat_tau_20, signaling = "WNT", color.heatmap = "Reds")
netAnalysis_contribution(cellchat_tau_20, signaling = "WNT")

netVisual_heatmap(cellchat_tau_20, signaling = "SOMATOSTATIN", color.heatmap = "Reds")
netAnalysis_contribution(cellchat_tau_20, signaling = "SOMATOSTATIN")

netVisual_heatmap(cellchat_tau_20, signaling = "PDGF", color.heatmap = "Reds")
netAnalysis_contribution(cellchat_tau_20, signaling = "PDGF")

netAnalysis_contribution(cellchat_tau_20, signaling = "CSF")

#---------Visualize interactions and pathways for wt 10mo----------------
cellchat_wt_10 <- cellchat_list$visium_wt_10mo
gg1 <- netVisual_heatmap(cellchat_wt_10)
gg2 <- netVisual_heatmap(cellchat_wt_10, measure = "weight")
gg1+gg2

cellchat_wt_10@netP$pathways
par(mfrow=c(1,1))

# Visualize individual signaling pathways in each sample
netVisual_heatmap(cellchat_wt_10, signaling = "PSAP", color.heatmap = "Reds")
netAnalysis_contribution(cellchat_wt_10, signaling = "PSAP")

netVisual_heatmap(cellchat_wt_10, signaling = "PTN", color.heatmap = "Reds")
netAnalysis_contribution(cellchat_wt_10, signaling = "PTN")

netVisual_heatmap(cellchat_wt_10, signaling = "SLITRK", color.heatmap = "Reds")
netAnalysis_contribution(cellchat_wt_10, signaling = "SLITRK")

netVisual_heatmap(cellchat_wt_10, signaling = "NT", color.heatmap = "Reds")
netAnalysis_contribution(cellchat_wt_10, signaling = "NT")

netVisual_heatmap(cellchat_wt_10, signaling = "WNT", color.heatmap = "Reds")
netAnalysis_contribution(cellchat_wt_10, signaling = "WNT")

netVisual_heatmap(cellchat_wt_10, signaling = "SOMATOSTATIN", color.heatmap = "Reds")
netAnalysis_contribution(cellchat_wt_10, signaling = "SOMATOSTATIN")

netAnalysis_contribution(cellchat_wt_10, signaling = "PMCH")

#---------Visualize interactions and pathways for wt 20mo----------------
cellchat_wt_20 <- cellchat_list$visium_wt_20mo
gg1 <- netVisual_heatmap(cellchat_wt_20)
gg2 <- netVisual_heatmap(cellchat_wt_20, measure = "weight")
gg1+gg2

cellchat_wt_20@netP$pathways
par(mfrow=c(1,1))

# Visualize individual signaling pathways in each sample
netVisual_heatmap(cellchat_wt_20, signaling = "PSAP", color.heatmap = "Reds")
netAnalysis_contribution(cellchat_wt_20, signaling = "PSAP")

netVisual_heatmap(cellchat_wt_20, signaling = "PTN", color.heatmap = "Reds")
netAnalysis_contribution(cellchat_wt_20, signaling = "PTN")

netVisual_heatmap(cellchat_wt_20, signaling = "SLITRK", color.heatmap = "Reds")
netAnalysis_contribution(cellchat_wt_20, signaling = "SLITRK")

netVisual_heatmap(cellchat_wt_20, signaling = "NT", color.heatmap = "Reds")
netAnalysis_contribution(cellchat_wt_20, signaling = "NT")

netVisual_heatmap(cellchat_wt_20, signaling = "WNT", color.heatmap = "Reds")
netAnalysis_contribution(cellchat_wt_20, signaling = "WNT")

netVisual_heatmap(cellchat_wt_20, signaling = "SOMATOSTATIN", color.heatmap = "Reds")
netAnalysis_contribution(cellchat_wt_20, signaling = "SOMATOSTATIN")

netVisual_heatmap(cellchat_wt_20, signaling = "FGF", color.heatmap = "Reds")
netAnalysis_contribution(cellchat_wt_20, signaling = "FGF")

netVisual_heatmap(cellchat_wt_20, signaling = "PACAP", color.heatmap = "Reds")
netAnalysis_contribution(cellchat_wt_20, signaling = "PACAP")

netAnalysis_contribution(cellchat_wt_20, signaling = "PMCH")


#---------------- Compare interactions for each sample in merged object----------------
cellchat_merged <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))

gg1 <- compareInteractions(cellchat_merged, show.legend = FALSE)
gg2 <- compareInteractions(cellchat_merged, show.legend = FALSE, measure = "weight")
gg1 + gg2

#---------------------SAVE OUTPUT FOR SUBMISSION-------------------------------------
save.image("workspace.RData")
sessionInfo()
writeLines(capture.output(sessionInfo()), "../session_info.txt")
