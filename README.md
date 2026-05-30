# Neuro Bioinformatics Final Project

The final project for my neuroscience bioinformatics course required the creation of a multi-modal study exploring a biological question of interest using publicly available transcriptomics datasets. My project used single-cell and spatial transcriptomics to analyze and compare gene expression differences between 1) rats with AD-like tauopathy vs. their wild-type littermates and 2) early stage tauopathy (10 months) vs. late stage tauopathy (20 months). 

# Datasets

The original dataset can be [downloaded here] (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE305314), where the single-cell samples are the files labeled `snRNA-seq Sample # Condition Age` and the spatial samples are the files labeled `Visium Condition Age Sample`. The paper I used as a reference for my methods and clustering can be [found here] (https://www.nature.com/articles/s42003-025-08959-z). 

# Methods

## Single-cell

### Data & Preprocessing

Single-cell samples were loaded in based on their condition and age group. Samples were individually preprocessed and had doublets removed using the package scDblFinder. Afterwards, all samples were merged into a Seurat object with their condition and age kept in the metadata before defining a subset based on nFeature_RNA and percent.mt. Then, the object was preprocessed and clustered using dimensions 1:19 and clustering resolution 0.3, resulting in 25 clusters. 

<p>
    <img src="/sc_img/sc_pre_harmony_clusters.png" width="600" height="300" />
</p>

### Harmony Integration

The Seurat object using Harmony for batch correction using the dimensions 1:19 and a clustering resolution of 0.25, resulting in 22 clusters.

### Cell-Type Annotation

Top gene markers for each cluster were extracted using FindAllMarkers and clusters were annotated based on their cell type using [Enrichr] (https://maayanlab.cloud/Enrichr/#find). After finding the cell types of each cluster, clusters were annotated with their cell type.

<p>
    <img src="/sc_img/sc_annotated_clusters.png" />
</p>

### Proportion Test

The proportion of cell types were compared between condition and age using the package scProportionTest with permutation-based testing. 

<p>
    <img src="/sc_img/sc_condition_proportion.png" />
</p>

<p>
    <img src="/sc_img/sc_age_proportion.png" />
</p>

### Differential Expression

Top markers for each cluster split on condition then age were identified then visualized to compare differential expression between the groups.

## Spatial

### Data & Preprocessing

Datasets were loaded in and converted to HDF5 file type using Load10X_Spatial. Afterwards, the individual samples were loaded in using Load10X_Spatial with their metadata containing the sample's condition and age. SCTransform was used to normalize individual samples before merging them into one Seurat object. The object was then clustered with dimensions 1:40 and clustering resolution 0.5, resulting in 19 clusters. 

### Cluster Annotation

Spatially variable features were identified using Moran's I. Top gene markers for each cluster were paired with a region/cell of interest using Enrichr, while clusters not deemed to be of interest were grouped under the label "Other". 11 clusters of interest were identified. 

<p>
    <img src="/spatial_img/spatial_clusters_annotated.png" />
</p>

### Differential Expression

Top markers for each cluster split on condition then age were identified and visualized to compare differential expression between the groups. 

### Cell-Cell Communication

CellChat using the mouse database was using to analyze cell-cell communication. Interactions and pathways were identified and visualized for each sample. Afterwards, all cellchat objects were merged and inferred interactions between samples were identified and visualized. 

<p>
    <img src="/spatial_img/spatial_interactions.png" />
</p>

# Limitations & Challenges

Although I tried to follow the original paper as closely as possible, there were a few issues that made my project not as accurate/informative. 

1) The original paper used SCTransform to normalize both the single-cell and spatial data; however, I ran into memory issues when trying to use it on the single-cell data, so I was only able to use it to normalize the spatial datasets.

2) The original datasets also included scATAC data; however, I did not have the knowlege to analyze that type of data at the time of this project and was unable to effectively learn it due to time constraints. 

3) Cluster labels were similar but not identical to the paper's clusters; I had messed with different dimensions and clustering resolutions for several days in an attempt to get as close to the paper as possible, but ultimately stuck with my final parameters due to time constraints.

4) The original paper cleared up ambient RNA contamination and doublets more effectively than my project; I had tried using DoubletFInder and SoupX but ran into memory issues and  not enough cells when trying to use them. I found ScDblFinder to be the next best thing to detect and remove doublets from the datasets.
