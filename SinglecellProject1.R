# Setting the seed for reproducibility
set.seed(42)

# Load required packages
library(Seurat)
library(DoubletFinder)
library(SingleR)
library(ggplot2)
library(dplyr)
library(cowplot)
library(tibble)
library(ggpubr)
library(enrichR)
library(monocle3)
library(CellChat)


#NOTE
#Merged_Datasets.integrated is the seurat_obj in our code.


# Define the base directory path where your .rds files are stored
base_dir <- "/home/sare00008/ssp/scbi_ds1"


# Load data from .rds files
bmmc_d1t1 <- readRDS(file.path(base_dir, "GSM4138872_scRNA_BMMC_D1T1.rds"))
bmmc_d1t2 <- readRDS(file.path(base_dir, "GSM4138873_scRNA_BMMC_D1T2.rds"))
cd34_d2t1 <- readRDS(file.path(base_dir, "GSM4138874_scRNA_CD34_D2T1.rds"))
cd34_d3t1 <- readRDS(file.path(base_dir, "GSM4138875_scRNA_CD34_D3T1.rds"))

# Creating Seurat objects with quality control filters
bmmc_d1t1_seurat <- CreateSeuratObject(counts = bmmc_d1t1, project = "BMMC_D1T1")
bmmc_d1t2_seurat <- CreateSeuratObject(counts = bmmc_d1t2, project = "BMMC_D1T2")
cd34_d2t1_seurat <- CreateSeuratObject(counts = cd34_d2t1, project = "CD34_D2T1")
cd34_d3t1_seurat <- CreateSeuratObject(counts = cd34_d3t1, project = "CD34_D3T1")

# Create a sample sheet based on Table 1 metadata
sample_sheet <- data.frame(
  SampleID = c("BMMC_D1T1", "BMMC_D1T2", "CD34_D2T1", "CD34_D3T1"),
  Donor = c("D1", "D1", "D2", "D3"),
  Replicate = c("T1", "T2", "T1", "T1"),
  Sex = c("F", "F", "M", "F"),
  stringsAsFactors = FALSE
)

# Display the sample sheet
print(sample_sheet)

# Add metadata from sample sheet to each Seurat object
bmmc_d1t1_seurat <- AddMetaData(bmmc_d1t1_seurat, metadata = list(
  Donor = rep("D1", ncol(bmmc_d1t1_seurat)),
  Replicate = rep("T1", ncol(bmmc_d1t1_seurat)),
  Sex = rep("F", ncol(bmmc_d1t1_seurat))
))

bmmc_d1t2_seurat <- AddMetaData(bmmc_d1t2_seurat, metadata = list(
  Donor = rep("D1", ncol(bmmc_d1t2_seurat)),
  Replicate = rep("T2", ncol(bmmc_d1t2_seurat)),
  Sex = rep("F", ncol(bmmc_d1t2_seurat))
))

cd34_d2t1_seurat <- AddMetaData(cd34_d2t1_seurat, metadata = list(
  Donor = rep("D2", ncol(cd34_d2t1_seurat)),
  Replicate = rep("T1", ncol(cd34_d2t1_seurat)),
  Sex = rep("M", ncol(cd34_d2t1_seurat))
))

cd34_d3t1_seurat <- AddMetaData(cd34_d3t1_seurat, metadata = list(
  Donor = rep("D3", ncol(cd34_d3t1_seurat)),
  Replicate = rep("T1", ncol(cd34_d3t1_seurat)),
  Sex = rep("F", ncol(cd34_d3t1_seurat))
))

# Printing Seurat objects to view the summary
print(bmmc_d1t1_seurat)
print(bmmc_d1t2_seurat)
print(cd34_d2t1_seurat)
print(cd34_d3t1_seurat)

# Display the first few rows of the metadata for the Seurat object
# Example for bmmc_d1t1_seurat
#print(head(bmmc_d1t1_seurat@meta.data))

###If you want to see for other cells then just replace it with bmmc_d1t1###


# Week 2 Code

print("Week 2 preprocessing started...")

## Preprocessing of data

# Filtering - Step 1: Calculation of mitochondrial gene percentage
print("Step 1: Calculating mitochondrial gene percentage for each Seurat object...")
bmmc_d1t1_seurat[["percent.mt"]] <- PercentageFeatureSet(bmmc_d1t1_seurat, pattern = "^MT-")
bmmc_d1t2_seurat[["percent.mt"]] <- PercentageFeatureSet(bmmc_d1t2_seurat, pattern = "^MT-")
cd34_d2t1_seurat[["percent.mt"]] <- PercentageFeatureSet(cd34_d2t1_seurat, pattern = "^MT-")
cd34_d3t1_seurat[["percent.mt"]] <- PercentageFeatureSet(cd34_d3t1_seurat, pattern = "^MT-")

print("Step 1 completed.")

# Normalizing each Seurat object to create the "data" layer
bmmc_d1t1_seurat <- NormalizeData(bmmc_d1t1_seurat)
bmmc_d1t2_seurat <- NormalizeData(bmmc_d1t2_seurat)
cd34_d2t1_seurat <- NormalizeData(cd34_d2t1_seurat)
cd34_d3t1_seurat <- NormalizeData(cd34_d3t1_seurat)

# Step 2: Violin plot for QC metrics
print("Step 2: Creating QC Violin plots...")
plt <- VlnPlot(bmmc_d1t1_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(filename = "BMMC_D1T1_QC.png", plot = plt, width = 7, height = 3.5)

plt_bmmc_d1t2 <- VlnPlot(bmmc_d1t2_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(filename = "BMMC_D1T2_QC.png", plot = plt_bmmc_d1t2, width = 7, height = 3.5)

plt_cd34_d2t1 <- VlnPlot(cd34_d2t1_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(filename = "CD34_D2T1_QC.png", plot = plt_cd34_d2t1, width = 7, height = 3.5)

plt_cd34_d3t1 <- VlnPlot(cd34_d3t1_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(filename = "CD34_D3T1_QC.png", plot = plt_cd34_d3t1, width = 7, height = 3.5)

print("Step 2 completed.")

# Step 3: Filtering based on nFeature_RNA and nCount_RNA
print("Step 3: Filtering cells based on QC metrics...")
bmmc_d1t1_seurat <- subset(bmmc_d1t1_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & nCount_RNA > 200 & nCount_RNA < 5500)
bmmc_d1t2_seurat <- subset(bmmc_d1t2_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & nCount_RNA > 200 & nCount_RNA < 5500)
cd34_d2t1_seurat <- subset(cd34_d2t1_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & nCount_RNA > 200 & nCount_RNA < 7500)
cd34_d3t1_seurat <- subset(cd34_d3t1_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & nCount_RNA > 200 & nCount_RNA < 7500)

print("Step 3 completed.")

# Function to run DoubletFinder on a Seurat object
run_doubletfinder <- function(seurat_obj, PCs = 1:10, pN = 0.25, doublet_rate = 0.075) {
  # Step 1: Find Variable Features
  print(paste("Finding variable features for:", seurat_obj@project.name))
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  
  # Step 2: Scaling and PCA for dimension reduction
  print(paste("Scaling data and running PCA for:", seurat_obj@project.name))
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, npcs = 20)
  
  # Step 3: Run parameter sweep to find optimal pK value
  print("Running paramSweep to determine optimal pK...")
  sweep.res.list <- paramSweep(seurat_obj, PCs = PCs, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  # Step 4: Determine optimal pK based on BCmetric
  optimal_pk <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  print(paste("Optimal pK determined:", optimal_pk))
  
  # Step 5: Visualize pK selection (BCmetric vs pK plot)
  pk_plot <- ggplot(bcmvn, aes(x = as.numeric(as.character(pK)), y = BCmetric)) +
    geom_point() +
    geom_line() +
    xlab("pK") +
    ylab("BCmetric") +
    ggtitle(paste0(seurat_obj@project.name, ": pK Optimization Plot"))
  ggsave(filename = paste0(seurat_obj@project.name, "_DoubletFinder_pK_Plot.png"), plot = pk_plot, width = 7, height = 3.5)
  
  # Step 6: Estimate the number of doublets (nExp)
  nExp_poi <- round(doublet_rate * ncol(seurat_obj))
  
  # Step 7: Run DoubletFinder
  print("Running DoubletFinder...")
  seurat_obj <- doubletFinder(seurat_obj, PCs = PCs, pN = pN, pK = optimal_pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  
  # Step 8: Identify the doublet classification column in metadata
  doublet_col <- colnames(seurat_obj@meta.data)[grep("DF.classifications", colnames(seurat_obj@meta.data))]
  
  # Step 9: Print the number of doublets and singlets
  doublet_summary <- table(seurat_obj@meta.data[[doublet_col]])
  print(doublet_summary)
  
  # Specifically print the number of doublets
  num_doublets <- doublet_summary["Doublet"]
  print(paste("Number of doublets identified:", num_doublets))
  
  
  # Step 10: Visualize singlet vs doublet classification using UMAP
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
  umap_plot <- DimPlot(seurat_obj, group.by = doublet_col) +
    ggtitle(paste0("UMAP: Singlet vs Doublet Classification - ", seurat_obj@project.name))
  ggsave(filename = paste0(seurat_obj@project.name, "_DoubletFinder_Singlet_Doublet_UMAP.png"), plot = umap_plot, width = 7, height = 3.5)
  
  # Step 11: Remove doublets from Seurat object
  print("Removing doublets...")
  seurat_obj <- subset(seurat_obj, cells = rownames(seurat_obj@meta.data)[seurat_obj@meta.data[, doublet_col] == "Singlet"])
  
  return(seurat_obj)
}

# Applying DoubletFinder to each Seurat object individually before merging
bmmc_d1t1_seurat <- run_doubletfinder(bmmc_d1t1_seurat)
bmmc_d1t2_seurat <- run_doubletfinder(bmmc_d1t2_seurat)
cd34_d2t1_seurat <- run_doubletfinder(cd34_d2t1_seurat)
cd34_d3t1_seurat <- run_doubletfinder(cd34_d3t1_seurat)



# Merging data
Merged_Datasets <- merge(bmmc_d1t1_seurat, y = c(bmmc_d1t2_seurat, cd34_d2t1_seurat, cd34_d3t1_seurat),
                         add.cell.ids = c("BMMC_D1T1", "BMMC_D1T2", "CD34_D2T1", "CD34_D3T1"), project = 'MD')

# Normalize and scale merged data before batch correction
Merged_Datasets <- NormalizeData(Merged_Datasets)
Merged_Datasets <- FindVariableFeatures(Merged_Datasets, selection.method = "vst", nfeatures = 2000)
Merged_Datasets <- ScaleData(Merged_Datasets, features = rownames(Merged_Datasets))

# PCA for merged data
Merged_Datasets <- RunPCA(Merged_Datasets)

# Elbow plot
Merged_Dataset_EP <- ElbowPlot(Merged_Datasets)
ggsave(filename = "Merged_Dataset_elbowplot.png", plot = Merged_Dataset_EP, width = 7, height = 3.5)

# Clustering and UMAP before batch correction
Merged_Datasets <- FindNeighbors(Merged_Datasets, dims = 1:18)
Merged_Datasets <- FindClusters(Merged_Datasets, resolution = 0.3)
Merged_Datasets <- RunUMAP(Merged_Datasets, dims = 1:18)    #selected 1 to 18

# Plotting UMAP for merged data before batch correction (seurat_clusters)
umap_before_bc_clusters <- DimPlot(Merged_Datasets, reduction = "umap") +
  ggtitle("UMAP Before Batch Correction: Clusters")
ggsave(filename = "Merged_Datasets_UMAP_Before_Batch_Correction_Clusters.png", plot = umap_before_bc_clusters, width = 7, height = 3.5)


# Plotting UMAP based on orig.ident before batch correction (BMMC vs CD34)
umap_before_bc_orig_ident <- DimPlot(Merged_Datasets, reduction = "umap", group.by = "orig.ident") +
  ggtitle("UMAP Before Batch Correction: BMMC vs CD34")
ggsave(filename = "Merged_Datasets_UMAP_Before_Batch_Correction_Orig_Ident.png", plot = umap_before_bc_orig_ident, width = 7, height = 3.5)


# Batch Correction - Integration

# Combine all datasets into a list for integration
dataset_list <- list(bmmc_d1t1_seurat, bmmc_d1t2_seurat, cd34_d2t1_seurat, cd34_d3t1_seurat)

# Step 1: Normalizing each dataset individually
dataset_list <- lapply(X = dataset_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  return(x)
})

# Step 2: Finding integration anchors
integration_anchors <- FindIntegrationAnchors(object.list = dataset_list, anchor.features = 2000)

# Step 3: Integrating the datasets
integrated_data <- IntegrateData(anchorset = integration_anchors)

# Set the integrated assay as the default to work with corrected data
DefaultAssay(integrated_data) <- "integrated"

# Step 4: Scaling the integrated data
integrated_data <- ScaleData(integrated_data, features = rownames(integrated_data))

# Step 5: Running PCA on the integrated data
integrated_data <- RunPCA(integrated_data)

# Step 6: Clustering and UMAP for integrated data
integrated_data <- FindNeighbors(integrated_data, dims = 1:18)
integrated_data <- FindClusters(integrated_data, resolution = 0.3)
integrated_data <- RunUMAP(integrated_data, dims = 1:18)

# Explicitly save clustering results to metadata
integrated_data$seurat_clusters <- Idents(integrated_data) # Save cluster assignments in metadata


# Plotting UMAP for merged data after batch correction (seurat_clusters)
umap_after_bc_clusters <- DimPlot(integrated_data, reduction = "umap") +
  ggtitle("UMAP After Batch Correction: Clusters")
ggsave(filename = "Integrated_Dataset_UMAP_After_Batch_Correction_Clusters.png", plot = umap_after_bc_clusters, width = 7, height = 3.5)

# Plotting UMAP based on orig.ident after batch correction (BMMC vs CD34)
umap_after_bc_orig_ident <- DimPlot(integrated_data, reduction = "umap", group.by = "orig.ident") +
  ggtitle("UMAP After Batch Correction: BMMC vs CD34")
ggsave(filename = "Integrated_Dataset_UMAP_After_Batch_Correction_Orig_Ident.png", plot = umap_after_bc_orig_ident, width = 7, height = 3.5)



###WEEK 3###


#6.1
# Loading the reference dataset for annotation using SingleR
reference_path <- "/home/sare00008/ssp/celldex_annot.RDS" 
hpca <- readRDS(reference_path)

# Extract normalized data matrix from the integrated Seurat object for annotation
integrated_data_matrix <- GetAssayData(integrated_data, assay = "integrated", slot = "data")

# Use SingleR to annotate cells using the custom reference
print("Running SingleR for automatic cell annotation...")
singleR_results <- SingleR(
  test = integrated_data_matrix,
  ref = hpca,
  labels = hpca$label.main)

# Add the SingleR labels to the metadata of the Seurat object
integrated_data$SingleR_label <- singleR_results$labels

# Plot UMAP with SingleR annotation
umap_annotation_plot <- DimPlot(integrated_data, reduction = "umap", group.by = "SingleR_label", label = TRUE, repel= TRUE) +
  ggtitle("UMAP Plot with Cell Type Annotation by SingleR")

# Save the UMAP plot
ggsave(filename = "UMAP_SingleR_CellType_Annotation.png", plot = umap_annotation_plot, width = 8, height = 5)

#6.2

# Rename integrated_data to Merged_Datasets.integrated to maintain consistency
Merged_Datasets.integrated <- integrated_data

# Perform Differential Expression Analysis
cluster_markers <- FindAllMarkers(Merged_Datasets.integrated, only.pos = TRUE)
View(cluster_markers)

# Visualizing Marker Expression for Cell-Type Assignment

# Hematopoietic Stem Cells (HSCs) - Clusters 11 and 13
HSC_plot <- FeaturePlot(
  Merged_Datasets.integrated, 
  features = c("CD38"), 
  label = TRUE
) + ggtitle("Hematopoietic Stem Cells (HSCs)")
ggsave(filename = "HSC_Markers_11_13.png", plot = HSC_plot, width = 7, height = 4)

# Lymphoid-Primed Multipotent Progenitors (LMPP) - Cluster 2
LMPP_plot <- FeaturePlot(
  Merged_Datasets.integrated, 
  features = c("CSF3R", "CD52"), 
  label = TRUE
) + ggtitle("LMPP Markers - Cluster 2")
ggsave(filename = "LMPP_Markers_Cluster_2.png", plot = LMPP_plot, width = 7, height = 4)

# B Cells - Cluster 6
BCells_plot <- FeaturePlot(
  Merged_Datasets.integrated, 
  features = c("CD19"), 
  label = TRUE
) + ggtitle("B Cell Markers - Cluster 6")
ggsave(filename = "BCells_Markers_Cluster_6.png", plot = BCells_plot, width = 7, height = 4)

# Plasma Cells - Cluster 9
Plasma_plot <- FeaturePlot(
  Merged_Datasets.integrated, 
  features = c("MZB1"), 
  label = TRUE
) + ggtitle("Plasma Cell Markers - Cluster 9")
ggsave(filename = "Plasma_Markers_Cluster_9.png", plot = Plasma_plot, width = 7, height = 4)

# CD8+ T Cells - Clusters 0, 3, and 8
CD8_T_plot <- FeaturePlot(
  Merged_Datasets.integrated, 
  features = c("CD3D", "CD3E", "CD8A", "CD8B"), 
  label = TRUE
) + ggtitle("CD8+ T Cell Markers - Clusters 0, 3, 8")
ggsave(filename = "CD8_T_Markers_Clusters.png", plot = CD8_T_plot, width = 7, height = 4)

# Natural Killer (NK) Cells - Clusters 4 and 12
NK_plot <- FeaturePlot(
  Merged_Datasets.integrated, 
  features = c("FCGR3A", "NKG7", "KLRB1"), 
  label = TRUE
) + ggtitle("Natural Killer Cell Markers - Clusters 4 and 12")
ggsave(filename = "NK_Markers_Clusters_4_12.png", plot = NK_plot, width = 7, height = 4)

# Plasmacytoid Dendritic Cells (pDCs) - Cluster 10
pDC_plot <- FeaturePlot(
  Merged_Datasets.integrated, 
  features = c("IRF8", "IRF4", "IRF7"), 
  label = TRUE
) + ggtitle("pDC Markers - Cluster 10")
ggsave(filename = "pDC_Markers_Cluster_10.png", plot = pDC_plot, width = 7, height = 4)

# Conventional Dendritic Cells (cDCs) - Cluster 7
cDC_plot <- FeaturePlot(
  Merged_Datasets.integrated, 
  features = c("ITGAM", "SIRPA"), 
  label = TRUE
) + ggtitle("cDC Markers - Cluster 7")
ggsave(filename = "cDC_Markers_Cluster_7.png", plot = cDC_plot, width = 7, height = 4)

# Monocytes - Cluster 1
Monocyte_plot <- FeaturePlot(
  Merged_Datasets.integrated, 
  features = c("CD68", "S100A12"), 
  label = TRUE
) + ggtitle("Monocyte Markers - Cluster 1")
ggsave(filename = "Monocyte_Markers_Cluster_1.png", plot = Monocyte_plot, width = 7, height = 4)

# Basophils - Cluster 5
Basophils_plot <- FeaturePlot(
  Merged_Datasets.integrated, 
  features = c("GATA2"), 
  label = TRUE
) + ggtitle("Basophil Markers - Cluster 5")
ggsave(filename = "Basophil_Markers_Cluster_5.png", plot = Basophils_plot, width = 7, height = 4)

# Assigning Cell-Type Names to Clusters
DefaultAssay(Merged_Datasets.integrated) <- "RNA"

Merged_Datasets.integrated <- RenameIdents(Merged_Datasets.integrated, 
                                           c(
                                             "0" = "CD8_T", 
                                             "1" = "Monocytes", 
                                             "2" = "LMPP", 
                                             "3" = "CD8_T", 
                                             "4" = "NK", 
                                             "5" = "Basophils", 
                                             "6" = "B_Cells", 
                                             "7" = "cDCs", 
                                             "8" = "CD8_T", 
                                             "9" = "Plasma", 
                                             "10" = "pDC", 
                                             "11" = "HSC", 
                                             "12" = "NK", 
                                             "13" = "HSC"
                                           )
)

# UMAP Plot with Manual Cell-Type Annotation
UMAP_Manual_Annotation <- DimPlot(
  Merged_Datasets.integrated, 
  reduction = "umap", 
  label = TRUE, 
  pt.size = 0.5
) + ggtitle("UMAP with Manual Cell-Type Annotations") + NoLegend()
ggsave(filename = "UMAP_Manual_Cell_Type_Annotation.png", plot = UMAP_Manual_Annotation, width = 7, height = 4)

## Combine the two plots side by side
combined_plot <- plot_grid(
  umap_annotation_plot, 
  UMAP_Manual_Annotation, 
  labels = c("A", "B"), # Label the plots as A and B
  label_size = 15,      # Adjust the label size
  ncol = 2              # Arrange the plots in a single row
)

# Save the combined plot
ggsave(filename = "UMAP_Annotation_Comparison.png", 
       plot = combined_plot, 
       width = 14, 
       height = 7)

# Display the combined plot
print(combined_plot)

## Set marker genes to visualize
marker_genes <- c("CD8A", "CD3D", "CD19")  # Example markers for T cells and B cells

# Violin Plot for Marker Genes
violin_plot <- VlnPlot(
  object = integrated_data, 
  features = marker_genes, 
  group.by = "seurat_clusters", 
  pt.size = 0.1, 
  ncol = 3
) 
ggsave(filename = "Violin_Plot_Marker_Genes.png", plot = violin_plot, width = 12, height = 5)

# UMAP Plot for Each Marker Gene
for (gene in marker_genes) {
  umap_plot <- FeaturePlot(
    object = integrated_data, 
    features = gene, 
    reduction = "umap", 
    cols = c("lightgray", "blue"),  # Gradient for expression levels
    label = TRUE
  ) + ggtitle(paste("UMAP Plot for", gene, "Expression"))
  ggsave(filename = paste0("UMAP_Plot_", gene, ".png"), plot = umap_plot, width = 8, height = 5)
}

#6,3 Bonus
# Compute cell-type proportions per sample
cell_type_proportions <- Merged_Datasets.integrated@meta.data %>%
  dplyr::group_by(orig.ident, SingleR_label) %>%
  dplyr::summarize(count = n()) %>%
  dplyr::mutate(proportion = count / sum(count))

# Print the computed proportions for verification
print(cell_type_proportions)

# Plot cell-type proportions as a bar plot
bar_plot <- ggplot(cell_type_proportions, aes(x = orig.ident, y = proportion, fill = SingleR_label)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Cell-Type Proportions Across Samples",
    x = "Sample",
    y = "Proportion",
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the bar plot
ggsave("Cell_Type_Proportions_Bar_Plot.png", plot = bar_plot, width = 8, height = 6, dpi = 300)

# Display the bar plot
print(bar_plot)

# Compute summary statistics to explain how samples vary in terms of cell-type proportions
proportion_summary <- cell_type_proportions %>%
  dplyr::group_by(SingleR_label) %>%
  dplyr::summarize(
    mean_proportion = mean(proportion),
    sd_proportion = sd(proportion)
  ) %>%
  arrange(desc(mean_proportion))

# Print summary statistics
print("Summary of cell-type proportions across samples:")
print(proportion_summary)


###7
# Set default assay to RNA for differential expression analysis
DefaultAssay(Merged_Datasets.integrated) <- "RNA"

# Join layers if required
Merged_Datasets.integrated <- JoinLayers(Merged_Datasets.integrated)

# 1. B Cells vs T Cells: Differential Expression Analysis**
b_vs_t_markers <- FindMarkers(
  Merged_Datasets.integrated, 
  ident.1 = "B_Cells", 
  ident.2 = "CD8_T", 
  logfc.threshold = 0.25, 
  min.pct = 0.1
)

# Save results to a CSV file
write.csv(b_vs_t_markers, "B_vs_T_Markers.csv")

# Volcano plot for B cells vs T cells
library(ggplot2)
volcano_b_vs_t <- ggplot(b_vs_t_markers, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = avg_log2FC > 0), alpha = 0.7) +
  scale_color_manual(values = c("red", "blue")) +
  labs(title = "B Cells vs T Cells", x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value") +
  theme_minimal()

# Save volcano plot
ggsave(filename = "Volcano_B_vs_T.png", plot = volcano_b_vs_t, width = 7, height = 4)

# 2. T Cells vs Monocytes: Differential Expression Analysis**
t_vs_mono_markers <- FindMarkers(
  Merged_Datasets.integrated, 
  ident.1 = "CD8_T", 
  ident.2 = "Monocytes", 
  logfc.threshold = 0.25, 
  min.pct = 0.1
)

# Save results to a CSV file
write.csv(t_vs_mono_markers, "T_vs_Monocytes_Markers.csv")

# Volcano plot for T cells vs Monocytes
volcano_t_vs_mono <- ggplot(t_vs_mono_markers, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = avg_log2FC > 0), alpha = 0.7) +
  scale_color_manual(values = c("red", "blue")) +
  labs(title = "T Cells vs Monocytes", x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value") +
  theme_minimal()

# Save volcano plot
ggsave(filename = "Volcano_T_vs_Monocytes.png", plot = volcano_t_vs_mono, width = 7, height = 4)


#7.2

# Top 5 markers for B Cells vs T Cells
top_b_vs_t <- b_vs_t_markers %>%
  rownames_to_column(var = "gene") %>%
  arrange(p_val_adj) %>%
  slice_head(n = 5)

# Add comparison group for later plotting
top_b_vs_t$comparison <- "B_Cells vs T_Cells"

# Top 5 markers for T Cells vs Monocytes
top_t_vs_mono <- t_vs_mono_markers %>%
  rownames_to_column(var = "gene") %>%
  arrange(p_val_adj) %>%
  slice_head(n = 5)

# Add comparison group for later plotting
top_t_vs_mono$comparison <- "T_Cells vs Monocytes"

# Combine both comparisons into a single dataframe
top_markers_combined <- bind_rows(top_b_vs_t, top_t_vs_mono)

# Calculate -log10(p-value) for significance
top_markers_combined <- top_markers_combined %>%
  mutate(log_pval = -log10(p_val_adj))

#review the contents of the top_markers_combined
print(top_markers_combined)

library(ggplot2)


# Dynamic midpoint for color gradient
midpoint_value <- mean(top_markers_combined$avg_log2FC, na.rm = TRUE)

# Dot Plot
dot_plot <- ggplot(top_markers_combined, aes(x = comparison, y = gene)) +
  geom_point(aes(size = log_pval, color = avg_log2FC)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(
    x = "Cell Type Comparison",
    y = "Top Differentially Expressed Genes",
    size = "-log10(p-value)",
    color = "log2(Fold Change)",
    title = "Top 5 Differentially Expressed Genes"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot
ggsave("Top_DEG_Dot_Plot.png", plot = dot_plot, width = 8, height = 6)

# Display the plot
print(dot_plot)

# 8.1.1 Differential Expression Analysis for BMMC vs CD34 (All Cell Types)

# Perform Differential Expression Analysis
bmmc_vs_cd34 <- FindMarkers(
  object = Merged_Datasets.integrated, 
  ident.1 = "BMMC", 
  ident.2 = "CD34", 
  group.by = "orig.ident", 
  logfc.threshold = 0.25, 
  min.pct = 0.1
)

# Add Gene Names as a Column
bmmc_vs_cd34$gene <- rownames(bmmc_vs_cd34)

# Extract the Top 5 Differentially Expressed Genes (DEGs)
top5_bmmc_vs_cd34 <- bmmc_vs_cd34 %>%
  arrange(p_val_adj) %>%
  slice_head(n = 5) %>%
  select(gene, avg_log2FC, p_val_adj)

# Print the Top 5 DEGs
print("Top 5 Differentially Expressed Genes (BMMC vs CD34):")
print(top5_bmmc_vs_cd34)


library(ggpubr)

# Create the table as a data frame
table_bmmc_vs_cd34 <- as.data.frame(top5_bmmc_vs_cd34)

# Convert the table to a ggtexttable object
table_plot_bmmc_vs_cd34 <- ggtexttable(
  table_bmmc_vs_cd34,
  rows = NULL, # Remove row numbers
  theme = ttheme("light")
)

# Save the table as an image
ggsave(
  filename = "Top5_BMMC_vs_CD34_Table.png",
  plot = table_plot_bmmc_vs_cd34,
  width = 5,
  height = 3
)

# Display the table plot
print(table_plot_bmmc_vs_cd34)

# *8.1.2 Differential Expression Analysis for Monocytes in BMMC vs CD34

# Subset for Monocyte Cells Only
monocytes <- subset(Merged_Datasets.integrated, subset = SingleR_label == "Monocyte")

# Perform Differential Expression Analysis for Monocytes
monocyte_bmmc_vs_cd34 <- FindMarkers(
  object = monocytes,
  ident.1 = "BMMC",
  ident.2 = "CD34",
  group.by = "orig.ident",
  logfc.threshold = 0.25,
  min.pct = 0.1
)

# Add Gene Names as a Column
monocyte_bmmc_vs_cd34$gene <- rownames(monocyte_bmmc_vs_cd34)

# Extract the Top 5 Differentially Expressed Genes for Monocytes
top5_monocyte_bmmc_vs_cd34 <- monocyte_bmmc_vs_cd34 %>%
  arrange(p_val_adj) %>%
  slice_head(n = 5) %>%
  select(gene, avg_log2FC, p_val_adj)

# Print the Top 5 DEGs for Monocytes
print("Top 5 Differentially Expressed Genes in Monocytes (BMMC vs CD34):")
print(top5_monocyte_bmmc_vs_cd34)



# Create the table as a data frame
table_monocyte_bmmc_vs_cd34 <- as.data.frame(top5_monocyte_bmmc_vs_cd34)

# Convert the table to a ggtexttable object
table_plot_monocyte_bmmc_vs_cd34 <- ggtexttable(
  table_monocyte_bmmc_vs_cd34,
  rows = NULL, # Remove row numbers
  theme = ttheme("light")
)

# Save the table as an image
ggsave(
  filename = "Top5_Monocyte_BMMC_vs_CD34_Table.png",
  plot = table_plot_monocyte_bmmc_vs_cd34,
  width = 5,
  height = 3
)

# Display the table plot
print(table_plot_monocyte_bmmc_vs_cd34)


###8,.2
##group assignment

# Step 1: Check unique values in orig.ident
print("Unique values in orig.ident:")
print(unique(Merged_Datasets.integrated$orig.ident))

# Step 2: Assign groups directly
Merged_Datasets.integrated$group <- NA  # Initialize with NA

# Assign BMMC group
Merged_Datasets.integrated$group[Merged_Datasets.integrated$orig.ident %in% c("BMMC", "BMMC_D1T1", "BMMC_D1T2")] <- "BMMC"

# Assign CD34 group
Merged_Datasets.integrated$group[Merged_Datasets.integrated$orig.ident %in% c("CD34", "CD34_D2T1", "CD34_D3T1")] <- "CD34"

# Step 3: Verify group assignment
print("Group Counts:")
print(table(Merged_Datasets.integrated$group))

# Step 4: Check unassigned cells, if any
unassigned_cells <- Merged_Datasets.integrated$orig.ident[is.na(Merged_Datasets.integrated$group)]
if (length(unassigned_cells) > 0) {
  print("Unassigned cells found:")
  print(unique(unassigned_cells))
} else {
  print("All cells are assigned to a group.")
}

# Step 5: Set identities
Idents(Merged_Datasets.integrated) <- Merged_Datasets.integrated$group
print("Unique identities set:")
print(unique(Idents(Merged_Datasets.integrated)))


# Perform pathway enrichment analysis
go_enrichment_plot <- DEenrichRPlot(
  object = Merged_Datasets.integrated,  # Seurat object
  ident.1 = "BMMC",                     # Group 1: BMMC
  ident.2 = "CD34",                     # Group 2: CD34
  enrich.database = "GO_Biological_Process_2021",  # Enrichment database
  logfc.threshold = 0.25,               # Minimum log fold change
  p.val.cutoff = 0.05,                  # P-value cutoff
  max.genes = 500,                      # Maximum genes for enrichment
  num.pathway = 10,                     # Number of pathways to display
  test.use = "wilcox",                  # Statistical test
  balanced = TRUE                       # Show both enriched and depleted pathways
)

# Save the enrichment plot
ggsave(
  filename = "GO_Enrichment_BMMC_vs_CD34.png", 
  plot = go_enrichment_plot, 
  width = 10, 
  height = 6
)

# Display the plot
#print(go_enrichment_plot)


###WEEK 4###
# Select Subset for Trajectory Analysis (Solution for 9.1)

# Ensure identities are set to the correct manual annotations
#Idents(Merged_Datasets.integrated) <- factor(Merged_Datasets.integrated$seurat_clusters, 
                                             #levels = names(c(
                                               #"0" = "CD8_T", 
                                               #"1" = "Monocytes", 
                                               #"2" = "LMPP", 
                                               #"3" = "CD8_T", 
                                               #"4" = "NK", 
                                               #"5" = "Basophils", 
                                               #"6" = "B_Cells", 
                                               #"7" = "cDCs", 
                                               #"8" = "CD8_T", 
                                               #"9" = "Plasma", 
                                               #"10" = "pDC", 
                                               #"11" = "HSC", 
                                               #"12" = "NK", 
                                               #"13" = "HSC"
#                                             #)))

# Define clusters of interest for trajectory analysis
#selected_clusters <- c("11", "6", "9")
#subset_seurat <- subset(Merged_Datasets.integrated, idents = selected_clusters)
# Verify the subset
#table(Idents(subset_seurat))



# Convert Seurat object to Monocle 3 CDS format
#cds <- as.cell_data_set(subset_seurat)

# Preprocess the data using PCA
#cds <- preprocess_cds(cds, num_dim = 50)

# Reduce dimensions using UMAP
#cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA")

# Perform clustering on the UMAP reduction
#cds <- cluster_cells(cds, reduction_method = "UMAP")

# Learn the trajectory graph
#cds <- learn_graph(cds)

# Plot UMAP with clusters
#trajectory_plot <- plot_cells(
  #cds, 
  #color_cells_by = "cluster", 
  #show_trajectory_graph = TRUE, 
  #label_groups_by_cluster = TRUE
#)
#ggsave("Trajectory_UMAP_Selected_Clusters.png", plot = trajectory_plot, width = 8, height = 6, dpi = 300)

# Define root cells (likely HSCs) as the starting point
#root_cells <- colnames(cds)[colData(cds)$seurat_clusters == "HSC"]
#cds <- order_cells(cds, root_pr_nodes = root_cells)

# Plot pseudotime trajectory
#pseudotime_plot <- plot_cells(
  #cds, 
  #color_cells_by = "pseudotime",
  #label_branch_points = TRUE
#)
#ggsave("Pseudotime_Trajectory.png", plot = pseudotime_plot, width = 8, height = 6, dpi = 300)

#print("Trajectory analysis and pseudotime plots completed and saved.")

######GETTING ERROR WITH MONOCLE3 VERSION AFTER SPENDING A LOT OF TIME ON INSTALLING IT AND ITS DEPENDENCIES

##10 Cell cell communication
# Define the two datasets
BMMC_data <- subset(Merged_Datasets.integrated, subset = orig.ident == "BMMC")
CD34_data <- subset(Merged_Datasets.integrated, subset = orig.ident == "CD34")

# Convert Seurat objects to CellChat objects
cellchat_bmmc <- createCellChat(object = BMMC_data, group.by = "SingleR_label")
cellchat_cd34 <- createCellChat(object = CD34_data, group.by = "SingleR_label")

# Load the CellChat database
CellChatDB <- CellChatDB.human  # Use the human database
cellchat_bmmc@DB <- CellChatDB
cellchat_cd34@DB <- CellChatDB

# Preprocess data for communication analysis
cellchat_bmmc <- subsetData(cellchat_bmmc)  # Subset the expression data
cellchat_cd34 <- subsetData(cellchat_cd34)

# Compute communication probabilities
cellchat_bmmc <- identifyOverExpressedGenes(cellchat_bmmc)
cellchat_bmmc <- identifyOverExpressedInteractions(cellchat_bmmc)
cellchat_bmmc <- computeCommunProb(cellchat_bmmc)

cellchat_cd34 <- identifyOverExpressedGenes(cellchat_cd34)
cellchat_cd34 <- identifyOverExpressedInteractions(cellchat_cd34)
cellchat_cd34 <- computeCommunProb(cellchat_cd34)

# Filter out lowly expressed interactions
cellchat_bmmc <- filterCommunication(cellchat_bmmc, min.cells = 10)
cellchat_cd34 <- filterCommunication(cellchat_cd34, min.cells = 10)

# Identify communication pathways
cellchat_bmmc <- computeCommunProbPathway(cellchat_bmmc)
cellchat_cd34 <- computeCommunProbPathway(cellchat_cd34)

# Aggregate communication networks
cellchat_bmmc <- aggregateNet(cellchat_bmmc)
cellchat_cd34 <- aggregateNet(cellchat_cd34)

# Visualize overall communication strength
groupSize_bmmc <- as.numeric(table(cellchat_bmmc@meta$labels))
groupSize_cd34 <- as.numeric(table(cellchat_cd34@meta$labels))

netVisual_circle(cellchat_bmmc@net$count, vertex.weight = groupSize_bmmc, title.name = "Number of Interactions - BMMC")
netVisual_circle(cellchat_cd34@net$count, vertex.weight = groupSize_cd34, title.name = "Number of Interactions - CD34")

# Compare signaling pathways between BMMC and CD34
pathways_bmmc <- cellchat_bmmc@netP$pathways
pathways_cd34 <- cellchat_cd34@netP$pathways
common_pathways <- intersect(pathways_bmmc, pathways_cd34)

print("Common pathways between BMMC and CD34:")
print(common_pathways)

# Visualize one common pathway using circle plot
pathway_of_interest <- common_pathways[1]  # Pick the first pathway
netVisual_aggregate(cellchat_bmmc, signaling = pathway_of_interest, layout = "circle")
ggsave(paste0("BMMC_", pathway_of_interest, "_CirclePlot.png"))

netVisual_aggregate(cellchat_cd34, signaling = pathway_of_interest, layout = "circle")
ggsave(paste0("CD34_", pathway_of_interest, "_CirclePlot.png"))

# Save CellChat objects for future use
saveRDS(cellchat_bmmc, file = "CellChat_BMMC.rds")
saveRDS(cellchat_cd34, file = "CellChat_CD34.rds")

print("Analysis completed and plots saved.")
