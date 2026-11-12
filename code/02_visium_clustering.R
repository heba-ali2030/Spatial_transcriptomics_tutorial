# Module 2: dimensionality reduction, clustering and visualization
# Tutorial section: 10x Visium > Dimensionality reduction, clustering, and visualization
source("code/_common.R")
brain <- load_checkpoint("01_brain_sct", "01_visium_qc_normalize.R")

# ---- 1. The standard single-cell workflow, run on spots ----------------------
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
ElbowPlot(brain, ndims = 50)  # how many PCs carry signal? The tutorial uses 30.
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)

# ---- 2. Clusters in UMAP space and in tissue space ---------------------------
p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
save_fig(p1 + p2, "02_umap_vs_tissue", width = 12)

# Highlight a few clusters one at a time, which is easier to read than 15 colors.
save_fig(
  SpatialDimPlot(brain,
                 cells.highlight = CellsByIdentities(object = brain, idents = c(2, 1, 4, 3, 5, 8)),
                 facet.highlight = TRUE, ncol = 3),
  "02_cluster_highlight", width = 12, height = 8
)

table(Idents(brain))

# ---- 3. Interactive plots (run these in RStudio, not with Rscript) -----------
if (interactive()) {
  SpatialDimPlot(brain, interactive = TRUE)                     # hover to see cluster per spot
  SpatialFeaturePlot(brain, features = "Ttr", interactive = TRUE) # change gene in the app
  LinkedDimPlot(brain)                                          # select in UMAP, see it on tissue
}

save_checkpoint(brain, "02_brain_clustered")

# ---- Exercises ----------------------------------------------------------------
# 1. Rerun FindClusters with resolution = 0.3, 0.8 and 1.5, saving each into
#    metadata (they are stored as SCT_snn_res.*). Plot them side by side with
#    SpatialDimPlot(brain, group.by = "SCT_snn_res.0.3") etc. Which looks most anatomical?
# 2. Using the Allen Brain Atlas (atlas.brain-map.org) as a guide, name the
#    anatomical region of five clusters (cortex layers, hippocampus, striatum,
#    ventricles, white matter...). Write the names as comments here.
# 3. Do clusters that are close in UMAP also sit close together in the tissue?
# 4. Try dims = 1:10 instead of 1:30. What changes?
