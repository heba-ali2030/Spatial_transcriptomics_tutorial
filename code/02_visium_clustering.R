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


# check specific cluster for specific marker 
p <- SpatialDimPlot(brain, cells.highlight = CellsByIdentities(brain, idents = 0)) +
  SpatialFeaturePlot(brain, features = "Hpca") +
  SpatialFeaturePlot(brain, features = "Ppp1r1b")
save_fig(p, "02_check_cluster0", width = 15, height = 5)

SpatialDimPlot(brain, cells.highlight = CellsByIdentities(brain, idents = c(1, 5, 8, 11)),
               facet.highlight = TRUE, ncol = 4)


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

LinkedDimPlot(brain)
save_checkpoint(brain, "02_brain_clustered")

# ---- Exercises ----------------------------------------------------------------
# 1. Rerun FindClusters with resolution = 0.3, 0.8 and 1.5, saving each into
#    metadata (they are stored as SCT_snn_res.*). Plot them side by side with
#    SpatialDimPlot(brain, group.by = "SCT_snn_res.0.3") etc. Which looks most anatomical?

brain <- FindClusters(brain, resolution = 0.3, verbose = FALSE)
brain <- FindClusters(brain, resolution = 1.5, verbose = FALSE)

p1 <- SpatialDimPlot(brain, group.by = "SCT_snn_res.0.3", label = TRUE, label.size = 3)
p2 <- SpatialDimPlot(brain, group.by = "SCT_snn_res.0.8", label = TRUE, label.size = 3)
p3 <- SpatialDimPlot(brain, group.by = "SCT_snn_res.1.5", label = TRUE, label.size = 3)
p1 + p2 + p3
# 2. Using the Allen Brain Atlas (atlas.brain-map.org) as a guide, name the
#    anatomical region of five clusters (cortex layers, hippocampus, striatum,
#    ventricles, white matter...). Write the names as comments here.

# 3. Do clusters that are close in UMAP also sit close together in the tissue?
# not really 
# 4. Try dims = 1:10 instead of 1:30. What changes?
test <- FindNeighbors(brain, dims = 1:10)
test <- FindClusters(test, verbose = FALSE)

Idents(brain) <- "SCT_snn_res.0.8"
save_checkpoint(brain, "02_brain_clustered")
SpatialDimPlot(test, label = TRUE, label.size = 3) + SpatialDimPlot(brain, label = TRUE, label.size = 3)
