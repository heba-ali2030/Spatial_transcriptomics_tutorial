# Module 4: subsetting anatomical regions, and integration with scRNA-seq
# Tutorial sections: Subset out anatomical regions; Integration with single-cell data
source("code/_common.R")
brain <- load_checkpoint("02_brain_clustered", "02_visium_clustering.R")

# ---- 1. Subset by cluster ------------------------------------------------------
# The tutorial keeps clusters 1, 2, 3, 4, 6, 7 as "cortex". Cluster numbers can
# change between Seurat versions, so CHECK figures/02_cluster_highlight.png and
# edit this vector to the clusters that lie in the cortex on your plot.
# (In a test run with Seurat 5.6 the dorsal cortex arc was clusters 2, 3, 4, 6;
#  cluster 1 was elsewhere. Deciding this yourself is part of the exercise.)
cortex.clusters <- c(1, 2, 3, 4, 6, 7)
brain_sub <- subset(brain, idents = cortex.clusters)
save_fig(SpatialDimPlot(brain_sub, label = TRUE, label.size = 3), "04_cortex_clusters")

# ---- 2. Crop the image to a rectangle ----------------------------------------
# Crop() only changes what is drawn; the spots stay in the object.
brain_crop <- brain_sub
brain_crop[["anterior1"]] <- Crop(brain_crop[["anterior1"]], x = c(150, 500), y = c(0, 425))
p1 <- SpatialDimPlot(brain_crop, crop = TRUE, label = TRUE)
p2 <- SpatialDimPlot(brain_crop, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
save_fig(p1 + p2, "04_crop", width = 12)

# ---- 3. Keep only spots above a line (the dorsal cortex) ---------------------
# Coordinates come from the low-resolution image. Plot y as -y so the picture
# matches the tissue orientation (image rows grow downward).
coords <- GetTissueCoordinates(brain_sub, scale = "lowres")
head(coords)
m <- 0.6   # slope of the cut line
b <- -550  # intercept
ggplot(coords, aes(x, -y)) + geom_point(size = 0.4) +
  geom_abline(slope = m, intercept = b, colour = "red") + coord_fixed()

cortex_coords <- subset(coords, -y >= (m * x + b))
cortex <- subset(brain_sub, cells = cortex_coords$cell)
# If the cut looks wrong on your data, change m and b and re-plot.
p1 <- SpatialDimPlot(cortex, crop = TRUE, label = TRUE)
p2 <- SpatialDimPlot(cortex, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
save_fig(p1 + p2, "04_cortex_line_cut", width = 12)

# ---- 4. Prepare the scRNA-seq reference --------------------------------------
# Allen Institute mouse cortex, ~14,000 cells, SMART-Seq2. Downloaded by 00_setup.R.
allen_reference <- load_reference("allen_cortex.rds")
allen_reference <- UpdateSeuratObject(allen_reference)
table(allen_reference$subclass)
# ncells = 3000 speeds up SCTransform by fitting the model on a subsample.
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)
DimPlot(allen_reference, group.by = "subclass", label = TRUE) + NoLegend()

# Re-normalize the cortex subset on its own before mapping.
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE) %>% RunPCA(verbose = FALSE)

# ---- 5. Label transfer -------------------------------------------------------
# For each spot, TransferData gives a score per reference cell type. A spot can
# score for several types at once, since it mixes several cells.
# recompute.residuals = FALSE keeps the query's own SCT model. Without it, Seurat
# re-normalizes the Visium spots with the SMART-Seq2 reference's model; the two
# technologies differ so much that only ~42 anchors are found (vs ~1150 here),
# and TransferData then fails with "Please set k.weight to be smaller...".
# (The tutorial page calls FindTransferAnchors(reference, query) without these
# arguments; on Seurat 5.6 that also finds only ~42 anchors.)
anchors <- FindTransferAnchors(reference = allen_reference, query = cortex,
                               normalization.method = "SCT",
                               recompute.residuals = FALSE)
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$subclass,
                                  prediction.assay = TRUE,
                                  weight.reduction = cortex[["pca"]], dims = 1:30)
cortex[["predictions"]] <- predictions.assay
rownames(cortex[["predictions"]])  # one "feature" per cell type, plus "max"

DefaultAssay(cortex) <- "predictions"
save_fig(SpatialFeaturePlot(cortex, features = c("L2/3 IT", "L4"),
                            pt.size.factor = 1.6, ncol = 2, crop = TRUE),
         "04_layer_predictions")

# ---- 6. Which cell types are spatially organized? ----------------------------
# The same Moran's I statistic, applied to prediction scores instead of genes.
cortex <- FindSpatiallyVariableFeatures(cortex, assay = "predictions", image = "anterior1",
                                        selection.method = "moransi",
                                        features = rownames(cortex), r.metric = 5,
                                        layer = "data")
top.clusters <- head(SpatiallyVariableFeatures(cortex, method = "moransi"), 4)
save_fig(SpatialPlot(object = cortex, features = top.clusters, ncol = 2),
         "04_top_spatial_celltypes", height = 8)
save_fig(SpatialFeaturePlot(cortex, features = c("Astro", "L2/3 IT", "L4", "L5 PT", "L5 IT",
                                                 "L6 CT", "L6 IT", "L6b", "Oligo"),
                            pt.size.factor = 1, crop = FALSE, alpha = c(0.1, 1)),
         "04_all_layers", width = 12, height = 10)

save_checkpoint(cortex, "04_cortex_predictions")

# ---- Exercises ----------------------------------------------------------------
# 1. Do the layer predictions appear in the right order from the brain surface
#    inward (L2/3 -> L4 -> L5 -> L6 -> L6b)? Describe what you see.
# 2. Plot "max" from the predictions assay. Where is mapping confidence lowest, and why?
# 3. Change m and b to cut out a different region. Can you keep only the deep layers?
# 4. Why might label transfer and deconvolution (module 7) give different answers?
