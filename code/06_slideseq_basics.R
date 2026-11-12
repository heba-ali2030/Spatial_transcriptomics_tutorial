# Module 6: Slide-seq QC, normalization and clustering
# Tutorial sections: Slide-seq > Dataset, Data preprocessing
source("code/_common.R")

# ---- 1. Load ---------------------------------------------------------------
# Slide-seq v2 of the mouse hippocampus. Beads are ~10 µm, so each is close to
# one cell, but each bead captures far fewer UMIs than a Visium spot.
slide.seq <- LoadData("ssHippo")
slide.seq <- UpdateSeuratObject(slide.seq)
slide.seq
Images(slide.seq)
summary(slide.seq$nCount_Spatial)  # compare with the Visium numbers from module 1

# ---- 2. QC -----------------------------------------------------------------
# Counts are so skewed that they are easier to read on a log scale.
p1 <- VlnPlot(slide.seq, features = "nCount_Spatial", pt.size = 0, log = TRUE) + NoLegend()
slide.seq$log_nCount_Spatial <- log(slide.seq$nCount_Spatial)
p2 <- SpatialFeaturePlot(slide.seq, features = "log_nCount_Spatial") +
  theme(legend.position = "right")
save_fig(wrap_plots(p1, p2), "06_slideseq_qc", width = 12)

# ---- 3. Normalize and cluster ---------------------------------------------
slide.seq <- SCTransform(slide.seq, assay = "Spatial", ncells = 3000, verbose = FALSE)
slide.seq <- RunPCA(slide.seq)
slide.seq <- RunUMAP(slide.seq, dims = 1:30)
slide.seq <- FindNeighbors(slide.seq, dims = 1:30)
slide.seq <- FindClusters(slide.seq, resolution = 0.3, verbose = FALSE)

# stroke = 0 removes bead outlines, which would otherwise hide the colors.
p1 <- DimPlot(slide.seq, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(slide.seq, stroke = 0)
save_fig(p1 + p2, "06_slideseq_clusters", width = 12)

save_fig(
  SpatialDimPlot(slide.seq,
                 cells.highlight = CellsByIdentities(object = slide.seq, idents = c(1, 6, 13)),
                 facet.highlight = TRUE),
  "06_slideseq_highlight", width = 12, height = 5
)

save_checkpoint(slide.seq, "06_slideseq_clustered")

# ---- Exercises ----------------------------------------------------------------
# 1. List three ways Slide-seq QC differs from Visium QC (UMIs per bead, number
#    of beads, image or no image). Write them as comments.
# 2. Should you filter low-count beads? Try subset(slide.seq, nCount_Spatial > 100)
#    and see how much of the tissue structure survives.
# 3. Find the dentate gyrus and CA1/CA3 in the cluster plot. Which clusters are they?
# 4. Run FindAllMarkers(slide.seq, only.pos = TRUE) and check the top genes of
#    the clusters you named.
