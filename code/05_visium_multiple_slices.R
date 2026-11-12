# Module 5: working with multiple slices
# Tutorial section: Working with multiple slices in Seurat
source("code/_common.R")
brain <- load_checkpoint("01_brain_sct", "01_visium_qc_normalize.R")

# ---- 1. Load and normalize the second slice ----------------------------------
# "posterior1" is the back half of the same sagittal section.
brain2 <- LoadData("stxBrain", type = "posterior1")
brain2 <- UpdateSeuratObject(brain2)
brain2 <- SCTransform(brain2, assay = "Spatial", verbose = FALSE)

# ---- 2. Merge ----------------------------------------------------------------
# merge() keeps both images; each spot remembers its slice in orig.ident.
brain.merge <- merge(brain, brain2)
Images(brain.merge)
table(brain.merge$orig.ident)

# Each slice was normalized separately; use the union of their variable genes.
DefaultAssay(brain.merge) <- "SCT"
VariableFeatures(brain.merge) <- c(VariableFeatures(brain), VariableFeatures(brain2))

# ---- 3. Joint clustering -----------------------------------------------------
brain.merge <- RunPCA(brain.merge, verbose = FALSE)
brain.merge <- FindNeighbors(brain.merge, dims = 1:30)
brain.merge <- FindClusters(brain.merge, verbose = FALSE)
brain.merge <- RunUMAP(brain.merge, dims = 1:30)

# If the two slices separate in UMAP by orig.ident, there is a batch effect.
save_fig(DimPlot(brain.merge, reduction = "umap", group.by = c("ident", "orig.ident")),
         "05_umap_merged", width = 12)
save_fig(SpatialDimPlot(brain.merge), "05_both_slices", width = 12)
save_fig(SpatialFeaturePlot(brain.merge, features = c("Hpca", "Plp1")),
         "05_hpca_plp1_both", width = 12, height = 8)

save_checkpoint(brain.merge, "05_brain_merged")

# ---- Exercises ----------------------------------------------------------------
# 1. Do the clusters line up across the two slices at the cut (for example the
#    cortex continuing from anterior to posterior)?
# 2. Count spots per cluster per slice: table(Idents(brain.merge), brain.merge$orig.ident).
#    Which clusters exist in only one slice? Is that anatomy or batch?
# 3. Advanced: read ?IntegrateLayers and the Seurat integration vignette. When
#    would you integrate rather than merge?
