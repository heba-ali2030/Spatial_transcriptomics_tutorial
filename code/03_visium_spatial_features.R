# Module 3: spatially variable features
# Tutorial section: 10x Visium > Identification of Spatially Variable Features
source("code/_common.R")
brain <- load_checkpoint("02_brain_clustered", "02_visium_clustering.R")

# ---- 1. Approach A: differential expression between clusters ----------------
# Uses the clusters as prior knowledge. Clusters 5 and 7 in the tutorial are
# neighboring regions; your numbers may differ slightly, so check module 2's plots.
de_markers <- FindMarkers(brain, ident.1 = 5, ident.2 = 7)
head(de_markers)  # avg_log2FC > 0 means higher in cluster 5
save_fig(
  SpatialFeaturePlot(brain, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3),
  "03_de_markers", width = 12, height = 5
)

# ---- 2. Approach B: spatial statistics, no clusters needed -------------------
# Moran's I measures whether neighboring spots have similar expression.
# The default method, "markvariogram", is slower; the tutorial uses "moransi".
# Testing the top 1000 variable genes takes a few minutes (faster with Rfast2).
brain <- FindSpatiallyVariableFeatures(
  brain, assay = "SCT",
  features = VariableFeatures(brain)[1:1000],
  selection.method = "moransi"
)

# Seurat 5 and 6 both store the result in the SCT assay's feature metadata.
top.features <- head(SpatiallyVariableFeatures(brain, method = "moransi"), 6)
top.features
save_fig(
  SpatialFeaturePlot(brain, features = top.features, ncol = 3, alpha = c(0.1, 1)),
  "03_moransi_top6", width = 12, height = 8
)

save_checkpoint(brain, "03_brain_svf")

# ---- Exercises ----------------------------------------------------------------
# 1. In your own words (as a comment): what does a high Moran's I mean, and why
#    can a gene be spatially variable without being a marker of any one cluster?
# 2. Look at the full Moran's I table:
#      mi <- brain[["SCT"]][[]] ; head(mi[order(mi$moransi.spatially.variable.rank), ])
#    (column names can differ by version: use colnames(mi) to find them).
# 3. Find one gene in the Moran's I top 50 that is NOT in the top 50 of
#    FindAllMarkers(brain, only.pos = TRUE). Plot it. Where is it expressed?
# 4. Plot a gene with a Moran's I near 0. What does its pattern look like?
