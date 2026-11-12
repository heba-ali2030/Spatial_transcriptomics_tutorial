# Module 8: capstone — analyze a Visium dataset on your own
# Goal: repeat modules 1–3 without looking back, then extend them.
source("code/_common.R")

# ---- 1. Get a dataset -------------------------------------------------------
# Browse https://www.10xgenomics.com/datasets (filter: Spatial Gene Expression).
# Download from one sample: the "Feature / barcode matrix HDF5 (filtered)" file
# and the "Spatial imaging data" archive. Put both in data/my_visium/ and unpack
# the imaging archive there, so you have:
#   data/my_visium/filtered_feature_bc_matrix.h5
#   data/my_visium/spatial/   (tissue_lowres_image.png, scalefactors_json.json, ...)
data.dir <- "data/my_visium"

if (dir.exists(data.dir)) {
  obj <- Load10X_Spatial(data.dir = data.dir, filename = "filtered_feature_bc_matrix.h5")
  obj

  # ---- 2. Your turn: fill in each step --------------------------------------
  # a. QC: VlnPlot and SpatialFeaturePlot of nCount_Spatial and nFeature_Spatial.
  #    Add percent mitochondrial reads:
  #      obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^(MT|mt)-")
  # b. Decide on a filter, if any, and justify it in a comment.
  # c. SCTransform(assay = "Spatial").
  # d. PCA, neighbors, clusters, UMAP. Choose dims from ElbowPlot.
  # e. SpatialDimPlot of clusters; name the regions using the tissue's anatomy.
  # f. FindAllMarkers and FindSpatiallyVariableFeatures(selection.method = "moransi").
  # g. Optional: find a matching scRNA-seq reference and run label transfer (module 4).

} else {
  message("Put a Visium sample in ", data.dir, " first (see the instructions above).")
}

# ---- 3. Write up ------------------------------------------------------------
# Write half a page answering:
#   * What tissue regions did clustering find?
#   * Which genes are most spatially variable, and do they match known biology?
#   * What would you do differently with Visium HD or an imaging platform like Xenium?
#
# Next tutorials:
#   https://satijalab.org/seurat/articles/spatial_vignette_2          (imaging: Xenium, MERFISH...)
#   https://satijalab.org/seurat/articles/visiumhd_analysis_vignette  (Visium HD)
