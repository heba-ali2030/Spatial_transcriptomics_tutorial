# Module 1: Visium data, QC and normalization
# Tutorial sections: 10x Visium > Dataset, Data preprocessing, Gene expression visualization
source("code/_common.R")

# ---- 1. Load the data ---------------------------------------------------------
# stxBrain = 10x Visium, sagittal mouse brain. "anterior1" is the front half.
# ?stxBrain 
brain <- LoadData("stxBrain", type = "anterior1")
brain <- UpdateSeuratObject(brain)
brain

# Look inside the object. Note the "Spatial" assay and the image "anterior1".
Assays(brain)
Images(brain)
dim(brain)                        # genes x spots
head(brain@meta.data)             # nCount_Spatial = UMIs per spot, nFeature_Spatial = genes per spot
head(GetTissueCoordinates(brain)) # spot positions on the image
brain[["anterior1"]]              # the image object (image, coordinates, scale factors)

# ---- 2. QC: sequencing depth per spot ----------------------------------------
# Depth varies a lot, and it follows the anatomy rather than being random noise.
p1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
p2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
save_fig(wrap_plots(p1, p2), "01_nCount_spatial")

summary(brain$nCount_Spatial)
summary(brain$nFeature_Spatial)

# ---- 3. Normalization with SCTransform ---------------------------------------
# Log-normalization forces every spot to the same total. Here, low counts in white
# matter reflect biology, so the tutorial uses SCTransform instead.
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
DefaultAssay(brain)  # now "SCT"
brain
# ---- 4. Visualize genes on the tissue ----------------------------------------
# Hpca = hippocampus marker; Ttr = choroid plexus marker.
save_fig(SpatialFeaturePlot(brain, features = c("Hpca", "Ttr")), "01_hpca_ttr")

# pt.size.factor changes the spot size; alpha = c(min, max) transparency
# lets low-expressing spots fade so the H&E image shows through.
p1 <- SpatialFeaturePlot(brain, features = "Ttr", pt.size.factor = 1)
p2 <- SpatialFeaturePlot(brain, features = "Ttr", alpha = c(0.1, 1))
save_fig(p1 + p2, "01_ttr_size_alpha")

save_checkpoint(brain, "01_brain_sct")

# ---- Exercises ----------------------------------------------------------------
# 1. Plot nFeature_Spatial on the tissue. Is it similar to nCount_Spatial? 
SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
SpatialFeaturePlot(brain, features = "nFeature_Spatial")

# check if they are correlates 
FeatureScatter(brain, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")

# 2. Make a log-normalized copy and compare it with SCT for one gene:
ln <- NormalizeData(brain, assay = "Spatial")

DefaultAssay(ln) <- "Spatial"
SpatialFeaturePlot(ln, "Hpca") + SpatialFeaturePlot(brain, "Hpca")

# 3. Pick two genes you know (for example Mbp for white matter, Gad1 for
#    interneurons) and plot them. Do they appear where you expect?
SpatialFeaturePlot(ln, "Mbp") + SpatialFeaturePlot(brain, "Mbp")
SpatialFeaturePlot(brain, features = c("Mbp", "Gad1"))

# 4. Try image.alpha = 0 in SpatialFeaturePlot. What does it change?
SpatialFeaturePlot(brain, features = "Hpca", image.alpha = 0)
