# Module 7: Slide-seq label transfer, spatial genes, and RCTD deconvolution
# Tutorial sections: Slide-seq > Integration with a scRNA-seq reference;
#                    Identification of spatially variable features; RCTD
source("code/_common.R")
slide.seq <- load_checkpoint("06_slideseq_clustered", "06_slideseq_basics.R")

# Mouse hippocampus scRNA-seq reference (Saunders et al. 2018). From 00_setup.R.
ref <- load_reference("mouse_hippocampus_reference.rds")
ref <- UpdateSeuratObject(ref)
table(ref$celltype)

# ==== Part A: label transfer =================================================
anchors <- FindTransferAnchors(reference = ref, query = slide.seq,
                               normalization.method = "SCT", npcs = 50)
predictions.assay <- TransferData(anchorset = anchors, refdata = ref$celltype,
                                  prediction.assay = TRUE,
                                  weight.reduction = slide.seq[["pca"]], dims = 1:50)
slide.seq[["predictions"]] <- predictions.assay

DefaultAssay(slide.seq) <- "predictions"
save_fig(SpatialFeaturePlot(slide.seq,
                            features = c("Dentate Principal cells", "CA3 Principal cells",
                                         "Entorhinal cortex", "Endothelial tip",
                                         "Ependymal", "Oligodendrocyte"),
                            alpha = c(0.1, 1)),
         "07_prediction_scores", width = 12, height = 8)

# Beads are near single-cell, so taking the top-scoring type as a label is reasonable.
slide.seq$predicted.id <- GetTransferPredictions(slide.seq)
Idents(slide.seq) <- "predicted.id"
save_fig(SpatialDimPlot(slide.seq,
                        cells.highlight = CellsByIdentities(object = slide.seq,
                          idents = c("CA3 Principal cells", "Dentate Principal cells",
                                     "Endothelial tip")),
                        facet.highlight = TRUE),
         "07_predicted_id", width = 12, height = 5)

# ==== Part B: spatially variable genes =======================================
# Slide-seq beads are not on a grid, so Moran's I first bins them into a
# 100 x 100 grid (x.cuts, y.cuts). This step takes several minutes.
DefaultAssay(slide.seq) <- "SCT"
slide.seq <- FindSpatiallyVariableFeatures(slide.seq, assay = "SCT", layer = "scale.data",
                                           features = VariableFeatures(slide.seq)[1:1000],
                                           selection.method = "moransi",
                                           x.cuts = 100, y.cuts = 100)
# max.cutoff = "q95" caps the color scale so a few extreme beads don't wash it out.
save_fig(SpatialFeaturePlot(slide.seq,
                            features = head(SpatiallyVariableFeatures(slide.seq, method = "moransi"), 6),
                            ncol = 3, alpha = c(0.1, 1), max.cutoff = "q95"),
         "07_slideseq_moransi", width = 12, height = 8)

save_checkpoint(slide.seq, "07_slideseq_predictions")

# ==== Part C: RCTD deconvolution (spacexr) ===================================
# RCTD models each bead as a mix of reference cell types. doublet_mode = "doublet"
# assumes at most two cell types per bead. The tutorial reloads ssHippo here,
# which throws away Parts A and B; we build the RCTD input from a fresh copy
# instead and keep slide.seq.
#
# Memory: RCTD needs a lot of RAM. On a 16 GB laptop, free everything from
# Parts A and B first (slide.seq was saved above and is reloaded at the end),
# subsample the reference, and use few cores (each core copies the data).
# You can start a fresh R session at this point: run source("code/_common.R") and the
# "ref <- ..." lines at the top, then continue from here.
library(spacexr)
rm(list = setdiff(ls(), c("ref", lsf.str())))
gc()

# Keep at most 1000 cells per cell type (Dentate Principal cells has 10,000).
set.seed(1)
keep <- unlist(lapply(split(colnames(ref), ref$celltype),
                      function(cells) sample(cells, min(length(cells), 1000))))
ref <- subset(ref, cells = keep)
gc()

raw <- UpdateSeuratObject(LoadData("ssHippo"))

# Query: bead coordinates and raw counts
counts <- raw[["Spatial"]]$counts
coords <- GetTissueCoordinates(raw)
if ("cell" %in% colnames(coords)) rownames(coords) <- coords$cell
coords <- coords[, c("x", "y")]
query <- SpatialRNA(coords, counts, colSums(counts))

# Reference: counts, cell-type labels and total UMIs per cell
ref.counts <- ref[["RNA"]]$counts
cluster <- as.factor(ref$celltype)
names(cluster) <- colnames(ref)
nUMI <- ref$nCount_RNA
names(nUMI) <- colnames(ref)
reference <- Reference(ref.counts, cluster, nUMI)
rm(raw, ref, ref.counts, counts); gc()

# The tutorial uses max_cores = 8. More cores is faster but uses more memory.
RCTD <- create.RCTD(query, reference, max_cores = 2)
rm(query, reference); gc()
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")   # ~2 h 40 min on 2 cores (41,000 beads)
saveRDS(RCTD, file.path("objects", "07_rctd.rds"))

slide.seq <- load_checkpoint("07_slideseq_predictions", "Parts A and B above")
slide.seq <- AddMetaData(slide.seq, metadata = RCTD@results$results_df)
p1 <- SpatialDimPlot(slide.seq, group.by = "first_type", stroke = 0)
p2 <- SpatialDimPlot(slide.seq, group.by = "second_type", stroke = 0)
save_fig(p1 | p2, "07_rctd_types", width = 14)

save_checkpoint(slide.seq, "07_slideseq_rctd")

# ---- Exercises ----------------------------------------------------------------
# 1. Compare label transfer with RCTD:
#      table(slide.seq$predicted.id, slide.seq$first_type)
#    Where do they agree, and where not?
# 2. table(slide.seq$spot_class): what share of beads are "singlet",
#    "doublet_certain" and "reject"? What does this say about bead resolution?
# 3. Try doublet_mode = "full" on a subset of beads. How do the weights differ?
# 4. Pick one gene from the Moran's I top list and find which cell type expresses it.
