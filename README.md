# Learning Seurat spatial analysis in R

A hands-on course built on the Seurat tutorial
[Analysis, visualization, and integration of spatial datasets](https://satijalab.org/seurat/articles/spatial_vignette).

**Start here:** open `Spatial_transcriptomics_tutorial.Rproj` in RStudio. That sets the
working directory to the project root, which every script expects. (From a terminal,
`cd` to this folder and run, for example, `Rscript code/01_visium_qc_normalize.R`.)

Each module is an R script in `code/`. Run it chunk by chunk in RStudio (Ctrl/Cmd + Enter),
read the comments, look at each plot, then do the **exercises** at the bottom.
Each module saves a checkpoint in `objects/`, so the next module can start without
recomputing everything.

The tutorial targets Seurat 6.0. These scripts were written to work on Seurat 5.x too

## Plan (about 2–3 weeks at 1 hour a day)

| Done | Module | Script | Time | Topic |
|---|---|---|---|---|
| [ ] | 0 | `code/00_setup.R` | 1 day | Install packages and download data |
| [ ] | 1 | `code/01_visium_qc_normalize.R` | 2 days | Seurat spatial object, QC, SCTransform |
| [ ] | 2 | `code/02_visium_clustering.R` | 2 days | PCA, clustering, UMAP, spatial plots |
| [ ] | 3 | `code/03_visium_spatial_features.R` | 2 days | Markers and spatially variable genes (Moran's I) |
| [ ] | 4 | `code/04_visium_subset_integrate.R` | 3 days | Cropping, and label transfer from scRNA-seq |
| [ ] | 5 | `code/05_visium_multiple_slices.R` | 1 day | Merging two slices |
| [ ] | 6 | `code/06_slideseq_basics.R` | 2 days | Slide-seq QC and clustering |
| [ ] | 7 | `code/07_slideseq_integration_rctd.R` | 3 days | Label transfer and RCTD deconvolution |
| [ ] | 8 | `code/08_capstone.R` | open | Analyze your own Visium dataset |

Modules 1–5 (Visium) and 6–7 (Slide-seq) are independent tracks. Module 0 is needed by both.

**Measured run times** (whole script with `Rscript`, 16 GB MacBook, Seurat 5.6):
01 ≈ 1 min · 02 ≈ 1 min · 03 ≈ 1.5 min · 04 ≈ 2 min · 05 ≈ 2 min · 06 ≈ 3.5 min ·
07 parts A–B ≈ 8 min, part C (RCTD, 2 cores) ≈ **2 h 40 min**. Start RCTD and do something else.

**Changes from the tutorial page** (each is explained in a comment in the script):
- Module 4: `FindTransferAnchors(..., normalization.method = "SCT", recompute.residuals = FALSE)`.
  The page's default call finds only 42 anchors on Seurat 5.6, and then `TransferData` fails.
- Module 7: RCTD keeps your earlier results, subsamples the reference to 1000 cells per
  type, and uses 2 cores. With the page's settings, R ran out of memory on 16 GB.

## Folder layout

```
Spatial_transcriptomics_tutorial.Rproj   open this in RStudio
code/      00_setup.R ... 08_capstone.R, plus _common.R (shared helpers)
data/      reference .rds files downloaded by 00_setup.R (not committed)
objects/   checkpoints saved by each module (not committed)
figures/   PNGs saved by each module (not committed)
```

## Key concepts before you start

- **Spot (Visium):** a 55 µm circle that captures RNA from about 1–10 cells. A Visium
  "cell" in Seurat is really a spot, so one column can be a mixture of cell types.
- **Bead (Slide-seq):** a 10 µm bead, close to single-cell size, but with far fewer UMIs per bead.
- **Seurat spatial object:** counts live in the `Spatial` assay. The tissue image,
  spot coordinates and scale factors live in `object@images` (for example `brain[["anterior1"]]`).
- **Scale factors:** convert spot coordinates from full-resolution pixels to the
  low-resolution image that is plotted.

## Glossary

| Term | Meaning |
|---|---|
| SCTransform (SCT) | Normalization with regularized negative binomial regression. It removes the effect of sequencing depth, which varies with tissue anatomy in spatial data. |
| Variable features | Genes with the most variation across spots. PCA uses them. |
| Anchors | Pairs of matching cells or spots between a reference and a query. Label transfer is built on them. |
| Label transfer | Predicting, for each spot, a score per reference cell type (`FindTransferAnchors` + `TransferData`). |
| Deconvolution | Estimating the proportions of cell types in each spot or bead (RCTD from `spacexr`). |
| Moran's I | A spatial autocorrelation statistic. High values mean a gene's expression is similar in neighboring spots. |
| Markvariogram | Another spatial statistic (the Seurat default). It is slower than Moran's I. |

## Common problems

| Problem | Fix |
|---|---|
| `readRDS("/brahms/shared/...")` fails | That path is on the authors' server. The scripts read from `data/` instead; run `00_setup.R` first. |
| Dropbox download gives a tiny HTML file | The link must end in `?dl=1`. Delete the file in `data/` and rerun `00_setup.R`. |
| `FindSpatiallyVariableFeatures` is very slow | Install `Rfast2` (`00_setup.R` tries to), or use fewer features. |
| RCTD section discards your Slide-seq results | The tutorial reloads `ssHippo` before RCTD. Module 7 keeps a separate object instead. |
| Out of memory | Close other R sessions, run `gc()`, and start each module in a fresh R session from its checkpoint. |
| An argument is "unused" | Seurat 5 and 6 differ slightly. Check `?FunctionName` and the comments in the script. |

## Self-check questions

1. Why does `nCount_Spatial` vary across a brain section, and why does that rule out plain log-normalization?
2. What is the difference between `SpatialFeaturePlot` and `SpatialDimPlot`?
3. What do Moran's I values near 0 and near 1 mean?
4. Why is label transfer onto Visium spots a *score* rather than one label?
5. What does `doublet_mode = "doublet"` assume about each Slide-seq bead?
6. When would you merge two slices, and when would you integrate them?
