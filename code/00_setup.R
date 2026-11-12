# Module 0: setup
# Installs the packages the tutorial needs and downloads the datasets.
# Run once. Downloads are about 3 GB in total (the two references are 1 GB and 1.7 GB).

# Run from the project root (the folder that contains code/), so data/ is created there:
#   Rscript code/00_setup.R      or open Spatial_transcriptomics_tutorial.Rproj in RStudio
if (!file.exists("code/00_setup.R")) {
  stop("Set the working directory to the project root (the folder that contains code/).",
       call. = FALSE)
}

options(timeout = 3600)  # large downloads; the default 60 s is too short

# ---- 1. CRAN / Bioconductor packages -----------------------------------------
cran <- c("Seurat", "ggplot2", "patchwork", "dplyr", "remotes", "hdf5r")
missing <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing)

# glmGamPoi makes SCTransform much faster (Bioconductor)
if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("glmGamPoi")
}

# SeuratData provides the example datasets (GitHub)
if (!requireNamespace("SeuratData", quietly = TRUE)) {
  remotes::install_github("satijalab/seurat-data")
}

# spacexr provides RCTD, used in module 7 (GitHub)
if (!requireNamespace("spacexr", quietly = TRUE)) {
  remotes::install_github("dmcable/spacexr", build_vignettes = FALSE)
}

# Optional: Rfast2 speeds up Moran's I in FindSpatiallyVariableFeatures
if (!requireNamespace("Rfast2", quietly = TRUE)) {
  try(install.packages("Rfast2"))
}

# ---- 2. Example datasets from SeuratData -------------------------------------
library(SeuratData)
installed <- rownames(InstalledData())
if (!"stxBrain" %in% installed) InstallData("stxBrain")  # 10x Visium, mouse brain
if (!"ssHippo"  %in% installed) InstallData("ssHippo")   # Slide-seq v2, mouse hippocampus

# ---- 3. scRNA-seq references (Dropbox links from the tutorial) ----------------
# The tutorial reads these from /brahms/shared/..., which only exists on the
# authors' server. We download them to data/ instead. "?dl=1" forces a direct download.
dir.create("data", showWarnings = FALSE)
refs <- c(
  "allen_cortex.rds" =
    "https://www.dropbox.com/s/cuowvm4vrf65pvq/allen_cortex.rds?dl=1",
  "mouse_hippocampus_reference.rds" =
    "https://www.dropbox.com/s/cs6pii5my4p3ke3/mouse_hippocampus_reference.rds?dl=1"
)
for (f in names(refs)) {
  dest <- file.path("data", f)
  if (!file.exists(dest)) {
    message("Downloading ", f, " ...")
    download.file(refs[[f]], dest, mode = "wb")
  }
  # A failed Dropbox download is a small HTML page, not an .rds file.
  if (file.size(dest) < 1e6) {
    warning(dest, " is only ", file.size(dest), " bytes; delete it and rerun.")
  }
}

# ---- 4. Check ----------------------------------------------------------------
for (p in c("Seurat", "SeuratData", "spacexr", "glmGamPoi", "Rfast2")) {
  ok <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("%-10s %s\n", p, if (ok) as.character(packageVersion(p)) else "NOT INSTALLED"))
}
print(InstalledData()[, c("Dataset", "Version")])
print(file.info(file.path("data", names(refs)))[, "size", drop = FALSE])

# ---- Concepts to read before module 1 ----------------------------------------
# * 10x Visium technology overview: spots, capture areas, H&E image.
# * Slide-seq: barcoded beads at ~10 µm resolution.
# * ?SeuratObject::VisiumV2 — how the image and coordinates are stored.
