# Shared helpers, sourced at the top of every module with source("code/_common.R").
# Scripts live in code/, but they run with the PROJECT ROOT as the working
# directory, so data/, objects/ and figures/ resolve next to code/.
# Easiest: open Spatial_transcriptomics_tutorial.Rproj in RStudio.
if (!dir.exists("code") || !file.exists("code/_common.R")) {
  stop("Working directory is ", getwd(), ".\n",
       "Set it to the project root (the folder that contains code/), e.g. setwd(\"..\").",
       call. = FALSE)
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratData)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

for (d in c("data", "objects", "figures")) dir.create(d, showWarnings = FALSE)

# Save a ggplot / patchwork figure to figures/ and also print it.
save_fig <- function(plot, name, width = 10, height = 6) {
  ggsave(file.path("figures", paste0(name, ".png")), plot,
         width = width, height = height, dpi = 120, bg = "white")
  print(plot)
  invisible(plot)
}

# Load a checkpoint from an earlier module, with a clear error if it is missing.
load_checkpoint <- function(name, made_by) {
  path <- file.path("objects", paste0(name, ".rds"))
  if (!file.exists(path)) {
    stop("Missing ", path, ". Run ", made_by, " first.", call. = FALSE)
  }
  readRDS(path)
}

save_checkpoint <- function(object, name) {
  saveRDS(object, file.path("objects", paste0(name, ".rds")))
  message("Saved objects/", name, ".rds")
}

# Read a reference .rds downloaded by 00_setup.R.
load_reference <- function(file) {
  path <- file.path("data", file)
  if (!file.exists(path)) stop("Missing ", path, ". Run 00_setup.R first.", call. = FALSE)
  readRDS(path)
}
