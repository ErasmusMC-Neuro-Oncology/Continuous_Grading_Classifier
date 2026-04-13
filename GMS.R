#!/usr/bin/env R

library(minfi)
library(readxl)
library(matrixStats)   # colMedians — already a minfi dependency

# ------------------------------------------------------------------
# Helper: download the 1,389 GMS probe IDs from the supplement
# ------------------------------------------------------------------
get_gms_probes <- function() {
  url <- paste0(
    "https://static-content.springer.com/esm/",
    "art%3A10.1038%2Fs43018-025-01023-z/MediaObjects/",
    "43018_2025_1023_MOESM2_ESM.xlsx"
  )
  tmp <- tempfile(fileext = ".xlsx")
  download.file(url, tmp, mode = "wb", quiet = TRUE)
  df  <- read_excel(tmp, sheet = 2, skip = 1)
  df$Probe_ID[df$GMS == TRUE]
}

# ------------------------------------------------------------------
# calculate_gms()
#
# Option A — Sentrix IDs + folder containing the IDAT files:
#   sentrix_ids : character vector, e.g. c("201496850071_R02C01")
#   idat_dir    : directory that holds <sentrix_id>_Grn.idat / _Red.idat
#
# Option B — explicit IDAT file paths:
#   idat_grn / idat_red : character vectors of equal length
#
# Returns a named numeric vector: one median M-value per sample.
# ------------------------------------------------------------------
calculate_gms <- function(sentrix_ids = NULL,
                           idat_dir    = ".",
                           idat_grn    = NULL,
                           idat_red    = NULL) {

  if (!is.null(sentrix_ids)) {
    basenames <- file.path(idat_dir, sentrix_ids)
  } else {
    stopifnot(!is.null(idat_grn), !is.null(idat_red),
              length(idat_grn) == length(idat_red))
    basenames <- sub("_Grn\\.idat$", "", idat_grn)
  }

  rg   <- read.metharray(basenames, extended = TRUE, force = TRUE)
  noob <- preprocessNoob(rg, dyeMethod = "single")
  m    <- getM(noob)

  probes <- intersect(get_gms_probes(), rownames(m))
  message(sprintf("GMS: using %d / 1389 probes present in data", length(probes)))

  scores <- colMedians(m[probes, , drop = FALSE], na.rm = TRUE)
  setNames(scores, colnames(m))
}

# ------------------------------------------------------------------
# Usage
# ------------------------------------------------------------------

# Option A: Sentrix IDs
gms <- calculate_gms(
  sentrix_ids = c("201496850071_R02C01",
                  "201496850071_R03C01"),
  idat_dir    = "path/to/idat/folder"
)

# Option B: explicit IDAT paths
gms <- calculate_gms(
  idat_grn = c("tmp/201496850071_R02C01_Grn.idat"),
  idat_red = c("tmp/201496850071_R02C01_Red.idat")
)

print(gms)

