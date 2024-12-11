#!/usr/bin/env R

# This scripts attempts to determine the CGC independent of the (astrocytoma) classification process


idat_grn <- 'tmp/201496850071_R02C01_Grn.idat'
idat_red <- 'tmp/201496850071_R02C01_Red.idat'
array_sentrix_id <- gsub("^.+/([^/]+)_(Grn|Red).idat$","\\1", idat_grn)



# 1. load M-values ----


tmp <- data.frame(
  array_channel_green = idat_grn,
  array_channel_red = idat_red
)  |> 
  dplyr::mutate(Basename = gsub("_(Grn|Red).idat$","", array_channel_green)) |> 
  dplyr::mutate(array_channel_green = NULL, array_channel_red = NULL) |> 
  dplyr::mutate(Sample_Name = paste0(array_sentrix_id)) |>
  dplyr::mutate(Array = gsub("^.+_","", array_sentrix_id)) |>
  dplyr::mutate(Slide = gsub("^([0-9]+)_.+$","\\1", array_sentrix_id))



RGSet <- minfi::read.metharray.exp(targets = tmp, force = T) #red/green channel together
proc <- minfi::preprocessNoob(RGSet, offset = 0, dyeCorr = T, verbose = TRUE, dyeMethod="single")  #dyeMethod="reference"


mvalue <- minfi::ratioConvert(proc, what = "M") |>
  assays() |>
  purrr::pluck('listData') |>
  purrr::pluck("M") |>
  data.table::as.data.table(keep.rownames = "probe_id")



# 2. load predictor ----


predictor <- readRDS("assets/LGC_predictor_probe_based_lm.Rds")




# 3. apply ----

# select target probes
data <- mvalue |> 
  tibble::column_to_rownames('probe_id') |> 
  t() |> 
  as.data.frame() |> 
  dplyr::select(rownames(predictor$beta)) |> 
  #t() |> 
  as.matrix()


# apply lm
out <- glmnet::predict.glmnet(predictor, data) |> 
  as.data.frame() |> 
  dplyr::rename(`CGC[Ac]` = 1)


out



