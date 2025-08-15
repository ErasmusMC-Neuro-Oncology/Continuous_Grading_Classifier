#!/usr/bin/env R

# This scripts attempts to determine the CGC independent of the (astrocytoma) classification process
# The predictors are glmnet models exported to Rds format

# 0. config and dependencies

# install:
# - minfi
# - IlluminaHumanMethylationEPICmanifest
# - IlluminaHumanMethylationEPICanno.ilm10b4.hg19
# - IlluminaHumanMethylationEPICv2manifest
# - IlluminaHumanMethylationEPICv2anno.20a1.hg38
# - glmnet


# only thing needed to configure:

idat_grn <- 'tmp/201496850071_R02C01_Grn.idat'
idat_red <- 'tmp/201496850071_R02C01_Red.idat'
array_sentrix_id <- gsub("^.+/([^/]+)_(Grn|Red).idat$","\\1", idat_grn)


# then run code below :)


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


predictor_epicv2 <- readRDS("assets/CGCy_predictor_probe_based_lm_epicv2.Rds") # a.k.a. 950k
predictor_epicv1 <- readRDS("assets/CGCy_predictor_probe_based_lm.Rds") # a.k.a. 850k
predictor_450k   <- readRDS("assets/CGCy_predictor_probe_based_lm_450k.Rds")


if(annotation(RGSet)[1] == "IlluminaHumanMethylationEPICv2") {
  predictor <- predictor_epicv2
  suffix <- "_epicv2"
} else if (annotation(RGSet)[1] == "IlluminaHumanMethylation450k") {
  predictor <- predictor_450k
  suffix <- "_450k"
} else if (annotation(RGSet)[1] == "IlluminaHumanMethylationEPIC") {
  predictor <- predictor_epicv1
  suffix <- ""
}



# 3. apply to array ----


# acquire the exact same m-values
data <- mvalue |> 
  tibble::column_to_rownames('probe_id') |> 
  t() |> 
  as.data.frame() |> 
  dplyr::select(rownames(predictor$beta)) |> 
  as.matrix()


# apply lm to the data
out <- glmnet::predict.glmnet(predictor, data) |> 
  as.data.frame() |> 
  dplyr::rename(`CGCψ` = 1) |> 
  dplyr::rename_with(.fn = ~ paste0(., suffix), .cols = c('CGCψ'))



out




