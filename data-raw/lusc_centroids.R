#' Wiklerson's LUSC predictor
#'
#' @return centroids data set for LUSC
library(dplyr)
library(tidyverse)
data_dir  <- "C:/Users/u1072932/OneDrive - IQVIA/Documents/Projects/Astellas_LUAD_subtyping_SOW/Data"

lusc_centroids <- read.csv(file.path(data_dir, "LUSC/predictor.centroids.csv"))%>%
  column_to_rownames(var = 'X')

usethis::use_data(lusc_centroids, compress = 'xz', overwrite = TRUE)
