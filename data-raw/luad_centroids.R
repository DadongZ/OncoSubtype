#' Wiklerson's LUAD predictor
#'
#' @return centroids data set for LUAD
library(dplyr)
library(tidyverse)
data_dir  <- "C:/Users/u1072932/OneDrive - IQVIA/Documents/Projects/Astellas_LUAD_subtyping_SOW/Data"

luad_centroids <- read.csv(file.path(data_dir, "LUAD/wilkerson.2012.LAD.predictor.centroids.csv")) %>%
  column_to_rownames(var = 'X')

usethis::use_data(luad_centroids, compress = 'xz', overwrite = TRUE)
