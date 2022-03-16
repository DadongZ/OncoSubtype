#' TCGA HNSC predictor
#'
#' @return centroids data set for HNSC
#' @references Nature volume 517, pages576–582(2015)
#'
library(dplyr)
data_dir  <- "Z:/resource/TCGA/HNSC"
hnsc_centroids <- read.delim(file.path(data_dir, "classification_centroid_728genes.txt"))
usethis::use_data(hnsc_centroids, compress = 'xz', overwrite = TRUE)
