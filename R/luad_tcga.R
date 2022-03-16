#' @title Preprocessed LUAD TCGA data as trainning set for mRNA subtyping
#'
#' @description LUAD TCGA Data preprocessed as 1) \eqn{log2(FPKM+1)} 2) medium centered 3) 471genes intersecting with luad_centroids;
#' @format A data.frame with 471 features and 245 patients with LUAD subtypes.
"luad_tcga"
