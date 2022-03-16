#' @title Preprocessed LUSC TCGA data as trainning set for mRNA subtyping
#'
#' @description LUSC TCGA Data preprocessed as 1) \eqn{log2(FPKM+1)} 2) medium centered 3) 195 genes intersecting with lusc_centroids;
#' @format A data.frame with 195 features and 179 patients with LUSC subtypes.
"lusc_tcga"
