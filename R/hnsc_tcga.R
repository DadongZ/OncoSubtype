#' @title Preprocessed HNSC TCGA data as trainning set for mRNA subtyping
#'
#' @description HNSC TCGA Data preprocessed as 1) \eqn{log2(FPKM+1)} 2) medium centered 3) 694 genes intersecting with the centroids predictor;
#' @format A data.frame with 694 features and 279 patients with HNSC subtypes.
"hnsc_tcga"
