#' @title Preprocessed BLCA TCGA data as trainning set for mRNA subtyping
#'
#' @description BLCA TCGA Data preprocessed as 1) \eqn{log2(FPKM+1)} 2) medium centered 3) top 4000 highly variable genes;
#' @format A data.frame with 4000 features and 408 patients with BLCA subtypes.
"blca_tcga"
