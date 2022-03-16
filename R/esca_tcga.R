#' @title Preprocessed ESCA TCGA data as trainning set for mRNA subtyping
#'
#' @description ESCA TCGA Data preprocessed as 1) \eqn{log2(FPKM+1)} 2) medium centered 3) Customized 202 genes based on the TCGA ESCA paper (2017);
#' @format A data.frame with 202 features and 148 patients with ESCA subtypes (ESCC and EAC).
"esca_tcga"
