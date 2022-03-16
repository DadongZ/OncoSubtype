#' @title Preprocessed ESCC TCGA data as trainning set for mRNA subtyping
#'
#' @description ESCC TCGA Data preprocessed as 1) \eqn{log2(FPKM+1)} 2) medium centered of ESCA 3) 844 genes from LUSC and HNSC according to the ESCA paper (2017) Extended Data Figure 6;
#' @format A data.frame with 844 features and 90 patients with ESCC subtypes.
"escc_tcga"
