#' select highly variable genes from a expression matrix
#'
#' @param data a (normalized) matrix with rownames of features and colnames of samples
#' @param top number of top highly variable genes to output
#' @return subset with top ranked genes by the variances
#' @export
#' @import SummarizedExperiment
#' @importFrom stats var
#' @examples
#' library(mlsubtyping)
#' data <- get_median_centered(example_fpkm)
#' data <- assays(data)$centered
#' get_hvg(data)
get_hvg <- function(data, top = 1000){
  varorder <- order(apply(data, 1, var, na.rm = TRUE), decreasing=T)[1:top]
  data[varorder,]
}
