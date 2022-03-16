#' convert expression matrix to median-centered
#'
#' @param data a numeric matrix or `S4` object
#' @param log2 logical, if `TRUE` \eqn{log2(x+1)} will be applied before median centering. Defaut `TRUE`.
#' @return median-centered express matrix or an object with new slot "centered"
#' @export
#' @import SummarizedExperiment
#' @examples
#' get_median_centered(example_fpkm)
#'
#'
get_median_centered <- function(data, log2 = TRUE) {
  if(is.matrix(data)){
    if(log2){
      data <- log2(assay(data) + 1)
    }
    return (sweep(data, 1,  apply(data, 1, median, na.rm = TRUE)))
  }
  if(typeof(data) == 'S4'){
    edata <- assays(data)[[1]]
    if(log2){
      edata <- log2(edata + 1)
    }
    assays(data)$centered <- sweep(edata, 1,  apply(edata, 1, median, na.rm = TRUE))
    return(data)
  }
}
