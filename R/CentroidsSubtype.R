#' Predict the subtypes of selected cancer type based published papers
#'
#' @param data data set to predict the subtypes which is a numeric matrix with row names of features and column names of samples
#' @param disease character string of the disease to predict subtypes, currently support 'LUSC', 'LUAD'
#' @return an object of class "SubtypeClass" with four slots: genes used for predictiong, predicted subtypes of samples, a matrix of predicting scores, and the method.
#' @export
#' @import SummarizedExperiment
#' @import e1071
#' @importFrom methods new
#' @importFrom stats cor median predict
#' @examples
#' library(mlsubtyping)
#' data <- get_median_centered(example_fpkm)
#' data <- assays(data)$centered
#' rownames(data) <- rowData(example_fpkm)$external_gene_name
#' CentroidsSubtype(data, disease = 'HNSC')

CentroidsSubtype <- function(data, disease = 'LUSC') {
  centroids_data <- get(paste0(tolower(disease), '_centroids'))
  output <- new("SubtypeClass", features = intersect(rownames(data), rownames(lusc_centroids)))
  test_data <- data[output@features, ]
  predict_data <- centroids_data[output@features, ]
  cor_scores <- do.call(rbind.data.frame, lapply(colnames(test_data), function(x) {
    cor(test_data[,x], predict_data, use = 'complete') }))
  rownames(cor_scores) <- colnames(test_data)
  subtypes <- colnames(cor_scores)[apply(cor_scores, 1, which.max)]
  names(subtypes) <- colnames(test_data)
  output@subtypes <- as.character(subtypes)
  output@score <- as.matrix(cor_scores)
  output@train_set <- data.frame(predict_data)
  output@test_set <- data.frame(test_data)
  output@method <- 'Nearest centroids'
  return(output)
}
