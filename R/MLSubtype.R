#' Predict the subtypes of selected cancer type using machine learning
#'
#' @param data data set to predict the subtypes which is a numeric matrix with row names of features and column names of samples
#' @param disease character string of the disease to predict subtypes, currently support 'LUSC', 'LUAD', and 'BLCA'.
#' @param method character string of the method to use currently support 'rf'.
#' @param seed integer seed to use.
#' @param removeBatch whether do batch effect correction using \code{limma::removeBatchEffect}, default TRUE.
#' @return an object of class "SubtypeClass" with four slots: genes used for predictiong, predicted subtypes of samples, a matrix of predicting scores, and the method.
#' @export
#' @import SummarizedExperiment
#' @import caret
#' @import randomForest
#' @import e1071
#' @importFrom methods new
#' @importFrom limma removeBatchEffect
#' @importFrom stats cor median predict
#' @examples
#' library(mlsubtyping)
#' data <- get_median_centered(example_fpkm)
#' data <- assays(data)$centered
#' rownames(data) <- rowData(example_fpkm)$external_gene_name
#' MLSubtype(data, disease = 'LUAD', method = 'rf', seed = 123)

MLSubtype <- function(data, disease = 'LUSC', method = 'rf', removeBatch = TRUE, seed = NULL) {

  train_data <- get(paste0(tolower(disease), '_tcga'))

  output <- new("SubtypeClass", features = intersect(rownames(data), colnames(train_data)[-1]))

  if (method == 'rf') {
    test_data <- data[output@features, ]
    predict_data <- train_data[, c('mRNA_subtype', output@features)]

    if (removeBatch) {
      batch <- c(rep('test', ncol(test_data)), rep('train', nrow(predict_data)))
      mat_batch <- cbind(test_data, t(predict_data[,-which(colnames(predict_data)=='mRNA_subtype')]))
      mat_rmbatch <- removeBatchEffect(mat_batch, batch)
      test_data_new <- mat_rmbatch[, batch == 'test']
      predict_data_new <- t(mat_rmbatch[, batch == 'train'])
      if (identical(rownames(predict_data_new), rownames(predict_data))) {
          predict_data_new <- cbind.data.frame(mRNA_subtype = predict_data[,'mRNA_subtype'], predict_data_new)
      }
      predict_data <- predict_data_new
      test_data <- test_data_new
    }

    out_pred <- get_rf_pred(train_set = predict_data, test_set = test_data, seed = seed)
    output@subtypes <- colnames(out_pred$results)[apply(out_pred$results, 1, which.max)]
    output@score <- out_pred$results
    output@train_set <- data.frame(predict_data)
    output@test_set <- data.frame(test_data)
    output@method <- out_pred$model
  }
  return(output)
}
