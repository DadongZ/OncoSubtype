#' Predict the subtypes of selected cancer type
#'
#' @param train_set training set with rownames of samples, first column named 'mRNA_subtype' and the rest of features and expression values.
#' @param test_set test set with rownames of features and colnames of samples.
#' @param method character string of the method to use currently support 'rf'.
#' @param seed integer seed to use.
#' @return a matrix with column names of subtypes and predicted probabilities.
#' @import caret
#' @import randomForest

get_rf_pred <- function(train_set, test_set, method = 'rf', seed = NULL){

  if (!is.null(seed)) set.seed(seed)

  trControl = trainControl(method = 'repeatedcv',
                           number = 4,
                           repeats = 5,
                           savePredictions = 'final')

  rf_model <- train(mRNA_subtype  ~ . ,
                    data = train_set,
                    trControl = trControl,
                    method = method,
                    tuneLength = 5)

  out_pred <- predict(rf_model, newdata = t(test_set), type = 'prob')

  return(list(model = rf_model, results = out_pred))
}
