#' Plot heatmap of the train set or test set
#'
#' @param object a SubtypeClass object
#' @param set options could be 'test', 'train' or 'both'. Default 'test'.
#' @param ...  Parameters passed to \code{pheatmap}.
#' @return a pheatmap object
#' @export
#' @import pheatmap
#' @import tibble
#' @importFrom rlang .data
#' @importFrom dplyr mutate filter select relocate arrange bind_rows
#' @examples
#' library(mlsubtyping)
#' data <- get_median_centered(example_fpkm)
#' data <- assays(data)$centered
#' rownames(data) <- rowData(example_fpkm)$external_gene_name
#' object <- MLSubtype(data, disease = 'LUSC')
#' PlotHeat(object, set = 'both', fontsize = 10, show_rownames = FALSE, show_colnames = FALSE)

PlotHeat <- function(object, set = 'test', ...) {
  test_set <- t(object@test_set) %>% data.frame %>%
    rownames_to_column('Samples') %>%
    mutate(mRNA_subtype = object@subtypes,
           set = 'test_set') %>%
    relocate(.data$Samples, .data$set, .data$mRNA_subtype, object@features)
  train_set <- object@train_set %>% data.frame %>%
    rownames_to_column('Samples') %>%
    mutate(set = 'train_set') %>%
    relocate(.data$Samples, .data$set, .data$mRNA_subtype, object@features)

  joint_set <- bind_rows(test_set, train_set) %>%
    arrange(.data$set, .data$mRNA_subtype)

  rownames(joint_set) <- make.unique(joint_set$Samples)

  if (set == 'both') {
    plotdat <- joint_set
  }

  if (set == 'test') {
    plotdat <- dplyr::filter(joint_set,  set == 'test_set')
  }

  if (set == 'train') {
    plotdat <- dplyr::filter(joint_set,  set == 'train_set')
  }

  anno <- plotdat[, c('set', 'mRNA_subtype')]

  pheatmap(t(plotdat[, object@features]),
           cluster_cols = FALSE,
           annotation_col = anno, ...)
}
