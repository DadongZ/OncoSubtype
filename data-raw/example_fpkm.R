#' example FPKM object
#'
#' @return example FPKM object
library(dplyr)
library(tidyverse)
Ydir <- "Y:/Group_Projects/Astellas_Subtyping/Pavese_LUSCSub_02DEC2020/Data"
load(file.path(Ydir, "LUSC/TCGA_LUSC_FPKM.rda"))

xgenes <- intersect(lusc_centroids$Symbol, rowData(lusc_edat)$external_gene_name)
add_genes <- sample(rowData(lusc_edat)$external_gene_name, 1000, replace = FALSE)
sub_genes <- unique(c(xgenes, add_genes))

colData(lusc_edat)
example_fpkm <- subset(lusc_edat, external_gene_name%in%sub_genes, !is.na(paper_Expression.Subtype))
