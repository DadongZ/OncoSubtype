#' Get TCGA data ready for ML
library(tidyverse)
library(TCGAbiolinks)
library(dplyr)
library(SummarizedExperiment)
library(readxl)
#LUSC
Ydir <- "Y:/Group_Projects/Astellas_Subtyping/Pavese_LUSCSub_02DEC2020/Data"
load(file.path(Ydir, "TCGA/LUSC/TCGA_LUSC_FPKM.rda"))
lusc_norm <- log2(assay(lusc_edat)+1)
expr_median_centered <- sweep(lusc_norm, 1,  apply(lusc_norm, 1, median, na.rm=T))
assays(lusc_edat)$normalized <- lusc_norm
assays(lusc_edat)$median_centered_log2 <- expr_median_centered

lusc_tcga <- lusc_edat[, !is.na(colData(lusc_edat)$paper_Expression.Subtype)]
genelist <- intersect(rownames(lusc_centroids), rowData(lusc_tcga)$external_gene_name)
#select highly variable genes
enslist <- filter(rowData(lusc_tcga)%>%data.frame, external_gene_name%in%genelist)%>%pull(ensembl_gene_id)
lusc_tcga <- lusc_tcga[enslist, ]

edat <- assays(lusc_tcga)[['median_centered_log2']]
rownames(edat) <- rowData(lusc_tcga)$external_gene_name
lusc_tcga<-bind_cols(mRNA_subtype=as.character(colData(lusc_tcga)$paper_Expression.Subtype),
             t(edat)%>%data.frame%>%rownames_to_column('samples'))%>%
  column_to_rownames('samples')
use_data(lusc_tcga, overwrite = TRUE)

#lUAD
load(file.path(Ydir, "TCGA/LUAD/TCGA_LUAD_FPKM.rda"))
luad_norm <- log2(assay(luad_edat)+1)
expr_median_centered <- sweep(luad_norm, 1,  apply(luad_norm, 1, median, na.rm=T))
assays(luad_edat)$normalized <- luad_norm
assays(luad_edat)$median_centered_log2 <- expr_median_centered

luad_tcga <- luad_edat[, !is.na(colData(luad_edat)$paper_expression_subtype)]
genelist <- intersect(rownames(luad_centroids), rowData(luad_tcga)$external_gene_name)
enslist <- filter(rowData(luad_tcga)%>%data.frame, external_gene_name%in%genelist)%>%pull(ensembl_gene_id)
luad_tcga <- luad_tcga[enslist, ]

edat <- assays(luad_tcga)[['median_centered_log2']]
rownames(edat) <- rowData(luad_tcga)$external_gene_name
luad_tcga<-bind_cols(mRNA_subtype=as.character(colData(luad_tcga)$paper_expression_subtype),
                     t(edat)%>%data.frame%>%rownames_to_column('samples'))%>%
  mutate(mRNA_subtype = factor(recode(mRNA_subtype, `prox.-inflam` = 'proximal-inflammatory',
                               `prox.-prolif.` = 'proximal-proliferative',
                               `TRU` = 'terminal-respiratory-unit')))%>%
  column_to_rownames('samples')
use_data(luad_tcga, overwrite = TRUE)

## BLCA
## BLCA data RSEM normalized
load(file.path(Ydir, 'TCGA/BLCA/blca_rsem_normalized_exp.rda'))
blca_edat <- data
blca_norm <- log2(assay(blca_edat)+1)
expr_median_centered <- sweep(blca_norm, 1,  apply(blca_norm, 1, median, na.rm=T))
assays(blca_edat)$normalized <- blca_norm
assays(blca_edat)$median_centered_log2 <- expr_median_centered

blca_tcga <- blca_edat[, !is.na(colData(blca_edat)[["paper_mRNA cluster"]])]

#
genelist <- c("ALDH1L2","ARSI","FAP","GFPT2","MSN","NRP2","TNC","ABCC9","C19orf45","CASP1",
              "DOK6","HIC2","MYCL1","PI4KAP2","SULT2A1","CPXM2","DCN","DOK5","EVI2A","FGF7",
              "GYPC","ITGB2","KLHDC7A","MRVI1","MSX2","PDCD1LG2","PDZRN3","RNASE6","SLA",
              "SLC9A4","SMOC2","UPK1A","UPK2","UPK3A","IL15RA","KSR2","LOC100188947","PPFIBP2",
              "SLC37A2","SV2A")

blca_tcga <- blca_tcga[rownames(blca_tcga)%in%genelist, ]

edat <- assays(blca_tcga)[['median_centered_log2']]
blca_tcga<-bind_cols(mRNA_subtype=as.character(colData(blca_tcga)[["paper_mRNA cluster"]]),
                     t(edat)%>%data.frame%>%rownames_to_column('samples'))%>%
  column_to_rownames('samples')
use_data(blca_tcga, overwrite = TRUE)

## HNSC
load(file.path(Ydir, "TCGA/HNSCC/HNSCtcga_normalized_exp.rda"))
hnsc_edat <- data
hnsc_norm <- log2(assay(data)+1)
expr_median_centered <- sweep(hnsc_norm, 1,  apply(hnsc_norm, 1, median, na.rm=T))
assays(hnsc_edat)$normalized <- hnsc_norm
assays(hnsc_edat)$median_centered_log2 <- expr_median_centered

hnsc_tcga <- hnsc_edat[, !is.na(colData(hnsc_edat)$paper_RNA)]
genelist <- intersect(rownames(hnsc_centroids), rownames(hnsc_tcga))
hnsc_tcga <- hnsc_tcga[genelist, ]

edat <- assays(hnsc_tcga)[['median_centered_log2']]
hnsc_tcga<-bind_cols(mRNA_subtype=as.character(colData(hnsc_tcga)$paper_RNA),
                     t(edat)%>%
                       data.frame%>%rownames_to_column('samples'))%>%
                       column_to_rownames('samples')
use_data(hnsc_tcga, overwrite = TRUE)

## STAD
load(file.path(Ydir, "TCGA/STAD/TCGA_STAD_normalized_exp.rda"))
stad_edat <- data
stad_norm <- log2(assay(data)+1)
expr_median_centered <- sweep(stad_norm, 1,  apply(stad_norm, 1, median, na.rm=T))
assays(stad_edat)$normalized <- stad_norm
assays(stad_edat)$median_centered_log2 <- expr_median_centered

stad_tcga <- stad_edat[, !is.na(colData(stad_edat)$paper_Molecular.Subtype)]
genelist <- rownames(get_hvg(assays(stad_tcga)$median_centered_log2, top = 1400))
stad_tcga <- stad_tcga[genelist, ]
edat <- assays(stad_tcga)[['median_centered_log2']]
stad_tcga<-bind_cols(mRNA_subtype=as.character(colData(stad_tcga)$paper_Molecular.Subtype),
                     t(edat)%>%
                       data.frame%>%rownames_to_column('samples'))%>%
  column_to_rownames('samples')
use_data(stad_tcga, overwrite = TRUE)

## ESCA
load(file.path(Ydir, "TCGA/ESCA/TCGA_ESCA_normalized_exp.rda"))
esca_edat <- data
esca_norm <- log2(assay(data)+1)
expr_median_centered <- sweep(esca_norm, 1,  apply(esca_norm, 1, median, na.rm=T))
assays(esca_edat)$normalized <- esca_norm
assays(esca_edat)$median_centered_log2 <- expr_median_centered

esca_tcga <- esca_edat[, !is.na(colData(esca_edat)$paper_patient)]

#USE genes from paper extended data figure 2 for ESCC and EAC
genelist <- unique(read_excel(file.path(Ydir, 'TCGA/ESCA/ESCA_paper_pathway_genes.xlsx')) %>% pull(gene.symbol))
esca_tcga <- esca_tcga[rownames(esca_tcga) %in% genelist, ]
edat <- assays(esca_tcga)[['median_centered_log2']]
esca_tcga<-bind_cols(mRNA_subtype=as.character(colData(esca_tcga)$`paper_Histological Type - Oesophagus`),
                     t(edat)%>%
                     data.frame%>%rownames_to_column('samples'))%>%
                     column_to_rownames('samples') %>%
                     filter(!is.na(mRNA_subtype))
use_data(esca_tcga, overwrite = TRUE)

#ESCC
## get gene lists from LUSC and HNSC (according to the ESCA paper extended data figure 6)
lusc_genes <- colnames(lusc_tcga)[-1]
hnsc_genes <- colnames(hnsc_tcga)[-1]
escc_genes <- unique(c(lusc_genes, hnsc_genes))
escc_tcga <- esca_tcga[rownames(esca_tcga) %in% escc_genes, ]
edat <- assays(escc_tcga)[['median_centered_log2']]
escc_tcga<-bind_cols(mRNA_subtype=as.character(colData(escc_tcga)$`paper_ESCC subtype`),
                     t(edat)%>%
                       data.frame%>%rownames_to_column('samples'))%>%
  column_to_rownames('samples') %>%
  filter(!is.na(mRNA_subtype))
use_data(escc_tcga, overwrite = TRUE)
