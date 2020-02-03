#################################join At and Lj########################
library('tidyverse')
library('ComplexHeatmap')
library('limma')
library('DESeq2')
library('RColorBrewer')

setwd('/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results_orthologs/')

load('degres_condi_RBH_rmfull_ath.RData')
rldDataAth <- rldData %>%
  as.data.frame %>%
  rownames_to_column('ID')

load('degres_condi_RBH_rmfull_lotus.RData')
rldDataLotus <- rldData %>%
  as.data.frame %>%
  rownames_to_column('ID')

rldData <- inner_join(rldDataAth, rldDataLotus) %>%
  column_to_rownames('ID')

kmeansRes <- read_csv('kmeans_12_RBH_rmfull.csv')

scaleCountAt <- rldData[, 1:12, drop = FALSE] %>%
  t %>%
  scale %>%
  t

scaleCountLj <- rldData[, 13:24, drop = FALSE] %>%
  t %>%
  scale %>%
  t

scaleC <- cbind(scaleCountAt, scaleCountLj) %>%
  as.data.frame %>%
  rownames_to_column('ID') %>%
  as_tibble %>%
  bind_cols(kmeansRes %>% select(ID, cl))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plot common RBH~~~~~~~~~~~~~~~~~~~~~~~
ht_list <- Heatmap(matrix = scaleC %>% select(contains('_')),
        name = 'Scaled Counts',
        row_order = order(scaleC$cl) %>% rev,
        row_split = scaleC$cl,
        row_gap = unit(2, "mm"),
        column_order = 1 : 24,
        column_split = rep(c('At', 'Lj'), each = 12),
        show_column_names = TRUE,
        col = colorRampPalette(rev(brewer.pal(n = 10, name = 'Spectral')))(100),
        use_raster = FALSE)

pdf('kmeans_12_RBH_rmfull_heatmap.pdf')
draw(ht_list)
dev.off()
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plot ath DEGs~~~~~~~~~~~~~~~~~~~~~~~
## AtSC vs. LjSC
degs <- read_csv('SynCom_vs_Mock_RBH_rmfull_ath_sva_k.csv') %>%
  filter(AtSC_At_vs_LjSC_At_padj < 0.1, abs(AtSC_At_vs_LjSC_At_log2FoldChange) > log2(1.5)) %>%
  select(ID)

scaleCdegs <- scaleC %>%
  filter(ID %in% degs$ID)

ht_list <- Heatmap(matrix = scaleCdegs %>%
                     select(contains('_')),
                   name = 'Scaled Counts',
                   ## row_order = order(scaleCdegs$cl) %>% rev,
                   row_split = scaleCdegs$cl,
                   row_gap = unit(2, "mm"),
                   column_order = 1 : 24,
                   column_split = rep(c('At', 'Lj'), each = 12),
                   show_column_names = TRUE,
                   col = colorRampPalette(rev(brewer.pal(n = 10, name = 'Spectral')))(100),
                   use_raster = FALSE)

pdf('kmeans_12_RBH_rmfull_heatmap_athdegs.pdf')
draw(ht_list)
dev.off()
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######################################################################
