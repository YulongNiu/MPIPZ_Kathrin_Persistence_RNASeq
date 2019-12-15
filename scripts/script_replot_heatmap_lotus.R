library('tidyverse')
library('ComplexHeatmap')
library('limma')
library('DESeq2')
library('RColorBrewer')

setwd('/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results_rmfull/')

load('degres_condi_Mock_lotus.RData')

##~~~~~~~~~~~~~~~~~~~~~~useful funcs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
meanFe <- function(v) {

  require('magrittr')

  res <- v %>%
    split(rep(1 : 8, each = 6)) %>%
    sapply(mean, na.rm = TRUE)

  return(res)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~heatmap~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
kmeansRes <- read_csv('kmeans_10_lotus_likeath.csv') %>%
  select(ID, cl)

## rlog transformed
rawC <- rldData %>%
  as.data.frame %>%
  .[, c(-1:-4, -9:-12)] %>%
  rownames_to_column('ID') %>%
  as_tibble %>%
  inner_join(kmeansRes %>% select(ID, cl))

scaleC <- rawC %>%
  select(contains('L_')) %>%
  t %>%
  scale %>%
  t %>%
  as_tibble %>%
  bind_cols(rawC %>% select(ID, cl))

cairo_pdf('kmeans_10_lotus_likeath_heatmap2.pdf')
syncom <- HeatmapAnnotation(SynCom = rep(c('AtSC', 'LjSC', 'Mock+LjNodule218'), each = 4),
                            col = list(SynCom = c('Mock+LjNodule218' = '#1b9e77', 'AtSC' = '#d95f02', 'LjSC' = '#7570b3')),
                            gp = gpar(col = 'black'))

Heatmap(matrix = scaleC %>% select(contains('L_')),
        name = 'Scaled Counts',
        ## row_order = order(scaleC$cl) %>% rev,
        row_split = scaleC$cl,
        row_gap = unit(2, "mm"),
        column_order = 1 : 12,
        column_split = rep(c('AtSC', 'LjSC', 'Mock+LjNodule218'), each = 4),
        show_column_names = FALSE,
        col = colorRampPalette(rev(brewer.pal(n = 10, name = 'Spectral'))[c(-3, -4, -7, -8)])(10),
        top_annotation = c(syncom))
dev.off()
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
