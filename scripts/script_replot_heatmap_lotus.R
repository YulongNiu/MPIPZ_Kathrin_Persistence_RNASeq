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
    split(rep(1 : 4, each = 4)) %>%
    sapply(mean, na.rm = TRUE)

  return(res)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~heatmap all transcripts~~~~~~~~~~~~~~~~
kmeansRes <- read_csv('kmeans_10_lotus_collaborator.csv') %>%
  select(ID, cl)

## rlog transformed
rawC <- rldData %>%
  as.data.frame %>%
  .[, c(-1:-4)] %>%
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

cairo_pdf('kmeans_10_lotus_heatmap.pdf', height = 8)
syncom <- HeatmapAnnotation(SynCom = rep(c('AtSC', 'AtSC+LjNodule218', 'LjSC', 'Mock+LjNodule218'), each = 4),
                            col = list(SynCom = c('Mock+LjNodule218' = '#1b9e77', 'LjSC' = '#7570b3', 'AtSC' = '#d95f02', 'AtSC+LjNodule218' = '#e7298a')),
                            gp = gpar(col = 'black'))

Heatmap(matrix = scaleC %>% select(contains('L_')),
        name = 'Scaled Counts',
        row_order = order(scaleC$cl) %>% rev,
        row_split = scaleC$cl,
        row_gap = unit(2, "mm"),
        column_order = 1 : 16,
        column_split = rep(c('AtSC', 'AtSC+LjNodule218', 'LjSC', 'Mock+LjNodule218'), each = 4),
        show_column_names = FALSE,
        col = colorRampPalette(rev(brewer.pal(n = 7, name = 'RdYlBu')))(100),
        top_annotation = c(syncom))
dev.off()
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


##~~~~~~~~~~~~~~~~~~~~~~~~~~box plot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sampleN <- c('AtSC', 'AtSC+LjNodule218', 'LjSC', 'Mock+LjNodule218')

boxplotData <- rldData %>%
  .[, c(-1:-4)] %>%
  t %>%
  scale %>%
  t %>%
  apply(1, meanFe) %>%
  t %>%
  set_colnames(sampleN) %>%
  as.data.frame %>%
  rownames_to_column('ID') %>%
  as_tibble %>%
  inner_join(kmeansRes)

for (i in 1:10) {
  boxplotData %>%
    filter(cl == i) %>%
    select(-ID, -cl) %>%
    gather(key = 'Conditions', value = 'ScaleCounts') %>%
    mutate(Conditions = Conditions %>% factor(levels = sampleN)) %>%
    ggplot(aes(x = Conditions, y = ScaleCounts)) +
    geom_boxplot() +
    ylab('Scaled counts') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90),
          plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
          legend.text.align = 0,
          axis.text = element_text(size = 13),
          axis.title = element_text(size = 14),
          legend.text=element_text(size= 13),
          legend.title = element_text(size = 14))

  ggsave(paste0('boxplot_lotus_rmfull/kmeans_10_lotus_boxplot', i, '.pdf'))
  ggsave(paste0('boxplot_lotus_rmfull/kmeans_10_lotus_boxplot', i, '.jpeg'))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~plot DEGs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
wholeDEG <- read_csv('SynCom_vs_Mock_lotus_sva_k.csv')
kmeansRes <- read_csv('kmeans_10_lotus_collaborator.csv') %>%
  select(ID, cl)

fcsig <- wholeDEG %>%
  select(ends_with('FoldChange')) %>%
  transmute_all(list(~ case_when(. > log2(1.5) ~ 1,
                                 . < -log2(1.5) ~ 1,
                                 TRUE ~ 0))) %>%
  select(-matches('fullSC_'))

padjsig <- wholeDEG %>%
  select(ends_with('padj')) %>%
  abs %>%
  `<`(0.05) %>%
  as_tibble %>%
  transmute_all(list(~ if_else(is.na(.), FALSE, .))) %>%
  select(-matches('fullSC_'))

heatsig <- (padjsig * fcsig) %>%
  as_tibble %>%
  rowSums %>%
  {. >= 1} %>%
  which %>%
  dplyr::slice(wholeDEG, .) %>%
  inner_join(kmeansRes)

rawC <- rldData %>%
  as.data.frame %>%
  .[, c(-1:-4)] %>%
  rownames_to_column('ID') %>%
  as_tibble %>%
  inner_join(heatsig %>% select(ID, cl))

scaleC <- rawC %>%
  select(contains('L_')) %>%
  t %>%
  scale %>%
  t %>%
  as_tibble %>%
  bind_cols(rawC %>% select(ID, cl))

cairo_pdf('kmeans_10_lotus_heatmap_sig2.pdf')
syncom <- HeatmapAnnotation(SynCom = rep(c('AtSC', 'AtSC+LjNodule218', 'LjSC', 'Mock+LjNodule218'), each = 4),
                            col = list(SynCom = c('Mock+LjNodule218' = '#1b9e77', 'LjSC' = '#7570b3', 'AtSC' = '#d95f02', 'AtSC+LjNodule218' = '#e7298a')),
                            gp = gpar(col = 'black'))

Heatmap(matrix = scaleC %>% select(contains('L_')),
        name = 'Scaled Counts',
        ## row_order = order(scaleC$cl) %>% rev,
        row_split = scaleC$cl,
        row_gap = unit(2, "mm"),
        column_order = 1 : 16,
        column_split = rep(c('AtSC', 'AtSC+LjNodule218', 'LjSC', 'Mock+LjNodule218'), each = 4),
        show_column_names = FALSE,
        col = colorRampPalette(rev(brewer.pal(n = 10, name = 'Spectral'))[c(-3, -4, -7, -8)])(10),
        top_annotation = c(syncom))

dev.off()
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
