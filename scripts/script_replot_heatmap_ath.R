library('tidyverse')
library('ComplexHeatmap')
library('limma')
library('DESeq2')
library('RColorBrewer')

setwd('/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results_rmfull/')

load('degres_condi_Mock_ath.RData')

##~~~~~~~~~~~~~~~~~~~~~~useful funcs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
meanFe <- function(v) {

  require('magrittr')

  res <- v %>%
    split(rep(1 : 3, each = 4)) %>%
    sapply(mean, na.rm = TRUE)

  return(res)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~heatmap~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
kmeansRes <- read_csv('kmeans_10_ath.csv') %>%
  select(ID, cl)

## rlog transformed
rawC <- rldData %>%
  as.data.frame %>%
  .[, c(5:16)] %>%
  rownames_to_column('ID') %>%
  as_tibble %>%
  inner_join(kmeansRes %>% select(ID, cl))

scaleC <- rawC %>%
  select(contains('C_')) %>%
  t %>%
  scale %>%
  t %>%
  as_tibble %>%
  bind_cols(rawC %>% select(ID, cl))

cairo_pdf('kmeans_10_ath_heatmap2.pdf', height = 8)
syncom <- HeatmapAnnotation(SynCom = rep(c('AtSC', 'LjSC', 'Mock'), each = 4),
                            col = list(SynCom = c('Mock' = '#1b9e77', 'AtSC' = '#d95f02', 'LjSC' = '#7570b3')),
                            gp = gpar(col = 'black'))

Heatmap(matrix = scaleC %>% select(contains('C_')),
        name = 'Scaled Counts',
        ## row_order = order(scaleC$cl) %>% rev,
        row_split = scaleC$cl,
        row_gap = unit(2, "mm"),
        column_order = 1 : 12,
        column_split = rep(c('AtSC', 'LjSC', 'Mock'), each = 4),
        show_column_names = FALSE,
        col = colorRampPalette(rev(brewer.pal(n = 7, name = 'RdYlBu')))(100),
        top_annotation = c(syncom))
dev.off()
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~box plot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sampleN <- c('AtSC', 'LjSC', 'Mock')

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

  ggsave(paste0('boxplot_ath/kmeans_10_ath_boxplot', i, '.pdf'))
  ggsave(paste0('boxplot_ath/kmeans_10_ath_boxplot', i, '.jpeg'))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


