
## originally by Yulong Niu
## yulong.niu@hotmail.com

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

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~heatmap all transcripts~~~~~~~~~~~~~~~~~
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

cairo_pdf('hierarch_6_ath_heatmap2.pdf', height = 8)
syncom <- HeatmapAnnotation(SynCom = rep(c('AtSC', 'LjSC', 'Mock'), each = 4),
                            col = list(SynCom = c('Mock' = '#1b9e77', 'AtSC' = '#d95f02', 'LjSC' = '#7570b3')),
                            gp = gpar(col = 'black'))

Heatmap(matrix = scaleC %>% select(contains('C_')),
        name = 'Scaled Counts',
        row_km = 6,
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

##~~~~~~~~~~~~~~~~~~~~~~~~~~plot DEGs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
wholeDEG <- read_csv('SynCom_vs_Mock_ath_sva_k.csv')
kmeansRes <- read_csv('kmeans_10_ath.csv') %>%
  select(ID, cl)

fcsig <- wholeDEG %>%
  select(ends_with('FoldChange')) %>%
  transmute_all(list(~ case_when(. > log2(1.5) ~ 1,
                                 . < -log2(1.5) ~ 1,
                                 TRUE ~ 0))) %>%
  select(-contains('fullSC'))

padjsig <- wholeDEG %>%
  select(ends_with('padj')) %>%
  abs %>%
  `<`(0.05) %>%
  as_tibble %>%
  transmute_all(list(~ if_else(is.na(.), FALSE, .))) %>%
  select(-contains('fullSC'))

heatsig <- (padjsig * fcsig) %>%
  as_tibble %>%
  rowSums %>%
  {. >= 1} %>%
  which %>%
  dplyr::slice(wholeDEG, .) %>%
  inner_join(kmeansRes)

rawC <- rldData %>%
  as.data.frame %>%
  .[, c(5:16)] %>%
  rownames_to_column('ID') %>%
  as_tibble %>%
  inner_join(heatsig %>% select(ID, cl))

scaleC <- rawC %>%
  select(contains('C_')) %>%
  t %>%
  scale %>%
  t %>%
  as_tibble %>%
  bind_cols(rawC %>% select(ID, cl))

cairo_pdf('kmeans_10_ath_heatmap_sig2.pdf')
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
        col = colorRampPalette(rev(brewer.pal(n = 10, name = 'Spectral'))[c(-3, -4, -7, -8)])(10),
        top_annotation = c(syncom))

dev.off()

## GO analysis
library('clusterProfiler')
library('org.At.tair.db')

kall <- lapply(kmeansRes$cl %>% unique, function(x) {

  eachG <- kmeansRes %>% filter(cl == x) %>% .$ID %>% strsplit(split = '.', fixed = TRUE) %>% sapply('[[', 1) %>% unlist %>% unique

  return(eachG)

}) %>%
  set_names(kmeansRes$cl %>% unique %>% paste0('cluster', .))

kallGOBP <- compareCluster(geneCluster = kall,
                           fun = 'enrichGO',
                           OrgDb = 'org.At.tair.db',
                           keyType= 'TAIR',
                           ont = 'BP',
                           universe = keys(org.At.tair.db),
                           pAdjustMethod = 'BH',
                           pvalueCutoff=0.05,
                           qvalueCutoff=0.1)

enrichGO(gene = kall[[5]],
         OrgDb = 'org.At.tair.db',
         keyType= 'TAIR',
         ont = 'BP',
         universe = keys(org.At.tair.db),
         pAdjustMethod = 'BH',
         pvalueCutoff=0.05,
         qvalueCutoff=0.1) %>%
  as.data.frame %>%
  head

dotplot(kallGOBP, showCategory = 15, font.size = 8)
ggsave('kmeans10_ath_cp_BP_dotplot_15_DEG.jpg', width = 13, height = 12)
ggsave('kmeans10_ath_cp_BP_dotplot_15_DEG.pdf', width = 13, height = 12)

write_csv(as.data.frame(kallGOBP), 'kmeans10_ath_cp_BP_allgene.csv')

## extract all GO
library('GO.db')
library('org.At.tair.db')
library('magrittr')
library('foreach')
library('doParallel')

xx <- as.list(org.At.tairGO) %>% {

    keepIdx <- lapply(., is.na) %>% sapply(sum) %>% is_greater_than(0) %>% not

    return(.[keepIdx])
  }

CollpaseEachGO <- function(eachGOList, geneID) {

  require('foreach')
  require('tidyverse')

  res <- foreach(i = seq_along(eachGOList), .combine = bind_rows) %do% {
    return(unlist(eachGOList[[i]]))
  } %>%
    bind_rows %>%
  mutate(GeneID = geneID)

  return(res)
}

allGOAnno <- Term(GOTERM) %>%
  unlist %>%
  as.data.frame %>%
  rownames_to_column('GOID') %>%
  set_colnames(c('GOID', 'Term')) %>%
  as_tibble

registerDoParallel(cores = 11)
athGO <- foreach(i = seq_along(xx), .combine = bind_rows) %dopar% {
  res <- CollpaseEachGO(xx[[i]], names(xx)[i])
} %>%
  dplyr::filter(Ontology %in% c('BP')) %>%
  dplyr::select(-Evidence, -Ontology) %>%
  dplyr::distinct(GOID, GeneID) %>%
  dplyr::inner_join(allGOAnno) %>%
  dplyr::group_by(GeneID) %>%
  dplyr::summarize_all(~paste(., collapse = ';'))
stopImplicitCluster()

write_csv(athGO, 'table2_ath_GO_anno.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
