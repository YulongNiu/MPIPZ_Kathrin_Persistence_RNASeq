library('tidyverse')
library('limma')
library('DESeq2')

library('tidyverse')
library('ComplexHeatmap')
library('limma')
library('DESeq2')
library('RColorBrewer')


setwd('/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results_rmfull/')

##~~~~~~~~~~~~~~~~~~~~~~~~~ath DEGs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load('degres_condi_Mock_ath.RData')

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

heatsigAth <- (padjsig * fcsig) %>%
  as_tibble %>%
  rowSums %>%
  {. >= 1} %>%
  which %>%
  dplyr::slice(wholeDEG, .) %>%
  inner_join(kmeansRes)

write_csv(heatsigAth, '../results_orthologs/heatsigAth.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~lotus DEGs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load('degres_condi_Mock_lotus.RData')

wholeDEG <- read_csv('SynCom_vs_Mock_lotus_sva_k.csv')
kmeansRes <- read_csv('kmeans_10_lotus_rmfull_rmAtSC.csv') %>%
  select(ID, cl)

fcsig <- wholeDEG %>%
  select(ends_with('FoldChange')) %>%
  transmute_all(list(~ case_when(. > log2(1.5) ~ 1,
                                 . < -log2(1.5) ~ 1,
                                 TRUE ~ 0))) %>%
  select(-matches('fullSC_|AtSC_'))

padjsig <- wholeDEG %>%
  select(ends_with('padj')) %>%
  abs %>%
  `<`(0.05) %>%
  as_tibble %>%
  transmute_all(list(~ if_else(is.na(.), FALSE, .))) %>%
  select(-matches('fullSC_|AtSC_'))

heatsigLotus <- (padjsig * fcsig) %>%
  as_tibble %>%
  rowSums %>%
  {. >= 1} %>%
  which %>%
  dplyr::slice(wholeDEG, .) %>%
  inner_join(kmeansRes)

write_csv(heatsigLotus, '../results_orthologs/heatsigLotus.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

