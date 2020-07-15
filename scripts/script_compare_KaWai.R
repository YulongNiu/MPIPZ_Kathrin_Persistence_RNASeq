########################compare Ka-Wai and Kathrin######################
library('tidyverse')
library('magrittr')
library('foreach')

setwd('/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results_rmfull/')
load('degres_condi_Mock_ath.RData')

##~~~~~~~~~~~~~~~~~~~~~~~Kathrin DEGs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

KathrinDEG <- scaleC %>%
  select(ID, cl)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Ka-Wai DEGs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
KaWaiDEG <- read_csv('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/removeZero/kmeans10_soil_sig.csv') %>%
  select(ID, cl)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~compare DEGs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
KathrinKaWaiJS <- foreach (i = 1:10, .combine = rbind) %do% {

  KathrinEach <- KathrinDEG %>%
    filter(cl == i) %>%
    .$ID %>%
    substr(start = 1, stop = nchar(.) - 2) %>%
    unique

  JacSim <- foreach(j = 1:10, .combine = c)  %do% {
    KaWaiEach <- KaWaiDEG %>%
      filter(cl == j) %>%
      .$ID %>%
      unique

    eachJacSim <- length(intersect(KathrinEach, KaWaiEach)) / length(union(KathrinEach, KaWaiEach))

    return(eachJacSim)

  } %>%
  set_names(paste0('Kathrin', 1:10))

  return(JacSim)
} %>%
  set_rownames(paste0('Ka-Wai', 1:10)) %>%
  round(digits = 3)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~compare all genes~~~~~~~~~~~~~~~~~~~~~~~~~~~~
KathrinAll <- read_csv('kmeans_10_ath.csv') %>%
  select(ID, cl)

KaWaiAll <- read_csv('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/removeZero/kmeans10_soil.csv') %>%
  select(ID, cl)

KathrinKaWaiJS <- foreach (i = 1:10, .combine = rbind) %do% {

  KathrinEach <- KathrinAll %>%
    filter(cl == i) %>%
    .$ID %>%
    substr(start = 1, stop = nchar(.) - 2) %>%
    unique

  JacSim <- foreach(j = 1:10, .combine = c)  %do% {
    KaWaiEach <- KaWaiAll %>%
      filter(cl == j) %>%
      .$ID %>%
      unique

    eachJacSim <- length(intersect(KathrinEach, KaWaiEach)) / length(union(KathrinEach, KaWaiEach))

    return(eachJacSim)

  } %>%
  set_names(paste0('Kathrin', 1:10))

  return(JacSim)
} %>%
  set_rownames(paste0('Ka-Wai', 1:10)) %>%
  round(digits = 3)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########################################################################
