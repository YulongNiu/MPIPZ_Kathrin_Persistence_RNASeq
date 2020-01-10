library('tidyverse')
library('limma')
library('DESeq2')

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

##~~~~~~~~~~~~~~~~~~~compare collaborator best hits~~~~~~~~~~~~~~~~~~~
library('foreach')

kmeansAth <- read_csv('kmeans_10_ath.csv') %>%
  select(ID, cl) %>%
  dplyr::rename(athID = ID, athcl = cl)

kmeansLotus <- read_csv('kmeans_10_lotus_rmfull_rmAtSC.csv')

## extract Best TAIR
bestTAIR <- kmeansLotus$Best_TAIR %>%
  strsplit(split = ' ', fixed = TRUE) %>%
  sapply('[[', 1) %>% {
    logIdx <- is.na(.)
    .[logIdx] <- kmeansLotus$ID[logIdx]
    return(.)
  }

kmeansLotus <- tibble(lotusID = kmeansLotus$ID,
                      BestTAIR = bestTAIR,
                      lotuscl = kmeansLotus$cl)

mergeKmeans <- inner_join(kmeansAth, kmeansLotus, by = c('athID' = 'BestTAIR'))

ath2lotus <- foreach(i = 1:10) %do% {
  mergeKmeans %>%
    filter(athcl %in% i) %>%
    .$lotuscl %>%
    table %>%
    sort(decreasing = TRUE)
}

names(ath2lotus) <- 1:10

lotus2ath <- foreach(i = 1:10) %do% {
  mergeKmeans %>%
    filter(lotuscl %in% i) %>%
    .$athcl %>%
    table %>%
    sort(decreasing = TRUE)
}

names(lotus2ath) <- 1:10
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~orthogroups~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library('magrittr')
library('stringr')
library('foreach')

## format orthogroups
orthoGFile <- '/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/orthology_prediction_AtLj/AtLCollaborator/Results_Nov11/Orthogroups.csv'


orthoAnno <- readLines(orthoGFile) %>%
  lapply(function(x) {
    eachElem <- strsplit(x, split = '\t', fixed = TRUE) %>%
      unlist %>%
      .[nchar(.) > 0] ## remove empty

    return(eachElem)
  }) %>%
  .[sapply(., length) > 2] %>% ## remove single species
  sapply(function(x) {
    eachElm <- strsplit(x, split = ', ', fixed = TRUE) %>%
      unlist %>%
      str_trim

    return(eachElm)
  }) %>%
  lapply(function(x) {
    tibble(ID = x[-1],
           Orthogroup = x[1])
  }) %>%
  bind_rows

kmeansAth <- read_csv('kmeans_10_ath.csv') %>%
  select(ID, cl) %>%
  dplyr::rename(athID = ID, athcl = cl) %>%
  inner_join(orthoAnno, by = c('athID' = 'ID'))

kmeansLotus <- read_csv('kmeans_10_lotus_rmfull_rmAtSC.csv') %>%
  select(ID, cl) %>%
  dplyr::rename(lotusID = ID, lotuscl = cl) %>%
  inner_join(orthoAnno, by = c('lotusID' = 'ID'))

mergeKmeans <- inner_join(kmeansAth, kmeansLotus)

ath2lotus <- foreach(i = 1:10) %do% {
  mergeKmeans %>%
    filter(athcl %in% i) %>%
    .$lotuscl %>%
    table %>%
    sort(decreasing = TRUE)
}

names(ath2lotus) <- 1:10

lotus2ath <- foreach(i = 1:10) %do% {
  mergeKmeans %>%
    filter(lotuscl %in% i) %>%
    .$athcl %>%
    table %>%
    sort(decreasing = TRUE)
}

names(lotus2ath) <- 1:10
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
