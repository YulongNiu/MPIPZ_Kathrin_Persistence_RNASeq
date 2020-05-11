
## originally by Yulong Niu
## yulong.niu@hotmail.com

################################orthologs###############################
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

interMat <- matrix(ncol = 10, nrow = 10, dimnames = list(paste0('At', 1:10), paste0('Lj', 1:10)))
for (i in seq_len(10)) {

  interNum <- mergeKmeans %>%
    filter(athcl %in% i) %>%
    .$lotuscl %>%
    table

  unionNum <- ((mergeKmeans$athcl %in% i) %>% length) + (mergeKmeans$lotuscl %>% table) - interNum

  interMat[i, ] <- interNum / unionNum
}

interMat %>% round(digits = 5)

## only DEGs
mergeKmeansDEG <- read_csv('../results_orthologs/heatsigAth.csv') %>%
  select(ID) %>%
  dplyr::rename(athID = ID) %>%
  inner_join(mergeKmeans)

kmeansLotusDEG <- read_csv('../results_orthologs/heatsigLotus.csv') %>%
  select(ID) %>%
  dplyr::rename(lotusID = ID) %>%
  inner_join(mergeKmeansDEG)
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

save(orthoAnno, file = '../results_orthologs/orthoAnno.RData')

## all transcripts
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

interMat <- matrix(ncol = 10, nrow = 10, dimnames = list(paste0('At', 1:10), paste0('Lj', 1:10)))
for (i in seq_len(10)) {

  interNum <- mergeKmeans %>%
    filter(athcl %in% i) %>%
    .$lotuscl %>%
    table

  unionNum <- ((mergeKmeans$athcl %in% i) %>% length) + (mergeKmeans$lotuscl %>% table) - interNum

  interMat[i, ] <- interNum / unionNum
}

interMat %>% round(digits = 5)

## only DEGs
kmeansAthDEG <- read_csv('../results_orthologs/heatsigAth.csv') %>%
  select(ID, cl) %>%
  dplyr::rename(athID = ID, athcl = cl) %>%
  inner_join(orthoAnno, by = c('athID' = 'ID'))

kmeansLotusDEG <- read_csv('../results_orthologs/heatsigLotus.csv') %>%
  select(ID, cl) %>%
  dplyr::rename(lotusID = ID, lotuscl = cl) %>%
  inner_join(orthoAnno, by = c('lotusID' = 'ID'))

mergeKmeansDEG <- inner_join(kmeansAthDEG, kmeansLotusDEG)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~orthogroup compare~~~~~~~~~~~~~~~~
kmeansAth <- read_csv('../results_orthologs/kmeans_10_ath_og_rmfull.csv') %>%
  rename(athcl = cl)

kmeansLotus <- read_csv('../results_orthologs/kmeans_10_lotus_og_rmfull.csv') %>%
  rename(lotuscl = cl)

mergeKmeans <- inner_join(kmeansAth, kmeansLotus)

## jaccard similarity
interMat <- matrix(ncol = 10, nrow = 10, dimnames = list(paste0('At', 1:10), paste0('Lj', 1:10)))
for (i in seq_len(10)) {

  interNum <- mergeKmeans %>%
    filter(athcl %in% i) %>%
    .$lotuscl %>%
    table

  unionNum <- ((mergeKmeans$athcl %in% i) %>% length) + (mergeKmeans$lotuscl %>% table) - interNum

  interMat[i, ] <- interNum / unionNum
}

interMat %>% round(digits = 4)

interMat <- matrix(ncol = 10, nrow = 10, dimnames = list(paste0('At', 1:10), paste0('Lj', 1:10)))
for (i in seq_len(10)) {

  interNum <- mergeKmeans %>%
    filter(athcl %in% i) %>%
    .$lotuscl %>%
    table

  interMat[i, ] <- mergeKmeans %>%
    filter(athcl %in% i) %>%
    nrow %>%
    {interNum / .}
}

interMat %>% round(digits = 4)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#######################################################################

################################RBH###############################
library('tidyverse')
library('limma')
library('DESeq2')

setwd('/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results_orthologs/')

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~separate~~~~~~~~~~~~~~~~~~~~~~~~~~~
## all transcripts
kmeansAth <- read_csv('kmeans_12_RBH_ath_rmfull.csv') %>%
  select(ID, cl) %>%
  dplyr::rename(athcl = cl)

kmeansLotus <- read_csv('kmeans_12_RBH_lotus_rmfull.csv') %>%
  select(ID, cl) %>%
  dplyr::rename(lotuscl = cl)

mergeKmeans <- inner_join(kmeansAth, kmeansLotus)

ath2lotus <- foreach(i = 1:12) %do% {
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

interMat <- matrix(ncol = 12, nrow = 12, dimnames = list(paste0('At', 1:12), paste0('Lj', 1:12)))
for (i in seq_len(12)) {

  interNum <- mergeKmeans %>%
    filter(athcl %in% i) %>%
    .$lotuscl %>%
    table

  unionNum <- ((mergeKmeans$athcl %in% i) %>% length) + (mergeKmeans$lotuscl %>% table) - interNum

  interMat[i, ] <- interNum / unionNum
}

interMat %>% round(digits = 5)

## only DEGs
kmeansAthDEG <- read_csv('../results_orthologs/heatsigAth.csv') %>%
  select(ID, cl) %>%
  dplyr::rename(athID = ID, athcl = cl) %>%
  inner_join(orthoAnno, by = c('athID' = 'ID'))

kmeansLotusDEG <- read_csv('../results_orthologs/heatsigLotus.csv') %>%
  select(ID, cl) %>%
  dplyr::rename(lotusID = ID, lotuscl = cl) %>%
  inner_join(orthoAnno, by = c('lotusID' = 'ID'))

mergeKmeansDEG <- inner_join(kmeansAthDEG, kmeansLotusDEG)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~map RBH to raw heatmap~~~~~~~~~~~~~~~~~
## raw cluster
rawPath <- '/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results_rmfull/'

## all transcripts
kmeansAth <- read_csv('kmeans_12_RBH_ath_rmfull.csv') %>%
  select(ID, ID_At) %>%
  dplyr::rename(RBHID = ID) %>% {
    rawCluster <- file.path(rawPath, 'kmeans_10_ath.csv') %>%
      read_csv

      inner_join(rawCluster, ., c('ID' = 'ID_At'))
  } %>%
  select(RBHID, cl) %>%
  dplyr::rename(athcl = cl)

kmeansLotus <- read_csv('kmeans_12_RBH_lotus_rmfull.csv') %>%
  select(ID, ID_Lj) %>%
  dplyr::rename(RBHID = ID) %>% {
    rawCluster <- file.path(rawPath, 'kmeans_10_lotus_rmfull_rmAtSC.csv') %>%
      read_csv

    inner_join(rawCluster, ., c('ID' = 'ID_Lj'))
  } %>%
  select(RBHID, cl) %>%
  dplyr::rename(lotuscl = cl)

mergeKmeans <- inner_join(kmeansAth, kmeansLotus)

interMat <- matrix(ncol = 10, nrow = 10, dimnames = list(paste0('At', 1:10), paste0('Lj', 1:10)))
for (i in seq_len(10)) {

  interNum <- mergeKmeans %>%
    filter(athcl %in% i) %>%
    .$lotuscl %>%
    table

  unionNum <- ((mergeKmeans$athcl %in% i) %>% length) + (mergeKmeans$lotuscl %>% table) - interNum

  interMat[i, ] <- interNum / unionNum
}

interMat %>% round(digits = 5)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################

########################Ath cluster to Lotus######################
library('tidyverse')
library('scales')
library('cowplot')

setwd('/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results_orthologs/')

clNum <- 12

athCluster <- read_csv('kmeans_12_RBH_ath_rmfull.csv') %>%
  select(ID, cl) %>%
  inner_join(read_csv('kmeans_12_RBH_ath_lotus_rmfull_clusterGene.csv')) %>%
  select(ID, AtSC_At : Mock_Lj, cl)

sampleNAth <- c('AtSC_At', 'LjSC_At', 'Mock_At')
sampleNLotus <- c('AtSCMloti_Lj', 'LjSC_Lj', 'Mock_Lj')

clusterCore <- athCluster %>%
  select(ID, AtSC_At : Mock_At, cl) %>%
  group_by(cl) %>%
  summarise_at(-1, mean, na.rm = TRUE) %>% ## mean of each cluster
  mutate(cl = paste0('cluster_', cl) %>%
           factor(levels = paste0('cluster_', cl))) %>%
  gather(Sample, NorExpress, -1) %>%
  mutate(host = Sample %>%
           strsplit(split = '_', fixed = TRUE) %>%
           sapply('[', 2) %>%
           paste(cl, ., sep = '_')) %>%
  mutate(Sample = Sample %>% factor(levels = sampleNAth, ordered = TRUE))

athClusterPlot <- ggplot(clusterCore, aes(Sample, NorExpress, col = cl, group = host)) +
  geom_point() +
  geom_line() +
  facet_wrap(. ~ cl, ncol = 2) +
  ylab('Scaled counts') +
  scale_color_manual(values = hue_pal()(clNum),
                     breaks = athCluster$cl %>%
                       table %>%
                       names %>%
                       paste0('cluster_', .),
                     labels = athCluster$cl %>%
                       table %>%
                       {paste0('cluster_', names(.), ' ', .)},
                     guide = guide_legend(title = paste0('kmeans (k = ',clNum, ')'))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

clusterCore <- athCluster %>%
  select(ID, AtSCMloti_Lj, LjSC_Lj, Mock_Lj, cl) %>%
  group_by(cl) %>%
  summarise_at(-1, mean, na.rm = TRUE) %>% ## mean of each cluster
  mutate(cl = paste0('cluster_', cl) %>%
           factor(levels = paste0('cluster_', cl))) %>%
  gather(Sample, NorExpress, -1) %>%
  mutate(host = Sample %>%
           strsplit(split = '_', fixed = TRUE) %>%
           sapply('[', 2) %>%
           paste(cl, ., sep = '_')) %>%
  mutate(Sample = Sample %>% factor(levels = sampleNLotus, ordered = TRUE))

lotusClusterPlot <- ggplot(clusterCore, aes(Sample, NorExpress, col = cl, group = host)) +
  geom_point() +
  geom_line() +
  facet_wrap(. ~ cl, ncol = 2) +
  ylab('Scaled counts') +
  scale_color_manual(values = hue_pal()(clNum),
                     breaks = athCluster$cl %>%
                       table %>%
                       names %>%
                       paste0('cluster_', .),
                     labels = athCluster$cl %>%
                       table %>%
                       {paste0('cluster_', names(.), ' ', .)},
                     guide = guide_legend(title = paste0('kmeans (k = ',clNum, ')'))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

plot_grid(athClusterPlot, lotusClusterPlot)
ggsave('kmeans_12_RBH_ath2lotus.jpg', width = 15)
##################################################################

########################Lotus cluster to Ath######################
library('tidyverse')
library('scales')
library('cowplot')

setwd('/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results_orthologs/')

clNum <- 12

lotusCluster <- read_csv('kmeans_12_RBH_lotus_rmfull.csv') %>%
  select(ID, cl) %>%
  inner_join(read_csv('kmeans_12_RBH_ath_lotus_rmfull_clusterGene.csv')) %>%
  select(ID, AtSC_At : Mock_Lj, cl)

sampleNAth <- c('AtSC_At', 'LjSC_At', 'Mock_At')
sampleNLotus <- c('AtSCMloti_Lj', 'LjSC_Lj', 'Mock_Lj')

clusterCore <- lotusCluster %>%
  select(ID, AtSC_At : Mock_At, cl) %>%
  group_by(cl) %>%
  summarise_at(-1, mean, na.rm = TRUE) %>% ## mean of each cluster
  mutate(cl = paste0('cluster_', cl) %>%
           factor(levels = paste0('cluster_', cl))) %>%
  gather(Sample, NorExpress, -1) %>%
  mutate(host = Sample %>%
           strsplit(split = '_', fixed = TRUE) %>%
           sapply('[', 2) %>%
           paste(cl, ., sep = '_')) %>%
  mutate(Sample = Sample %>% factor(levels = sampleNAth, ordered = TRUE))

athClusterPlot <- ggplot(clusterCore, aes(Sample, NorExpress, col = cl, group = host)) +
  geom_point() +
  geom_line() +
  facet_wrap(. ~ cl, ncol = 2) +
  ylab('Scaled counts') +
  scale_color_manual(values = hue_pal()(clNum),
                     breaks = lotusCluster$cl %>%
                       table %>%
                       names %>%
                       paste0('cluster_', .),
                     labels = lotusCluster$cl %>%
                       table %>%
                       {paste0('cluster_', names(.), ' ', .)},
                     guide = guide_legend(title = paste0('kmeans (k = ',clNum, ')'))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

clusterCore <- lotusCluster %>%
  select(ID, AtSCMloti_Lj, LjSC_Lj, Mock_Lj, cl) %>%
  group_by(cl) %>%
  summarise_at(-1, mean, na.rm = TRUE) %>% ## mean of each cluster
  mutate(cl = paste0('cluster_', cl) %>%
           factor(levels = paste0('cluster_', cl))) %>%
  gather(Sample, NorExpress, -1) %>%
  mutate(host = Sample %>%
           strsplit(split = '_', fixed = TRUE) %>%
           sapply('[', 2) %>%
           paste(cl, ., sep = '_')) %>%
  mutate(Sample = Sample %>% factor(levels = sampleNLotus, ordered = TRUE))

lotusClusterPlot <- ggplot(clusterCore, aes(Sample, NorExpress, col = cl, group = host)) +
  geom_point() +
  geom_line() +
  facet_wrap(. ~ cl, ncol = 2) +
  ylab('Scaled counts') +
  scale_color_manual(values = hue_pal()(clNum),
                     breaks = lotusCluster$cl %>%
                       table %>%
                       names %>%
                       paste0('cluster_', .),
                     labels = lotusCluster$cl %>%
                       table %>%
                       {paste0('cluster_', names(.), ' ', .)},
                     guide = guide_legend(title = paste0('kmeans (k = ',clNum, ')'))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

plot_grid(lotusClusterPlot, athClusterPlot)
ggsave('kmeans_12_RBH_lotus2ath.jpg', width = 15)
##################################################################
