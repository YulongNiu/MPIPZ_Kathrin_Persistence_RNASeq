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


###########################At GO --> Lotus##############################
## AtSC LjSC Mock
## all gene clusters:
## c2: up-regulated in both AtSC and LjSC
## c4: down-regulated in both AtSC and LjSC, LjSC more
## c3 and c5: up-regulated in only AtSC

meanLotus <- function(v) {

  require('magrittr')

  res <- v %>%
    split(rep(1 : 4, each = 4)) %>%
    sapply(mean, na.rm = TRUE)

  return(res)
}

library('tidyverse')
library('magrittr')
library('ggplot2')

setwd('/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results_rmfull/')

AtCluster <- read_csv('kmeans_10_ath.csv') %>%
  select(ID, cl) %>%
  mutate(TAIRGene =  ID %>% strsplit(split = '.', fixed = TRUE) %>% sapply('[[', 1)) %>%
  rename(TAIR = ID)

##~~~~~~~~~~~~~~~~~~~~~~~~load Lotus~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load('degres_condi_Mock_lotus.RData')

rawC <- rldData %>%
  as.data.frame %>%
  .[, c(-1 : -4)] %>%
  rownames_to_column('ID') %>%
  as_tibble

LjScale <- rawC %>%
  select(-ID) %>%
  t %>%
  scale %>%
  t %>%
  as_tibble %>%
  bind_cols(rawC %>% select(ID))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lotusAnno <- read_tsv('/extDisk1/Biotools/RefData/lotus_gifu_collaborator_v1.2/LjGifu_1.2_functional_annotation.txt', comment = '#') %>%
  dplyr::select(`Protein-Accession`, `Best BlastHit against 'tair'`, `Human-Readable-Description`) %>%
  set_colnames(c('GID', 'TAIR', 'GENENAME')) %>%
  dplyr::mutate(TAIR = TAIR %>% strsplit(split = ' ', fixed = TRUE) %>% sapply('[[', 1)) %>%
  dplyr::mutate(TAIRGene = TAIR %>% strsplit(split = '.', fixed = TRUE) %>% sapply('[[', 1))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~select cluster~~~~~~~~~~~~~~~~~~~~~~~~
clNum <- c(3, 5)

lotusBH <- AtCluster %>%
  filter(cl %in% clNum) %>%
  select(TAIRGene) %>%
  distinct %>%
  inner_join(lotusAnno) %>%
  select(GID, GENENAME) %>%
  distinct %>%
  mutate(GGene = GID %>% strsplit(split = '.', fixed = TRUE) %>% sapply('[[', 1))

## library('clusterProfiler')
## library('org.Ljaponicus.eg.db')

## enrichGO(gene = lotusBH$GGene,
##          OrgDb = 'org.Ljaponicus.eg.db',
##          keyType= 'GID',
##          ont = 'BP',
##          universe = keys(org.Ljaponicus.eg.db),
##          pAdjustMethod = 'BH',
##          pvalueCutoff=0.05,
##          qvalueCutoff=0.1) %>%
##   as.data.frame %>%
##   .[, 2]

lotusBH %>%
  select(GID) %>%
  inner_join(LjScale, c('GID' = 'ID')) %>%
  select(-GID) %>%
  apply(1, meanLotus) %>%
  t %>%
  as_tibble %>%
  set_colnames(c('L_AtSC', 'L_AtSCMloti', 'L_LjSC', 'L_mock')) %>%
  gather(key = 'Sample', value = 'ScaleC') %>%
  ggplot(aes(x = Sample, y = ScaleC)) +
  geom_boxplot()

ggsave('At_cl35_lotus_scale.jpg')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########################################################################
