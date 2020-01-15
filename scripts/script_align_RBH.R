###########################DEGs##################################
##~~~~~~~~~~~~~~~~~~~~~~~~~~~useful funcs~~~~~~~~~~~~~~~~~~~~~~~~
checkZeros <- function(v, threshold) {
  res <- ifelse(sum(v == 0) > threshold, FALSE, TRUE)
  return(res)
}

checkPersis <- function(v, threshold) {

  require('magrittr')

  res <- v %>%
    split(rep(1 : 4, each = 4)) %>%
    sapply(checkZeros, threshold) %>%
    all

  return(res)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library('tximport')
library('rhdf5')
library('magrittr')
library('DESeq2')
library('tidyverse')
library('foreach')
library('ParaMisc')

load('/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results_orthologs/rbhAnno.RData')

annoAth <- read_csv('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/Ensembl_ath_Anno.csv',
                    col_types = cols(Chromosome = col_character())) %>%
  mutate_all(list(~replace(., is.na(.), ''))) %>%
  inner_join(rbhAnno) %>%
  setNames(paste0(names(.), '_At'))
c
annoLotus <- read_csv('/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results/lotus_gifu_collaborator_v1p2_Anno.csv',
                 col_types = cols(Chromosome = col_character())) %>%
  mutate_all(list(~replace(., is.na(.), ''))) %>%
  inner_join(rbhAnno) %>%
  setNames(paste0(names(.), '_Lj'))

## check collaborator annotation and Othofinder
annoMerge <- inner_join(annoAth, annoLotus, by = c('RBH_At' = 'RBH_Lj')) %>%
  dplyr::rename(RBH = RBH_At) %>%
  select(RBH, everything())

##~~~~~~~~~~~~~~~~~~~~load ath alignments~~~~~~~~~~~~~~~~~~~~~~~~~~
wd <- '/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/align_data/ath'
setwd(wd)

labelanno <- read_delim('../../results/mapping.txt', delim = '\t') %>%
  dplyr::rename(ID = `library_number`, SampleAnno = `library_name`) %>%
  mutate(ID = ID %>% str_replace('\\.', '_')) %>%
  filter(species %>% str_detect('Ath'))

slabel <- labelanno$SampleAnno %>%
  paste0('_ath_kallisto')

kres <- file.path(wd, slabel, 'abundance.h5') %>%
  set_names(labelanno$SampleAnno) %>%
  tximport(type = 'kallisto', txOut = TRUE)

condi <- c('fullSC_At', 'AtSC_At', 'LjSC_At', 'Mock_At')

sampleTable <- data.frame(condition = factor(rep(condi, each = 4), levels = condi)) %>%
  set_rownames(colnames(kres$counts))

kresAth <- DESeqDataSetFromTximport(kres, sampleTable, ~ condition)
kresAth %<>%
  {.[match(annoMerge$ID_At, rownames(.))]} %>%
  set_rownames(annoMerge$RBH)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~load lotus alignments~~~~~~~~~~~~~~~~~~~~~~~~~~
wd <- '/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/align_data/lotus_collaborator'
setwd(wd)

labelanno <- read_delim('../../results/mapping.txt', delim = '\t') %>%
  dplyr::rename(ID = `library_number`, SampleAnno = `library_name`) %>%
  mutate(ID = ID %>% str_replace('\\.', '_')) %>%
  filter(species %>% str_detect('Lj'))

slabel <- labelanno$SampleAnno %>%
  paste0('_lotus_kallisto')

kres <- file.path(wd, slabel, 'abundance.h5') %>%
  set_names(labelanno$SampleAnno) %>%
  tximport(type = 'kallisto', txOut = TRUE)

condi <- c('fullSC_Lj', 'AtSC_Lj', 'AtSCMloti_Lj', 'LjSC_Lj', 'Mock_Lj')

sampleTable <- data.frame(condition = factor(rep(condi, each = 4), levels = condi))
rownames(sampleTable) <- colnames(kres$counts)

kresLotus <- DESeqDataSetFromTximport(kres, sampleTable, ~ condition)
kresLotus %<>%
  {.[match(annoMerge$ID_Lj, rownames(kresLotus))]} %>%
  set_rownames(annoMerge$RBH)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

save(kresAth, kresLotus, file = '/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results_orthologs/kresRBH.RData')
#################################################################
