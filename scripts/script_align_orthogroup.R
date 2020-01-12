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

load('/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results_orthologs/orthoAnno.RData')

annoAth <- read_csv('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/Ensembl_ath_Anno.csv',
                    col_types = cols(Chromosome = col_character())) %>%
  mutate_all(list(~replace(., is.na(.), ''))) %>%
  inner_join(orthoAnno) %>%
  setNames(paste0(names(.), '_At'))

annoLotus <- read_csv('/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results/lotus_gifu_collaborator_v1p2_Anno.csv',
                 col_types = cols(Chromosome = col_character())) %>%
  mutate_all(list(~replace(., is.na(.), ''))) %>%
  inner_join(orthoAnno) %>%
  setNames(paste0(names(.), '_Lj'))

## check collaborator annotation and Othofinder
annoMerge <- inner_join(annoAth, annoLotus, by = c('Orthogroup_At' = 'Orthogroup_Lj')) %>%
  dplyr::rename(Orthogroup = Orthogroup_At) %>%
  select(Orthogroup, everything())

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
kresAthCounts <- assay(kresAth)

kresAthOrthoCounts <- split(annoAth$ID_At, annoAth$Orthogroup_At) %>%
  lapply(function(x) {
    eachCount <- match(x, rownames(kresAthCounts)) %>%
      kresAthCounts[., , drop = FALSE] %>%
      colSums

    return(eachCount)
  }) %>% do.call(rbind, .)

kresAthOrtho  <- DESeqDataSetFromMatrix(kresAthOrthoCounts, sampleTable, ~ condition)
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

sampleTable$condition %<>% relevel(ref = 'Mock')

kresLotus <- DESeqDataSetFromTximport(kres, sampleTable, ~ condition)
kresLotusCounts <- assay(kresLotus)

kresLotusOrthoCounts <- split(annoLotus$ID_Lj, annoLotus$Orthogroup_Lj) %>%
  lapply(function(x) {
    eachCount <- match(x, rownames(kresLotusCounts)) %>%
      kresLotusCounts[., , drop = FALSE] %>%
      colSums

    return(eachCount)
  }) %>% do.call(rbind, .)

kresLotusOrtho <- DESeqDataSetFromMatrix(kresLotusOrthoCounts, sampleTable, ~ condition)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

save(kresAthOrtho, kresLotusOrtho, file = '/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results_orthologs/kresOrtho.RData')
#################################################################
