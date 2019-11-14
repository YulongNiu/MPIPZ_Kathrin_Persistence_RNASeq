#######################merge plants and SynCom########################
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~useful funcs~~~~~~~~~~~~~~~~~~~~~~~~~~~
BacSavePath <- function(strainName) {
  ## INPUT: `strainName`: a character vector containing strain names.
  ## OUTPUT: a character vector containing save path of strain.

  majorPath <- '/netscratch/dep_psl/grp_rgo/yniu/ref'
  rootPath <- file.path(majorPath, 'RootRef')
  ljPath <- file.path(majorPath, 'LjRef')

  res <- character(length = length(strainName))

  res[str_detect(strainName, '^Root')] <- rootPath
  res[str_detect(strainName, '^Lj')] <- ljPath

  return(res)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library('tidyverse')
library('magrittr')
library('foreach')
library('doParallel')

ZLESS_PATH <- '/bin/zless'
GZIP_PATH <- '/bin/gzip'
REF_PATH <- '/netscratch/dep_psl/grp_rgo/yniu/ref'
SAVE_PATH <- file.path(REF_PATH, 'Kathrin_SynCom')
CAT_PATH <- '/bin/cat'

setwd('/netscratch/dep_psl/grp_rgo/yniu/KathrinPersistence/results/')

syncom <- read_csv('ath_lj_syncom.csv')

cond <- c('C_fSC', 'C_AtSC', 'C_LjSC',
          'L_fSC', 'L_AtSC', 'L_AtSCMloti', 'L_LjSC', 'L_mock')

## syncom component
syncomAt <- syncom$At_strains
syncomLj <- syncom$Lj_strains
syncomList <- list(c(syncomAt, syncomLj), syncomAt, syncomLj,
                   c(syncomAt, syncomLj), syncomAt, c(syncomAt, 'LjNodule218'), syncomLj, 'LjNodule218') %>%
  set_names(cond)


registerDoParallel(cores = 7)
foreach (i = seq_along(cond), .combine = c) %dopar% {

  fn <- cond[i] %>%
    file.path(SAVE_PATH, .)

  if(!file.exists(fn)) {
    dir.create(fn)
  } else {}

  if (str_detect(cond[i], '^C')) {
    plantGenome <- file.path(REF_PATH, 'ath/Arabidopsis_thaliana.TAIR10.latentvirus1.dna.toplevel.fa')
    plantcDNA <- file.path(REF_PATH, 'ath/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz')
  } else if (str_detect(cond[i], '^L')) {
    plantGenome <- file.path(REF_PATH,'lotus_gifu_collaborator_v1.2/LjGifu1.1_pseudomol.fa')
    plantcDNA <- file.path(REF_PATH, 'lotus_gifu_collaborator_v1.2/LjGifu1.2_cds.fasta')
  } else {}

  ## merge genome
  outfile <- cond[i] %>%
    paste0('_genome.fasta') %>%
    file.path(fn, .)

  commandL <- syncomList[[i]] %>%
    paste0('.fna') %>%
    {
      eachfp <- BacSavePath(.)
      file.path(eachfp, .)
    } %>%
    paste(collapse = ' ') %>%
    paste(ZLESS_PATH, plantGenome, ., '>', outfile, collapse = ' ')

  system(commandL)

  ## merge cDNA
  outfile <- cond[i] %>%
    paste0('_cDNA.fasta') %>%
    file.path(fn, .)

  commandL <- syncomList[[i]] %>%
    paste0('.ffn') %>%
    {
      eachfp <- BacSavePath(.)
      file.path(eachfp, .)
    } %>%
    paste(collapse = ' ') %>%
    paste(ZLESS_PATH, plantcDNA, ., '>', outfile, collapse = ' ')

  system(commandL)

  return(NULL)
}
stopImplicitCluster()
######################################################################
