
## originally by Yulong Niu
## yulong.niu@hotmail.com

##########################Best hit######################################
##~~~~~~~~~~~~~~~~~~~~~~~~~Useful funcs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RBH_ <- function(eachf, allt, threshold = 1.0e-6) {
  ## INPUT: `eachf` is tibble of from genes.
  ##        `allt` is tibble of all to genes.
  ##        `threshold` is the BLAST evalue threshold
  ## OUTPUT: A tibble of reciprocal best hits.

  require('tidyverse')

  eachf %<>%
    filter(evalue < threshold) %>%
    arrange(bitscore %>% desc())

  ## check genef exists
  if (nrow(eachf) == 0) {
    return(eachf)
  } else {}

  genef1 <- eachf$qseqid[1]
  genet1 <- eachf$sseqid[1]

  eacht <- allt %>%
    filter(qseqid %in% genet1, evalue < threshold) %>%
    arrange(bitscore %>% desc())
  genet2 <- eacht$sseqid[1]

  ## check genet exists
  if (nrow(eacht) == 0) {
    return(eacht)
  } else {}

  if (genef1 == genet2) {
    return(bind_rows(eachf[1, , drop = FALSE],
              eacht[1, , drop = FALSE]))
  } else {
    return(eachf %>% filter(FALSE))
  }
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library('tidyverse')

setwd('/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/orthology_prediction_AtLj/AtLCollaborator/Results_Nov11/WorkingDirectory/')

blastf <- read_delim('Blast0_1.txt',
                     delim = '\t',
                     col_names = c('qseqid', 'sseqid',
                                   'pident', 'length',
                                   'mismatch', 'gapopen',
                                   'qstart', 'qend',
                                   'sstart', 'send',
                                   'evalue', 'bitscore'))

blastt <- read_delim('Blast1_0.txt',
                     delim = '\t',
                     col_names = c('qseqid', 'sseqid',
                                   'pident', 'length',
                                   'mismatch', 'gapopen',
                                   'qstart', 'qend',
                                   'sstart', 'send',
                                   'evalue', 'bitscore'))

rbhMat <- blastf %>%
  group_by(qseqid) %>%
  do(RBH_(., blastt)) %>%
  ungroup

save(rbhMat, file = '/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results_orthologs/rbhMat.RData')
#######################################################################


################################RBM anno###############################
library('tidyverse')
library('magrittr')

setwd('/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results_orthologs')

load('rbhMat.RData')

## preprocess anno
anno <- read_delim('/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/orthology_prediction_AtLj/AtLCollaborator/Results_Nov11/WorkingDirectory/SequenceIDs.txt', delim = '\t', col_names = c('Annotation')) %>%
  .$Annotation %>%
  strsplit(split = ':', fixed = TRUE) %>%
  lapply(function(x) {
    x %>%
      str_trim %>%
      {
        eachBlastID <- .[1]
        eachID <- .[2] %>%
          strsplit(split = ' ', fixed = TRUE) %>%
          sapply('[', 1)
        return(tibble(BlastID = eachBlastID, ID = eachID))
      }
  }) %>%
  bind_rows

rbhAnno <- rbhMat %>%
  as_tibble %>%
  slice(seq(1, nrow(.), 2)) %>%
  inner_join(anno, by = c('qseqid' = 'BlastID')) %>%
  inner_join(anno, by = c('sseqid' = 'BlastID')) %>%
  rename(ID_ath = ID.x, ID_lotus = ID.y) %>%
  select(ID_ath, ID_lotus, everything()) %>% {
    tibble(ID = c(.$ID_ath, .$ID_lotus),
          RBH = rep(1:nrow(.), 2) %>% paste0('RBH', .))
  }

save(rbhAnno, file = 'rbhAnno.RData')
#######################################################################
