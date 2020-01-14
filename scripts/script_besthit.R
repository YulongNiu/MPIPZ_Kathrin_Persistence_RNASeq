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

  genef <- eachf$qseqid[1]

  eacht <- allt %>%
    filter(sseqid %in% genef, evalue < threshold) %>%
    arrange(bitscore %>% desc())
  genet <- eacht$sseqid[1]

  ## check genet exists
  if (nrow(eacht) == 0) {
    return(eacht)
  } else {}

  if (genef == genet) {
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
  do(RBH_(., blastt))
#######################################################################
