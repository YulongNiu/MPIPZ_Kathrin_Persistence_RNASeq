########################prepare Lotus GO################################
## ## no Lotus japonicus Gifu
## library('AnnotationHub')
## hub <- AnnotationHub()

library('tidyverse')
library('magrittr')
library('AnnotationForge')

## from NCBI
## not Gifu !!!
## makeOrgPackageFromNCBI(version = '0.1',
##                                author = 'Some One <so@someplace.org>',
##                                maintainer = 'Some One <so@someplace.org>',
##                                outputDir = '.',
##                                tax_id = '34305',
##                                genus = 'Lotus',
##                                species = 'corniculatus')

setwd('/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results')

lotusGO <- read_tsv('/extDisk1/Biotools/RefData/lotus_gifu_collaborator_v1.2/Lotusjaponicus_Gifu_v1.2_geneOntology.gaf',
                      comment = '!',
                    col_names = c('DB', 'DB Object ID', 'DB Object Symbol', 'Qualifier', 'GO ID', 'DB:Reference', 'Evidence Code', 'With From', 'Aspect', 'DB Object Name', 'DB Object Synonym', 'DB Object Type', 'Taxon', 'Date', 'Assigned By'))

## GO
fGO <- lotusGO %>%
  dplyr::select('DB Object ID', 'GO ID', 'Evidence Code') %>%
  set_colnames(c('GID', 'GO', 'EVIDENCE'))

## Gene name
fSym <- lotusGO %>%
  dplyr::select(c('DB Object ID', 'DB Object Name')) %>%
  set_colnames(c('GID', 'GENENAME')) %>%
  dplyr::mutate(SYMBOL = GID %>% strsplit(split = '.', fixed = TRUE) %>% sapply('[[', 1)) %>%
  dplyr::select(GID, SYMBOL, everything())

## chromosome
fChr <- read_csv('lotus_gifu_collaborator_Anno_July15.csv') %>%
  dplyr::select('ID', 'Chromosome') %>%
  set_colnames(c('GID', 'CHROMOSOME'))

left_join(fSym, fChr) %>%
  write_csv('tmp1.csv')


########################################################################
