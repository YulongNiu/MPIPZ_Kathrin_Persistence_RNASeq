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
  dplyr::select('DB Object ID', 'GO ID', 'Evidence Code', 'Aspect') %>%
  set_colnames(c('GID', 'GO', 'EVIDENCE', 'Aspect')) %>%
  distinct

## Gene name
fSym <- lotusGO %>%
  dplyr::select('DB Object ID', 'DB Object Name') %>%
  set_colnames(c('GID', 'GENENAME')) %>%
  dplyr::mutate(SYMBOL = GID %>% strsplit(split = '.', fixed = TRUE) %>% sapply('[[', 1)) %>%
  dplyr::select(GID, SYMBOL, everything()) %>%
  distinct

## chromosome
fChr <- read_csv('lotus_gifu_collaborator_Anno_July15.csv') %>%
  dplyr::select(ID, Chromosome) %>%
  set_colnames(c('GID', 'CHROMOSOME')) %>%
  {left_join(fSym, .)} %>%
  dplyr::select(GID, CHROMOSOME) %>%
  distinct

makeOrgPackage(gene_info = fSym,
               chromosome = fChr,
               go = fGO,
               version = '0.1',
               maintainer = 'Yulong Niu <yulong.niu@hotmail.com>',
               author = 'Yulong Niu <yulong.niu@hotmail.com>',
               outputDir = '.',
               tax_id = 'Gifu',
               genus = 'Lotus',
               species = 'japonicus',
               goTable = 'go')

write_csv(fGO %>% dplyr::filter(Aspect == 'P'), 'lotus_GO_BP.csv')
########################################################################



########################prepare Lotus GO gene level#################
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
  dplyr::select('DB Object ID', 'GO ID', 'Evidence Code', 'Aspect') %>%
  set_colnames(c('GID', 'GO', 'EVIDENCE', 'Aspect')) %>%
  mutate(GENE = GID %>% strsplit(split = '.', fixed = TRUE) %>% sapply('[[', 1)) %>%
  dplyr::select(-GID) %>%
  dplyr::select(GENE, everything()) %>%
  distinct

## Gene name
fSym <- lotusGO %>%
  dplyr::select('DB Object ID', 'DB Object Name') %>%
  set_colnames(c('TRANSCRIPT', 'GENENAME')) %>%
  dplyr::mutate(GENE = TRANSCRIPT %>% strsplit(split = '.', fixed = TRUE) %>% sapply('[[', 1)) %>%
  dplyr::mutate(GENENAME = GENENAME %>% strsplit(split = ';', fixed = TRUE) %>% sapply('[[', 1)) %>%
  dplyr::select(GENE, GENENAME) %>%
  distinct %>%
  dplyr::slice(.$GENE %>% duplicated %>% not %>% which) %>%
  distinct

## chromosome
fChr <- read_csv('lotus_gifu_collaborator_Anno_July15.csv') %>%
  dplyr::select(ID, Chromosome) %>%
  set_colnames(c('TRANSCRIPT', 'CHROMOSOME')) %>%
  dplyr::mutate(GENE = TRANSCRIPT %>% strsplit(split = '.', fixed = TRUE) %>% sapply('[[', 1)) %>%
  dplyr::select(-TRANSCRIPT) %>%
  distinct %>%
  {left_join(fSym, .)} %>%
  dplyr::select(GENE, CHROMOSOME) %>%
  distinct

makeOrgPackage(gene_info = fSym %>% dplyr::rename(GID = GENE),
               chromosome = fChr %>% dplyr::rename(GID = GENE),
               go = fGO %>% dplyr::rename(GID = GENE) %>% dplyr::select(-Aspect),
               version = '0.1',
               maintainer = 'Yulong Niu <yulong.niu@hotmail.com>',
               author = 'Yulong Niu <yulong.niu@hotmail.com>',
               outputDir = '.',
               tax_id = 'Gifu',
               genus = 'Lotus',
               species = 'japonicus',
               goTable = 'go')

write_csv(fGO %>% dplyr::filter(Aspect == 'P'), 'lotus_GO_BP.csv')
########################################################################
