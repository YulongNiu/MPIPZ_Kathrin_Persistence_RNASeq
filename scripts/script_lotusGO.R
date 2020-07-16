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
library('tidyverse')
library('magrittr')
library('AnnotationForge')

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


############################Lotus from besthit##########################
library('tidyverse')
library('doParallel')
library('magrittr')
library('GO.db')
library('org.At.tair.db')
library('AnnotationForge')

setwd('/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results')

lotusAnno <- read_tsv('/extDisk1/Biotools/RefData/lotus_gifu_collaborator_v1.2/LjGifu_1.2_functional_annotation.txt', comment = '#') %>%
  dplyr::select(`Protein-Accession`, `Best BlastHit against 'tair'`, `Human-Readable-Description`) %>%
  set_colnames(c('GID', 'TAIR', 'GENENAME')) %>%
  dplyr::mutate(TAIR = TAIR %>% strsplit(split = ' ', fixed = TRUE) %>% sapply('[[', 1)) %>%
  dplyr::mutate(TAIRGene = TAIR %>% strsplit(split = '.', fixed = TRUE) %>% sapply('[[', 1))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~TAIR GO~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
xx <- org.At.tairGO %>%
  as.list %>% {

    keepIdx <- lapply(., is.na) %>% sapply(sum) %>% is_greater_than(0) %>% not

    return(.[keepIdx])
  }

xx <- lotusAnno$TAIRGene %>%
  unique %>%
  {.[!is.na(.)]} %>%
  {xx[match(., names(xx))]}

CollpaseEachGO <- function(eachGOList, geneID) {

  require('foreach')
  require('tidyverse')

  res <- foreach(i = seq_along(eachGOList), .combine = bind_rows) %do% {
    return(unlist(eachGOList[[i]]))
  } %>%
  bind_rows %>%
  mutate(GeneID = geneID)

  return(res)
}

## Time consuming
registerDoParallel(cores = 12)

## foreach(i = 1:500, .combine = bind_rows) %dopar% {
##   res <- CollpaseEachGO(xx[[i]], names(xx)[i])
## }

TAIRGO <- foreach(i = seq_along(xx), .combine = bind_rows) %dopar% {
  res <- CollpaseEachGO(xx[[i]], names(xx)[i])
}

TAIRGO %<>%
  dplyr::select(GOID, GeneID, Evidence) %>%
  distinct
stopImplicitCluster()
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fGO <- inner_join(lotusAnno, TAIRGO, by = c('TAIRGene' = 'GeneID')) %>%
  dplyr::mutate(GID = GID %>% strsplit(split = '.', fixed = TRUE) %>% sapply('[[', 1)) %>%
  dplyr::mutate(Evidence = 'ISS') %>%
  dplyr::rename(GO = GOID, EVIDENCE = Evidence) %>%
  dplyr::select(GID, GO, EVIDENCE) %>%
  distinct

fSym <- lotusAnno %>%
  dplyr::mutate(GID = GID %>% strsplit(split = '.', fixed = TRUE) %>% sapply('[[', 1)) %>%
  dplyr::select(GID, GENENAME) %>%
  dplyr::slice(.$GID %>% duplicated %>% not %>% which) %>%
  distinct

## chromosome
fChr <- read_csv('lotus_gifu_collaborator_Anno_July15.csv') %>%
  dplyr::select(ID, Chromosome) %>%
  set_colnames(c('TRANSCRIPT', 'CHROMOSOME')) %>%
  dplyr::mutate(GID = TRANSCRIPT %>% strsplit(split = '.', fixed = TRUE) %>% sapply('[[', 1)) %>%
  dplyr::select(-TRANSCRIPT) %>%
  distinct %>%
  {inner_join(fSym, .)} %>%
  dplyr::select(GID, CHROMOSOME) %>%
  distinct

inner_join(fChr, fGO)

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
########################################################################
