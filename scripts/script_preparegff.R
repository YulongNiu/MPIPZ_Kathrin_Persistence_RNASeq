########################deal with gff3########################
library('stringr')
library('utils')
library('readr')
library('dplyr')
library('magrittr')

gffPath <- '/extDisk1/Biotools/RefData/lotus_gifu_collaborator_v1.2/20190809_Lj_Gifu_v1.2_predictedGenes.gff3'
gffAnno <- read_tsv(gffPath,
                    col_names = c('chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'),
                    col_types = cols(chromosome = col_character()),
                    comment = '#')

##~~~~~~~~~~~~~~~~~~~~~~~~~~gene and lncgene table~~~~~~~~~~~~~~
geneAnno <- gffAnno %>%
  filter(feature %in% c('mRNA'))

## "ID=LotjaGi0g1v0000400.1;Parent=LotjaGi0g1v0000400;human_readable_description=protein kinase family protein;interpro_id=IPR011009 (Protein kinase-like domain), IPR011992 (EF-hand domain pair);gene_ontology_id=N.A.;best_blast_hit_tair=AT5G19450.1 calcium-dependent protein kinase 19;best_blast_hit_sprot=sp|Q84SL0|CDPKK_ORYSJ Calcium-dependent protein kinase 20;best_blast_hit_trembl_plants=tr|A0A151TU37|A0A151TU37_CAJCA Calcium-dependent protein kinase 32"
noteAnno <- geneAnno %>%
  select(attribute) %>%
  unlist %>%
  unname %>%
  str_trim

geneTable <- tibble(ID = str_extract(noteAnno, '(?<=ID=).*?(?=;)'),
                    cDNA = str_extract(noteAnno, '(?<=Parent=).*?(?=;)'),
                    Chromosome = geneAnno$chromosome,
                    Start = geneAnno$start,
                    End = geneAnno$end,
                    Strand = geneAnno$strand,
                    Interpro_ID = str_extract(noteAnno, '(?<=interpro_id=).*?(?=;)') %>% {if_else(is.na(.), '', .)} %>% sapply(URLdecode),
                    Ontology_ID = str_extract(noteAnno, '(?<=gene_ontology_id=).*?(?=;)') %>% {if_else(is.na(.), '', .)} %>% sapply(URLdecode),
                    Best_TAIR = str_extract(noteAnno, '(?<=best_blast_hit_tair=).*?(?=;)') %>% {if_else(is.na(.), '', .)} %>% sapply(URLdecode),
                    Best_sport = str_extract(noteAnno, '(?<=best_blast_hit_sprot=).*?(?=;)') %>% {if_else(is.na(.), '', .)} %>% sapply(URLdecode),
                    Best_trembl_plants = str_extract(noteAnno, '(?<=best_blast_hit_trembl_plants=).*?(?=;)') %>% {if_else(is.na(.), '', .)} %>% sapply(URLdecode),
                    Description = str_extract(noteAnno, '(?<=human_readable_description=).*?(?=;)') %>% {if_else(is.na(.), '', .)} %>% sapply(URLdecode)) %>%
  mutate_all(list(~str_replace(., 'N.A.', '')))


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write_csv(geneTable, '/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results/lotus_gifu_collaborator_v1p2_Anno.csv')

##~~~~~~~~~~~~~~~~~~check cDNA in k~~~~~~~~~~~~~~~~~~~~~~~~~~~
kcDNA <- read_tsv('/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/align_data/lotus_collaborator/L_mock_4_lotus_kallisto//abundance.tsv')

## kcDNA 48359
## athAnno 54013

## all kcDNA in athAnno
sum(kcDNA$target_id %in% geneTable$ID)

## lncRNA  miRNA  ncRNA   rRNA snoRNA  snRNA   tRNA
##   3879    325    377     15    287     82    689
anti_join(geneTable, kcDNA, by = c('ID' = 'target_id')) %>%
  write_csv('tmp1.csv')
  .$BioType %>%
  table
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##############################################################
