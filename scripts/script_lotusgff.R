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
  filter(feature %in% c('gene', 'ncRNA_gene'))

## "ID=LotjaGi0g1v0000100"
noteAnno <- geneAnno %>%
  select(attribute) %>%
  unlist %>%
  unname %>%
  str_trim

geneTable <- tibble(ID = str_extract(noteAnno, '(?<=ID=).*'),
                    Gene = str_extract(noteAnno, '(?<=Name=).*?(?=;)') %>% {if_else(is.na(.), '', .)} %>% sapply(URLdecode),
                    BioType = str_extract(noteAnno, '(?<=biotype=).*?(?=;)'),
                    Description = str_extract(noteAnno, '(?<=description=).*?(?=;)') %>% {if_else(is.na(.), '', .)} %>% sapply(URLdecode))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~cDNA table~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cDNAanno <- gffAnno %>%
  filter(feature %in% c('lnc_RNA', 'miRNA', 'mRNA', 'ncRNA',
                        'rRNA', 'snoRNA', 'snRNA', 'tRNA'))

## "ID=LotjaGi0g1v0000100.1;Parent=LotjaGi0g1v0000100;human_readable_description=Molybdenum cofactor sulfurase;interpro_id=IPR015424 (Pyridoxal phosphate-dependent transferase);gene_ontology_id=GO:0003824;best_blast_hit_tair=AT5G51920.1 Pyridoxal phosphate (PLP)-dependent transferases superfamily protein;best_blast_hit_sprot=sp|Q9C5X8|MOCOS_ARATH Molybdenum cofactor sulfurase;best_blast_hit_trembl_plants=tr|A0A0S3TBS8|A0A0S3TBS8_PHAAN Uncharacterized protein"
noteAnno <- cDNAanno %>%
  select(attribute) %>%
  unlist %>%
  unname %>%
  str_trim

cDNATable <- tibble(cDNA = str_extract(noteAnno, '(?<=ID=).*?(?=;)'),
                    ID = str_extract(noteAnno, '(?<=Parent=).*?(?=;)'),
                    Chromosome = cDNAanno$chromosome,
                    Start = cDNAanno$start,
                    End = cDNAanno$end,
                    Strand = cDNAanno$strand)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## merge cDNA table and gene table
creAnno <- inner_join(geneTable, cDNATable, by = 'ID') %>%
  select(-ID) %>%
  select(cDNA, everything()) %>%
  rename(ID = cDNA) %>%
  mutate(Length = abs(End - Start)) %>%
  select(ID, Gene, Chromosome:Length, BioType, Description)

write_csv(creAnno, '/extDisk1/RESEARCH/MPIPZ_Chlamy_RNASeq/results/Ensembl_cre_Anno.csv')

##~~~~~~~~~~~~~~~~~~check cDNA in k~~~~~~~~~~~~~~~~~~~~~~~~~~~
kcDNA <- read_tsv('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/align_data/Mock_Flg22_1_ath_kallisto/abundance.tsv')

## kcDNA 48359
## athAnno 54013

## all kcDNA in athAnno
sum(kcDNA$target_id %in% athAnno$cDNA)

## lncRNA  miRNA  ncRNA   rRNA snoRNA  snRNA   tRNA
##   3879    325    377     15    287     82    689
anti_join(athAnno, kcDNA, by = c('cDNA' = 'target_id')) %>%
  .$BioType %>%
  table
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##############################################################
