########################compare Ka-Wai and Kathrin######################
library('tidyverse')
library('magrittr')
library('foreach')

setwd('/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results_rmfull/')
load('degres_condi_Mock_ath.RData')

##~~~~~~~~~~~~~~~~~~~~~~~Kathrin DEGs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
wholeDEG <- read_csv('SynCom_vs_Mock_ath_sva_k.csv')
kmeansRes <- read_csv('kmeans_10_ath.csv') %>%
  select(ID, cl)

fcsig <- wholeDEG %>%
  select(ends_with('FoldChange')) %>%
  transmute_all(list(~ case_when(. > log2(1.5) ~ 1,
                                 . < -log2(1.5) ~ 1,
                                 TRUE ~ 0))) %>%
  select(-contains('fullSC'))

padjsig <- wholeDEG %>%
  select(ends_with('padj')) %>%
  abs %>%
  `<`(0.05) %>%
  as_tibble %>%
  transmute_all(list(~ if_else(is.na(.), FALSE, .))) %>%
  select(-contains('fullSC'))

heatsig <- (padjsig * fcsig) %>%
  as_tibble %>%
  rowSums %>%
  {. >= 1} %>%
  which %>%
  dplyr::slice(wholeDEG, .) %>%
  inner_join(kmeansRes)

rawC <- rldData %>%
  as.data.frame %>%
  .[, c(5:16)] %>%
  rownames_to_column('ID') %>%
  as_tibble %>%
  inner_join(heatsig %>% select(ID, cl))

scaleC <- rawC %>%
  select(contains('C_')) %>%
  t %>%
  scale %>%
  t %>%
  as_tibble %>%
  bind_cols(rawC %>% select(ID, cl))

KathrinDEG <- scaleC %>%
  select(ID, cl)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Ka-Wai DEGs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
KaWaiDEG <- read_csv('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/removeZero/kmeans10_soil_sig.csv') %>%
  select(ID, cl)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~compare DEGs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
KathrinKaWaiJS <- foreach (i = 1:10, .combine = rbind) %do% {

  KathrinEach <- KathrinDEG %>%
    filter(cl == i) %>%
    .$ID %>%
    substr(start = 1, stop = nchar(.) - 2) %>%
    unique

  JacSim <- foreach(j = 1:10, .combine = c)  %do% {
    KaWaiEach <- KaWaiDEG %>%
      filter(cl == j) %>%
      .$ID %>%
      unique

    eachJacSim <- length(intersect(KathrinEach, KaWaiEach)) / length(union(KathrinEach, KaWaiEach))

    return(eachJacSim)

  } %>%
  set_names(paste0('Kathrin', 1:10))

  return(JacSim)
} %>%
  set_rownames(paste0('Ka-Wai', 1:10)) %>%
  round(digits = 3)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~compare all genes~~~~~~~~~~~~~~~~~~~~~~~~~~~~
KathrinAll <- read_csv('kmeans_10_ath.csv') %>%
  select(ID, cl)

KaWaiAll <- read_csv('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/removeZero/kmeans10_soil.csv') %>%
  select(ID, cl)

KathrinKaWaiJS <- foreach (i = 1:10, .combine = rbind) %do% {

  KathrinEach <- KathrinAll %>%
    filter(cl == i) %>%
    .$ID %>%
    substr(start = 1, stop = nchar(.) - 2) %>%
    unique

  JacSim <- foreach(j = 1:10, .combine = c)  %do% {
    KaWaiEach <- KaWaiAll %>%
      filter(cl == j) %>%
      .$ID %>%
      unique

    eachJacSim <- length(intersect(KathrinEach, KaWaiEach)) / length(union(KathrinEach, KaWaiEach))

    return(eachJacSim)

  } %>%
  set_names(paste0('Kathrin', 1:10))

  return(JacSim)
} %>%
  set_rownames(paste0('Ka-Wai', 1:10)) %>%
  round(digits = 3)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########################################################################


###########################At GO --> Lotus##############################
## AtSC LjSC Mock
## all gene clusters:
## c2: up-regulated in both AtSC and LjSC
## c4: down-regulated in both AtSC and LjSC, LjSC more
## c3 and c5: up-regulated in only AtSC

meanLotus <- function(v) {

  require('magrittr')

  res <- v %>%
    split(rep(1 : 4, each = 4)) %>%
    sapply(mean, na.rm = TRUE)

  return(res)
}

library('tidyverse')
library('magrittr')
library('ggplot2')

setwd('/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results_rmfull/')

AtCluster <- read_csv('kmeans_10_ath.csv') %>%
  select(ID, cl) %>%
  mutate(TAIRGene =  ID %>% strsplit(split = '.', fixed = TRUE) %>% sapply('[[', 1)) %>%
  rename(TAIR = ID)

LjCluster <- read_csv('kmeans_10_lotus_collaborator.csv') %>%
  dplyr::select(ID, cl)

##~~~~~~~~~~~~~~~~~~~~~~~~load Lotus~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load('degres_condi_Mock_lotus.RData')

rawC <- rldData %>%
  as.data.frame %>%
  .[, c(-1 : -4)] %>%
  rownames_to_column('ID') %>%
  as_tibble

LjScale <- rawC %>%
  select(-ID) %>%
  t %>%
  scale %>%
  t %>%
  as_tibble %>%
  bind_cols(rawC %>% select(ID))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lotusAnno <- read_tsv('/extDisk1/Biotools/RefData/lotus_gifu_collaborator_v1.2/LjGifu_1.2_functional_annotation.txt', comment = '#') %>%
  dplyr::select(`Protein-Accession`, `Best BlastHit against 'tair'`, `Human-Readable-Description`) %>%
  set_colnames(c('GID', 'TAIR', 'GENENAME')) %>%
  dplyr::mutate(TAIR = TAIR %>% strsplit(split = ' ', fixed = TRUE) %>% sapply('[[', 1)) %>%
  dplyr::mutate(TAIRGene = TAIR %>% strsplit(split = '.', fixed = TRUE) %>% sapply('[[', 1))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~select cluster~~~~~~~~~~~~~~~~~~~~~~~~
clNum <- c(3, 5)

lotusBH <- AtCluster %>%
  filter(cl %in% clNum) %>%
  select(TAIRGene) %>%
  distinct %>%
  inner_join(lotusAnno) %>%
  select(GID, GENENAME) %>%
  distinct %>%
  mutate(GGene = GID %>% strsplit(split = '.', fixed = TRUE) %>% sapply('[[', 1))

## library('clusterProfiler')
## library('org.Ljaponicus.eg.db')

## enrichGO(gene = lotusBH$GGene,
##          OrgDb = 'org.Ljaponicus.eg.db',
##          keyType= 'GID',
##          ont = 'BP',
##          universe = keys(org.Ljaponicus.eg.db),
##          pAdjustMethod = 'BH',
##          pvalueCutoff=0.05,
##          qvalueCutoff=0.1) %>%
##   as.data.frame %>%
##   .[, 2]

lotusBH %>%
  select(GID) %>%
  inner_join(LjScale, c('GID' = 'ID')) %>%
  select(-GID) %>%
  apply(1, meanLotus) %>%
  t %>%
  as_tibble %>%
  set_colnames(c('L_AtSC', 'L_AtSCMloti', 'L_LjSC', 'L_mock')) %>%
  gather(key = 'Sample', value = 'ScaleC') %>%
  ggplot(aes(x = Sample, y = ScaleC)) +
  geom_boxplot()

ggsave('At_cl35_lotus_scale.jpg')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~defense At --> Lj~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library('ComplexHeatmap')
library('RColorBrewer')

GOAtall <- read_csv('kmeans10_ath_cp_BP_allgene.csv')

## At defense related GO BP terms
defenseAtGO <- c('GO:0010200', 'GO:0034050', 'GO:0002682',
                 'GO:0034050', 'GO:0009626', 'GO:0045088',
                 'GO:0050776', 'GO:0012501', 'GO:0010363',
                 'GO:0043067', 'GO:0080135', 'GO:0060548',
                 'GO:0002679', 'GO:0045730', 'GO:0009611')

defenseAtGenes <- GOAtall %>%
  dplyr::filter(ID %in% defenseAtGO) %>%
  .$geneID %>%
  strsplit(split = '/', fixed = TRUE) %>%
  unlist %>%
  unique %>%
  bind_cols %>%
  set_colnames(c('TAIRGene'))


## Lj cluster number
defenseLj <- inner_join(lotusAnno, defenseAtGenes) %>%
  inner_join(inner_join(LjCluster, LjScale), by = c('GID' = 'ID')) %>%
  select(-TAIR, -GENENAME, -TAIRGene)


cairo_pdf('defense_At2Lj_heatmap.pdf')
Heatmap(matrix = defenseLj %>% select(contains('_')),
        name = 'Scaled Counts',
        ## row_order = order(scaleC$cl) %>% rev,
        row_split = defenseLj$cl,
        row_gap = unit(2, "mm"),
        column_order = 1 : 16,
        column_split = rep(c('AtSC', 'AtSCMloti', 'LjSC', 'Mock'), each = 4),
        show_column_names = FALSE,
        col = colorRampPalette(rev(brewer.pal(n = 10, name = 'Spectral'))[c(-3, -4, -7, -8)])(10))
dev.off()

## Lj with At cluster number
defenseLj <- inner_join(lotusAnno, defenseAtGenes) %>%
  inner_join(inner_join(AtCluster, defenseAtGenes)) %>%
  inner_join(LjScale, by = c('GID' = 'ID')) %>%
  select(-TAIR, -GENENAME, -TAIRGene)

cairo_pdf('defense_At2Lj_AtclusterID_heatmap.pdf')
Heatmap(matrix = defenseLj %>% select(contains('_')),
        name = 'Scaled Counts',
        ## row_order = order(scaleC$cl) %>% rev,
        row_split = defenseLj$cl,
        row_gap = unit(2, "mm"),
        column_order = 1 : 16,
        column_split = rep(c('AtSC', 'AtSCMloti', 'LjSC', 'Mock'), each = 4),
        show_column_names = FALSE,
        col = colorRampPalette(rev(brewer.pal(n = 10, name = 'Spectral'))[c(-3, -4, -7, -8)])(10))
dev.off()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########################################################################


###########################select interesting gene#######################
library('tidyverse')
library('DESeq2')

setwd('/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results_rmfull/')

##~~~~~~~~~~~~~~~~~~~~~~~~~select ath~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load('degres_condi_Mock_ath.RData')

kmeansRes <- read_csv('kmeans_10_ath.csv') %>%
  mutate_at(c('Gene', 'Description'), list(~replace(., is.na(.), ''))) %>%
  select(ID, cl)

## rlog transformed
rawC <- rldData %>%
  as.data.frame %>%
  .[, c(5:16)] %>%
  rownames_to_column('ID') %>%
  as_tibble %>%
  inner_join(kmeansRes %>% select(ID, cl))

scaleCAt <- rawC %>%
  select(contains('C_')) %>%
  t %>%
  scale %>%
  t %>%
  as_tibble %>%
  bind_cols(rawC %>% select(ID, cl))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~select lotus~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load('degres_condi_Mock_lotus.RData')

kmeansRes <- read_csv('kmeans_10_lotus_rmfull_rmAtSC.csv') %>%
  mutate_at(c('Gene', 'Description'), list(~replace(., is.na(.), ''))) %>%
  select(ID, cl)

## rlog transformed
rawC <- rldData %>%
  as.data.frame %>%
  .[, c(-1:-8)] %>%
  rownames_to_column('ID') %>%
  as_tibble %>%
  inner_join(kmeansRes %>% select(ID, cl))

scaleCLj3 <- rawC %>%
  select(contains('L_')) %>%
  t %>%
  scale %>%
  t %>%
  as_tibble %>%
  bind_cols(rawC %>% select(ID, cl))

kmeansRes <- read_csv('kmeans_10_lotus_collaborator.csv') %>%
  mutate_at(c('Gene', 'Description'), list(~replace(., is.na(.), ''))) %>%
  select(ID, cl)

## rlog transformed
rawC <- rldData %>%
  as.data.frame %>%
  .[, c(-1:-4)] %>%
  rownames_to_column('ID') %>%
  as_tibble %>%
  inner_join(kmeansRes %>% select(ID, cl))

scaleCLj4 <- rawC %>%
  select(contains('L_')) %>%
  t %>%
  scale %>%
  t %>%
  as_tibble %>%
  bind_cols(rawC %>% select(ID, cl))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~select genes~~~~~~~~~~~~~~~~~~~~~~~~~~~
## select genes from At
AtG <- read_csv('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/Ensembl_ath_Anno.csv',
                col_types = cols(Chromosome = col_character())) %>%
  mutate_at(c('Gene', 'Description'), list(~replace(., is.na(.), ''))) %>% {
    selectIdx <- (str_detect(.$Gene, regex('MYB|WRKY|(^LYK\\d+$)|(^ERF\\d+$)', ignore_case = TRUE)) |
                  str_detect(.$Description, regex('((^| |//)MYB)|((^| |//)WRKY)|((^| |//)LYK\\d+)|((^| |//)ERF\\d+)', ignore_case = TRUE))) %>%
      which

    dplyr::slice(., selectIdx)
  } %>%
  select(ID, Gene, Description) %>%
  dplyr::rename(IDAt = ID, GeneAt = Gene, DescriptionAt = Description)

## select genes from At --> Lj
LjG1 <- read_csv('/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results/lotus_gifu_collaborator_v1p2_Anno.csv',
                 col_types = cols(Chromosome = col_character())) %>%
  mutate_at(c('Best_TAIR', 'Best_sport', 'Best_trembl_plants', 'Description'), list(~replace(., is.na(.), ''))) %>%
  dplyr::filter(nchar(Best_TAIR) > 0) %>%
  select(ID, Best_TAIR : Description) %>%
  mutate(Best_TAIR = Best_TAIR %>% strsplit(' ', fixed = TRUE) %>% sapply('[[', 1)) %>%
  dplyr::rename(IDLj = ID, DescriptionLj = Description) %>%
  inner_join(AtG, c('Best_TAIR' = 'IDAt')) %>%
  dplyr::select(IDLj : DescriptionLj)

## select genes from Lj anno
LjG2 <- read_csv('/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results/lotus_gifu_collaborator_v1p2_Anno.csv',
                 col_types = cols(Chromosome = col_character())) %>%
  mutate_at(c('Best_TAIR', 'Best_sport', 'Best_trembl_plants', 'Description'), list(~replace(., is.na(.), ''))) %>%
  dplyr::filter(nchar(Best_TAIR) > 0) %>%
  {
    selectIdx <- (str_detect(.$Best_sport, regex('((^| |//)MYB)|((^| |//)WRKY)|((^| |//)LYK\\d+)|((^| |//)ERF\\d+)', ignore_case = TRUE)) |
                  str_detect(.$Best_trembl_plants, regex('((^| |//)MYB)|((^| |//)WRKY)|((^| |//)LYK\\d+)|((^| |//)ERF\\d+)', ignore_case = TRUE)) |
      str_detect(.$Description, regex('((^| |//)MYB)|((^| |//)WRKY)|((^| |//)LYK\\d+)|((^| |//)ERF\\d+)', ignore_case = TRUE))) %>%
      which

    dplyr::slice(., selectIdx)
  } %>%
  dplyr::filter(!is.na(Best_TAIR)) %>%
  select(ID, Best_TAIR : Description) %>%
  mutate(Best_TAIR = Best_TAIR %>% strsplit(' ', fixed = TRUE) %>% sapply('[[', 1)) %>%
  dplyr::rename(IDLj = ID, DescriptionLj = Description)

LjG <- bind_rows(LjG1, LjG2) %>%
  distinct

AtLjG <- inner_join(AtG, LjG, by = c('IDAt' = 'Best_TAIR'))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#########################################################################
