
## originally by Yulong Niu
## yulong.niu@hotmail.com

#######################GO analysis############################
setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/geneset')

library('goseq')
library('GO.db')
library('foreach')
library('doMC')
library('KEGGAPI')
library('BioCycAPI')
library('magrittr')
library('dplyr')
library('tibble')
library('readr')
library('stringr')

registerDoMC(12)

load('athGO.RData')
load('athKEGG.RData')
load('athBioCyc.RData')

load('/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results_orthologs/rbhAnno.RData')
rbhAnno %<>%
  filter(ID %>% str_detect('^AT'))

anno <- read_csv('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/Ensembl_ath_Anno.csv',
                 col_types = cols(Chromosome = col_character())) %>%
  mutate(Gene = Gene %>% {if_else(is.na(.), '', .)}) %>%
  mutate(Description = Description %>% {if_else(is.na(.), '', .)}) %>%
  dplyr::select(ID, Length)

kmeansBkg <- read_csv('/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results_orthologs/kmeans_12_RBH_rmfull.csv',
                      col_types = cols(Chromosome = col_character())) %>%
  inner_join(rbhAnno, c('ID' = 'RBH')) %>%
  dplyr::rename(RBH = ID, ID = ID.y) %>%
  inner_join(anno)

savepath <- '/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results_orthologs/geneset/fullbc'

##~~~~~~~~~~~~~~~select genesets~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
athGO %<>%
  lapply(function(x) {x[x %in% kmeansBkg$ID]}) %>%
  {
    l <- sapply(., length) > 0
    .[l]
  }
GOMat <- foreach(i = 1:length(athGO), .combine = rbind) %dopar% {
  eachMat <- cbind(athGO[[i]], names(athGO)[i])
  return(eachMat)
} %>% as.data.frame

athKEGG %<>%
  lapply(function(x) {x[x %in% kmeansBkg$ID]}) %>%
  {
    l <- sapply(., length) > 0
    .[l]
  }
KEGGMat <- foreach(i = 1:length(athKEGG), .combine = rbind) %dopar% {
  eachMat <- cbind(athKEGG[[i]], names(athKEGG)[i])
  return(eachMat)
} %>% as.data.frame


athBioCyc %<>%
  lapply(function(x) {x[x %in% kmeansBkg$ID]}) %>%
  {
    l <- sapply(., length) > 0
    .[l]
  }
BioCycMat <- foreach(i = 1:length(athBioCyc), .combine = rbind) %dopar% {
  eachMat <- cbind(athBioCyc[[i]], names(athBioCyc)[i])
  return(eachMat)
} %>% as.data.frame
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~whole cluster gene-set~~~~~~~~~~~~~~~~~
for (i in kmeansBkg$cl %>% unique) {

  prefix <- 'kmeans10_1stadd'

  degVec <- (kmeansBkg$cl == i) %>%
    as.integer %>%
    set_names(kmeansBkg$ID)

  pwf <- nullp(degVec, bias.data = kmeansBkg$Length)

  ## GO
  GOTestWithCat <- goseq(pwf, gene2cat = GOMat, use_genes_without_cat = FALSE) %>%
    as_tibble %>%
    filter(!is.na(ontology))

  termCat <- c('BP', 'MF', 'CC')
  for (j in termCat) {
    write.csv(GOTestWithCat %>% filter(ontology == j),
              paste0(prefix, '_cluster', i, '_', j, '.csv') %>% file.path(savepath, .))
  }

  ## KEGG
  pathAnno <- getKEGGPathAnno('ath') %>%
    as_tibble %>%
    mutate(Annotation = Annotation %>% substr(., 1, nchar(.) - 37))

  KEGGTestWithCat <- goseq(pwf, gene2cat = KEGGMat, use_genes_without_cat = FALSE) %>%
    as_tibble %>%
    inner_join(., pathAnno, by = c('category' = 'pathID')) %>%
    mutate(ontology = 'KEGG')

  write.csv(KEGGTestWithCat,
            paste0(prefix, '_cluster', i, '_KEGG.csv') %>% file.path(savepath, .))

  ## BioCyc
  pathAnno <- getCycPathway('ARA') %>%
    rename(Annotation = pathAnno) %>%
    mutate(Annotation = Annotation %>% str_replace_all('<.*?>', ''))

  BioCycTestWithCat <- goseq(pwf, gene2cat = BioCycMat, use_genes_without_cat = FALSE) %>%
    as_tibble %>%
    inner_join(., pathAnno, by = c('category' = 'pathID')) %>%
    mutate(ontology = 'BioCyc')

  write.csv(BioCycTestWithCat,
            paste0(prefix, '_cluster', i, '_BioCyc.csv') %>% file.path(savepath, .))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
