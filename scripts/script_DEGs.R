###########################DEGs##################################
##~~~~~~~~~~~~~~~~~~~~~~~~~~~useful funcs~~~~~~~~~~~~~~~~~~~~~~~~
checkZeros <- function(v, threshold) {
  res <- ifelse(sum(v == 0) > threshold, FALSE, TRUE)
  return(res)
}

checkPersis <- function(v, threshold) {

  require('magrittr')

  res <- v %>%
    split(rep(1 : 5, each = 4)) %>%
    sapply(checkZeros, threshold) %>%
    all

  return(res)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


library('tximport')
library('rhdf5')
library('magrittr')
library('DESeq2')
library('tibble')
library('readr')
library('dplyr')
library('stringr')
library('foreach')

anno <- read_csv('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/Ensembl_ath_Anno.csv',
                 col_types = cols(Chromosome = col_character())) %>%
  mutate(Gene = Gene %>% {if_else(is.na(.), '', .)}) %>%
  mutate(Description = Description %>% {if_else(is.na(.), '', .)})


##~~~~~~~~~~~~~~~~~~~~load k alignments~~~~~~~~~~~~~~~~~~~~~~~~~~

wd <- '/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/align_data'
setwd(wd)

labelanno <- read_delim('../results/mapping.txt', delim = '\t') %>%
  dplyr::rename(ID = `library_number`, SampleAnno = `library_name`) %>%
  mutate(ID = ID %>% str_replace('\\.', '_')) %>%
  filter(species %>% str_detect('Lj'))

slabel <- labelanno$SampleAnno %>%
  paste0('_lotus_kallisto')

kres <- file.path(wd, slabel, 'abundance.h5') %>%
  set_names(labelanno$SampleAnno) %>%
  tximport(type = 'kallisto', txOut = TRUE)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~normalization~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
condi <- c('fullSC', 'AtSC', 'LjSC', 'Mock')

condi <- c('fullSC', 'AtSC', 'AtSCMloti', 'LjSC', 'Mock')

sampleTable <- data.frame(condition = factor(rep(condi, each = 4), levels = condi))
rownames(sampleTable) <- colnames(kres$counts)

sampleTable$condition %<>% relevel(ref = 'Mock')

degres <- DESeqDataSetFromTximport(kres, sampleTable, ~condition)

## remove 0|0|0|x and |0|0|0|0
degres %<>%
  estimateSizeFactors %>%
  counts(normalized = TRUE) %>%
  apply(1, checkPersis, 1) %>%
  degres[., ] %>%
  DESeq

## count transformation
rld <- rlog(degres)
ntd <- normTransform(degres)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PCA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library('directlabels')
library('ggplot2')
library('RColorBrewer')
library('limma')
library('sva')

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## remove low count
thres <- 0
rldData <- assay(rld)

rl <- apply(rldData, 1, function(x){
  return(sum(x > thres) == length(x))
})
rldData %<>% .[rl, ]

## batch correction limma
## rldData %<>% removeBatchEffect(rep(1 : 4, 4) %>% factor)

## ## batch correction sva
## modcombat <- model.matrix(~1, data = sampleTable)
## rldData %<>% ComBat(dat = ., batch = rep(rep(1 : 4, 10) %>% factor) %>% factor, mod = modcombat, par.prior = TRUE, prior.plots = FALSE)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cols <- colData(rld)[, 1] %>% factor(., labels = brewer.pal(5, name = 'Set1'))

## 1 - 2 C
pca <- prcomp(t(rldData))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
percentVar <- round(100 * percentVar)
pca1 <- pca$x[,1]
pca2 <- pca$x[,2]
pcaData <- data.frame(PC1 = pca1, PC2 = pca2, Group = colData(rld)[, 1], ID = rownames(colData(rld)))
ggplot(pcaData, aes(x = PC1, y = PC2, colour = Group)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_dl(aes(label = ID, color = Group), method = 'smart.grid') +
  scale_colour_manual(values = levels(cols))
ggsave('PCA_1stadd_sva.pdf', width = 15, height = 12)
ggsave('PCA_1stadd_sva.jpg', width = 15, height = 12)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

######################################################################
