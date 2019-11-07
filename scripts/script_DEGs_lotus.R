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
library('ParaMisc')
library('GUniFrac')

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
  degres[., ]

save(degres, file = '../results/degres_condi_Mock_lotus.RData')

degres %<>% DESeq

## count transformation
rld <- rlog(degres)
ntd <- normTransform(degres)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~subsample~~~~~~~~~~~~~~~~~~~~~~~~~~
tmp1 <- kres$counts %>%
  t %>%
  Rarefy(., 3000000) %>%
  .$otu.tab.rff
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~DEGs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cond <- degres %>%
  resultsNames %>%
  str_extract('(?<=condition_).*') %>%
  .[!is.na(.)]

resRaw <- lapply(cond,
                 function(x) {
                   degres %>%
                     results(name = paste0('condition_', x)) %T>%
                     summary %>%
                     as_tibble %>%
                     select(pvalue, padj, log2FoldChange) %>%
                     rename_all(.funs = list(~paste0(x, '_', .)))
                 }) %>%
  bind_cols

res <- cbind.data.frame(as.matrix(mcols(degres)[, 1:10]), assay(ntd), stringsAsFactors = FALSE) %>%
  rownames_to_column(., var = 'ID') %>%
  as_tibble %>%
  bind_cols(resRaw) %>%
  inner_join(anno, by = 'ID') %>%
  select(ID, Gene : Description, C_fSC_1 : LjSC_vs_Mock_log2FoldChange) %>%
  arrange(fullSC_vs_Mock_padj)

write_csv(res, '../results/SynCom_vs_Mock_ath_k.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~heatmap~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library('pheatmap')

pheatmap(assay(ntd),
         show_rownames=FALSE)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PCA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
## batch correction limma - lotus
cutMat <- CutSeqEqu(ncol(rld), 4)
for (i in seq_len(ncol(cutMat))) {
  eachCols <- cutMat[1, i] : cutMat[2, i]
  rldData[, eachCols] %<>% removeBatchEffect(c(1, 1, 2, 2) %>% factor)
}

## ## batch correction sva
## modcombat <- model.matrix(~1, data = sampleTable)
## rldData %<>% ComBat(dat = ., batch = rep(rep(1 : 4, 5) %>% factor) %>% factor, mod = modcombat, par.prior = TRUE, prior.plots = FALSE)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cols <- c('#E41A1C', '#377EB8', '#4DAF4A', '#FF7F00', '#984EA3')

## full
sampleIdx <- 1:20
colorIdx <- 1:5

## without AtSC
sampleIdx <- (1:20)[-5:-8]
colorIdx <- (1:5)[-3]

## without AtSCMloti
sampleIdx <- (1:20)[-9:-12]
colorIdx <- (1:5)[-4]

## 1 - 2 C
pca <- prcomp(t(rldData[, sampleIdx]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
percentVar <- round(100 * percentVar)
pca1 <- pca$x[,1]
pca2 <- pca$x[,2]
pcaData <- data.frame(PC1 = pca1, PC2 = pca2, Group = colData(rld)[sampleIdx, 1], ID = rownames(colData(rld))[sampleIdx])
ggplot(pcaData, aes(x = PC1, y = PC2, colour = Group)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_dl(aes(label = ID, color = Group), method = 'smart.grid') +
  scale_colour_manual(values = cols[colorIdx])
ggsave('../results/PCA_ath_limma_noAtSCMloti.pdf', width = 15, height = 12)
ggsave('../results/PCA_ath_limma_noAtSCMloti.jpg', width = 15, height = 12)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

######################################################################
