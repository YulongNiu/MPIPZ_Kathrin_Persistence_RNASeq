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
library('tidyverse')
library('foreach')
library('ParaMisc')
library('GUniFrac')

anno <- read_csv('/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results/lotus_gifu_collaborator_v1p2_Anno.csv',
                 col_types = cols(Chromosome = col_character())) %>%
  mutate_all(list(~replace(., is.na(.), '')))

##~~~~~~~~~~~~~~~~~~~~load k alignments~~~~~~~~~~~~~~~~~~~~~~~~~~
wd <- '/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/align_data/lotus_collaborator'
setwd(wd)

labelanno <- read_delim('../../results/mapping.txt', delim = '\t') %>%
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
setwd('/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results_rmfull/')

condi <- c('fullSC', 'AtSC', 'AtSCMloti', 'LjSC', 'Mock')

sampleTable <- data.frame(condition = factor(rep(condi, each = 4), levels = condi))
rownames(sampleTable) <- colnames(kres$counts)

sampleTable$condition %<>% relevel(ref = 'Mock')

degres <- DESeqDataSetFromTximport(kres, sampleTable, ~condition)

## remove 0|0|x|x, 0|0|0|x, 0|0|0|0
degres %<>%
  estimateSizeFactors %>%
  counts(normalized = TRUE) %>%
  apply(1, checkPersis, 1) %>%
  degres[., ]

degres %<>% DESeq

## count transformation
rld <- rlog(degres)
ntd <- normTransform(degres)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~subsample~~~~~~~~~~~~~~~~~~~~~~~~~~
tmp1 <- kres$counts %>%
  t %>%
  Rarefy(., 3000000) %>%
  .$otu.tab.rff
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~hidden batch effect~~~~~~~~~~~~~~~~~~~~~
library('sva')
library('ggplot2')

dat <- rld %>%
  assay %>%
  {.[rowMeans(.) > 1, ]}
mod <- model.matrix(~ condition, colData(degres))
mod0 <- model.matrix(~ 1, colData(degres))

## manual detect surrogate variance
svnum <- 4
svseq <- svaseq(dat, mod, mod0, n.sv = svnum)

## auto detect sv
svobj <- sva(dat, mod, mod0)
svnum <- svobj$sv %>% ncol

svobj$sv %>%
  set_colnames(paste0('sv', seq_len(svnum))) %>%
  as_tibble %>%
  gather(key = 'sv', value = 'value') %>%
  mutate(condition = colData(degres) %>%
           .$condition %>%
           rep(svnum) %>%
           as.character,
         sample = rep(colnames(degres), svnum)) %>%
  mutate(group = paste(sv, condition, sep = '_')) %>%
  ggplot(aes(sample, value, colour = sv, group = group)) +
  geom_point() +
  geom_line() +
  theme(axis.text.x = element_text(angle = 90))
ggsave('auto_lotus_sv.jpg')
ggsave('auto_lotus_sv.pdf')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~DEGs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
degres$sv1 <- svobj$sv[, 1]
degres$sv2 <- svobj$sv[, 2]
degres$sv3 <- svobj$sv[, 3]
design(degres) <- ~sv1 + sv2 + sv3 + condition

degres <- DESeq(degres)

cond <- list(c('fullSC', 'Mock'),
             c('AtSC', 'Mock'),
             c('AtSCMloti', 'Mock'),
             c('LjSC', 'Mock'),
             c('fullSC', 'AtSC'),
             c('fullSC', 'AtSCMloti'),
             c('fullSC', 'LjSC'),
             c('AtSC', 'AtSCMloti'),
             c('AtSC', 'LjSC'))

resRaw <- lapply(cond,
                 function(x) {
                   degres %>%
                     results(contrast = c('condition', x)) %T>%
                     summary %>%
                     as_tibble %>%
                     select(pvalue, padj, log2FoldChange) %>%
                     rename_all(.funs = list(~paste0(paste(x, collapse = '_vs_'), '_', .)))
                 }) %>%
  bind_cols

res <- cbind.data.frame(as.matrix(mcols(degres)[, 1:10]), assay(rld), stringsAsFactors = FALSE) %>%
  rownames_to_column(., var = 'ID') %>%
  as_tibble %>%
  bind_cols(resRaw) %>%
  inner_join(anno, by = 'ID') %>%
  select(ID, Gene : Description, L_fSC_1 : AtSC_vs_LjSC_log2FoldChange) %>%
  arrange(fullSC_vs_Mock_padj)

write_csv(res, 'SynCom_vs_Mock_lotus_sva_k.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~heatmap~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library('pheatmap')

pheatmap(assay(ntd),
         show_rownames=FALSE)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PCA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library('ggrepel')
library('ggplot2')
library('RColorBrewer')
library('limma')
library('sva')

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- rld %>%
  assay %>%
  {.[rowMeans(.) > 1, ]}

group <- sampleTable$condition
design <- model.matrix(~ group)
rldData <- dat %>%
  removeBatchEffect(covariates = svobj$sv,
                    design = design)

## ## batch correction limma - lotus
## cutMat <- CutSeqEqu(ncol(rld), 4)
## for (i in seq_len(ncol(cutMat))) {
##   eachCols <- cutMat[1, i] : cutMat[2, i]
##   rldData[, eachCols] %<>% removeBatchEffect(c(1, 1, 2, 2) %>% factor)
## }
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cols <-  brewer.pal(5, name = 'Dark2')

## full
sampleIdx <- 1:20
colorIdx <- 1:5

## ## without AtSC
## sampleIdx <- (1:20)[-5:-8]
## colorIdx <- (1:5)[-3]

## without fullSC
sampleIdx <- (1:20)[-1:-4]
colorIdx <- (1:5)[c(1, 2, 4, 3)]

## without fullSC and AtSCMloti
sampleIdx <- (1:20)[c(-1:-4, -9:-12)]
colorIdx <- (1:5)[c(1, 2, 3)]

## 1 - 2 C
pca <- prcomp(t(rldData[, sampleIdx]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
percentVar <- round(100 * percentVar)
pca1 <- pca$x[,1]
pca2 <- pca$x[,2]
pcaData <- tibble(PC1 = pca1, PC2 = pca2, Group = colData(rld)[sampleIdx, 1], ID = rownames(colData(rld))[sampleIdx]) %>%
  mutate(SynCom = rep(c('AtSC', 'LjSC', 'Mock+LjNodule218'), each = 4) %>% factor(levels = c('Mock+LjNodule218', 'AtSC', 'LjSC')))

ggplot(pcaData, aes(x = PC1, y = PC2, colour = SynCom)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_colour_manual(name = 'SynCom',
                      values = cols[colorIdx]) +
  stat_ellipse(aes(x = PC1, y = PC2, group = SynCom), type = 't', linetype = 2, level = 0.8) +
  coord_fixed(1) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
        legend.text.align = 0,
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 14),
        legend.text=element_text(size= 13),
        legend.title = element_text(size = 14))

ggsave('PCA_lotus_likeath_sva.pdf', width = 10)
ggsave('PCA_lotus_likeath_sva.jpg', width = 10)

save(degres, rldData, file = 'degres_condi_Mock_lotus.RData')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######################################################################
