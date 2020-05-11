
## originally by Yulong Niu
## yulong.niu@hotmail.com

library('tximport')
library('rhdf5')
library('magrittr')
library('DESeq2')
library('tidyverse')
library('foreach')
library('ParaMisc')

setwd('/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results_orthologs')

load('kresRBH.RData')

load('/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results_orthologs/rbhAnno.RData')

annoLotus <- read_csv('/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results/lotus_gifu_collaborator_v1p2_Anno.csv',
                      col_types = cols(Chromosome = col_character())) %>%
  mutate_all(list(~replace(., is.na(.), ''))) %>%
  inner_join(rbhAnno) %>%
  setNames(paste0(names(.), '_Lj'))

##~~~~~~~~~~~~~~~~~~~normalization~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sampleTable <- colData(kresLotus)[9:20, , drop = FALSE]
sampleTable$condition %<>% droplevels

degres <- kresLotus %>%
  assay %>%
  .[, 9:20, drop = FALSE] %>%
  DESeqDataSetFromMatrix(., sampleTable, ~ condition)

degres %<>% DESeq

## count transformation
rld <- rlog(degres)
ntd <- normTransform(degres)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~hidden batch effect~~~~~~~~~~~~~~~~~~~~~
library('sva')
library('ggplot2')

dat <- rld %>%
  assay %>%
  {.[rowMeans(.) > 1, ]}
mod <- model.matrix(~condition, colData(degres))
mod0 <- model.matrix(~1, colData(degres))

## manual detect surrogate variance
svnum <- 4
svseq <- svaseq(dat, mod, mod0, n.sv = svnum)

## ## auto detect sv
## svobj <- sva(dat, mod, mod0)
## svnum <- svobj$sv %>% ncol

svseq$sv %>%
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
ggsave('auto_ath_sv_3.jpg')
ggsave('auto_ath_sv_3.pdf')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~DEGs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
degres$sv1 <- svseq$sv[, 1]
degres$sv2 <- svseq$sv[, 2]
degres$sv3 <- svseq$sv[, 3]
degres$sv4 <- svseq$sv[, 4]
design(degres) <- ~sv1 + sv2 + sv3 + sv4 + condition

degres <- DESeq(degres)

cond <- list(c('AtSCMloti_Lj', 'Mock_Lj'),
             c('LjSC_Lj', 'Mock_Lj'),
             c('AtSCMloti_Lj', 'LjSC_Lj'))

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
  inner_join(annoLotus, by = c('ID' = 'RBH_Lj')) %>%
  select(ID, ID_Lj : Description_Lj, L_AtSCMloti_1 : AtSCMloti_Lj_vs_LjSC_Lj_log2FoldChange) %>%
  arrange(AtSCMloti_Lj_vs_Mock_Lj_padj)

write_csv(res, 'SynCom_vs_Mock_RBH_rmfull_lotus_sva_k.csv')
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
  removeBatchEffect(covariates = svseq$sv,
                    design = design)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cols <- brewer.pal(3, name = 'Dark2')

sampleIdx <- 1:12
colorIdx <- 1:3

## 1 - 2 C
pca <- prcomp(t(rldData[, sampleIdx]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
percentVar <- round(100 * percentVar)
pca1 <- pca$x[,1]
pca2 <- pca$x[,2]
pcaData <- data.frame(PC1 = pca1, PC2 = pca2, Group = colData(rld)[sampleIdx, 1], ID = rownames(colData(rld))[sampleIdx])
ggplot(pcaData, aes(x = PC1, y = PC2, colour = Group, label = ID)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_colour_manual(name = 'SynCom',
                      values = cols[colorIdx]) +
  stat_ellipse(aes(x = PC1, y = PC2, group = Group), type = 't', linetype = 2, level = 0.8) +
  coord_fixed(1) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
        legend.text.align = 0,
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 14),
        legend.text=element_text(size= 13),
        legend.title = element_text(size = 14))

ggsave('PCA_RBH_rmfull_lotus_sva.pdf', width = 10)
ggsave('PCA_RBH_rmfull_lotus_sva.jpg', width = 10)

save(degres, rldData, file = 'degres_condi_RBH_rmfull_lotus.RData')
########################################################################
