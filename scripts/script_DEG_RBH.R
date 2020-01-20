################################normalization#########################
library('tximport')
library('rhdf5')
library('magrittr')
library('DESeq2')
library('tidyverse')
library('foreach')
library('ParaMisc')

setwd('/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results_orthologs')

##~~~~~~~~~~~~~~~~~~~~~~~~~~try normal within independent~~~~~~~~~~~~~~
load('degres_condi_RBH_rmfull_ath.RData')
rldDataAth <- rldData %>%
  as.data.frame %>%
  rownames_to_column('ID')

load('degres_condi_RBH_rmfull_lotus.RData')
rldDataLotus <- rldData %>%
  as.data.frame %>%
  rownames_to_column('ID')

rldData <- inner_join(rldDataAth, rldDataLotus) %>%
  column_to_rownames('ID')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

load('kresRBH.RData')

load('/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results_orthologs/rbhAnno.RData')

annoAth <- read_csv('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/Ensembl_ath_Anno.csv',
                    col_types = cols(Chromosome = col_character())) %>%
  mutate_all(list(~replace(., is.na(.), ''))) %>%
  inner_join(rbhAnno) %>%
  setNames(paste0(names(.), '_At'))

annoLotus <- read_csv('/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results/lotus_gifu_collaborator_v1p2_Anno.csv',
                      col_types = cols(Chromosome = col_character())) %>%
  mutate_all(list(~replace(., is.na(.), ''))) %>%
  inner_join(rbhAnno) %>%
  setNames(paste0(names(.), '_Lj'))

## check collaborator annotation and Othofinder
annoMerge <- inner_join(annoAth, annoLotus, by = c('RBH_At' = 'RBH_Lj')) %>%
  dplyr::rename(RBH = RBH_At) %>%
  select(RBH, everything())

## check
(rownames(kresAth) == rownames(kresLotus)) %>%
  sum(.) == nrow(kresAth)

## ## full
## sampleIdx <- 1:36

## without full
sampleIdx <- c(5:16, 25:36)

sampleTable <- rbind(colData(kresAth), colData(kresLotus)) %>%
  .[sampleIdx, , drop = FALSE]

sampleTable$condition %<>% droplevels

kresRBH <- cbind(assay(kresAth), assay(kresLotus)) %>%
  .[, sampleIdx, drop = FALSE] %>%
  DESeqDataSetFromMatrix(., sampleTable, ~ condition)

degres <- DESeq(kresRBH)

## count transformation
rld <- rlog(degres)
ntd <- normTransform(degres)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


##~~~~~~~~~~~~~~~~~~~~~~~~~~hidden batch effect between group~~~~~~~~~~
library('sva')
library('ggplot2')
library('limma')

dat <- rld %>%
  assay %>%
  {.[rowMeans(.) > 1, ]}

groupCond <- rep(c('Ath', 'Lotus'), each = 12) %>%
  as.factor %>%
  {data.frame(condition = .)} %>%
  set_rownames(rownames(colData(degres)))

mod <- model.matrix(~condition, groupCond)
mod0 <- model.matrix(~1, groupCond)

## ## manual detect surrogate variance
## svnum <- 4
## svseq <- svaseq(dat, mod, mod0, n.sv = svnum)

## auto detect sv
svobj <- sva(dat, mod, mod0)
svnum <- svobj$sv %>% ncol

group <- groupCond$condition
design <- model.matrix(~ group)
dat %<>%
  removeBatchEffect(covariates = svobj$sv,
                    design = design)

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
ggsave('auto_og_sv_6.jpg')
ggsave('auto_og_sv_6.pdf')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~hidden batch effect in group~~~~~~~~~~~~~~~~
library('sva')
library('ggplot2')
library('limma')

## ath norm
athCond <- colData(degres)[1:12, ]
athCond$condition %<>% droplevels

mod <- model.matrix(~condition, athCond)
mod0 <- model.matrix(~1, athCond)

## svobjAth <- svaseq(dat[, 1:12], mod, mod0, n.sv = 3)
svobjAth <- sva(dat[, 1:12], mod, mod0)
svnumAth <- svobjAth$sv %>% ncol

groupAth <- sampleTable$condition[1:12] %<>% droplevels
designAth <- model.matrix(~ groupAth)
dat[, 1:12] %<>%
  removeBatchEffect(covariates = svobjAth$sv,
                    design = designAth)

## lotus norm
lotusCond <- colData(degres)[13:24, ]
lotusCond$condition %<>% droplevels

mod <- model.matrix(~condition, lotusCond)
mod0 <- model.matrix(~1, lotusCond)

svobjLotus <- sva(dat[, 13:24], mod, mod0)
svnumLotus <- svobjLotus$sv %>% ncol

groupLotus <- sampleTable$condition[13:24] %<>% droplevels
designLotus <- model.matrix(~ groupLotus)
dat[, 13:24] %<>%
  removeBatchEffect(covariates = svobjLotus$sv,
                    design = designLotus)

svobjLotus$sv %>%
  set_colnames(paste0('sv', seq_len(svnumLotus))) %>%
  as_tibble %>%
  gather(key = 'sv', value = 'value') %>%
  mutate(condition = lotusCond %>%
           .$condition %>%
           rep(svnumLotus) %>%
           as.character,
         sample = rep(colnames(dat[, 13:24]), svnumLotus)) %>%
  mutate(group = paste(sv, condition, sep = '_')) %>%
  ggplot(aes(sample, value, colour = sv, group = group)) +
  geom_point() +
  geom_line() +
  theme(axis.text.x = element_text(angle = 90))
ggsave('auto_og_sv_6.jpg')
ggsave('auto_og_sv_6.pdf')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~DEGs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
degres$sv1 <- svobj$sv[, 1]
degres$sv2 <- svobj$sv[, 2]
degres$sv3 <- svobj$sv[, 3]
degres$sv4 <- svobj$sv[, 4]
degres$sv5 <- svobj$sv[, 5]
degres$sv6 <- svobj$sv[, 6]
design(degres) <- ~sv1 + sv2 + sv3 + sv4 + sv5 + sv6 + condition

degres <- DESeq(degres)

cond <- list(c('AtSC_At', 'Mock_At'),
             c('LjSC_At', 'Mock_At'),
             c('AtSC_At', 'LjSC_At'),
             c('AtSCMloti_Lj', 'Mock_Lj'),
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
  inner_join(annoMerge, by = c('ID' = 'Orthogroup')) %>%
  select(ID, Gene_At : Description_Lj, C_fSC_1 : AtSCMloti_Lj_vs_LjSC_Lj_log2FoldChange) %>%
  arrange(AtSC_At_vs_Mock_At_padj)

write_csv(res, 'SynCom_vs_Mock_og_rmfull_sva_k.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PCA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library('ggrepel')
library('ggplot2')
library('RColorBrewer')

cols <- brewer.pal(9, name = 'Set1')

## ## full
## colorIdx <- 1:9

## without full
colorIdx <- c(2:4, 7:9)

## 1 - 2 C
rldData <- dat
pca <- prcomp(t(rldData))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
percentVar <- round(100 * percentVar)
pca1 <- pca$x[,1]
pca2 <- pca$x[,2]
## pcaData <- data.frame(PC1 = pca1, PC2 = pca2, Group = colData(rld)[, 1], ID = rownames(colData(rld)))

pcaData <- data.frame(PC1 = pca1, PC2 = pca2, Group = rep(c('C_AtSC', 'C_LjSC', 'C_mock', 'L_AtSCMloti', 'L_LjSC', 'L_mock'), each = 4), ID = colnames(rldData))
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

ggsave('PCA_RBH_rmfull_sva.pdf', width = 15)
ggsave('PCA_RBH_rmfull_sva.jpg', width = 15)

save(degres, rldData, file = 'degres_condi_RBH_rmfull.RData')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######################################################################


