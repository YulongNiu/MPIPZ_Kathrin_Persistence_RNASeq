######################hierarchical clustering####################
setwd('/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results_orthologs/')

library('readr')
library('magrittr')
library('tibble')
library('gplots')
library('dendextend')
library('dynamicTreeCut')
library('ggplot2')
library('tidyr')
library('DESeq2')
library('dplyr')
library('RColorBrewer')
library('gridExtra')
library('cluster')

load('degres_condi_og_rmfull.RData')
deganno <- read_csv('SynCom_vs_Mock_og_rmfull_sva_k.csv')

##~~~~~~~~~~~~~~~~~~~~~~useful funcs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
meanFlg22 <- function(v) {

  require('magrittr')

  res <- v %>%
    split(rep(1 : 3, each = 4)) %>%
    sapply(mean, na.rm = TRUE)

  return(res)
}

##p value calculation from WGCNA
corPvalueStudent <- function(cor, nSamples) {

  ## ref: https://courses.lumenlearning.com/introstats1/chapter/testing-the-significance-of-the-correlation-coefficient/
  T <- sqrt(nSamples - 2) * cor / sqrt(1 - cor^2)

  p <- apply(T, 1:2, function(x) {
    if (x < 0) {
      eachp <- 1 -  pt(x, nSamples - 2, lower.tail = FALSE)
    } else {
      eachp <- pt(x, nSamples - 2, lower.tail = FALSE)
    }
  })

  return(p)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~k-means cluster~~~~~~~~~~~~~~~~~~~~~~~~~
rmfull_all <- c(5:16, 25:36)
sampleN <- c('AtSC_At', 'LjSC_At', 'Mock_At', 'AtSCMloti_Lj', 'LjSC_Lj', 'Mock_Lj')

rmfull_ath <- c(5:16)
sampleN <- c('AtSC_At', 'LjSC_At', 'Mock_At')

rmfull_lotus <- c(25:36)
sampleN <- c('AtSCMloti_Lj', 'LjSC_Lj', 'Mock_Lj')

rawCount <- rldData[, rmfull_ath]

## mean value of normalized count
meanCount <- rawCount %>%
  apply(1, meanFlg22) %>%
  t
colnames(meanCount) <- sampleN

## scale
scaleCount <- meanCount %>%
  t %>%
  scale %>%
  t
scaleCount %<>% .[complete.cases(.), ]

## Cluster rows by Pearson correlation
hr <- scaleCount %>%
  t %>%
  cor(method = 'pearson') %>%
  {1 - .} %>%
  as.dist %>%
  hclust(method = 'complete')

## Clusters columns by Spearman correlation
hc <- scaleCount %>%
  cor(method = 'spearman') %>%
  {1 - .} %>%
  as.dist %>%
  hclust(method = 'complete')

cairo_pdf('hierarch_ath.pdf')
heatmap.2(meanCount,
          Rowv = as.dendrogram(hr),
          Colv = as.dendrogram(hc),
          col = greenred(100),
          scale = 'row',
          margins = c(7, 7),
          cexCol = 0.7,
          labRow = F,
          main = 'Heatmap.2',
          trace = 'none')
dev.off()

hc %>%
  as.dendrogram(method = 'average') %>%
  plot(main = 'Sample Clustering',
       ylab = 'Height')

hr %>%
  as.dendrogram(method = 'average') %>%
  plot(leaflab = 'none',
       main = 'Gene Clustering',
       ylab = 'Height')


z_var <- apply(meanCount, 1, var)
z_mean <- apply(meanCount, 1, mean)
plot(log2(z_mean), log2(z_var), pch = '.')
abline(h = 1, col='red')
abline(v = 1, col='red')
text(x = 13,
     y = 23,
     labels = 'variance > 1 &\n mean > 1',
     col = 'red')

## filter
## meanCount %<>% .[which(z_var > 0 & z_mean > 0), ]

## choose groups
## 1. sum of squared error
wss <- (nrow(scaleCount) - 1) * sum(apply(scaleCount, 2, var))

for (i in 2:20) {
  wss[i] <- sum(kmeans(scaleCount,
                       centers=i,
                       algorithm = 'MacQueen')$withinss)
}

ggplot(tibble(k = 1:20, wss = wss), aes(k, wss)) +
  geom_point(colour = '#D55E00', size = 3) +
  geom_line(linetype = 'dashed') +
  xlab('Number of clusters') +
  ylab('Sum of squared error')
ggsave('kmeans_sse_og_rmfull.pdf')
ggsave('kmeans_sse_og_rmfull.jpg')


## 2. Akaike information criterion
kmeansAIC = function(fit){
  m = ncol(fit$centers)
  n = length(fit$cluster)
  k = nrow(fit$centers)
  D = fit$tot.withinss
  return(D + 2*m*k)
}

aic <- numeric(20)
for (i in 1:20) {
  fit <- kmeans(x = scaleCount, centers = i, algorithm = 'MacQueen')
  aic[i] <- kmeansAIC(fit)
}

ggplot(tibble(k = 1:20, aic = aic), aes(k, wss)) +
  geom_point(colour = '#009E73', size = 3) +
  geom_line(linetype = 'dashed') +
  xlab('Number of clusters') +
  ylab('Akaike information criterion')
ggsave('kmeans_AIC_og_rmfull.pdf')
ggsave('kmeans_AIC_og_rmfull.jpg')

## execute
kClust10 <- kmeans(scaleCount, centers = 10, algorithm = 'MacQueen', nstart = 1000, iter.max = 20)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


##~~~~~~~~~~~~~~~~~~~~~~~~~~plot patterns~~~~~~~~~~~~~~~~~~~~~~~~
cl <- kClust10$cluster
prefix <- 'kmeans_10'

clusterGene <- scaleCount %>%
  as.data.frame %>%
  rownames_to_column(var = 'ID') %>%
  as_tibble %>%
  {
    cl <- as.data.frame(cl) %>%
      rownames_to_column(var = 'ID')
    inner_join(., cl)
  }

## plot core cluster
clusterCore <- clusterGene %>%
  group_by(cl) %>%
  summarise_at(-1, mean, na.rm = TRUE) %>% ## mean of each cluster
  mutate(cl = cl %>% paste0('cluster_', .)) %>%
  gather(Sample, NorExpress, -1)
clusterCore$Sample %<>% factor(levels = sampleN, ordered = TRUE)

ggplot(clusterCore, aes(Sample, NorExpress, col = cl, group = cl)) +
  geom_point() +
  geom_line() +
  facet_wrap(. ~ cl, ncol = 2) +
  ylab('Scaled counts') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  guides(colour = guide_legend(title = 'kmeans (k=10)'))
ggsave(paste0(prefix, '_ath_og_rmfull.pdf'), height = 10)
ggsave(paste0(prefix, '_ath_og_rmfull.jpg'), height = 10)

## plot all genes
clusterGenePlot <- clusterGene %>%
  gather(Sample, NorExpress, -ID, -cl) %>%
  mutate(cl = cl %>% paste0('cluster_', .))
clusterGenePlot$Sample %<>% factor(levels = sampleN, ordered = TRUE)

clusterCorePlot <- clusterCore %>% dplyr::mutate(ID = 1 : nrow(clusterCore))
ggplot(clusterGenePlot, aes(Sample, NorExpress, group = ID)) +
  geom_line(color = 'grey30', alpha = 0.01) +
  facet_wrap(. ~ cl, ncol = 2) +
  geom_point(data = clusterCorePlot, aes(Sample, NorExpress, col = cl, group = ID)) +
  geom_line(data = clusterCorePlot, aes(Sample, NorExpress, group = cl, col = cl)) +
  ylab('Scaled counts') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  guides(colour = guide_legend(title = 'kmeans (k=16)'))
ggsave(paste0(prefix, '_genes_lotus_rmfull_rmAtSC.pdf'), width = 10, dpi = 320)
ggsave(paste0(prefix, '_genes_lotus_rmfull_rmAtSC.jpg'), width = 10, dpi = 320)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write_csv(clusterGene, '_ath_og_rmfull.csv')

