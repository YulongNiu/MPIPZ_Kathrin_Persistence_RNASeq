# RNA-Seq data (remove full SynCom) from Kathrin Wippel #

<!-- content start -->

**Table of Contents**

- [1. Cluster](#1-cluster)
    - [1.1 Arabidopsis](#11-arabidopsis)
    - [1.2 Lotus](#12-lotus)
- [2. Orthology mapping](#2-orthology-mapping)
    - [2.1 Orthogroup](#21-orthogroup)
    - [2.2 TAIR Best hit](#22-tair-best-hit)
- [References](#references)
    
<!-- content end -->
 
## 1 Cluster

### 1.1 Arabidopsis

*Arabidopsis thaliana* Col-0

* PCA plot

PCA plot of auto sva corrected data

![PCA_ath_sva](results_rmfull/PCA_ath_sva.jpg)

* K-means cluster

![kmeans_10_ath](results_rmfull/kmeans_10_ath.jpg)

* Heatmap with whole transcriptome

![kmeans_10_ath_heatmap2](results_rmfull/kmeans_10_ath_heatmap2.jpg)

* Heatmap with DEGs

![kmeans_10_ath_heatmap_sig2](results_rmfull/kmeans_10_ath_heatmap_sig2.jpg)

### 1.2 Lotus

*Lotus japonicus* Gifu

* PCA plot

PCA plot of auto sva corrected data

![PCA_lotus_sva](results_rmfull/PCA_lotus_sva.jpg)

![PCA_lotus_likeath_sva](results_rmfull/PCA_lotus_likeath_sva.jpg)

![PCA_lotus_rmfull_rmAtSC_sva](results_rmfull/PCA_lotus_rmfull_rmAtSC_sva.jpg)

* K-means cluster

![kmeans_10_lotus](results_rmfull/kmeans_10_lotus.jpg)

![kmeans_10_lotus_likeath](results_rmfull/kmeans_10_lotus_likeath.jpg)

![kmeans_10_lotus_likeath](results_rmfull/kmeans_10_lotus_rmfull_rmAtSC.jpg)

* Heatmap with whole transcriptome

![kmeans_10_lotus_heatmap2](results_rmfull/kmeans_10_lotus_heatmap2.jpg)

![kmeans_10_lotus_likeath_heatmap2](results_rmfull/kmeans_10_lotus_likeath_heatmap2.jpg)

![kmeans_10_lotus_rmfull_rmAtSC_heatmap2](results_rmfull/kmeans_10_lotus_rmfull_rmAtSC_heatmap2.jpg)

* Heatmap with DEGs

![kmeans_10_lotus_heatmap_sig2](results_rmfull/kmeans_10_lotus_heatmap_sig2.jpg)

![kmeans_10_lotus_likeath_heatmap_sig2](results_rmfull/kmeans_10_lotus_heatmap_likeath_sig2.jpg)

![kmeans_10_lotus_heatmap_rmfull_rmAtSC_sig2](results_rmfull/kmeans_10_lotus_heatmap_rmfull_rmAtSC_sig2.jpg)

## 2. Orthology mapping

Row At-transcriptome

```
## number of transcripts for k-means cluster
25725

## number of transcripts in each cluster

   1    2    3    4    5    6    7    8    9   10 
1981 3022 1865 2659 1830 2478 1882 3711 3431 2866 
```

Row Lj-transcriptome (collaborator's annotation)

```
## number of transcripts for k-means cluster
24840

## number of transcripts in each cluster

   1    2    3    4    5    6    7    8    9   10 
2476 2246 2207 2689 2051 2718 2821 2379 3001 2252
```

### 2.1 Orthogroup



### 2.2 TAIR Best hit


