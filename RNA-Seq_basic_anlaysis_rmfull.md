# RNA-Seq data (remove full SynCom) from Kathrin Wippel #

<!-- content start -->

**Table of Contents**

- [1. Cluster](#3-cluster)
    - [1.1 Arabidopsis](#31-arabidopsis)
    - [1.2 Lotus](#32-lotus)
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

* Heatmap

![kmeans_10_heatmap_ath](results/kmeans_10_heatmap_ath.jpg)

### 3.2 Lotus

*Lotus japonicus* Gifu

* PCA plot

PCA plot of raw data

![PCA_lotus_raw](results/PCA_lotus_raw.jpg)

PCA plot of auto sva corrected data

![auto_lotus_sv](results/auto_lotus_sv_collaborator.jpg)

![PCA_lotus_sva](results/PCA_lotus_sva_collaborator.jpg)

PCA plot of limma correction data (`1 2` and `3 4`)

![PCA_lotus_limma](results/PCA_lotus_limma.jpg)

![PCA_lotus_noAtSC_limma](results/PCA_lotus_limma_noAtSC.jpg)

![PCA_lotus_noAtSCMloti_limma](results/PCA_lotus_limma_noAtSCMloti.jpg)

* K-means cluster

![kmeans_sse](results/kmeans_sse_lotus_collaborator.jpg)

![kmeans_AIC](results/kmeans_AIC_lotus_collaborator.jpg)

![kmeans_10_lotus_collaborator](results/kmeans_10_lotus_collaborator.jpg)
  
![kmeans_10_gene_lotus_collaborator](results/kmeans_10_genes_lotus_collaborator.jpg)

* Trait

```
     fullSC AtSC AtSCMloti LjSC LjNodule218
        1    0         0    0           1
        0    1         0    0           0
        0    0         1    0           1
        0    0         0    1           1
        0    0         0    0           1
```

![kmeans_10_trait_lotus_collaborator](results/kmeans_10_trait_lotus_collaborator.jpg)

* Heatmap

![kmeans_10_heatmap_lotus_collaborator](results/kmeans_10_heatmap_lotus_collaborator.jpg)
