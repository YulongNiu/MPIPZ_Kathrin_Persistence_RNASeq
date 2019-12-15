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

![kmeans_10_ath_heatmap](results_rmfull/kmeans_10_ath_heatmap.jpg)

![kmeans_10_ath_heatmap_2](results_rmfull/kmeans_10_ath_heatmap_2.jpg)

### 3.2 Lotus

*Lotus japonicus* Gifu

* PCA plot

PCA plot of auto sva corrected data

![PCA_lotus_sva](results/PCA_lotus_sva_collaborator.jpg)

* K-means cluster

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
