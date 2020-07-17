# Compare clusters (remove full SynCom) with Ka-Wai's flowpot #

<!-- content start -->

**Table of Contents**

- [1. ClusterAt](#1-clusterat)
    - [1.1 DEGs](#11-degs)
    - [1.2 All genes](#12-all-genes)
- [2. ClusterLj](#1-clusterlj)
    - [2.1 Remove fullSynCom](#21-remove-fullsyncom)
    - [2.2 Remove fullSynCom AtSynCom](#22-remove-fullsyncom-atsyncom)
<!-- content end -->

## 1. ClusterAt

### 1.1 DEGs

* Kathrin

![kmeans_10_ath_heatmap_sig2](results_rmfull/kmeans_10_ath_heatmap_sig2.jpg)

* Ka-Wai

![kmeans10_heatmap_soil_sig2](results_rmfull/kmeans10_heatmap_soil_sig2.jpg)

* GO BP enrichment

![kmeans10_ath_cp_BP_dotplot_15_DEG](results_rmfull/kmeans10_ath_cp_BP_dotplot_15_DEG.jpg)

```
         Kathrin1 Kathrin2 Kathrin3 Kathrin4 Kathrin5 Kathrin6 Kathrin7 Kathrin8 Kathrin9 Kathrin10
Ka-Wai1       0      0      0      0      0      0      0      0      0       0
Ka-Wai2       0      0      0      0      0      0      0      0      0       0
Ka-Wai3       0      0      0      0      0      0      0      0      0       0
Ka-Wai4       0      0      0      0      0      0      0      0      0       0
Ka-Wai5       0      0      0      0      0      0      0      0      0       0
Ka-Wai6       0      0      0      0      0      0      0      0      0       0
Ka-Wai7       0      0      0      0      0      0      0      0      0       0
Ka-Wai8       0      0      0      0      0      0      0      0      0       0
Ka-Wai9       0      0      0      0      0      0      0      0      0       0
Ka-Wai10      0      0      0      0      0      0      0      0      0       0
```

### 1.2 All genes

* GO BP enrichment

![kmeans10_ath_cp_BP_dotplot_15_allgene](results_rmfull/kmeans10_ath_cp_BP_dotplot_15_allgene.jpg)

```
         Kathrin1 Kathrin2 Kathrin3 Kathrin4 Kathrin5 Kathrin6 Kathrin7 Kathrin8 Kathrin9 Kathrin10
Ka-Wai1       0      0      0      0      0      0      0      0      0       0
Ka-Wai2       0      0      0      0      0      0      0      0      0       0
Ka-Wai3       0      0      0      0      0      0      0      0      0       0
Ka-Wai4       0      0      0      0      0      0      0      0      0       0
Ka-Wai5       0      0      0      0      0      0      0      0      0       0
Ka-Wai6       0      0      0      0      0      0      0      0      0       0
Ka-Wai7       0      0      0      0      0      0      0      0      0       0
Ka-Wai8       0      0      0      0      0      0      0      0      0       0
Ka-Wai9       0      0      0      0      0      0      0      0      0       0
Ka-Wai10      0      0      0      0      0      0      0      0      0       0
```

## 2. ClusterLj

### 2.1 Remove fullSynCom

* Heatmap with DEGs

![kmeans_10_lotus_heatmap_sig2](results_rmfull/kmeans_10_lotus_heatmap_sig2.jpg)

* GO BP enrichment with raw GO anno

![kmeans10_rmfullSC_cp_BP_dotplot_15_allgene](results_rmfull/kmeans10_rmfullSC_cp_BP_dotplot_15_allgene.jpg)


* GO BP enrichment with besthit GO anno

![kmeans10_rmfullSC_cp_BP_dotplot_15_allgene](results_rmfull/kmeans10_rmfullSC_cp_BP_besthit_dotplot_15_allgene.jpg)

### 2.2 Remove fullSynCom AtSynCom

* Heatmap with DEGs

![kmeans_10_lotus_heatmap_rmfull_rmAtSC_sig2](results_rmfull/kmeans_10_lotus_heatmap_rmfull_rmAtSC_sig2.jpg)

* GO BP enrichment with GO anno

![kmeans10_rmfull_rmAtSC_cp_BP_dotplot_15_allgene](results_rmfull/kmeans10_rmfull_rmAtSC_cp_BP_dotplot_15_allgene.jpg)


* GO BP enrichment with besthit GO anno

![kmeans10_rmfull_rmAtSC_cp_BP_dotplot_15_allgene](results_rmfull/kmeans10_rmfull_rmAtSC_cp_BP_besthit_dotplot_15_allgene.jpg)
