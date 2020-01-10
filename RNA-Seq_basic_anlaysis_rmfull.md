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

Raw At-transcriptome

```
## number of transcripts for k-means cluster
25725

## number of transcripts in each cluster

   1    2    3    4    5    6    7    8    9   10 
1981 3022 1865 2659 1830 2478 1882 3711 3431 2866 
```

Raw Lj-transcriptome (collaborator's annotation)

```
## number of transcripts for k-means cluster
24840

## number of transcripts in each cluster

   1    2    3    4    5    6    7    8    9   10 
2476 2246 2207 2689 2051 2718 2821 2379 3001 2252
```

### 2.1 Orthogroup

A total of `12498` orthogroups were detected from Eik's analysis. Each orthogroup contains at least two transcripts, one is from At and the other is from Lj. 

```
# A tibble: 20,152 x 5
   athID       athcl Orthogroup lotusID              lotuscl
   <chr>       <dbl> <chr>      <chr>                  <dbl>
 1 AT2G25730.2     2 OG0001627  LotjaGi1g1v0391000.1      10
 2 AT3G62150.3     8 OG0000909  LotjaGi5g1v0238400.1       9
 3 AT3G62150.3     8 OG0000909  LotjaGi5g1v0238400.2       9
 4 AT3G48060.3     8 OG0001635  LotjaGi1g1v0183300.1       7
 5 AT1G06410.2     8 OG0003173  LotjaGi3g1v0253100.1       9
 6 AT4G31970.1     9 OG0006747  LotjaGi1g1v0349800.1       3
 7 AT3G25820.3     4 OG0000069  LotjaGi4g1v0043400.1       7
 8 AT5G49960.1     8 OG0006150  LotjaGi6g1v0363700.3       1
 9 AT2G25730.3    10 OG0001627  LotjaGi1g1v0391000.1      10
10 AT5G13800.3     8 OG0001729  LotjaGi3g1v0117700.1      10
# … with 20,142 more rows
```

Combine At-transcriptome with orthogroups

```
## number of transcripts for k-means cluster
20152

## number of transcripts in each cluster

   1    2    3    4    5    6    7    8    9   10 
1564 2452 1470 1943 1414 2032 1513 3011 2465 2288 
```

combine Lj-transcriptome (collaborator's annotation) with orthogroups

```
## number of transcripts for k-means cluster
20152

## number of transcripts in each cluster

   1    2    3    4    5    6    7    8    9   10 
2131 1698 1731 2251 1545 2241 2453 1837 2471 1794 
```

Intersection of At and Lj clusters

```
     Lj1 Lj2 Lj3 Lj4 Lj5 Lj6 Lj7 Lj8 Lj9 Lj10
At1  163 125 149 204 124 193 153 147 185  121
At2  288 199 230 245 175 237 295 221 307  255
At3  158 110 146 160 125 171 196 110 170  124
At4  191 171 166 222 139 207 263 170 243  171
At5  145 130 131 153  95 133 222 110 172  123
At6  220 177 160 243 179 224 232 192 253  152
At7  190 119 128 165 110 171 154 154 195  127
At8  314 279 240 324 225 340 375 261 390  263
At9  227 197 213 278 156 283 295 262 288  266
At10 235 191 168 257 217 282 268 210 268  192
```

### 2.2 TAIR Best hit


