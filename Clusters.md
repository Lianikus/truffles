Clusters from STRUCTURE
================
Lia Baumann
2024-09-10

# Cluster analysis

``` r
# Convert genind to a matrix if needed
geno_matrix <- tab(myData_genind_allMarkersOnly, NA.method = "mean")

# Set K to the same number as the optimal STRUCTURE K
optimal_k <- 6 # Replace with your determined K

# Run k-means
set.seed(123) # For reproducibility
kmeans_result <- kmeans(geno_matrix, centers = optimal_k)

# Identify constant columns
constant_cols <- apply(geno_matrix, 2, function(x) var(x, na.rm = TRUE) == 0)

# Remove constant columns
geno_matrix_filtered <- geno_matrix[, !constant_cols]

# Run PCA on the filtered matrix
pca_result <- prcomp(geno_matrix_filtered, scale = TRUE)

# Extract PCA coordinates
pca_coords <- as.data.frame(pca_result$x[, 1:2]) # First two principal components

# Add cluster assignments
pca_coords$STRUCTURE_cluster <- apply(qmatrix_k6, 1, which.max)  # Max probability per row
# Assuming `STRUCTURE_cluster` contains cluster labels as integers from 2 to 7
pca_coords$STRUCTURE_cluster <- as.factor(pca_coords$STRUCTURE_cluster - 1)

pca_coords$KMeans_cluster <- kmeans_result$cluster


# Plot STRUCTURE assignments
p1 <- ggplot(pca_coords, aes(x = PC1, y = PC2, color = as.factor(STRUCTURE_cluster))) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "STRUCTURE Clustering", color = "Cluster") +
  theme_minimal()

# Plot k-means assignments
p2 <- ggplot(pca_coords, aes(x = PC1, y = PC2, color = as.factor(KMeans_cluster))) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_brewer(palette = "Set2") +
  labs(title = "k-means Clustering", color = "Cluster") +
  theme_minimal()

# Arrange plots side by side
grid.arrange(p1, p2, ncol = 2)
```

![](Clusters_files/figure-gfm/compare%20structure%20with%20k-means-1.png)<!-- -->

``` r
# Compute Adjusted Rand Index
ari_score <- adjustedRandIndex(pca_coords$STRUCTURE_cluster, pca_coords$KMeans_cluster)
print(paste("Adjusted Rand Index between STRUCTURE and k-means:", ari_score))
```

    ## [1] "Adjusted Rand Index between STRUCTURE and k-means: 0.0941534752897646"

``` r
n_per_cluster <- T_all_withMLGs_withClusters_allMarkersOnly %>%
  group_by(`clusterssign(q>0.5)`) %>%
  summarise(n=n()) %>%
  arrange() %>%
  kable()
n_per_cluster
```

| clusterssign(q\>0.5) |   n |
|:---------------------|----:|
| cluster1             | 408 |
| cluster2             | 253 |
| cluster3             | 568 |
| cluster4             | 174 |
| cluster5             | 233 |
| cluster6             | 849 |
| mixed                | 154 |

``` r
T_all_withMLGs_withClusters_allMarkersOnly$Sampling_date <- as_date(T_all_withMLGs_withClusters_allMarkersOnly$Sampling_date)
T_all_withMLGs_withClusters_allMarkersOnly$`clusterssign(q>0.5)` <- factor(T_all_withMLGs_withClusters_allMarkersOnly$`clusterssign(q>0.5)`, ordered=TRUE, levels = c("cluster1","cluster2","cluster3","cluster4","cluster5","cluster6","mixed"))

ggplot(T_all_withMLGs_withClusters_allMarkersOnly, aes(x=Sampling_date, y= fct_rev(Site_1_abrev))) +
  geom_point(aes(colour=factor('clusterssign(q<0.5)'))) +
  theme_publish(base_family = "Times") +
  scale_colour_brewer(palette= "Greens") +
  scale_x_date(date_labels="%Y", date_breaks = "5 years") +
  labs(title = "Overview of samples in truffle monitoring dataset", y = "Sampling site abbreviation", x ="Sampling date", colour="Genetic cluster")
```

    ## Warning: Removed 3 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    ## not found in Windows font database

    ## Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    ## not found in Windows font database

    ## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    ## family not found in Windows font database

    ## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    ## family not found in Windows font database

    ## Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    ## not found in Windows font database

    ## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    ## family not found in Windows font database

    ## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    ## family not found in Windows font database

    ## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    ## family not found in Windows font database

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## font family not found in Windows font database

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## font family not found in Windows font database

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## font family not found in Windows font database

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## font family not found in Windows font database

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## font family not found in Windows font database

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## font family not found in Windows font database

![](Clusters_files/figure-gfm/how%20many%20clusters-1.png)<!-- -->

``` r
setPop(myData_genind_withClusters) <- ~cluster
myData_genind_allMarkersOnly_withoutCluster3_SY <- popsub(myData_genind_withClusters,sublist=c("1","2","4","5","6"))
poppr(myData_genind_allMarkersOnly_withoutCluster3_SY)
```

    ##     Pop    N MLG eMLG   SE    H     G lambda   E.5  Hexp    Ia  rbarD
    ## 1     5  233  47 40.7 2.06 2.92 10.28  0.903 0.530 0.264 1.510 0.1484
    ## 2     6  849  94 37.5 3.54 2.87  7.28  0.863 0.378 0.408 1.745 0.1503
    ## 3     2  253  73 55.0 3.08 3.16 11.05  0.909 0.446 0.340 1.374 0.1391
    ## 4     4  174  81 81.0 0.00 3.98 33.42  0.970 0.618 0.361 0.440 0.0376
    ## 5     1  408  94 51.9 3.88 3.25 11.64  0.914 0.431 0.291 1.807 0.1714
    ## 6 Total 1917 389 79.9 5.29 4.53 29.34  0.966 0.310 0.494 0.779 0.0653
    ##                                              File
    ## 1 myData_genind_allMarkersOnly_withoutCluster3_SY
    ## 2 myData_genind_allMarkersOnly_withoutCluster3_SY
    ## 3 myData_genind_allMarkersOnly_withoutCluster3_SY
    ## 4 myData_genind_allMarkersOnly_withoutCluster3_SY
    ## 5 myData_genind_allMarkersOnly_withoutCluster3_SY
    ## 6 myData_genind_allMarkersOnly_withoutCluster3_SY

``` r
setPop(myData_genind_allMarkersOnly_withoutCluster3_SY) <- ~Pop/SamplingYear
n_samples_allMarkersOnly_withoutCluster3_SY <- poppr(myData_genind_allMarkersOnly_withoutCluster3_SY) %>%
  separate(Pop, sep="_", into=c("Pop","SamplingYear"))
```

    ## Warning: Expected 2 pieces. Missing pieces filled with `NA` in 1 rows [106].

``` r
n_samples_allMarkersOnly_withoutCluster3_SY %>%
  select(Pop,SamplingYear,N) %>%
  pivot_wider(names_from=SamplingYear, values_from=N) %>%
  select(.,-"NA") %>%
  replace(is.na(.),0) %>%
  kable()
```

| Pop   | 2011 | 2012 | 2013 | 2014 | 2015 | 2016 | 2017 | 2018 | 2019 | 2020 | 2021 | 2022 | 2023 |
|:------|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|
| ALD   |   24 |   12 |   15 |    7 |    3 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |
| RIE   |   14 |    9 |   11 |   12 |    3 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |
| TRO   |   16 |   10 |    8 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |
| SCG   |    9 |    3 |    1 |    2 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |
| BOH   |   11 |   13 |    6 |   10 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |
| BOB   |    6 |   35 |   50 |   29 |   15 |   52 |  142 |   63 |   27 |    0 |    0 |    0 |    0 |
| FRB   |    4 |    1 |    1 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |
| UEB   |   27 |    6 |   54 |   36 |   27 |    0 |    0 |    0 |    3 |    0 |    0 |    0 |    0 |
| FRE   |    0 |    1 |    2 |    0 |    0 |    0 |    2 |    0 |    0 |    0 |    0 |    0 |    0 |
| SCL   |    0 |   51 |    7 |   24 |    1 |    9 |    5 |    0 |    0 |    0 |    0 |    0 |    0 |
| SCD   |    0 |   23 |    4 |    6 |   15 |   43 |    3 |    0 |    0 |    0 |    0 |    0 |    0 |
| WSL   |    0 |    1 |    0 |    0 |    0 |   58 |   61 |   14 |   72 |   72 |    3 |    0 |    0 |
| BUR   |    0 |    0 |    2 |    2 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |
| SCS   |    0 |    0 |    7 |   15 |    0 |   19 |    2 |    0 |    0 |    0 |    0 |    0 |    0 |
| NEU   |    0 |    0 |    4 |    3 |    0 |   13 |    2 |    0 |    0 |    4 |    1 |    0 |    0 |
| UST   |    0 |    0 |   23 |    8 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |
| KON   |    0 |    0 |    3 |   12 |    2 |   25 |   17 |   41 |   26 |   25 |   27 |   11 |   23 |
| FRI   |    0 |    0 |   14 |    4 |    0 |    8 |   36 |    1 |    0 |    0 |    0 |    0 |    0 |
| BAR   |    0 |    0 |    0 |    1 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |
| LIM   |    0 |    0 |    0 |    0 |    0 |   52 |   66 |   53 |    0 |    0 |    0 |    0 |    0 |
| BRU   |    0 |    0 |    0 |    0 |    0 |    8 |    3 |   35 |    0 |    0 |    0 |    0 |    0 |
| HAN   |    0 |    0 |    0 |    0 |    0 |    0 |   35 |   17 |    3 |    0 |    0 |    0 |    0 |
| GEN   |    0 |    0 |    0 |    0 |    0 |    0 |    0 |    1 |    0 |    3 |    1 |    5 |    0 |
| Total |    0 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |

``` r
grp.1pca <- find.clusters(myData_genind_allMarkersOnly, n.pca=1, n.clust=6)
grp.3pca <- find.clusters(myData_genind_allMarkersOnly, n.pca=3, n.clust=6)
# 3 PC axis retained and 6 clusters
grp.5pca <- find.clusters(myData_genind_allMarkersOnly, n.pca=5, n.clust=6)
grp.21pca <- find.clusters(myData_genind_allMarkersOnly, n.pca=21, n.clust=6)

options(max.print=4000)
grp.1pca_df <- as_tibble(data.frame(grp.1pca$grp), rownames="Sample")
grp.3pca_df <- as_tibble(data.frame(grp.3pca$grp), rownames="Sample")
grp.5pca_df <- as_tibble(data.frame(grp.5pca$grp), rownames="Sample")
grp.21pca_df <- as_tibble(data.frame(grp.21pca$grp), rownames="Sample")
#write.csv(grp_df, file="kmeans_clusters_allMarkersOnly.csv",sep=";")
compare_clusters <- left_join(T_all_withMLGs_withClusters_allMarkersOnly,grp.3pca_df,by=c("Code"="Sample"))
compare_clusters <- left_join(compare_clusters,grp.5pca_df,by=c("Code"="Sample"))
compare_clusters <- left_join(compare_clusters,grp.1pca_df,by=c("Code"="Sample"))
compare_clusters <- left_join(compare_clusters,grp.21pca_df,by=c("Code"="Sample"))
compare_clusters
```

    ## # A tibble: 2,639 × 11
    ##    Code   MLG    Site_1_abrev truffle_year Sampling_year Sampling_date
    ##    <chr>  <chr>  <chr>               <dbl>         <dbl> <date>       
    ##  1 TROSS3 MLG13  TRO                  2010          2011 2011-03-27   
    ##  2 TROSS4 MLG13  TRO                  2010          2011 2011-03-27   
    ##  3 BUR315 MLG7   BUR                  2022          2022 2022-08-13   
    ##  4 FREI1  MLG354 FRE                  2010          2011 2011-01-11   
    ##  5 ALD1   MLG363 ALD                  2010          2011 2011-03-05   
    ##  6 ALD2   MLG363 ALD                  2010          2011 2011-03-05   
    ##  7 FREI2  MLG354 FRE                  2010          2011 2011-03-05   
    ##  8 FREI3  MLG354 FRE                  2010          2011 2011-03-05   
    ##  9 RIETH1 MLG285 RIE                  2010          2011 2011-03-05   
    ## 10 RIETH2 MLG185 RIE                  2010          2011 2011-03-05   
    ## # ℹ 2,629 more rows
    ## # ℹ 5 more variables: `clusterssign(q>0.5)` <ord>, grp.3pca.grp <fct>,
    ## #   grp.5pca.grp <fct>, grp.1pca.grp <fct>, grp.21pca.grp <fct>

``` r
#write.csv(compare_clusters,file="compare_clusters.csv")

#if I compare clusters from k-means and from STRUCTURE, they don't seem very similar and consistent
```

``` r
setPop(myData_genind_withClusters) <- ~cluster
(truffleClusters.matFst <- mat_pw_fst(myData_genind_withClusters))
```

    ##                3          5          6      mixed          2           
    ## 3     0.00000000 0.18852779 0.16042288 0.14556082 0.15681700 0.08580212
    ## 5     0.18852779 0.00000000 0.06682615 0.05585950 0.10448772 0.08325061
    ## 6     0.16042288 0.06682615 0.00000000 0.04422653 0.06141397 0.04752407
    ## mixed 0.14556082 0.05585950 0.04422653 0.00000000 0.08439333 0.05846984
    ## 2     0.15681700 0.10448772 0.06141397 0.08439333 0.00000000 0.02798734
    ##       0.08580212 0.08325061 0.04752407 0.05846984 0.02798734 0.00000000
    ## 4     0.14764463 0.07984067 0.07383716 0.04164793 0.07656007 0.06515767
    ## 1     0.17925558 0.10568646 0.07142598 0.05915624 0.11459039 0.08143417
    ##                4          1
    ## 3     0.14764463 0.17925558
    ## 5     0.07984067 0.10568646
    ## 6     0.07383716 0.07142598
    ## mixed 0.04164793 0.05915624
    ## 2     0.07656007 0.11459039
    ##       0.06515767 0.08143417
    ## 4     0.00000000 0.07432676
    ## 1     0.07432676 0.00000000

``` r
boxplot(truffleClusters.matFst, col=funky(nPop(myData_genind_withClusters)), las=3,
xlab="Cluster", ylab="Fst")
```

![](Clusters_files/figure-gfm/FST%20with%20clusters-1.png)<!-- -->

``` r
#only clusters 1,2,4,5,6 allMarkersOnly
setPop(myData_genind_allMarkersOnly_withoutCluster3_SY) <- ~cluster
(myData_genind_allMarkersOnly_withoutCluster3_SY.matFst <- mat_pw_fst(myData_genind_allMarkersOnly_withoutCluster3_SY))
```

    ##            5          6          2          4          1
    ## 5 0.00000000 0.06682615 0.10448772 0.07984067 0.10568646
    ## 6 0.06682615 0.00000000 0.06141397 0.07383716 0.07142598
    ## 2 0.10448772 0.06141397 0.00000000 0.07656007 0.11459039
    ## 4 0.07984067 0.07383716 0.07656007 0.00000000 0.07432676
    ## 1 0.10568646 0.07142598 0.11459039 0.07432676 0.00000000

``` r
boxplot(myData_genind_allMarkersOnly_withoutCluster3_SY.matFst, col=funky(nPop(myData_genind_allMarkersOnly_withoutCluster3_SY)), las=3,
xlab="Cluster", ylab="Fst")
```

![](Clusters_files/figure-gfm/FST%20with%20clusters-2.png)<!-- -->
