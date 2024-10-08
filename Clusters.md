Clusters from STRUCTURE
================
Lia Baumann
2024-09-10

# Cluster analysis

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
grp <- find.clusters(myData_genind_allMarkersOnly, n.pca=3)
```

![](Clusters_files/figure-gfm/k-means%20clusters-1.png)<!-- -->

    ## Choose the number of clusters (>=2):

``` r
# 1 PC axis retained and 6 clusters

grp$size
```

    ## [1] 2639

``` r
options(max.print=4000)
grp_df <- as_tibble(data.frame(grp$grp), rownames="Sample")
#write.csv(grp_df, file="kmeans_clusters_allMarkersOnly.csv",sep=";")
compare_clusters <- left_join(T_all_withMLGs_withClusters_allMarkersOnly,grp_df,by=c("Code"="Sample"))

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
