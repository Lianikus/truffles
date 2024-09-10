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

| clusterssign(q\>0.5) |    n |
|:---------------------|-----:|
| cluster1             |  544 |
| cluster2             |  304 |
| cluster3             |  702 |
| cluster4             |  179 |
| cluster5             |  310 |
| cluster6             | 1278 |
| mixed                |  154 |

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

    ## Warning: Removed 4 rows containing missing values or values outside the scale range
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
