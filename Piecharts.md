Pie Charts
================
Lia Baumann
2024-10-22

``` r
map1 <- mapmixture(anc.mat,pops,crs=2056)
map1
```

![](Piecharts_files/figure-gfm/calculate%20Piechart-1.png)<!-- -->

``` r
map2 <- mapmixture(
  admixture_df = anc.mat,
  coords_df = pops,
  cluster_cols = c("#f1a340","#998ec3","salmon2","plum","springgreen3","turquoise3"),
  cluster_names = c("Cluster 1","Cluster 2","Cluster 3","Cluster 4", "Cluster 5","Cluster 6"),
  crs = 2056,
  basemap = rnaturalearthhires::countries10[, c("geometry")],
  boundary = c(xmin=5.5, xmax=10.5, ymin=45, ymax=53),
  pie_size = 0.3,
  pie_border = 0.3,
  pie_border_col = "white",
  pie_opacity = 1,
  land_colour = "#d9d9d9",
  sea_colour = "#deebf7",
  expand = TRUE,
  arrow = TRUE,
  arrow_size = 1.5,
  arrow_position = "bl",
  scalebar = TRUE,
  scalebar_size = 1.5,
  scalebar_position = "tl",
  plot_title = "Ancestry Map",
  plot_title_size = 12,
  axis_title_size = 10,
  axis_text_size = 8
)
map2
```

![](Piecharts_files/figure-gfm/calculate%20Piechart-2.png)<!-- -->

``` r
map3_ohneHannover <- mapmixture(
  admixture_df = anc.mat,
  coords_df = pops,
  cluster_cols = c("#f1a340","#998ec3","salmon2","plum","springgreen3","turquoise3"),
  cluster_names = c("Cluster 1","Cluster 2","Cluster 3","Cluster 4", "Cluster 5","Cluster 6"),
  crs = 2056,
  basemap = rnaturalearthhires::countries10[, c("geometry")],
  boundary = c(xmin=5.5, xmax=10.5, ymin=45, ymax=48.5),
  pie_size = 0.3,
  pie_border = 0.3,
  pie_border_col = "white",
  pie_opacity = 1,
  land_colour = "#d9d9d9",
  sea_colour = "#deebf7",
  expand = TRUE,
  arrow = TRUE,
  arrow_size = 1.5,
  arrow_position = "bl",
  scalebar = TRUE,
  scalebar_size = 1.5,
  scalebar_position = "tl",
  plot_title = "Ancestry Map",
  plot_title_size = 12,
  axis_title_size = 10,
  axis_text_size = 8
)
map3_ohneHannover
```

![](Piecharts_files/figure-gfm/calculate%20Piechart-3.png)<!-- -->

``` r
map4 <- mapmixture(
  admixture_df = anc.mat,
  coords_df = pops,
  cluster_cols = c("#f1a340","#998ec3","salmon2","plum","springgreen3","turquoise3"),
  cluster_names = c("Cluster 1","Cluster 2","Cluster 3","Cluster 4", "Cluster 5","Cluster 6"),
  crs = 2056,
  boundary = c(xmin=5.5, xmax=10.5, ymin=45, ymax=48.5),
  pie_size = 0.3
)+
  # Adjust theme options
  theme(
    legend.position = "top",
    plot.margin = margin(l = 10, r = 10),
  )+
  # Adjust the size of the legend keys
  guides(fill = guide_legend(override.aes = list(size = 5, alpha = 1)))

# Traditional structure barplot
structure_barplot <- structure_plot(
  admixture_df = anc.mat,
  type = "structure",
  cluster_cols = c("#f1a340","#998ec3","salmon2","plum","springgreen3","turquoise"),
  site_dividers = TRUE,
  site_labels_x=1,
  divider_width = 0.4,
  labels = "Population",
  flip_axis = FALSE,
  site_ticks_size = -0.05,
  site_labels_y = -0.35,
  site_labels_size = 2.2
)+
  # Adjust theme options
  theme(
    axis.title.y = element_text(size = 8, hjust = 1),
    axis.text.y = element_text(size = 5),
  )
structure_barplot
```

![](Piecharts_files/figure-gfm/calculate%20Piechart-4.png)<!-- -->

``` r
facet_barplot <- structure_plot(anc.mat,
  type = "facet",
  cluster_cols = c("#f1a340","#998ec3","salmon2","plum","springgreen3","turquoise"),
  facet_col = 2,
  ylabel = "Ancestry proportions",
)+
  theme(
    axis.title.y = element_text(size = 10),
    axis.text.y = element_text(size = 5),
    strip.text = element_text(size = 6, vjust = 1, margin = margin(t=1.5, r=0, b=1.5, l=0)),
  )

# Arrange plots
grid.arrange(map4, structure_barplot, nrow = 2, heights = c(4,1))
```

![](Piecharts_files/figure-gfm/calculate%20Piechart-5.png)<!-- -->

``` r
grid.arrange(map4, facet_barplot, ncol = 2, widths = c(3,2))
```

![](Piecharts_files/figure-gfm/calculate%20Piechart-6.png)<!-- -->
