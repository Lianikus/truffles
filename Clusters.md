Clusters from STRUCTURE
================
Lia Baumann
2024-09-10

# Cluster analysis

``` r
setPop(myData_genind_aMO_k5) <- ~cluster
pca_result_k5 <- dudi.pca(tab(myData_genind_aMO_k5), center = TRUE, scale = TRUE, scannf = FALSE, nf = 2)

# Extract scores for plotting
pca_scores5 <- as.data.frame(pca_result_k5$li)
pca_scores5$Group <- pop(myData_genind_aMO_k5)  # Add population or strata information for coloring

setPop(myData_genind_aMO_k5) <- ~Pop
population_vector5 <- pop(myData_genind_aMO_k5)  # Create a vector with population names
# Convert pop5ulation vector to a data frame
population_df5 <- data.frame(Pop = population_vector5)
pca_scores5$Pop <- population_df5

ggplot(pca_scores5, aes(x = Axis1, y = Axis2, color = Group)) +
  geom_point(size = 3) +
  theme_minimal() +
  geom_text(aes(label = Pop$Pop), vjust = -0.5, hjust = 0.5, size = 3) + 
  labs(title = "PCA mit 5 clusters",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  scale_color_manual(values = color_palette_k5, name = "Cluster")
```

![](Clusters_files/figure-gfm/compute%20PCA%20with%205%20clusters-1.png)<!-- -->

``` r
setPop(myData_genind_aMO_k6) <- ~cluster
pca_result_k6 <- dudi.pca(tab(myData_genind_aMO_k6), center = TRUE, scale = TRUE, scannf = FALSE, nf = 2)

# Extract scores for plotting
pca_scores6 <- as.data.frame(pca_result_k6$li)
pca_scores6$Group <- pop(myData_genind_aMO_k6)  # Add population or strata information for coloring

setPop(myData_genind_aMO_k6) <- ~Pop
population_vector6 <- pop(myData_genind_aMO_k6)  # Create a vector with population names
# Convert population vector to a data frame
population_df6 <- data.frame(Pop = population_vector6)
pca_scores6$Pop <- population_df6

ggplot(pca_scores6, aes(x = Axis1, y = Axis2, color = Group)) +
  geom_point(size = 3) +
  theme_minimal() +
  geom_text(aes(label = Pop$Pop), vjust = -0.5, hjust = 0.5, size = 3) + 
  labs(title = "PCA mit 6 clusters",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  scale_color_manual(values = color_palette_k6, name = "Cluster")
```

![](Clusters_files/figure-gfm/compute%20PCA%20with%206%20clusters-1.png)<!-- -->

``` r
setPop(myData_genind_aMO_k5) <- ~cluster
myData_genind_aMO_k5_minus3 <- popsub(myData_genind_aMO_k5,exclude="3")
pca_result_k5_minus3 <- dudi.pca(tab(myData_genind_aMO_k5_minus3), center = TRUE, scale = TRUE, scannf = FALSE, nf = 2)

# Extract scores for plotting
pca_scores5_minus3 <- as.data.frame(pca_result_k5_minus3$li)
pca_scores5_minus3$Group <- pop(myData_genind_aMO_k5_minus3)  # Add population or strata information for coloring

setPop(myData_genind_aMO_k5_minus3) <- ~Pop
population_vector5_minus3 <- pop(myData_genind_aMO_k5_minus3)  # Create a vector with population names
# Convert pop5ulation vector to a data frame
population_df5_minus3 <- data.frame(Pop = population_vector5_minus3)
pca_scores5_minus3$Pop <- population_df5_minus3

ggplot(pca_scores5_minus3, aes(x = Axis1, y = Axis2, color = Group)) +
  geom_point(size = 3) +
  theme_minimal() +
  geom_text(aes(label = Pop$Pop), vjust = -0.5, hjust = 0.5, size = 3) + 
  labs(title = "PCA mit 4 clusters",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  scale_color_manual(values = color_palette_k5, name = "Group")
```

![](Clusters_files/figure-gfm/pca%20k5%20without%20cluster%203-1.png)<!-- -->

``` r
setPop(myData_genind_aMO_k6) <- ~cluster
myData_genind_aMO_k6_minus3 <- popsub(myData_genind_aMO_k6,exclude="3")
pca_result_k6_minus3 <- dudi.pca(tab(myData_genind_aMO_k6_minus3), center = TRUE, scale = TRUE, scannf = FALSE, nf = 2)

# Extract scores for plotting
pca_scores6_minus3 <- as.data.frame(pca_result_k6_minus3$li)
pca_scores6_minus3$Group <- pop(myData_genind_aMO_k6_minus3)  # Add population or strata information for coloring

setPop(myData_genind_aMO_k6_minus3) <- ~Pop
population_vector6_minus3 <- pop(myData_genind_aMO_k6_minus3)  # Create a vector with population names
# Convert pop5ulation vector to a data frame
population_df6_minus3 <- data.frame(Pop = population_vector6_minus3)
pca_scores6_minus3$Pop <- population_df6_minus3

ggplot(pca_scores6_minus3, aes(x = Axis1, y = Axis2, color = Group)) +
  geom_point(size = 3) +
  theme_minimal() +
  geom_text(aes(label = Pop$Pop), vjust = -0.5, hjust = 0.5, size = 3) + 
  labs(title = "PCA mit 5 clusters",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  scale_color_manual(values = color_palette_k6, name = "Cluster")
```

![](Clusters_files/figure-gfm/pca%20k6%20without%20cluster%203-1.png)<!-- -->

``` r
# Assuming your other table is named `sample_info` and the column is named `sample_codes`
sample_codes_to_keep <- qmatrix_k6_ind$IndID2
# Find indices of the individuals in the genind object that match the sample codes
indices_to_keep <- which(indNames(myData_genind_allMarkersOnly) %in% sample_codes_to_keep)
# Subset the genind object by these indices
filtered_genind <- myData_genind_allMarkersOnly[indices_to_keep, ]

# Convert genind to a matrix
cc_geno_matrix <- tab(filtered_genind, NA.method = "mean")

# Set K to the same number as the optimal STRUCTURE K
optimal_k <- 6 # Replace with your determined K

# Run k-means
set.seed(123) # For reproducibility
kmeans_result <- kmeans(cc_geno_matrix, centers = optimal_k)

# Identify constant columns
constant_cols <- apply(cc_geno_matrix, 2, function(x) var(x, na.rm = TRUE) == 0)

# Remove constant columns
cc_geno_matrix_filtered <- cc_geno_matrix[, !constant_cols]

# Run PCA on the filtered matrix
cc_pca_result <- prcomp(cc_geno_matrix_filtered, scale = TRUE)

# Extract PCA coordinates
cc_pca_coords <- as.data.frame(cc_pca_result$x[, 1:2]) # First two principal components

# Add cluster assignments
cc_pca_coords$STRUCTURE_cluster <- apply(qmatrix_k6, 1, which.max)  # Max probability per row
# Assuming `STRUCTURE_cluster` contains cluster labels as integers from 2 to 7
cc_pca_coords$STRUCTURE_cluster <- as.factor(cc_pca_coords$STRUCTURE_cluster - 1)

cc_pca_coords$KMeans_cluster <- kmeans_result$cluster


# Plot STRUCTURE assignments
p1 <- ggplot(cc_pca_coords, aes(x = PC1, y = PC2, color = as.factor(STRUCTURE_cluster))) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "STRUCTURE Clustering", color = "Cluster") +
  theme_minimal()

# Plot k-means assignments
p2 <- ggplot(cc_pca_coords, aes(x = PC1, y = PC2, color = as.factor(KMeans_cluster))) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_brewer(palette = "Set2") +
  labs(title = "k-means Clustering", color = "Cluster") +
  theme_minimal()

# Arrange plots side by side
grid.arrange(p1, p2, ncol = 2)
```

![](Clusters_files/figure-gfm/compare%20structure%20with%20k-means%20clonecorrected-1.png)<!-- -->

``` r
# Compute Adjusted Rand Index
ari_score <- adjustedRandIndex(cc_pca_coords$STRUCTURE_cluster, cc_pca_coords$KMeans_cluster)
print(paste("Adjusted Rand Index between STRUCTURE and k-means:", ari_score))
```

    ## [1] "Adjusted Rand Index between STRUCTURE and k-means: 0.00501927925404147"

``` r
# Assuming your other table is named `sample_info` and the column is named `sample_codes`
sample_codes_to_keep <- qmatrix_k5_ind$IndID2
# Find indices of the individuals in the genind object that match the sample codes
indices_to_keep <- which(indNames(myData_genind_allMarkersOnly) %in% sample_codes_to_keep)
# Subset the genind object by these indices
filtered_genind <- myData_genind_allMarkersOnly[indices_to_keep, ]

# Convert genind to a matrix
cc_geno_matrix <- tab(filtered_genind, NA.method = "mean")

# Set K to the same number as the optimal STRUCTURE K
optimal_k <- 5 # Replace with your determined K

# Run k-means
set.seed(123) # For reproducibility
kmeans_result <- kmeans(cc_geno_matrix, centers = optimal_k)

# Identify constant columns
constant_cols <- apply(cc_geno_matrix, 2, function(x) var(x, na.rm = TRUE) == 0)

# Remove constant columns
cc_geno_matrix_filtered <- cc_geno_matrix[, !constant_cols]

# Run PCA on the filtered matrix
cc_pca_result <- prcomp(cc_geno_matrix_filtered, scale = TRUE)

# Extract PCA coordinates
cc_pca_coords <- as.data.frame(cc_pca_result$x[, 1:2]) # First two principal components

# Add cluster assignments
cc_pca_coords$STRUCTURE_cluster <- apply(qmatrix_k5_ind, 1, which.max)  # Max probability per row

cc_pca_coords$KMeans_cluster <- kmeans_result$cluster


# Plot STRUCTURE assignments
p1 <- ggplot(cc_pca_coords, aes(x = PC1, y = PC2, color = as.factor(STRUCTURE_cluster))) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "STRUCTURE Clustering", color = "Cluster") +
  theme_minimal()

# Plot k-means assignments
p2 <- ggplot(cc_pca_coords, aes(x = PC1, y = PC2, color = as.factor(KMeans_cluster))) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_brewer(palette = "Set2") +
  labs(title = "k-means Clustering", color = "Cluster") +
  theme_minimal()

# Arrange plots side by side
grid.arrange(p1, p2, ncol = 2)
```

![](Clusters_files/figure-gfm/compare%20structure%20with%20k-means%20clonecorrected%20k%205%20locPrior-1.png)<!-- -->

``` r
# Compute Adjusted Rand Index
ari_score <- adjustedRandIndex(cc_pca_coords$STRUCTURE_cluster, cc_pca_coords$KMeans_cluster)
print(paste("Adjusted Rand Index between STRUCTURE and k-means:", ari_score))
```

    ## [1] "Adjusted Rand Index between STRUCTURE and k-means: 0.00817191282176695"

``` r
kmeans_on_q <- kmeans(cc_geno_matrix, centers = optimal_k)
# Create STRUCTURE-based assignments from the Q-matrix
structure_assignments <- apply(cc_geno_matrix, 1, which.max)  # STRUCTURE cluster assignment based on max probability

# k-means assignments
kmeans_assignments <- kmeans_on_q$cluster

# Calculate ARI to compare STRUCTURE and k-means on Q-matrix
ari_score <- adjustedRandIndex(structure_assignments, kmeans_assignments)
print(paste("Adjusted Rand Index between STRUCTURE and k-means on Q-matrix:", ari_score))
```

    ## [1] "Adjusted Rand Index between STRUCTURE and k-means on Q-matrix: 0.151471899835961"

``` r
# Run PCA on Q-matrix
pca_q <- prcomp(cc_geno_matrix_filtered, scale = TRUE)
pca_coords <- as.data.frame(pca_q$x[, 1:2])  # Get the first two principal components

# Add STRUCTURE and k-means assignments
pca_coords$STRUCTURE_cluster <- as.factor(structure_assignments)
pca_coords$KMeans_cluster <- as.factor(kmeans_assignments)

# Plot STRUCTURE and k-means on Q-matrix side by side
p1 <- ggplot(pca_coords, aes(x = PC1, y = PC2, color = STRUCTURE_cluster)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(title = "STRUCTURE Clustering") +
  theme_minimal()

p2 <- ggplot(pca_coords, aes(x = PC1, y = PC2, color = KMeans_cluster)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(title = "k-means on STRUCTURE Q-matrix") +
  theme_minimal()

grid.arrange(p1, p2, ncol = 2)
```

![](Clusters_files/figure-gfm/running%20k-means%20on%20cluster%20results-1.png)<!-- -->

``` r
# Assuming your other table is named `sample_info` and the column is named `sample_codes`
sample_codes_to_keep <- qmatrix_k5_ind$IndID2
# Find indices of the individuals in the genind object that match the sample codes
indices_to_keep <- which(indNames(myData_genind_allMarkersOnly) %in% sample_codes_to_keep)
# Subset the genind object by these indices
filtered_genind <- myData_genind_allMarkersOnly[indices_to_keep, ]

# Convert genind to a matrix
cc_geno_matrix <- tab(filtered_genind, NA.method = "mean")
# Compute distance matrix from the Q-matrix
dist_matrix <- dist(qmatrix_k5_ind, method = "euclidean")
```

    ## Warning in dist(qmatrix_k5_ind, method = "euclidean"): NAs introduced by
    ## coercion

``` r
# Perform hierarchical clustering
hclust_result <- hclust(dist_matrix, method = "ward.D2")
hierarchical_assignments <- cutree(hclust_result, k = 5)

# Create STRUCTURE-based assignments from the Q-matrix
structure_assignments <- apply(cc_geno_matrix, 1, which.max)  # STRUCTURE cluster assignment based on max probability

# Calculate ARI between STRUCTURE and hierarchical clustering
ari_score_hierarchical <- adjustedRandIndex(structure_assignments, hierarchical_assignments)
print(paste("Adjusted Rand Index between STRUCTURE and hierarchical clustering:", ari_score_hierarchical))
```

    ## [1] "Adjusted Rand Index between STRUCTURE and hierarchical clustering: 0.00703088974830359"

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
setPop(myData_genind_aMO_k5) <- ~cluster
truffleClusters_k5.matFst <- mat_pw_fst(myData_genind_aMO_k5)
truffleClusters_k5.matFst
```

    ##               3          2      mixed          4          1          5
    ## 3     0.0000000 0.15651725 0.14552757 0.15140046 0.16905531 0.20140461
    ## 2     0.1565172 0.00000000 0.04330260 0.04196748 0.09258745 0.10010754
    ## mixed 0.1455276 0.04330260 0.00000000 0.03443488 0.04090685 0.06172819
    ## 4     0.1514005 0.04196748 0.03443488 0.00000000 0.05783851 0.09922283
    ## 1     0.1690553 0.09258745 0.04090685 0.05783851 0.00000000 0.09476120
    ## 5     0.2014046 0.10010754 0.06172819 0.09922283 0.09476120 0.00000000

``` r
boxplot(truffleClusters_k5.matFst, col=color_palette_k5, las=3,
xlab="Cluster", ylab="Fst")
```

![](Clusters_files/figure-gfm/FST%20with%20clusters-1.png)<!-- -->

``` r
#only clusters 1,2,4,5,6 allMarkersOnly
setPop(myData_genind_aMO_k5_minus3) <- ~cluster
(myData_genind_k5_minus3.matFst <- mat_pw_fst(myData_genind_aMO_k5_minus3))
```

    ##                2      mixed          4          1          5
    ## 2     0.00000000 0.04330260 0.04196748 0.09258745 0.10010754
    ## mixed 0.04330260 0.00000000 0.03443488 0.04090685 0.06172819
    ## 4     0.04196748 0.03443488 0.00000000 0.05783851 0.09922283
    ## 1     0.09258745 0.04090685 0.05783851 0.00000000 0.09476120
    ## 5     0.10010754 0.06172819 0.09922283 0.09476120 0.00000000

``` r
boxplot(myData_genind_k5_minus3.matFst, col=color_palette_k5, las=3,
xlab="Cluster", ylab="Fst")
```

![](Clusters_files/figure-gfm/FST%20with%20clusters-2.png)<!-- -->

``` r
#colors are wrong
# Calculate pairwise Fst values
fst_values_k5 <- pairwise.neifst(myData_genind_aMO_k5)

# Convert the Fst matrix to a data frame
fst_k5_df <- as.data.frame(fst_values_k5)
fst_k5_long <- melt(fst_k5_df, varnames = c("1","2","3","4","5","mixed"), value.name = "Fst")
```

    ## No id variables; using all as measure variables

``` r
fst_long <- na.omit(fst_k5_long)

ggplot(fst_k5_long, aes(x = variable, y = Fst, fill=variable)) +
  geom_boxplot() +
  scale_fill_manual(values = color_palette_k5) + 
  theme_minimal() +
  labs(title = "Comparison of Fst Values Across Clusters",
       x = "Cluster",
       y = "Fst Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

    ## Warning: Removed 6 rows containing non-finite outside the scale range
    ## (`stat_boxplot()`).

![](Clusters_files/figure-gfm/FST%20with%20clusters-3.png)<!-- -->

``` r
setPop(myData_genind_aMO_k5) <- ~cluster
truffleClusters_k5.matFst <- mat_pw_fst(myData_genind_aMO_k5)
truffleClusters_k5.matFst
```

    ##               3          2      mixed          4          1          5
    ## 3     0.0000000 0.15651725 0.14552757 0.15140046 0.16905531 0.20140461
    ## 2     0.1565172 0.00000000 0.04330260 0.04196748 0.09258745 0.10010754
    ## mixed 0.1455276 0.04330260 0.00000000 0.03443488 0.04090685 0.06172819
    ## 4     0.1514005 0.04196748 0.03443488 0.00000000 0.05783851 0.09922283
    ## 1     0.1690553 0.09258745 0.04090685 0.05783851 0.00000000 0.09476120
    ## 5     0.2014046 0.10010754 0.06172819 0.09922283 0.09476120 0.00000000

``` r
boxplot(truffleClusters_k5.matFst, col=color_palette_k5, las=3,
xlab="Cluster", ylab="Fst")
```

![](Clusters_files/figure-gfm/FST%20with%20clusters-4.png)<!-- -->

``` r
#same with k=6 now
setPop(myData_genind_aMO_k6_minus3) <- ~cluster
(myData_genind_k6_minus3.matFst <- mat_pw_fst(myData_genind_aMO_k6_minus3))
```

    ##                5          6      mixed          2          4          1
    ## 5     0.00000000 0.06682615 0.05585950 0.10448772 0.07984067 0.10568646
    ## 6     0.06682615 0.00000000 0.04422653 0.06141397 0.07383716 0.07142598
    ## mixed 0.05585950 0.04422653 0.00000000 0.08439333 0.04164793 0.05915624
    ## 2     0.10448772 0.06141397 0.08439333 0.00000000 0.07656007 0.11459039
    ## 4     0.07984067 0.07383716 0.04164793 0.07656007 0.00000000 0.07432676
    ## 1     0.10568646 0.07142598 0.05915624 0.11459039 0.07432676 0.00000000

``` r
# Calculate pairwise Fst values
fst_values_k6 <- pairwise.neifst(myData_genind_aMO_k6)

# Convert the Fst matrix to a data frame
fst_k6_df <- as.data.frame(fst_values_k6)
fst_k6_long <- melt(fst_k6_df, varnames = c("1","2","3","4","5","6","mixed"), value.name = "Fst")
```

    ## No id variables; using all as measure variables

``` r
fst_k6_long <- na.omit(fst_k6_long)

ggplot(fst_k6_long, aes(x = variable, y = Fst, fill=variable)) +
  geom_boxplot() +
  scale_fill_manual(values = color_palette_k6) + 
  theme_minimal() +
  labs(title = "Comparison of Fst Values Across Clusters",
       x = "Cluster",
       y = "Fst Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

![](Clusters_files/figure-gfm/FST%20with%20clusters-5.png)<!-- -->
