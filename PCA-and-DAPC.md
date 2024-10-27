PCA and DAPC
================
Lia Baumann
2024-06-04

## 1. PCA

``` r
#uncorrected PCA with allMarkersOnly dataset
setPop(myData_genind_allMarkersOnly) <- ~Pop
x.pops_allMarkersOnly <- tab(myData_genind_allMarkersOnly,
                             freq=TRUE, NA.method="mean")

pca.pops.allMarkersOnly <- dudi.pca(x.pops_allMarkersOnly,center=TRUE,scale=FALSE, nf=3, scannf=FALSE)
s.label(pca.pops.allMarkersOnly$li)
```

![](PCA-and-DAPC_files/figure-gfm/pca-1.png)<!-- -->

``` r
# Calculate explained variance for each PCA axis
explained_variance <- pca.pops.allMarkersOnly$eig / sum(pca.pops.allMarkersOnly$eig) * 100
barplot(explained_variance, main = "Explained Variance per PCA Axis", xlab = "Principal Component", ylab = "Variance Explained (%)")
```

![](PCA-and-DAPC_files/figure-gfm/pca-2.png)<!-- -->

``` r
#looking at first two axes:
s.class(pca.pops.allMarkersOnly$li,fac=pop(myData_genind_allMarkersOnly),col=transp(funky(15),0.6),axesel=FALSE,cstar=0,cpoint=3)
add.scatter.eig(pca.pops.allMarkersOnly$eig[1:50],3,1,2, ratio=.17)
```

![](PCA-and-DAPC_files/figure-gfm/pca-3.png)<!-- -->

``` r
#looking at 2 and 3 axis:
s.class(pca.pops.allMarkersOnly$li, fac=pop(myData_genind_allMarkersOnly),
xax=2, yax=3, col=transp(funky(15),.6),
axesel=FALSE, cstar=0, cpoint=3)
add.scatter.eig(pca.pops.allMarkersOnly$eig[1:50],3,2,3, ratio=.2)
```

![](PCA-and-DAPC_files/figure-gfm/pca-4.png)<!-- -->

``` r
#emf(file="pca.allMarkersOnly.emf")
#s.class(pca.pops.allMarkersOnly$li, fac=pop(myData_genind_allMarkersOnly),
 #       col=funky(15))
#dev.off()


#clonecorrected PCA  with samplingYear
setPop(cc_myData_genind_allMarkersOnly_SY) <- ~Pop
x.pops_cc.SY <- tab(cc_myData_genind_allMarkersOnly_SY,
                    freq=TRUE, NA.method="mean")
pca.pops_cc.SY <- dudi.pca(df = x.pops_cc.SY,
                           center = TRUE, scale = FALSE, scannf = FALSE, nf = 3)

# Calculate explained variance for each PCA axis
cc.explained_variance <- pca.pops_cc.SY$eig / sum(pca.pops_cc.SY$eig) * 100
barplot(cc.explained_variance, main = "Explained Variance per PCA Axis", xlab = "Principal Component", ylab = "Variance Explained (%)")
```

![](PCA-and-DAPC_files/figure-gfm/pca-5.png)<!-- -->

``` r
#looking at first two axes:
s.class(pca.pops_cc.SY$li,
        fac=pop(cc_myData_genind_allMarkersOnly_SY),
        col=transp(funky(15),0.6),axesel=FALSE,cstar=0,cpoint=3)
add.scatter.eig(pca.pops_cc.SY$eig[1:50],3,1,2, ratio=.17)
```

![](PCA-and-DAPC_files/figure-gfm/pca-6.png)<!-- -->

``` r
#looking at 2 and 3 axis:
s.class(pca.pops_cc.SY$li, fac=pop(cc_myData_genind_allMarkersOnly_SY),
xax=2, yax=3, col=transp(funky(15),.6),
axesel=FALSE, cstar=0, cpoint=3)
add.scatter.eig(pca.pops_cc.SY$eig[1:50],3,2,3, ratio=.2)
```

![](PCA-and-DAPC_files/figure-gfm/pca-7.png)<!-- -->

``` r
#emf(file="pca.cc.allMarkersOnly_SY.emf")
#s.class(pca.pops_cc.SY$li,
#        fac=pop(cc_myData_genind_allMarkersOnly_SY),
#        col=funky(15))
#dev.off()


#clonecorrected PCA  with TruffleYear
setPop(cc_myData_genind_allMarkersOnly_TY) <- ~Pop
x.pops_cc.TY <- tab(cc_myData_genind_allMarkersOnly_TY,
                    freq=TRUE, NA.method="mean")
pca.pops_cc.TY <- dudi.pca(df = x.pops_cc.TY,
                           center = TRUE, scale = FALSE, scannf = FALSE, nf = 3)
#looking at first two axes:
s.class(pca.pops_cc.TY$li,
        fac=pop(cc_myData_genind_allMarkersOnly_TY),
        col=transp(funky(15),0.6),axesel=FALSE,cstar=0,cpoint=3)
add.scatter.eig(pca.pops_cc.TY$eig[1:50],3,1,2, ratio=.17)
```

![](PCA-and-DAPC_files/figure-gfm/pca-8.png)<!-- -->

``` r
#looking at 2 and 3 axis:
s.class(pca.pops_cc.TY$li, fac=pop(cc_myData_genind_allMarkersOnly_TY),
xax=2, yax=3, col=transp(funky(15),.6),
axesel=FALSE, cstar=0, cpoint=3)
add.scatter.eig(pca.pops_cc.TY$eig[1:50],3,2,3, ratio=.2)
```

![](PCA-and-DAPC_files/figure-gfm/pca-9.png)<!-- -->

``` r
#emf(file="pca.cc.allMarkersOnly_TY.emf")
#s.class(pca.pops_cc.TY_2$li,
#        fac=pop(cc_myData_genind_allMarkersOnly_TY),
#        col=funky(15))
#dev.off()

#plot(pca.pops_cc.TY$li, col=
```

# 2. DAPC

DAPC was pioneered by Jombart and colleagues (Jombart et al., 2010) and
can be used to infer the number of clusters of genetically related
individuals. In this multivariate statistical approach variance in the
sample is partitioned into a between-group and within- group component,
in an effort to maximize discrimination between groups. In DAPC, data is
first transformed using a principal components analysis (PCA) and
subsequently clusters are identified using discriminant analysis (DA).
This tutorial is based on the vignette written by Thibaut Jombart. We
encourage the user to explore this vignette further. The vignette can
also be opened within R by executing adegenetTutorial or
<https://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf>

Usual approaches such as Principal Component Analysis (PCA) or Principal
Coordinates Analysis (PCoA / MDS) focus on V AR(X) (–\> total variance =
(variance between groups) + (variance within groups)). That is, they
only describe the global diversity, possibly overlooking differences
between groups. On the contrary, DAPC optimizes variances between groups
while minimizing variances within groups: it seeks synthetic variables,
the discriminant functions, which show differences between groups as
best as possible while minimizing variation within clusters.

The trade-off between power of discrimination and over-fitting can be
measured by the ascore, which is simply the difference between the
proportion of successful reassignment of the analysis (observed
discrimination) and values obtained using random groups (random
discrimination). It can be seen as the proportion of successful
reassignment corrected for the number of retained PCs. It is implemented
by a.score, which relies on repeating the DAPC analysis using randomized
groups, and computing a-scores for each group, as well as the average
a-score

``` r
setPop(myData_genind_allMarkersOnly) <- ~Pop
dapc2 <- dapc(myData_genind_allMarkersOnly, n.da=100, n.pca=10)
temp <- a.score(dapc2)
names(temp)
```

    ## [1] "tab"       "pop.score" "mean"

``` r
temp$pop.score
```

    ##         FRE         ALD         RIE         TRO         SCG         BOH 
    ##  0.95700935  0.87096774  0.56451613  0.97058824  1.00000000  0.67500000 
    ##         BOB         FRB         UEB         SCL         SCD         WSL 
    ## -0.03146853  0.46153846  0.95424837  0.86868687  0.81000000  0.95618375 
    ##         BUR         SCS         NEU         UST         KON         FRI 
    ##  0.96143498  0.97674419  0.78125000  0.68571429  0.85141509  0.96666667 
    ##         BAR         BRU         LIM         HAN         GEN 
    ##  0.80000000  0.71910112  0.73055556  0.44615385  0.90000000

``` r
temp$mean
```

    ## [1] 0.7772307

``` r
dapc2 <- dapc(myData_genind_allMarkersOnly, n.da=100, n.pca=50)
temp <- optim.a.score(dapc2)
```

![](PCA-and-DAPC_files/figure-gfm/define%20a-score-1.png)<!-- -->

``` r
#optimal score is 21 PCAs

#now the same corrected
setPop(cc_myData_genind_allMarkersOnly_SY) <- ~Pop
cc.dapc2 <- dapc(cc_myData_genind_allMarkersOnly_SY, n.da=100, n.pca=10)
temp <- a.score(cc.dapc2)
names(temp)
```

    ## [1] "tab"       "pop.score" "mean"

``` r
temp$pop.score
```

    ##        FRE        ALD        RIE        TRO        SCG        BOH        BOB 
    ##  0.7857143  0.6571429  0.6240741  0.9000000  0.8666667  0.7515152 -0.0840796 
    ##        FRB        UEB        SCL        SCD        WSL        BUR        SCS 
    ##  0.0000000  0.9358491  0.8500000  0.7437500  0.8461538  0.8687500  0.7666667 
    ##        NEU        UST        KON        FRI        BAR        BRU        LIM 
    ##  0.8428571  0.2608696  0.9450000  0.7500000  0.5000000  0.0000000  0.1875000 
    ##        HAN        GEN 
    ##  0.3076923  0.7800000

``` r
temp$mean
```

    ## [1] 0.6124401

``` r
cc.dapc2 <- dapc(cc_myData_genind_allMarkersOnly_SY, n.da=100, n.pca=50)
temp <- optim.a.score(cc.dapc2)
```

![](PCA-and-DAPC_files/figure-gfm/define%20a-score-2.png)<!-- -->

``` r
temp
```

    ## $pop.score
    ## $pop.score$`1`
    ##        FRE        ALD        RIE        TRO        SCG        BOH        BOB 
    ##  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000 -0.1691542 
    ##        FRB        UEB        SCL        SCD        WSL        BUR        SCS 
    ##  0.0000000  0.3962264  0.7500000  0.2916667  0.0000000  0.0000000  0.0000000 
    ##        NEU        UST        KON        FRI        BAR        BRU        LIM 
    ##  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000 
    ##        HAN        GEN 
    ##  0.0000000  0.0000000 
    ## 
    ## $pop.score$`5`
    ##         FRE         ALD         RIE         TRO         SCG         BOH 
    ##  0.64285714  0.66666667  0.55555556  0.40000000  0.22222222  0.60606061 
    ##         BOB         FRB         UEB         SCL         SCD         WSL 
    ## -0.16567164  0.00000000  0.90566038  0.78958333  0.52083333  0.84615385 
    ##         BUR         SCS         NEU         UST         KON         FRI 
    ##  0.87500000  0.88888889  0.00000000  0.04347826  0.88333333  0.30000000 
    ##         BAR         BRU         LIM         HAN         GEN 
    ## -0.20000000  0.00000000  0.12500000  0.00000000  0.00000000 
    ## 
    ## $pop.score$`10`
    ##         FRE         ALD         RIE         TRO         SCG         BOH 
    ##  0.77142857  0.64761905  0.62777778  0.90000000  0.86666667  0.75757576 
    ##         BOB         FRB         UEB         SCL         SCD         WSL 
    ## -0.07661692 -0.07000000  0.94339623  0.83750000  0.73750000  0.84615385 
    ##         BUR         SCS         NEU         UST         KON         FRI 
    ##  0.87187500  0.77777778  0.85714286  0.26086957  0.95000000  0.75000000 
    ##         BAR         BRU         LIM         HAN         GEN 
    ##  0.60000000  0.00000000  0.20833333  0.26923077  0.80000000 
    ## 
    ## $pop.score$`15`
    ##         FRE         ALD         RIE         TRO         SCG         BOH 
    ##  0.76428571  0.66190476  0.66296296  0.86000000  0.74444444  0.75151515 
    ##         BOB         FRB         UEB         SCL         SCD         WSL 
    ## -0.02835821  0.25000000  0.95849057  0.75625000  0.80416667  0.79230769 
    ##         BUR         SCS         NEU         UST         KON         FRI 
    ##  0.85312500  0.77777778  0.85000000  0.67391304  0.87000000  0.76500000 
    ##         BAR         BRU         LIM         HAN         GEN 
    ##  0.00000000 -0.01666667  0.42083333  0.63846154  0.76000000 
    ## 
    ## $pop.score$`20`
    ##         FRE         ALD         RIE         TRO         SCG         BOH 
    ## 0.592857143 0.766666667 0.735185185 0.920000000 0.711111111 0.727272727 
    ##         BOB         FRB         UEB         SCL         SCD         WSL 
    ## 0.009452736 0.360000000 0.924528302 0.708333333 0.845833333 0.792307692 
    ##         BUR         SCS         NEU         UST         KON         FRI 
    ## 0.837500000 0.744444444 0.921428571 0.734782609 0.870000000 0.795000000 
    ##         BAR         BRU         LIM         HAN         GEN 
    ## 0.000000000 0.205555556 0.316666667 0.930769231 0.620000000 
    ## 
    ## $pop.score$`25`
    ##        FRE        ALD        RIE        TRO        SCG        BOH        BOB 
    ## 0.58571429 0.80000000 0.66481481 0.90000000 0.61111111 0.71818182 0.09353234 
    ##        FRB        UEB        SCL        SCD        WSL        BUR        SCS 
    ## 0.79000000 0.90943396 0.71875000 0.81041667 0.71538462 0.81250000 0.77777778 
    ##        NEU        UST        KON        FRI        BAR        BRU        LIM 
    ## 0.77142857 0.77391304 0.87500000 0.78000000 0.00000000 0.37222222 0.37083333 
    ##        HAN        GEN 
    ## 0.81538462 0.54000000 
    ## 
    ## $pop.score$`30`
    ##       FRE       ALD       RIE       TRO       SCG       BOH       BOB       FRB 
    ## 0.5285714 0.8523810 0.6592593 0.7700000 0.7555556 0.7030303 0.1402985 0.7100000 
    ##       UEB       SCL       SCD       WSL       BUR       SCS       NEU       UST 
    ## 0.8924528 0.6916667 0.7875000 0.6538462 0.7781250 0.7444444 0.7500000 0.7565217 
    ##       KON       FRI       BAR       BRU       LIM       HAN       GEN 
    ## 0.8550000 0.7600000 0.0000000 0.5222222 0.4500000 0.9076923 0.4800000 
    ## 
    ## $pop.score$`35`
    ##       FRE       ALD       RIE       TRO       SCG       BOH       BOB       FRB 
    ## 0.5285714 0.8333333 0.6388889 0.8000000 0.7444444 0.6969697 0.2004975 0.5900000 
    ##       UEB       SCL       SCD       WSL       BUR       SCS       NEU       UST 
    ## 0.8792453 0.6229167 0.7895833 0.6615385 0.7843750 0.6333333 0.7642857 0.6434783 
    ##       KON       FRI       BAR       BRU       LIM       HAN       GEN 
    ## 0.8700000 0.7100000 0.0000000 0.4888889 0.4250000 0.8538462 0.3200000 
    ## 
    ## $pop.score$`40`
    ##       FRE       ALD       RIE       TRO       SCG       BOH       BOB       FRB 
    ## 0.6071429 0.5714286 0.6351852 0.8100000 0.5888889 0.6787879 0.2577114 0.5200000 
    ##       UEB       SCL       SCD       WSL       BUR       SCS       NEU       UST 
    ## 0.8622642 0.6729167 0.7937500 0.6846154 0.7687500 0.5222222 0.7642857 0.6260870 
    ##       KON       FRI       BAR       BRU       LIM       HAN       GEN 
    ## 0.8550000 0.5750000 0.0000000 0.5055556 0.4333333 0.7846154 0.2600000 
    ## 
    ## $pop.score$`45`
    ##       FRE       ALD       RIE       TRO       SCG       BOH       BOB       FRB 
    ## 0.6928571 0.4714286 0.6574074 0.7600000 0.5888889 0.5969697 0.2368159 0.4900000 
    ##       UEB       SCL       SCD       WSL       BUR       SCS       NEU       UST 
    ## 0.8528302 0.6854167 0.7354167 0.6153846 0.7843750 0.4888889 0.6142857 0.6347826 
    ##       KON       FRI       BAR       BRU       LIM       HAN       GEN 
    ## 0.8283333 0.5800000 0.0000000 0.4888889 0.4291667 0.5153846 0.5800000 
    ## 
    ## $pop.score$`50`
    ##       FRE       ALD       RIE       TRO       SCG       BOH       BOB       FRB 
    ## 0.6428571 0.4571429 0.5870370 0.6400000 0.6000000 0.6151515 0.2990050 0.5400000 
    ##       UEB       SCL       SCD       WSL       BUR       SCS       NEU       UST 
    ## 0.8358491 0.6875000 0.7645833 0.6000000 0.7312500 0.4444444 0.7285714 0.6695652 
    ##       KON       FRI       BAR       BRU       LIM       HAN       GEN 
    ## 0.8216667 0.6250000 0.0000000 0.4555556 0.3833333 0.5000000 0.5400000 
    ## 
    ## 
    ## $mean
    ##          1          5         10         15         20         25         30 
    ## 0.05516256 0.38720095 0.61453175 0.63349625 0.65520414 0.66114779 0.65863336 
    ##         35         40         45         50 
    ## 0.62953028 0.59902349 0.57945746 0.57254402 
    ## 
    ## $pred
    ## $pred$x
    ##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
    ## [26] 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50
    ## 
    ## $pred$y
    ##  [1] 0.05516256 0.14408033 0.23063859 0.31246843 0.38720095 0.45281156
    ##  [7] 0.50865294 0.55442209 0.58981603 0.61453175 0.62879558 0.63495107
    ## [13] 0.63587110 0.63442853 0.63349625 0.63529244 0.63941656 0.64481336
    ## [19] 0.65042764 0.65520414 0.65834014 0.66004282 0.66077187 0.66098696
    ## [25] 0.66114779 0.66156528 0.66195539 0.66188531 0.66092223 0.65863336
    ## [31] 0.65472491 0.64945912 0.64323725 0.63646055 0.62953028 0.62278692
    ## [37] 0.61632791 0.61018990 0.60440954 0.59902349 0.59406974 0.58959161
    ## [43] 0.58563375 0.58224082 0.57945746 0.57729870 0.57566107 0.57441148
    ## [49] 0.57341683 0.57254402
    ## 
    ## 
    ## $best
    ## [1] 27

``` r
#optimal score is 23 PCAs
```

``` r
# Idea: to find the optimal number of clusters, we use k-means, which will maximise the variance between groups (B(X)).
# Therefore we run k-means sequentially with rising numbers of k and compare it with the Bayesian Information Criterion BIC where we will choose the clustering solution corresponding to the lowest BIC.
# First, a PCA is performed but no information is lost due to keeping the principal components and therefore all variation in the original data. The analysis will go faster if we reduce the number of PCs to be retained, though.
# Clusters can be identified using find.clusters() which will first run the PCA and then the kmeans()

# Route 1: uncorrected AllMarkersOnly
setPop(myData_genind_allMarkersOnly) <- ~Pop
grp <- find.clusters(myData_genind_allMarkersOnly, max.n.clust=6, max.n.pca=40, n.pca=21)
```

![](PCA-and-DAPC_files/figure-gfm/dapc%20define%20clusters-1.png)<!-- -->

    ## Choose the number of clusters (>=2):

``` r
#no clear model exists, the more the merrier but it's not very useful / more than 30 PCAs retained don't give a lot more information so we keep it at 30 for the moment
table(pop(myData_genind_allMarkersOnly), grp$grp)
```

    ##      
    ##         1
    ##   FRE 214
    ##   ALD  62
    ##   RIE  62
    ##   TRO  34
    ##   SCG  15
    ##   BOH  40
    ##   BOB 429
    ##   FRB  13
    ##   UEB 153
    ##   SCL  99
    ##   SCD 100
    ##   WSL 283
    ##   BUR 223
    ##   SCS  43
    ##   NEU  32
    ##   UST  35
    ##   KON 212
    ##   FRI 120
    ##   BAR   1
    ##   BRU 178
    ##   LIM 216
    ##   HAN  65
    ##   GEN  10

``` r
table.value(table(pop(myData_genind_allMarkersOnly), grp$grp), col.lab=paste("inf", 1:6),
row.lab=paste("ori", 1:6))
```

![](PCA-and-DAPC_files/figure-gfm/dapc%20define%20clusters-2.png)<!-- -->

``` r
#this looks horrible and I don't understand it at all

# Route 2: corrected AllMarkersOnly
setPop(cc_myData_genind_allMarkersOnly_SY) <- ~Pop
cc_grp <- find.clusters(cc_myData_genind_allMarkersOnly_SY, max.n.clust=6, max.n.pca=40, n.pca=23)
```

![](PCA-and-DAPC_files/figure-gfm/dapc%20define%20clusters-3.png)<!-- -->

    ## Choose the number of clusters (>=2):

``` r
#no clear model exists, the more the merrier but it's not very useful
table(pop(cc_myData_genind_allMarkersOnly_SY), cc_grp$grp)
```

    ##      
    ##         1
    ##   FRE  14
    ##   ALD  21
    ##   RIE  54
    ##   TRO  10
    ##   SCG   9
    ##   BOH  33
    ##   BOB 201
    ##   FRB  10
    ##   UEB  53
    ##   SCL  48
    ##   SCD  48
    ##   WSL  13
    ##   BUR  32
    ##   SCS   9
    ##   NEU  14
    ##   UST  23
    ##   KON  60
    ##   FRI  20
    ##   BAR   1
    ##   BRU  18
    ##   LIM  24
    ##   HAN  13
    ##   GEN   5

``` r
table.value(table(pop(cc_myData_genind_allMarkersOnly_SY), cc_grp$grp), col.lab=paste("inf", 1:6),
row.lab=paste("ori", 1:6))
```

![](PCA-and-DAPC_files/figure-gfm/dapc%20define%20clusters-4.png)<!-- -->

DAPC aims to provide an efficient description of genetic clusters using
a few synthetic variables. These are constructed as linear combinations
of the original variables (alleles) which have the largest between-group
variance and the smallest within-group variance. Coefficients of the
alleles used in the linear combination are called loadings, while the
synthetic variables are themselves referred to as discriminant
functions. While these are different from the admixture coefficients of
software like STRUCTURE, they can still be interpreted as proximities of
individuals to the different clusters. Membership probabilities also
provide indications of how clear-cut genetic clusters are. Loose
clusters will result in fairly flat distributions of membership
probabilities of individuals across clusters, pointing to possible
admixture. Lastly, using the allele loadings, it is possible to
represent new individuals (which have not participated to the analysis)
onto the factorial planes, and derive membership probabilities as welll.
Such individuals are referred to as supplementary individuals.

``` r
#uncorrected AllMarkersOnly, only visualising the DAPC only preset clusters
setPop(myData_genind_allMarkersOnly) <- ~Pop
dapc.allMarkers <- dapc(myData_genind_allMarkersOnly,var.contrib = TRUE, scale = FALSE, n.pca = 21, n.da = nPop(myData_genind_allMarkersOnly) - 1)
scatter(dapc.allMarkers, cell = 0, pch = 18:23, cstar = 0, mstree = TRUE, lwd = 2, lty = 2, scree.pca=TRUE,posi.pca="bottomleft", scree.da=TRUE, posi.da="top")
```

![](PCA-and-DAPC_files/figure-gfm/dapc-1.png)<!-- -->

``` r
#clonecorrected  AllMarkersOnly
setPop(cc_myData_genind_allMarkersOnly_SY) <- ~Pop
cc_dapc.allMarkers <- dapc(cc_myData_genind_allMarkersOnly_SY, var.contrib = TRUE, scale = FALSE, n.pca = 23, n.da = nPop(cc_myData_genind_allMarkersOnly_SY) - 1)
scatter(cc_dapc.allMarkers, cell = 0, pch = 18:23, cstar = 0, mstree = TRUE, lwd = 2, lty = 2, scree.pca=TRUE,posi.pca="bottomleft", scree.da=TRUE, posi.da="top")
```

![](PCA-and-DAPC_files/figure-gfm/dapc3-1.png)<!-- -->
