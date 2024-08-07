PCA
================
Lia Baumann
2024-06-04

## PCA

from
<https://www.molecularecologist.com/2015/07/30/pca-of-multilocus-genotypes-in-r/>

``` r
library(scatterplot3d)
library(rgl)
library(tidyverse)
library(readxl)
library(knitr)
T_all <- read_excel("Tuaest_Master_genMonit_haploid_MP_Mai2024.xlsx")
cc_FullTable_Tuaest_allMarkersOnly <- read.csv("C:/Users/liaba/OneDrive - Eidg. Forschungsanstalt WSL/R/truffles/cc_TUAEST_ALLMARKERSONLY_HWE.csv")
T_all <- separate_wider_delim(T_all,Code_Analyses2, delim="_", names = c("Code","Yr"))
T_all_withMLGs_allMarkersOnly <- inner_join(T_all,cc_FullTable_Tuaest_allMarkersOnly,by=c("Code"="Sample")) %>%
  select(Code, MLG, Site_1_abrev, truffle_year, Sampling_year, Sampling_date, aest06_1.x, aest07_1.x, aest15_1.x, aest26_1.x, aest28_1.x, aest35_1.x, aest36_1.x, aest01_1.x, aest10_1.x, aest18_1.x, aest24_1.x, aest25_1.x, aest29_1.x, aest31_1.x)

T_all_allMarkersOnly.pca<-prcomp(T_all_withMLGs_allMarkersOnly[,7:20], center=TRUE,scale.=TRUE)
pcs123<-T_all_allMarkersOnly.pca$x[,1:3]
WSL<-grep("WSL",T_all_withMLGs_allMarkersOnly$Site_1_abrev)
BUR<-grep("BUR",T_all_withMLGs_allMarkersOnly$Site_1_abrev)
BRU<-grep("BRU",T_all_withMLGs_allMarkersOnly$Site_1_abrev)
SCG<-grep("SCG",T_all_withMLGs_allMarkersOnly$Site_1_abrev)

p<-scatterplot3d(pcs123[WSL,1],pcs123[WSL,2],pcs123[WSL,3],color=c("red"),
                 pch=20,xlab="PC1",ylab="PC2",zlab="PC3",
                 xlim=c(min(pcs123[,1]),max(pcs123[,1])),
                 ylim=c(min(pcs123[,2]),max(pcs123[,2])),
                 zlim=c(min(pcs123[,3]),max(pcs123[,3])))

p$points3d(pcs123[BUR,1],pcs123[BUR,2], pcs123[BUR,3],col=c("blue"),pch=20)
p$points3d(pcs123[BRU,1],pcs123[BRU,2], pcs123[BRU,3],col=c("maroon"),pch=20)
p$points3d(pcs123[SCG,1],pcs123[SCG,2], pcs123[SCG,3],col=c("green"),pch=20)
legend(p$xyz.convert(-80, 20, -10), col= c("red","blue","maroon","green"), bg="white", pch=c(20,20,20,20), yjust=0,
  legend = c("WSL", "BUR", "BRU", "SCG"), cex = 1.1)
```

![](PCA_files/figure-gfm/pca%20test%201-1.png)<!-- -->
