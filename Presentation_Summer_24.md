Presentation2024
================
Lia Baumann
2024-07-10

Readme:

SY = Sampling Year TY = Truffle Year

Contents:

1.  Nei’s genetic diversity Hexp

- allMarkersOnly (SY and TY)
- rar = rarefied (SY and TY)
- cc = clone-corrected (SY and TY)

2.  Simpson’s diversity

3.  PCA

- allMarkersOnly
- cc = clone-corrected (SY and TY)

## 1. Nei’s genetic diversity

#### myData_genind_allMarkersOnly locus table mean over all data

``` r
setPop(myData_genind_allMarkersOnly) <- ~Pop/SamplingYear
locus_table(myData_genind_allMarkersOnly, lev="genotype")
```

    ##           summary
    ## locus      genotype   1-D  Hexp Evenness
    ##   aest06_1     7.00  0.67  0.67     0.66
    ##   aest07_1     6.00  0.64  0.64     0.74
    ##   aest15_1     3.00  0.36  0.36     0.60
    ##   aest26_1    11.00  0.77  0.77     0.86
    ##   aest28_1    21.00  0.82  0.82     0.79
    ##   aest35_1     7.00  0.47  0.47     0.64
    ##   aest36_1     7.00  0.57  0.57     0.80
    ##   aest01_1     7.00  0.74  0.74     0.79
    ##   aest10_1     7.00  0.77  0.77     0.87
    ##   aest18_1     5.00  0.73  0.73     0.89
    ##   aest24_1     7.00  0.55  0.55     0.53
    ##   aest25_1     5.00  0.57  0.57     0.73
    ##   aest29_1     9.00  0.69  0.69     0.78
    ##   aest31_1     7.00  0.67  0.67     0.80
    ##   mean         7.79  0.64  0.64     0.75

#### myData_genind_allMarkersOnly locus table per Population, BAR and SCHIF removed (Strata: Population / SamplingYear)

``` r
setPop(myData_genind_allMarkersOnly) <- ~Pop
myData_genind_allMarkersOnly_allButBAR_SCHIF <- popsub(myData_genind_allMarkersOnly, exclude=c("BAR","SCHIF"))
setPop(myData_genind_allMarkersOnly_allButBAR_SCHIF) <- ~Pop/SamplingYear

# create a data frame and transpose the result;
# use sapply to iterate over all populations calculated by seppop.
# then use the locus_table function on the level Genotype

locus_table_SY_allButBAR_SCHIF <- data.frame(t(
  sapply(seppop(myData_genind_allMarkersOnly_allButBAR_SCHIF),
    function(ls) poppr::locus_table(ls,lev="genotype"))))

# choose only the columns corresponding to Hexp values and rename them with their loci

locus_table_SY_hexp_allButBAR_SCHIF <- 
  select(locus_table_SY_allButBAR_SCHIF, X31:X45) %>%
  rename(aest06_1=X31,aest07_1=X32, aest15_1=X33, aest26_1=X34,aest28_1=X35, aest35_1=X36, aest36_1=X37, aest01_1=X38, aest10_1=X39, aest18_1=X40, aest24_1=X41, aest25_1=X42, aest29_1=X43, aest31_1=X44, mean=X45) %>%
  tibble::rownames_to_column("Pop") %>%
  separate(.,Pop,sep="_", into=c("Pop","SamplingYear"))
                                      
#pivot to have Hexp values in columns as value to use in ggplot later
locus_table_SY_hexp_allButBAR_SCHIF_pivot <- pivot_longer(locus_table_SY_hexp_allButBAR_SCHIF,cols=3:16,names_to="locus",values_to="value")

ggplot(locus_table_SY_hexp_allButBAR_SCHIF_pivot,aes(x=as.numeric(SamplingYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  facet_wrap(vars(Pop)) +
  labs(y="Nei's genetic diversity Hexp", x="Sampling year") +
  ggtitle("Nei's genetic diversity per locus, sampling year and population") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines"))
```

![](Presentation_Summer_24_files/figure-gfm/Hexp%20table%20Pop%20SamplingYear%20without%20empty%20facets-1.png)<!-- -->

``` r
#create high resolution graph
library(devEMF)
emf(file= "locus_table_SY_hexp_allButBAR_SCHIF_pivot.emf") # Opens a device
ggplot(locus_table_SY_hexp_allButBAR_SCHIF_pivot,aes(x=as.numeric(SamplingYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  facet_wrap(vars(Pop)) +
  labs(y="Nei's genetic diversity Hexp", x="Sampling year") +
  ggtitle("Nei's genetic diversity per locus, sampling year and population") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines"))
dev.off()  # Close the device
```

    ## png 
    ##   2

``` r
#create zoom-in to BOB
locus_table_SY_hexp_all_BOB <- locus_table_SY_hexp_allButBAR_SCHIF_pivot %>%
  filter(Pop=="BOB") %>%
  ggplot(aes(x=as.double(SamplingYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  labs(y="Nei's genetic diversity Hexp", x="Sampling year") +
  ggtitle("Nei's genetic diversity per locus, sampling year and population BOB") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines")) +
  scale_x_continuous(breaks=~axisTicks(., log = FALSE)) +
  ylim(0,1)
locus_table_SY_hexp_all_BOB
```

![](Presentation_Summer_24_files/figure-gfm/Hexp%20table%20Pop%20SamplingYear%20without%20empty%20facets-2.png)<!-- -->

``` r
emf(file="locus_table_SY_hexp_all_BOB.emf")
locus_table_SY_hexp_all_BOB
dev.off()
```

    ## png 
    ##   2

``` r
#create zoom-in to BUR
locus_table_SY_hexp_all_BUR <- locus_table_SY_hexp_allButBAR_SCHIF_pivot %>%
  filter(Pop=="BUR") %>%
  ggplot(aes(x=as.double(SamplingYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  labs(y="Nei's genetic diversity Hexp", x="Sampling year") +
  ggtitle("Nei's genetic diversity per locus, sampling year and population BUR") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines")) +
  scale_x_continuous(breaks=~axisTicks(., log = FALSE))+
  ylim(0,1)
locus_table_SY_hexp_all_BUR
```

![](Presentation_Summer_24_files/figure-gfm/Hexp%20table%20Pop%20SamplingYear%20without%20empty%20facets-3.png)<!-- -->

``` r
emf(file="locus_table_SY_hexp_all_BUR.emf")
locus_table_SY_hexp_all_BUR
dev.off()
```

    ## png 
    ##   2

``` r
#create zoom-in to FRE
locus_table_SY_hexp_all_FRE <- locus_table_SY_hexp_allButBAR_SCHIF_pivot %>%
  filter(Pop=="FRE") %>%
  ggplot(aes(x=as.double(SamplingYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  labs(y="Nei's genetic diversity Hexp", x="Sampling year") +
  ggtitle("Nei's genetic diversity per locus, sampling year and population FRE") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines")) +
  scale_x_continuous(breaks=~axisTicks(., log = FALSE))+
  ylim(0,1)
locus_table_SY_hexp_all_FRE
```

![](Presentation_Summer_24_files/figure-gfm/Hexp%20table%20Pop%20SamplingYear%20without%20empty%20facets-4.png)<!-- -->

``` r
emf(file="locus_table_SY_hexp_all_FRE.emf")
locus_table_SY_hexp_all_FRE
dev.off()
```

    ## png 
    ##   2

``` r
#create zoom-in to KON
locus_table_SY_hexp_all_KON <- locus_table_SY_hexp_allButBAR_SCHIF_pivot %>%
  filter(Pop=="KON") %>%
  ggplot(aes(x=as.double(SamplingYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  labs(y="Nei's genetic diversity Hexp", x="Sampling year") +
  ggtitle("Nei's genetic diversity per locus, sampling year and population KON") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines")) +
  scale_x_continuous(breaks=~axisTicks(., log = FALSE))+
  ylim(0,1)
locus_table_SY_hexp_all_KON
```

![](Presentation_Summer_24_files/figure-gfm/Hexp%20table%20Pop%20SamplingYear%20without%20empty%20facets-5.png)<!-- -->

``` r
emf(file="locus_table_SY_hexp_all_KON.emf")
locus_table_SY_hexp_all_KON
dev.off()
```

    ## png 
    ##   2

``` r
#create zoom-in to UEB
locus_table_SY_hexp_all_UEB <- locus_table_SY_hexp_allButBAR_SCHIF_pivot %>%
  filter(Pop=="UEB") %>%
  ggplot(aes(x=as.double(SamplingYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  labs(y="Nei's genetic diversity Hexp", x="Sampling year") +
  ggtitle("Nei's genetic diversity per locus, sampling year and population UEB") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines")) +
  scale_x_continuous(breaks=~axisTicks(., log = FALSE))+
  ylim(0,1)
locus_table_SY_hexp_all_UEB
```

![](Presentation_Summer_24_files/figure-gfm/Hexp%20table%20Pop%20SamplingYear%20without%20empty%20facets-6.png)<!-- -->

``` r
emf(file="locus_table_SY_hexp_all_UEB.emf")
locus_table_SY_hexp_all_UEB
dev.off()
```

    ## png 
    ##   2

``` r
n_samples_allMarkersOnly_allButBAR_SCHIF_SY <- poppr(myData_genind_allMarkersOnly_allButBAR_SCHIF) %>%
  separate(Pop, sep="_", into=c("Pop","SamplingYear")) %>%
  select(1:4) %>%
  pivot_longer(cols=c("N","MLG"),names_to="Number",values_to="Count")

n_samples_allMarkersOnly_allButBAR_SCHIF_SY %>%
  filter(Pop != "Total") %>%
  ggplot(aes(x=as.numeric(SamplingYear), y=Count, fill=Number)) +
  geom_col(stat="identity", position="dodge") +
  facet_wrap(vars(Pop)) +
  labs(y="", x="Sampling year", subtitle="Dataset: only samples with all markers") +
  ggtitle("Number of samples vs. number of MLGs") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines"))
```

![](Presentation_Summer_24_files/figure-gfm/Hexp%20table%20Pop%20SamplingYear%20without%20empty%20facets-7.png)<!-- -->

``` r
emf(file="n_samples_allMarkersOnly_allButBAR_SCHIF_SY.emf")
n_samples_allMarkersOnly_allButBAR_SCHIF_SY %>%
  filter(Pop != "Total") %>%
  ggplot(aes(x=as.numeric(SamplingYear), y=Count, fill=Number)) +
  geom_col(stat="identity", position="dodge") +
  facet_wrap(vars(Pop)) +
  labs(y="", x="Sampling year", subtitle="Dataset: only samples with all markers") +
  ggtitle("Number of samples vs. number of MLGs") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines"))
dev.off()
```

    ## png 
    ##   2

#### myData_genind_allMarkersOnly locus table per Population, BAR and SCHIF removed (Strata: Population / TruffleYear)

``` r
setPop(myData_genind_allMarkersOnly_allButBAR_SCHIF) <- ~Pop/TruffleYear

# create a data frame and transpose the result;
# use sapply to iterate over all populations calculated by seppop.
# then use the locus_table function on the level Genotype

locus_table_TY_allButBAR_SCHIF <- data.frame(t(
  sapply(seppop(myData_genind_allMarkersOnly_allButBAR_SCHIF),
    function(ls) poppr::locus_table(ls,lev="genotype"))))

# choose only the columns corresponding to Hexp values and rename them with their loci

locus_table_TY_hexp_allButBAR_SCHIF <- 
  select(locus_table_TY_allButBAR_SCHIF, X31:X45) %>%
  rename(aest06_1=X31,aest07_1=X32, aest15_1=X33, aest26_1=X34,aest28_1=X35, aest35_1=X36, aest36_1=X37, aest01_1=X38, aest10_1=X39, aest18_1=X40, aest24_1=X41, aest25_1=X42, aest29_1=X43, aest31_1=X44, mean=X45) %>%
  tibble::rownames_to_column("Pop") %>%
  separate(.,Pop,sep="_", into=c("Pop","TruffleYear"))
                                      
#pivot to have Hexp values in columns as value to use in ggplot later
locus_table_TY_hexp_allButBAR_SCHIF_pivot <- pivot_longer(locus_table_TY_hexp_allButBAR_SCHIF,cols=3:16,names_to="locus",values_to="value")

ggplot(locus_table_TY_hexp_allButBAR_SCHIF_pivot,aes(x=as.numeric(TruffleYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  facet_wrap(vars(Pop)) +
  labs(y="Nei's genetic diversity Hexp", x="Truffle year") +
  ggtitle("Nei's genetic diversity per locus, truffle year and population") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines"))
```

![](Presentation_Summer_24_files/figure-gfm/Hexp%20table%20Pop%20TruffleYear%20without%20empty%20facets-1.png)<!-- -->

``` r
#create high resolution graph
library(devEMF)
emf(file= "locus_table_TY_hexp_allButBAR_SCHIF_pivot.emf") # Opens a device
ggplot(locus_table_TY_hexp_allButBAR_SCHIF_pivot,aes(x=as.numeric(TruffleYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  facet_wrap(vars(Pop)) +
  labs(y="Nei's genetic diversity Hexp", x="Truffle year") +
  ggtitle("Nei's genetic diversity per locus, truffle year and population") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines"))
dev.off()  # Close the device
```

    ## png 
    ##   2

``` r
n_samples_allMarkersOnly_allButBAR_SCHIF_TY <- poppr(myData_genind_allMarkersOnly_allButBAR_SCHIF) %>%
  separate(Pop, sep="_", into=c("Pop","TruffleYear")) %>%
  select(1:4) %>%
  pivot_longer(cols=c("N","MLG"),names_to="Number",values_to="Count")

n_samples_allMarkersOnly_allButBAR_SCHIF_TY %>%
  filter(Pop != "Total") %>%
  ggplot(aes(x=as.numeric(TruffleYear), y=Count, fill=Number)) +
  geom_col(stat="identity", position="dodge") +
  facet_wrap(vars(Pop)) +
  labs(y="", x="Truffle year", subtitle="Dataset: only samples with all markers") +
  ggtitle("Number of samples vs. number of MLGs") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines"))
```

![](Presentation_Summer_24_files/figure-gfm/Hexp%20table%20Pop%20TruffleYear%20without%20empty%20facets-2.png)<!-- -->

``` r
emf(file="n_samples_allMarkersOnly_allButBAR_SCHIF_TY.emf")
n_samples_allMarkersOnly_allButBAR_SCHIF_TY %>%
  filter(Pop != "Total") %>%
  ggplot(aes(x=as.numeric(TruffleYear), y=Count, fill=Number)) +
  geom_col(stat="identity", position="dodge") +
  facet_wrap(vars(Pop)) +
  labs(y="", x="Truffle year", subtitle="Dataset: only samples with all markers") +
  ggtitle("Number of samples vs. number of MLGs") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines"))
dev.off()
```

    ## png 
    ##   2

``` r
#create zoom-in to BOB
locus_table_TY_hexp_all_BOB <- locus_table_TY_hexp_allButBAR_SCHIF_pivot %>%
  filter(Pop=="BOB") %>%
  ggplot(aes(x=as.double(TruffleYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  labs(y="Nei's genetic diversity Hexp", x="Truffle year") +
  ggtitle("Nei's genetic diversity per locus, Truffle year and population BOB") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines")) +
  scale_x_continuous(breaks=~axisTicks(., log = FALSE)) +
  ylim(0,1)
locus_table_TY_hexp_all_BOB
```

![](Presentation_Summer_24_files/figure-gfm/Hexp%20table%20Pop%20TruffleYear%20without%20empty%20facets-3.png)<!-- -->

``` r
emf(file="locus_table_TY_hexp_all_BOB.emf")
locus_table_TY_hexp_all_BOB
dev.off()
```

    ## png 
    ##   2

``` r
#create zoom-in to BUR
locus_table_TY_hexp_all_BUR <- locus_table_TY_hexp_allButBAR_SCHIF_pivot %>%
  filter(Pop=="BUR") %>%
  ggplot(aes(x=as.double(TruffleYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  labs(y="Nei's genetic diversity Hexp", x="Truffle year") +
  ggtitle("Nei's genetic diversity per locus, truffle year and population BUR") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines")) +
  scale_x_continuous(breaks=~axisTicks(., log = FALSE))+
  ylim(0,1)
locus_table_TY_hexp_all_BUR
```

![](Presentation_Summer_24_files/figure-gfm/Hexp%20table%20Pop%20TruffleYear%20without%20empty%20facets-4.png)<!-- -->

``` r
emf(file="locus_table_TY_hexp_all_BUR.emf")
locus_table_TY_hexp_all_BUR
dev.off()
```

    ## png 
    ##   2

``` r
#create zoom-in to FRE
locus_table_TY_hexp_all_FRE <- locus_table_TY_hexp_allButBAR_SCHIF_pivot %>%
  filter(Pop=="FRE") %>%
  ggplot(aes(x=as.double(TruffleYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  labs(y="Nei's genetic diversity Hexp", x="Truffle year") +
  ggtitle("Nei's genetic diversity per locus, truffle year and population FRE") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines")) +
  scale_x_continuous(breaks=~axisTicks(., log = FALSE))+
  ylim(0,1)
locus_table_TY_hexp_all_FRE
```

![](Presentation_Summer_24_files/figure-gfm/Hexp%20table%20Pop%20TruffleYear%20without%20empty%20facets-5.png)<!-- -->

``` r
emf(file="locus_table_TY_hexp_all_FRE.emf")
locus_table_TY_hexp_all_FRE
dev.off()
```

    ## png 
    ##   2

``` r
#create zoom-in to KON
locus_table_TY_hexp_all_KON <- locus_table_TY_hexp_allButBAR_SCHIF_pivot %>%
  filter(Pop=="KON") %>%
  ggplot(aes(x=as.double(TruffleYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  labs(y="Nei's genetic diversity Hexp", x="Truffle year") +
  ggtitle("Nei's genetic diversity per locus, truffle year and population KON") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines")) +
  scale_x_continuous(breaks=~axisTicks(., log = FALSE))+
  ylim(0,1)
locus_table_TY_hexp_all_KON
```

![](Presentation_Summer_24_files/figure-gfm/Hexp%20table%20Pop%20TruffleYear%20without%20empty%20facets-6.png)<!-- -->

``` r
emf(file="locus_table_TY_hexp_all_KON.emf")
locus_table_TY_hexp_all_KON
dev.off()
```

    ## png 
    ##   2

``` r
#create zoom-in to UEB
locus_table_TY_hexp_all_UEB <- locus_table_TY_hexp_allButBAR_SCHIF_pivot %>%
  filter(Pop=="UEB") %>%
  ggplot(aes(x=as.double(TruffleYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  labs(y="Nei's genetic diversity Hexp", x="Truffle year") +
  ggtitle("Nei's genetic diversity per locus, truffle year and population UEB") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines")) +
  scale_x_continuous(breaks=~axisTicks(., log = FALSE))+
  ylim(0,1)
locus_table_TY_hexp_all_UEB
```

![](Presentation_Summer_24_files/figure-gfm/Hexp%20table%20Pop%20TruffleYear%20without%20empty%20facets-7.png)<!-- -->

``` r
emf(file="locus_table_TY_hexp_all_UEB.emf")
locus_table_TY_hexp_all_UEB
dev.off()
```

    ## png 
    ##   2

### Rarefaction of Neis’ diversity index

Rarefied: Hexp \* (N/N-1). N= Anzahl samples (nicht Allele)

``` r
setPop(myData_genind_allMarkersOnly_allButBAR_SCHIF) <- ~Pop/SamplingYear
popdata_pop_year_rareNei_SY <- poppr(myData_genind_allMarkersOnly_allButBAR_SCHIF) %>%
  select(Pop, N)

# choose only the columns corresponding to Hexp values and rename them with their loci
locus_table_SY_hexp_N_allButBAR_SCHIF <- 
  select(locus_table_SY_allButBAR_SCHIF, X31:X45) %>%
  rename(aest06_1=X31,aest07_1=X32, aest15_1=X33, aest26_1=X34,aest28_1=X35, aest35_1=X36, aest36_1=X37, aest01_1=X38, aest10_1=X39, aest18_1=X40, aest24_1=X41, aest25_1=X42, aest29_1=X43, aest31_1=X44, mean=X45) %>%
    tibble::rownames_to_column("Pop") %>%
  left_join(.,popdata_pop_year_rareNei_SY,by="Pop") %>%
  separate(.,Pop,sep="_", into=c("Pop","SamplingYear"))

N_minus_SY <- locus_table_SY_hexp_N_allButBAR_SCHIF$N #number of samples
#calculate rarefaction over all the values (columns 3:17)
rarHexpSY_minus <- (N_minus_SY/(N_minus_SY-1))*select(locus_table_SY_hexp_N_allButBAR_SCHIF,3:17)
#add the pop and year again
rarHexpSY_minus <- bind_cols(rarHexpSY_minus,Pop = locus_table_SY_hexp_N_allButBAR_SCHIF$Pop, SamplingYear = locus_table_SY_hexp_N_allButBAR_SCHIF$SamplingYear)

rarHexpSY_minus_pivot <- pivot_longer(locus_table_SY_hexp_N_allButBAR_SCHIF,cols=3:16,names_to="locus",values_to="value")

ggplot(rarHexpSY_minus_pivot,aes(x=as.numeric(SamplingYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  facet_wrap(vars(Pop)) +
  labs(y="Nei's genetic diversity rarefied Hexp", x="Sampling year") +
  ggtitle("Nei's rarefied genetic diversity per locus, sampling year and population") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines"))
```

![](Presentation_Summer_24_files/figure-gfm/rarefy%20Nei%20SamplingYear%20without%20BAR%20and%20SCHIF-1.png)<!-- -->

``` r
#create high resolution graph
library(devEMF)
emf(file= "rar_locus_table_SY_hexp_allButBAR_SCHIF.emf") # Opens a device
ggplot(rarHexpSY_minus_pivot,aes(x=as.numeric(SamplingYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  facet_wrap(vars(Pop)) +
  labs(y="Nei's genetic diversity rarefied Hexp", x="Sampling year") +
  ggtitle("Nei's rarefied genetic diversity per locus, sampling year and population") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines"))
dev.off()  # Close the device
```

    ## png 
    ##   2

``` r
#create zoom-in to BOB
rarlocus_table_SY_hexp_all_BOB <- rarHexpSY_minus_pivot %>%
  filter(Pop=="BOB") %>%
  ggplot(aes(x=as.double(SamplingYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  labs(y="Nei's rarefied genetic diversity Hexp", x="Sampling year") +
  ggtitle("Nei's rarefied genetic diversity per locus, sampling year and population BOB") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines")) +
  scale_x_continuous(breaks=~axisTicks(., log = FALSE)) +
  ylim(0,1)
rarlocus_table_SY_hexp_all_BOB
```

![](Presentation_Summer_24_files/figure-gfm/rarefy%20Nei%20SamplingYear%20without%20BAR%20and%20SCHIF-2.png)<!-- -->

``` r
emf(file="rarlocus_table_SY_hexp_all_BOB.emf")
rarlocus_table_SY_hexp_all_BOB
dev.off()
```

    ## png 
    ##   2

``` r
#create zoom-in to BUR
rarlocus_table_SY_hexp_all_BUR <- rarHexpSY_minus_pivot %>%
  filter(Pop=="BUR") %>%
  ggplot(aes(x=as.double(SamplingYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  labs(y="Nei's rarefied genetic diversity Hexp", x="Sampling year") +
  ggtitle("Nei's rarefied genetic diversity per locus, sampling year and population BUR") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines")) +
  scale_x_continuous(breaks=~axisTicks(., log = FALSE))+
  ylim(0,1)
rarlocus_table_SY_hexp_all_BUR
```

![](Presentation_Summer_24_files/figure-gfm/rarefy%20Nei%20SamplingYear%20without%20BAR%20and%20SCHIF-3.png)<!-- -->

``` r
emf(file="rarlocus_table_SY_hexp_all_BUR.emf")
rarlocus_table_SY_hexp_all_BUR
dev.off()
```

    ## png 
    ##   2

``` r
#create zoom-in to FRE
rarlocus_table_SY_hexp_all_FRE <- rarHexpSY_minus_pivot %>%
  filter(Pop=="FRE") %>%
  ggplot(aes(x=as.double(SamplingYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  labs(y="Nei's rarefied genetic diversity Hexp", x="Sampling year") +
  ggtitle("Nei's rarefied genetic diversity per locus, sampling year and population FRE") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines")) +
  scale_x_continuous(breaks=~axisTicks(., log = FALSE))+
  ylim(0,1)
rarlocus_table_SY_hexp_all_FRE
```

![](Presentation_Summer_24_files/figure-gfm/rarefy%20Nei%20SamplingYear%20without%20BAR%20and%20SCHIF-4.png)<!-- -->

``` r
emf(file="rarlocus_table_SY_hexp_all_FRE.emf")
rarlocus_table_SY_hexp_all_FRE
dev.off()
```

    ## png 
    ##   2

``` r
#create zoom-in to KON
rarlocus_table_SY_hexp_all_KON <- rarHexpSY_minus_pivot %>%
  filter(Pop=="KON") %>%
  ggplot(aes(x=as.double(SamplingYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  labs(y="Nei's rarefied genetic diversity Hexp", x="Sampling year") +
  ggtitle("Nei's rarefied genetic diversity per locus, sampling year and population KON") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines")) +
  scale_x_continuous(breaks=~axisTicks(., log = FALSE))+
  ylim(0,1)
rarlocus_table_SY_hexp_all_KON
```

![](Presentation_Summer_24_files/figure-gfm/rarefy%20Nei%20SamplingYear%20without%20BAR%20and%20SCHIF-5.png)<!-- -->

``` r
emf(file="rarlocus_table_SY_hexp_all_KON.emf")
rarlocus_table_SY_hexp_all_KON
dev.off()
```

    ## png 
    ##   2

``` r
#create zoom-in to UEB
rarlocus_table_SY_hexp_all_UEB <- rarHexpSY_minus_pivot %>%
  filter(Pop=="UEB") %>%
  ggplot(aes(x=as.double(SamplingYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  labs(y="Nei's rarefied genetic diversity Hexp", x="Sampling year") +
  ggtitle("Nei's rarefied genetic diversity per locus, sampling year and population UEB") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines")) +
  scale_x_continuous(breaks=~axisTicks(., log = FALSE))+
  ylim(0,1)
rarlocus_table_SY_hexp_all_UEB
```

![](Presentation_Summer_24_files/figure-gfm/rarefy%20Nei%20SamplingYear%20without%20BAR%20and%20SCHIF-6.png)<!-- -->

``` r
emf(file="rarlocus_table_SY_hexp_all_UEB.emf")
rarlocus_table_SY_hexp_all_UEB
dev.off()
```

    ## png 
    ##   2

``` r
setPop(myData_genind_allMarkersOnly_allButBAR_SCHIF) <- ~Pop/TruffleYear
popdata_pop_year_rareNei_TY <- poppr(myData_genind_allMarkersOnly_allButBAR_SCHIF) %>%
  select(Pop, N)

# choose only the columns corresponding to Hexp values and rename them with their loci
locus_table_TY_hexp_N_allButBAR_SCHIF <- 
  select(locus_table_TY_allButBAR_SCHIF, X31:X45) %>%
  rename(aest06_1=X31,aest07_1=X32, aest15_1=X33, aest26_1=X34,aest28_1=X35, aest35_1=X36, aest36_1=X37, aest01_1=X38, aest10_1=X39, aest18_1=X40, aest24_1=X41, aest25_1=X42, aest29_1=X43, aest31_1=X44, mean=X45) %>%
    tibble::rownames_to_column("Pop") %>%
  left_join(.,popdata_pop_year_rareNei_TY,by="Pop") %>%
  separate(.,Pop,sep="_", into=c("Pop","TruffleYear"))

N_minus_TY <- locus_table_TY_hexp_N_allButBAR_SCHIF$N #number of samples
#calculate rarefaction over all the values (columns 3:17)
rarHexpTY_minus <- (N_minus_TY/(N_minus_TY-1))*select(locus_table_TY_hexp_N_allButBAR_SCHIF,3:17)
#add the pop and year again
rarHexpTY_minus <- bind_cols(Pop = locus_table_TY_hexp_N_allButBAR_SCHIF$Pop, TruffleYear = locus_table_TY_hexp_N_allButBAR_SCHIF$TruffleYear, rarHexpTY_minus)

rarHexpTY_minus_pivot <- pivot_longer(rarHexpTY_minus,cols=3:16,names_to="locus",values_to="value")

ggplot(rarHexpTY_minus_pivot,aes(x=as.numeric(TruffleYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  facet_wrap(vars(Pop)) +
  labs(y="Nei's genetic diversity rarefied Hexp", x="Truffle year") +
  ggtitle("Nei's rarefied genetic diversity per locus, truffle year and population") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines"))
```

![](Presentation_Summer_24_files/figure-gfm/rarefy%20Nei%20TruffleYear%20without%20BAR%20and%20SCHIF-1.png)<!-- -->

``` r
#create high resolution graph
library(devEMF)
emf(file= "rar_locus_table_TY_hexp_pivot_allButBAR_SCHIF.emf") # Opens a device
ggplot(rarHexpTY_minus_pivot,aes(x=as.numeric(TruffleYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  facet_wrap(vars(Pop)) +
  labs(y="Nei's genetic diversity rarefied Hexp", x="Truffle year") +
  ggtitle("Nei's rarefied genetic diversity per locus, truffle year and population") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines"))
dev.off()  # Close the device
```

    ## png 
    ##   2

``` r
#create zoom-in to BOB
rarlocus_table_TY_hexp_all_BOB <- rarHexpTY_minus_pivot %>%
  filter(Pop=="BOB") %>%
  ggplot(aes(x=as.double(TruffleYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  labs(y="Nei's rarefied genetic diversity Hexp", x="Truffle year") +
  ggtitle("Nei's rarefied genetic diversity per locus, truffle year and population BOB") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines")) +
  scale_x_continuous(breaks=~axisTicks(., log = FALSE)) +
  ylim(0,1)
rarlocus_table_TY_hexp_all_BOB
```

![](Presentation_Summer_24_files/figure-gfm/rarefy%20Nei%20TruffleYear%20without%20BAR%20and%20SCHIF-2.png)<!-- -->

``` r
emf(file="rarlocus_table_TY_hexp_all_BOB.emf")
rarlocus_table_TY_hexp_all_BOB
dev.off()
```

    ## png 
    ##   2

``` r
#create zoom-in to BUR
rarlocus_table_TY_hexp_all_BUR <- rarHexpTY_minus_pivot %>%
  filter(Pop=="BUR") %>%
  ggplot(aes(x=as.double(TruffleYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  labs(y="Nei's rarefied genetic diversity Hexp", x="Truffle year") +
  ggtitle("Nei's rarefied genetic diversity per locus, truffle year and population BUR") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines")) +
  scale_x_continuous(breaks=~axisTicks(., log = FALSE))+
  ylim(0,1)
rarlocus_table_TY_hexp_all_BUR
```

![](Presentation_Summer_24_files/figure-gfm/rarefy%20Nei%20TruffleYear%20without%20BAR%20and%20SCHIF-3.png)<!-- -->

``` r
emf(file="rarlocus_table_TY_hexp_all_BUR.emf")
rarlocus_table_TY_hexp_all_BUR
dev.off()
```

    ## png 
    ##   2

``` r
#create zoom-in to FRE
rarlocus_table_TY_hexp_all_FRE <- rarHexpTY_minus_pivot %>%
  filter(Pop=="FRE") %>%
  ggplot(aes(x=as.double(TruffleYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  labs(y="Nei's rarefied genetic diversity Hexp", x="Truffle year") +
  ggtitle("Nei's rarefied genetic diversity per locus, truffle year and population FRE") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines")) +
  scale_x_continuous(breaks=~axisTicks(., log = FALSE))+
  ylim(0,1)
rarlocus_table_TY_hexp_all_FRE
```

![](Presentation_Summer_24_files/figure-gfm/rarefy%20Nei%20TruffleYear%20without%20BAR%20and%20SCHIF-4.png)<!-- -->

``` r
emf(file="rarlocus_table_TY_hexp_all_FRE.emf")
rarlocus_table_TY_hexp_all_FRE
dev.off()
```

    ## png 
    ##   2

``` r
#create zoom-in to KON
rarlocus_table_TY_hexp_all_KON <- rarHexpTY_minus_pivot %>%
  filter(Pop=="KON") %>%
  ggplot(aes(x=as.double(TruffleYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  labs(y="Nei's rarefied genetic diversity Hexp", x="Truffle year") +
  ggtitle("Nei's rarefied genetic diversity per locus, truffle year and population KON") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines")) +
  scale_x_continuous(breaks=~axisTicks(., log = FALSE))+
  ylim(0,1)
rarlocus_table_TY_hexp_all_KON
```

![](Presentation_Summer_24_files/figure-gfm/rarefy%20Nei%20TruffleYear%20without%20BAR%20and%20SCHIF-5.png)<!-- -->

``` r
emf(file="rarlocus_table_TY_hexp_all_KON.emf")
rarlocus_table_TY_hexp_all_KON
dev.off()
```

    ## png 
    ##   2

``` r
#create zoom-in to UEB
rarlocus_table_TY_hexp_all_UEB <- rarHexpTY_minus_pivot %>%
  filter(Pop=="UEB") %>%
  ggplot(aes(x=as.double(TruffleYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  labs(y="Nei's rarefied genetic diversity Hexp", x="Truffle year") +
  ggtitle("Nei's rarefied genetic diversity per locus, truffle year and population UEB") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines")) +
  scale_x_continuous(breaks=~axisTicks(., log = FALSE))+
  ylim(0,1)
rarlocus_table_TY_hexp_all_UEB
```

![](Presentation_Summer_24_files/figure-gfm/rarefy%20Nei%20TruffleYear%20without%20BAR%20and%20SCHIF-6.png)<!-- -->

``` r
emf(file="rarlocus_table_TY_hexp_all_UEB.emf")
rarlocus_table_TY_hexp_all_UEB
dev.off()
```

    ## png 
    ##   2

#### clone-corrected myData_genind_allMarkersOnly locus table per Population (Strata: Population / SamplingYear), removed BAR and SCHIF

``` r
setPop(cc_myData_genind_allMarkersOnly_SY) <- ~Pop
cc_myData_genind_allMarkersOnly_allButBAR_SCHIF_SY <- popsub(cc_myData_genind_allMarkersOnly_SY, exclude=c("BAR","SCHIF"))
setPop(cc_myData_genind_allMarkersOnly_allButBAR_SCHIF_SY) <- ~Pop/SamplingYear

# create a data frame and transpose the result;
# use sapply to iterate over all populations calculated by seppop.
# then use the locus_table function on the level Genotype

cc_locus_table_allButBAR_SCHIF_SY <- data.frame(t(
  sapply(seppop(cc_myData_genind_allMarkersOnly_allButBAR_SCHIF_SY),
    function(ls) poppr::locus_table(ls,lev="genotype"))))

# choose only the columns corresponding to Hexp values and rename them with their loci

cc_locus_table_SY_hexp_allButBAR_SCHIF <- 
  select(cc_locus_table_allButBAR_SCHIF_SY, X31:X45) %>%
  rename(aest06_1=X31,aest07_1=X32, aest15_1=X33, aest26_1=X34,aest28_1=X35, aest35_1=X36, aest36_1=X37, aest01_1=X38, aest10_1=X39, aest18_1=X40, aest24_1=X41, aest25_1=X42, aest29_1=X43, aest31_1=X44, mean=X45) %>%
  tibble::rownames_to_column("Pop") %>%
  separate(.,Pop,sep="_", into=c("Pop","SamplingYear"))
                                      
cc_locus_table_SY_hexp_allButBAR_SCHIF_pivot <- pivot_longer(cc_locus_table_SY_hexp_allButBAR_SCHIF,cols=3:16,names_to="locus",values_to="value")

ggplot(cc_locus_table_SY_hexp_allButBAR_SCHIF_pivot,aes(x=as.numeric(SamplingYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  facet_wrap(vars(Pop)) +
  labs(y="Nei's genetic diversity Hexp", x="Sampling year") +
  ggtitle("Nei's genetic diversity per locus, sampling year and population (clone-corrected)") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines"))
```

![](Presentation_Summer_24_files/figure-gfm/cc%20Hexp%20table%20Pop%20SamplingYear%20removed%20BAR%20and%20SCHIF-1.png)<!-- -->

``` r
emf(file= "cc_locus_table_SY_hexp_allButBAR_SCHIF_pivot.emf") # Opens a device
ggplot(cc_locus_table_SY_hexp_allButBAR_SCHIF_pivot,aes(x=as.numeric(SamplingYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  facet_wrap(vars(Pop)) +
  labs(y="Nei's genetic diversity Hexp", x="Sampling year") +
  ggtitle("Nei's genetic diversity per locus, sampling year and population (clone-corrected)") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines"))
dev.off()
```

    ## png 
    ##   2

``` r
#create zoom-in to BOB
cc_locus_table_SY_hexp_all_BOB <- cc_locus_table_SY_hexp_allButBAR_SCHIF_pivot %>%
  filter(Pop=="BOB") %>%
  ggplot(aes(x=as.double(SamplingYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  labs(y="Nei's genetic diversity Hexp", x="Sampling year") +
  ggtitle("Nei's genetic diversity per locus, Sampling year and population BOB (clone-corrected)") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines")) +
  scale_x_continuous(breaks=~axisTicks(., log = FALSE)) +
  ylim(0,1)
rarlocus_table_SY_hexp_all_BOB
```

![](Presentation_Summer_24_files/figure-gfm/cc%20Hexp%20table%20Pop%20SamplingYear%20removed%20BAR%20and%20SCHIF-2.png)<!-- -->

``` r
emf(file="cc_locus_table_SY_hexp_all_BOB.emf")
cc_locus_table_SY_hexp_all_BOB
dev.off()
```

    ## png 
    ##   2

``` r
#create zoom-in to BUR
cc_locus_table_SY_hexp_all_BUR <- cc_locus_table_SY_hexp_allButBAR_SCHIF_pivot %>%
  filter(Pop=="BUR") %>%
  ggplot(aes(x=as.double(SamplingYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  labs(y="Nei's genetic diversity Hexp", x="Sampling year") +
  ggtitle("Nei's genetic diversity per locus, Sampling year and population BUR (clone-corrected)") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines")) +
  scale_x_continuous(breaks=~axisTicks(., log = FALSE)) +
  ylim(0,1)
rarlocus_table_SY_hexp_all_BUR
```

![](Presentation_Summer_24_files/figure-gfm/cc%20Hexp%20table%20Pop%20SamplingYear%20removed%20BAR%20and%20SCHIF-3.png)<!-- -->

``` r
emf(file="cc_locus_table_SY_hexp_all_BUR.emf")
cc_locus_table_SY_hexp_all_BUR
dev.off()
```

    ## png 
    ##   2

``` r
#create zoom-in to FRE
cc_locus_table_SY_hexp_all_FRE <- cc_locus_table_SY_hexp_allButBAR_SCHIF_pivot %>%
  filter(Pop=="FRE") %>%
  ggplot(aes(x=as.double(SamplingYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  labs(y="Nei's genetic diversity Hexp", x="Sampling year") +
  ggtitle("Nei's genetic diversity per locus, Sampling year and population FRE (clone-corrected)") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines")) +
  scale_x_continuous(breaks=~axisTicks(., log = FALSE)) +
  ylim(0,1)
rarlocus_table_SY_hexp_all_FRE
```

![](Presentation_Summer_24_files/figure-gfm/cc%20Hexp%20table%20Pop%20SamplingYear%20removed%20BAR%20and%20SCHIF-4.png)<!-- -->

``` r
emf(file="cc_locus_table_SY_hexp_all_FRE.emf")
cc_locus_table_SY_hexp_all_FRE
dev.off()
```

    ## png 
    ##   2

``` r
#create zoom-in to KON
cc_locus_table_SY_hexp_all_KON <- cc_locus_table_SY_hexp_allButBAR_SCHIF_pivot %>%
  filter(Pop=="KON") %>%
  ggplot(aes(x=as.double(SamplingYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  labs(y="Nei's genetic diversity Hexp", x="Sampling year") +
  ggtitle("Nei's genetic diversity per locus, Sampling year and population KON (clone-corrected)") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines")) +
  scale_x_continuous(breaks=~axisTicks(., log = FALSE)) +
  ylim(0,1)
rarlocus_table_SY_hexp_all_KON
```

![](Presentation_Summer_24_files/figure-gfm/cc%20Hexp%20table%20Pop%20SamplingYear%20removed%20BAR%20and%20SCHIF-5.png)<!-- -->

``` r
emf(file="cc_locus_table_SY_hexp_all_KON.emf")
cc_locus_table_SY_hexp_all_KON
dev.off()
```

    ## png 
    ##   2

``` r
#create zoom-in to UEB
cc_locus_table_SY_hexp_all_UEB <- cc_locus_table_SY_hexp_allButBAR_SCHIF_pivot %>%
  filter(Pop=="UEB") %>%
  ggplot(aes(x=as.double(SamplingYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  labs(y="Nei's genetic diversity Hexp", x="Sampling year") +
  ggtitle("Nei's genetic diversity per locus, Sampling year and population UEB (clone-corrected)") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines")) +
  scale_x_continuous(breaks=~axisTicks(., log = FALSE)) +
  ylim(0,1)
rarlocus_table_SY_hexp_all_UEB
```

![](Presentation_Summer_24_files/figure-gfm/cc%20Hexp%20table%20Pop%20SamplingYear%20removed%20BAR%20and%20SCHIF-6.png)<!-- -->

``` r
emf(file="cc_locus_table_SY_hexp_all_UEB.emf")
cc_locus_table_SY_hexp_all_UEB
dev.off()
```

    ## png 
    ##   2

#### clone-corrected myData_genind_allMarkersOnly locus table per Population (Strata: Population / TruffleYear), removed BAR, SCHIF and GEN

``` r
setPop(cc_myData_genind_allMarkersOnly_SY) <- ~Pop
cc_myData_genind_allMarkersOnly_allButBAR_SCHIF_GEN_TY <- popsub(cc_myData_genind_allMarkersOnly_TY, exclude=c("BAR","SCHIF","GEN"))
setPop(cc_myData_genind_allMarkersOnly_allButBAR_SCHIF_GEN_TY) <- ~Pop/TruffleYear

# create a data frame and transpose the result;
# use sapply to iterate over all populations calculated by seppop.
# then use the locus_table function on the level Genotype

cc_locus_table_allButBAR_SCHIF_GEN_TY <- data.frame(t(
  sapply(seppop(cc_myData_genind_allMarkersOnly_allButBAR_SCHIF_GEN_TY),
    function(ls) poppr::locus_table(ls,lev="genotype"))))

# choose only the columns corresponding to Hexp values and rename them with their loci

cc_locus_table_TY_hexp_allButBAR_SCHIF_GEN <- 
  select(cc_locus_table_allButBAR_SCHIF_GEN_TY, X31:X45) %>%
  rename(aest06_1=X31,aest07_1=X32, aest15_1=X33, aest26_1=X34,aest28_1=X35, aest35_1=X36, aest36_1=X37, aest01_1=X38, aest10_1=X39, aest18_1=X40, aest24_1=X41, aest25_1=X42, aest29_1=X43, aest31_1=X44, mean=X45) %>%
  tibble::rownames_to_column("Pop") %>%
  separate(.,Pop,sep="_", into=c("Pop","TruffleYear"))
                                      
cc_locus_table_TY_hexp_allButBAR_SCHIF_GEN_pivot <- pivot_longer(cc_locus_table_TY_hexp_allButBAR_SCHIF_GEN,cols=3:16,names_to="locus",values_to="value")

ggplot(cc_locus_table_TY_hexp_allButBAR_SCHIF_GEN_pivot,aes(x=as.numeric(TruffleYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  facet_wrap(vars(Pop)) +
  labs(y="Nei's genetic diversity Hexp", x="Truffle year") +
  ggtitle("Nei's genetic diversity per locus, truffle year and population (clone-corrected)") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines"))
```

![](Presentation_Summer_24_files/figure-gfm/cc%20Hexp%20table%20Pop%20TruffleYear%20removed%20BAR%20SCHIF%20GEN-1.png)<!-- -->

``` r
emf(file= "cc_locus_table_TY_hexp_allButBAR_SCHIF_GEN_pivot.emf") # Opens a device
ggplot(cc_locus_table_TY_hexp_allButBAR_SCHIF_GEN_pivot,aes(x=as.numeric(TruffleYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  facet_wrap(vars(Pop)) +
  labs(y="Nei's genetic diversity Hexp", x="Truffle year") +
  ggtitle("Nei's genetic diversity per locus, truffle year and population (clone-corrected)") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines"))
dev.off()
```

    ## png 
    ##   2

``` r
#create zoom-in to BOB
cc_locus_table_TY_hexp_all_BOB <- cc_locus_table_TY_hexp_allButBAR_SCHIF_GEN_pivot %>%
  filter(Pop=="BOB") %>%
  ggplot(aes(x=as.double(TruffleYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  labs(y="Nei's genetic diversity Hexp", x="Truffle year") +
  ggtitle("Nei's genetic diversity per locus, truffle year and population BOB (clone-corrected)") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines")) +
  scale_x_continuous(breaks=~axisTicks(., log = FALSE)) +
  ylim(0,1)
rarlocus_table_TY_hexp_all_BOB
```

![](Presentation_Summer_24_files/figure-gfm/cc%20Hexp%20table%20Pop%20TruffleYear%20removed%20BAR%20SCHIF%20GEN-2.png)<!-- -->

``` r
emf(file="cc_locus_table_TY_hexp_all_BOB.emf")
cc_locus_table_TY_hexp_all_BOB
dev.off()
```

    ## png 
    ##   2

``` r
#create zoom-in to BUR
cc_locus_table_TY_hexp_all_BUR <- cc_locus_table_TY_hexp_allButBAR_SCHIF_GEN_pivot %>%
  filter(Pop=="BUR") %>%
  ggplot(aes(x=as.double(TruffleYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  labs(y="Nei's genetic diversity Hexp", x="Truffle year") +
  ggtitle("Nei's genetic diversity per locus, truffle year and population BUR (clone-corrected)") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines")) +
  scale_x_continuous(breaks=~axisTicks(., log = FALSE)) +
  ylim(0,1)
rarlocus_table_TY_hexp_all_BUR
```

![](Presentation_Summer_24_files/figure-gfm/cc%20Hexp%20table%20Pop%20TruffleYear%20removed%20BAR%20SCHIF%20GEN-3.png)<!-- -->

``` r
emf(file="cc_locus_table_TY_hexp_all_BUR.emf")
cc_locus_table_TY_hexp_all_BUR
dev.off()
```

    ## png 
    ##   2

``` r
#create zoom-in to FRE
cc_locus_table_TY_hexp_all_FRE <- cc_locus_table_TY_hexp_allButBAR_SCHIF_GEN_pivot %>%
  filter(Pop=="FRE") %>%
  ggplot(aes(x=as.double(TruffleYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  labs(y="Nei's genetic diversity Hexp", x="Truffle year") +
  ggtitle("Nei's genetic diversity per locus, truffle year and population FRE (clone-corrected)") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines")) +
  scale_x_continuous(breaks=~axisTicks(., log = FALSE)) +
  ylim(0,1)
rarlocus_table_TY_hexp_all_FRE
```

![](Presentation_Summer_24_files/figure-gfm/cc%20Hexp%20table%20Pop%20TruffleYear%20removed%20BAR%20SCHIF%20GEN-4.png)<!-- -->

``` r
emf(file="cc_locus_table_TY_hexp_all_FRE.emf")
cc_locus_table_TY_hexp_all_FRE
dev.off()
```

    ## png 
    ##   2

``` r
#create zoom-in to KON
cc_locus_table_TY_hexp_all_KON <- cc_locus_table_TY_hexp_allButBAR_SCHIF_GEN_pivot %>%
  filter(Pop=="KON") %>%
  ggplot(aes(x=as.double(TruffleYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  labs(y="Nei's genetic diversity Hexp", x="Truffle year") +
  ggtitle("Nei's genetic diversity per locus, truffle year and population KON (clone-corrected)") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines")) +
  scale_x_continuous(breaks=~axisTicks(., log = FALSE)) +
  ylim(0,1)
rarlocus_table_TY_hexp_all_KON
```

![](Presentation_Summer_24_files/figure-gfm/cc%20Hexp%20table%20Pop%20TruffleYear%20removed%20BAR%20SCHIF%20GEN-5.png)<!-- -->

``` r
emf(file="cc_locus_table_TY_hexp_all_KON.emf")
cc_locus_table_TY_hexp_all_KON
dev.off()
```

    ## png 
    ##   2

``` r
#create zoom-in to UEB
cc_locus_table_TY_hexp_all_UEB <- cc_locus_table_TY_hexp_allButBAR_SCHIF_GEN_pivot %>%
  filter(Pop=="UEB") %>%
  ggplot(aes(x=as.double(TruffleYear), y=value)) +
  geom_point(pch=20, alpha=0.3, size=0.8) +
  geom_smooth(method="gam", formula=y~s(x,k=2), colour="darkgreen", linewidth=0.6) +
  labs(y="Nei's genetic diversity Hexp", x="Truffle year") +
  ggtitle("Nei's genetic diversity per locus, truffle year and population UEB (clone-corrected)") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines")) +
  scale_x_continuous(breaks=~axisTicks(., log = FALSE)) +
  ylim(0,1)
rarlocus_table_TY_hexp_all_UEB
```

![](Presentation_Summer_24_files/figure-gfm/cc%20Hexp%20table%20Pop%20TruffleYear%20removed%20BAR%20SCHIF%20GEN-6.png)<!-- -->

``` r
emf(file="cc_locus_table_TY_hexp_all_UEB.emf")
cc_locus_table_TY_hexp_all_UEB
dev.off()
```

    ## png 
    ##   2

## 2. Simpson’s diversity

``` r
setPop(myData_genind_allMarkersOnly) <- ~Pop/SamplingYear
popdata_all_pop_year <- poppr(myData_genind_allMarkersOnly) %>%
  separate(Pop,c("Pop","SamplingYear"))
```

    ## Warning: Expected 2 pieces. Missing pieces filled with `NA` in 1 rows [122].

``` r
popdata_all_pop_year_graph <- filter(popdata_all_pop_year,Pop!="Total") %>%
ggplot(aes(as.numeric(SamplingYear),lambda)) +
  geom_point(pch=20, size=0.9) +
  facet_wrap(vars(Pop)) +
  labs(y="Simpson's index", x="Sampling year") +
  ggtitle("Simpson's index over the sampling years and all populations") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines"))
popdata_all_pop_year_graph
```

![](Presentation_Summer_24_files/figure-gfm/plot%20diversity%20Simpson-1.png)<!-- -->

``` r
emf(file="simpson_all_SY.emf")
popdata_all_pop_year_graph
dev.off()
```

    ## png 
    ##   2

``` r
#rarefied:
N_all      <- popdata_all_pop_year$N      # number of samples
lambda_all <- popdata_all_pop_year$lambda # Simpson's index
popdata_all_pop_year$rarLambda_all <- (N_all/(N_all - 1)) * lambda_all             # Corrected Simpson's index

rar_popdata_all_pop_year_graph <- filter(popdata_all_pop_year,Pop!="Total") %>%
ggplot(aes(as.numeric(SamplingYear),rarLambda_all)) +
  geom_point(pch=20, size=0.9) +
  facet_wrap(vars(Pop)) +
  labs(y="Rarefied Simpson's index", x="Sampling year") +
  ggtitle("Rarefied Simpson's index over the sampling years and all populations") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines"))
rar_popdata_all_pop_year_graph
```

    ## Warning: Removed 11 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](Presentation_Summer_24_files/figure-gfm/plot%20diversity%20Simpson-2.png)<!-- -->

``` r
#remove BAR and SCHIF
rar_popdata_all_pop_year_graph1 <- filter(popdata_all_pop_year,! Pop %in% c("Total", "BAR","SCHIF")) %>%
ggplot(aes(as.numeric(SamplingYear),rarLambda_all)) +
  geom_point(pch=20, size=0.9) +
  facet_wrap(vars(Pop)) +
  labs(y="Rarefied Simpson's index", x="Sampling year", subtitle="BAR and SCHIF: sample size too small") +
  ggtitle("Rarefied Simpson's index over the sampling years") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines"))
rar_popdata_all_pop_year_graph1
```

    ## Warning: Removed 9 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](Presentation_Summer_24_files/figure-gfm/plot%20diversity%20Simpson-3.png)<!-- -->

``` r
emf(file="rarSimpson_all_SY_minusBARSCHIF.emf")
rar_popdata_all_pop_year_graph1
```

    ## Warning: Removed 9 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

``` r
dev.off()
```

    ## png 
    ##   2

``` r
# TruffleYears

setPop(myData_genind_allMarkersOnly) <- ~Pop/TruffleYear
popdata_all_pop_year <- poppr(myData_genind_allMarkersOnly) %>%
  separate(Pop,c("Pop","TruffleYear"))
```

    ## Warning: Expected 2 pieces. Missing pieces filled with `NA` in 1 rows [123].

``` r
popdata_all_pop_year_graph_TY <- filter(popdata_all_pop_year,Pop!="Total") %>%
ggplot(aes(as.numeric(TruffleYear),lambda)) +
  geom_point(pch=20, size=0.9) +
  facet_wrap(vars(Pop)) +
  labs(y="Simpson's index", x="Truffle year") +
  ggtitle("Simpson's index over the truffle years and all populations") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines"))
popdata_all_pop_year_graph_TY
```

![](Presentation_Summer_24_files/figure-gfm/plot%20diversity%20Simpson-4.png)<!-- -->

``` r
emf(file="simpson_all_TY.emf")
popdata_all_pop_year_graph_TY
dev.off()
```

    ## png 
    ##   2

``` r
#rarefied:
N_all      <- popdata_all_pop_year$N      # number of samples
lambda_all <- popdata_all_pop_year$lambda # Simpson's index
popdata_all_pop_year$rarLambda_all <- (N_all/(N_all - 1)) * lambda_all             # Corrected Simpson's index

rar_popdata_all_pop_year_graph_TY <- filter(popdata_all_pop_year,Pop!="Total") %>%
ggplot(aes(as.numeric(TruffleYear),rarLambda_all)) +
  geom_point(pch=20, size=0.9) +
  facet_wrap(vars(Pop)) +
  labs(y="Rarefied Simpson's index", x="Truffle year") +
  ggtitle("Rarefied Simpson's index over the truffle years") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines"))
rar_popdata_all_pop_year_graph_TY
```

    ## Warning: Removed 9 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](Presentation_Summer_24_files/figure-gfm/plot%20diversity%20Simpson-5.png)<!-- -->

``` r
emf(file="rarSimpson_all_TY.emf")
rar_popdata_all_pop_year_graph_TY
```

    ## Warning: Removed 9 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

``` r
dev.off()
```

    ## png 
    ##   2

``` r
#remove BAR and SCHIF
rar_popdata_all_pop_year_TY_graph1 <- filter(popdata_all_pop_year,! Pop %in% c("Total", "BAR","SCHIF")) %>%
ggplot(aes(as.numeric(TruffleYear),rarLambda_all)) +
  geom_point(pch=20, size=0.9) +
  facet_wrap(vars(Pop)) +
  labs(y="Rarefied Simpson's index", x="Truffle year", subtitle="BAR and SCHIF: sample size too small") +
  ggtitle("Rarefied Simpson's index over the truffle years") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines"))
rar_popdata_all_pop_year_TY_graph1
```

    ## Warning: Removed 7 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](Presentation_Summer_24_files/figure-gfm/plot%20diversity%20Simpson-6.png)<!-- -->

``` r
emf(file="rarSimpson_all_TY_minusBARSCHIF.emf")
rar_popdata_all_pop_year_TY_graph1
```

    ## Warning: Removed 7 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

``` r
dev.off()
```

    ## png 
    ##   2

``` r
#clone corrected, Sampling Year:

setPop(cc_myData_genind_allMarkersOnly_SY) <- ~Pop/SamplingYear
cc_popdata_all_pop_year <- poppr(cc_myData_genind_allMarkersOnly_SY) %>%
  separate(Pop,c("Pop","SamplingYear"))
```

    ## Warning: Expected 2 pieces. Missing pieces filled with `NA` in 1 rows [122].

``` r
cc_popdata_all_pop_year_graph <- filter(cc_popdata_all_pop_year,Pop!="Total") %>%
ggplot(aes(as.numeric(SamplingYear),lambda)) +
  geom_point(pch=20, size=0.9) +
  facet_wrap(vars(Pop)) +
  labs(y="Simpson's index", x="Sampling year") +
  ggtitle("Simpson's index over the sampling years (clone-corrected)") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines"))
cc_popdata_all_pop_year_graph
```

![](Presentation_Summer_24_files/figure-gfm/plot%20diversity%20Simpson-7.png)<!-- -->

``` r
emf(file="cc_simpson_all_SY.emf")
cc_popdata_all_pop_year_graph
dev.off()
```

    ## png 
    ##   2

``` r
#rarefied doesn't make sense

# TruffleYears

setPop(cc_myData_genind_allMarkersOnly_TY) <- ~Pop/TruffleYear
cc_popdata_all_pop_year <- poppr(cc_myData_genind_allMarkersOnly_TY) %>%
  separate(Pop,c("Pop","TruffleYear"))
```

    ## Warning: Expected 2 pieces. Missing pieces filled with `NA` in 1 rows [123].

``` r
cc_popdata_all_pop_year_TY_graph <- filter(cc_popdata_all_pop_year,Pop!="Total") %>%
ggplot(aes(as.numeric(TruffleYear),lambda)) +
  geom_point(pch=20, size=0.9) +
  facet_wrap(vars(Pop)) +
  labs(y="Simpson's index", x="Truffle year") +
  ggtitle("Simpson's index over the truffle years (clone-corrected)") +
  theme(aspect.ratio=0.4,
    strip.background = element_blank(),
    strip.text=element_text(size=7,hjust=0.1, vjust=0.5),
    panel.grid = element_blank(),
    panel.spacing=unit(0.1,"lines"))
cc_popdata_all_pop_year_TY_graph
```

![](Presentation_Summer_24_files/figure-gfm/plot%20diversity%20Simpson-8.png)<!-- -->

``` r
emf(file="cc_simpson_all_TY.emf")
cc_popdata_all_pop_year_TY_graph
dev.off()
```

    ## png 
    ##   2

## 3. PCA

``` r
#uncorrected PCA with allMarkersOnly dataset
setPop(myData_genind_allMarkersOnly) <- ~Pop
x.pops_allMarkersOnly <- tab(myData_genind_allMarkersOnly,
                             freq=TRUE, NA.method="mean")
#x.pops_allMarkersOnly[grep("FRE", rownames(x.pops_allMarkersOnly)), ]

pca.pops.allMarkersOnly <- dudi.pca(df = x.pops_allMarkersOnly,
                                    center = TRUE, scale = FALSE, scannf = FALSE, nf = 2)

s.class(pca.pops.allMarkersOnly$li, fac=pop(myData_genind_allMarkersOnly),
        col=funky(15))
```

![](Presentation_Summer_24_files/figure-gfm/pca-1.png)<!-- -->

``` r
emf(file="pca.allMarkersOnly.emf")
s.class(pca.pops.allMarkersOnly$li, fac=pop(myData_genind_allMarkersOnly),
        col=funky(15))
dev.off()
```

    ## png 
    ##   2

``` r
#clonecorrected PCA  with samplingYear
setPop(cc_myData_genind_allMarkersOnly_SY) <- ~Pop
x.pops_cc.SY <- tab(cc_myData_genind_allMarkersOnly_SY,
                    freq=TRUE, NA.method="mean")
pca.pops_cc.SY <- dudi.pca(df = x.pops_cc.SY,
                           center = TRUE, scale = FALSE, scannf = FALSE, nf = 2)
s.class(pca.pops_cc.SY$li,
        fac=pop(cc_myData_genind_allMarkersOnly_SY),
        col=funky(15))
```

![](Presentation_Summer_24_files/figure-gfm/pca-2.png)<!-- -->

``` r
emf(file="pca.cc.allMarkersOnly_SY.emf")
s.class(pca.pops_cc.SY$li,
        fac=pop(cc_myData_genind_allMarkersOnly_SY),
        col=funky(15))
dev.off()
```

    ## png 
    ##   2

``` r
#clonecorrected PCA  with TruffleYear
setPop(cc_myData_genind_allMarkersOnly_TY) <- ~Pop
x.pops_cc.TY <- tab(cc_myData_genind_allMarkersOnly_TY,
                    freq=TRUE, NA.method="mean")
pca.pops_cc.TY_2 <- dudi.pca(df = x.pops_cc.TY,
                           center = TRUE, scale = FALSE, scannf = FALSE, nf = 2)
pca.pops_cc.TY_3 <- dudi.pca(df = x.pops_cc.TY,
                           center = TRUE, scale = FALSE, scannf = FALSE, nf = 3)
s.class(pca.pops_cc.TY_2$li,
        fac=pop(cc_myData_genind_allMarkersOnly_TY),
        col=funky(15))
s.class(pca.pops_cc.TY_3$li,
        fac=pop(cc_myData_genind_allMarkersOnly_TY),
        col=funky(15))
```

![](Presentation_Summer_24_files/figure-gfm/pca-3.png)<!-- -->

``` r
emf(file="pca.cc.allMarkersOnly_TY.emf")
s.class(pca.pops_cc.TY_2$li,
        fac=pop(cc_myData_genind_allMarkersOnly_TY),
        col=funky(15))
dev.off()
```

    ## png 
    ##   2

``` r
#plot(pca.pops_cc.TY$li, col=
```
