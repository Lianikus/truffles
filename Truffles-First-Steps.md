Truffle Monitoring
================
Lia Baumann
2024-03-05

## Basis

The basis for the research about Truffle Monitoring data is described in
detail in these publications: Virginie Molinier et al. (2013) Virgine
Molinier et al. (2016) Virginie Molinier et al. (2016) Staubli et al.
(2022) Steidinger et al. (2022) Legendre and Fortin (2010) Kamvar,
Tabima, and Grünwald (2014)

## Monitoring Sites

![](Truffles-First-Steps_files/figure-gfm/Maps%20laden-1.png)<!-- -->

The dataset was already corrected in the following way: - Removal of
samples with less than 10 markers - Removal of samples with two mating
types.

# Number of samples per Site

| Site_1_abrev | Site                  |   n |
|:-------------|:----------------------|----:|
| ALD          | Aldingen              |  71 |
| BAR          | Barbengo              |   2 |
| BOB          | Bohlingen Buche       | 473 |
| BOH          | Bohlingen Hasel       |  54 |
| BRU          | Zürich Brücke         | 229 |
| BUR          | Bursins               | 295 |
| FRB          | Freiburg-Brumaleweg   |  20 |
| FRE          | Freiburg-Wittnau      | 268 |
| FRI          | Frick                 | 125 |
| GEN          | Genolier              |  11 |
| HAN          | Hannover              |  73 |
| KON          | Köniz                 | 272 |
| LIM          | Zürich Limmat         | 231 |
| NEU          | Neuchatel             |  38 |
| RIE          | Rietheim              |  67 |
| SCD          | Schlieren Dreieck     | 123 |
| SCG          | Schlieren Graben      |  24 |
| SCL          | Schlieren Grillplatz  | 298 |
| SCS          | Schiffenensee         |  50 |
| TRO          | Trossingen            |  37 |
| UEB          | Bohlingen Ueberlingen | 204 |
| UST          | Uster                 |  43 |
| WSL          | WSL                   | 283 |

Number of observations per site

![](Truffles-First-Steps_files/figure-gfm/Data%20structure%20per%20site-1.png)<!-- -->

\#Idee: alle Observationen auf einer Zeitachse und nach Standort
unterteilt anzeigen

![](Truffles-First-Steps_files/figure-gfm/all%20observations-1.png)<!-- -->

Show missing data

``` r
setPop(myData) <- ~Pop
(microsats_lt <- locus_table(myData))
```

    ## 
    ## allele = Number of observed alleles
    ## 1-D = Simpson index
    ## Hexp = Nei's 1978 gene diversity
    ## ------------------------------------------

    ##            summary
    ## locus       allele   1-D  Hexp Evenness
    ##   aest06_1    7.00  0.65  0.65     0.63
    ##   aest07_1    7.00  0.64  0.64     0.73
    ##   aest15_1    6.00  0.37  0.37     0.60
    ##   aest26_1   12.00  0.76  0.76     0.85
    ##   aest28_1   23.00  0.82  0.82     0.78
    ##   aest35_1    7.00  0.53  0.53     0.71
    ##   aest36_1    8.00  0.57  0.57     0.80
    ##   aest01_1    8.00  0.75  0.75     0.81
    ##   aeste10_1   8.00  0.78  0.78     0.89
    ##   aest18_1    5.00  0.72  0.72     0.87
    ##   aest24_1    9.00  0.55  0.55     0.52
    ##   aest25_1    5.00  0.60  0.60     0.79
    ##   aest29_1    9.00  0.69  0.69     0.78
    ##   aest31_1    8.00  0.68  0.68     0.83
    ##   mean        8.71  0.65  0.65     0.76

``` r
info_table(myData, type = "missing", plot = TRUE)
```

![](Truffles-First-Steps_files/figure-gfm/allele%20frequencies-1.png)<!-- -->

    ##           Locus
    ## Population aest06_1 aest07_1 aest15_1 aest26_1 aest28_1 aest35_1 aest36_1
    ##      FRE          .        .        .        .  0.01866        .        .
    ##      ALD          .        .        .        .        .        .        .
    ##      RIE          .        .        .        .        .        .        .
    ##      TRO          .        .        .        .        .        .        .
    ##      SCG          .        .        .        .  0.20833        .        .
    ##      BOH          .        .        .        .  0.07407        .        .
    ##      BOB    0.00211        .        .        .  0.01903        .  0.00211
    ##      FRB          .        .        .        .  0.20000        .        .
    ##      UEB          .        .  0.12255        .  0.05882        .        .
    ##      SCL          .        .  0.00336        .  0.63423        .  0.00336
    ##      SCD          .        .        .        .  0.05691        .        .
    ##      WSL          .        .        .        .        .        .        .
    ##      BUR          .        .        .  0.02712  0.00678  0.00678  0.00678
    ##      NEU          .        .        .        .  0.05263        .        .
    ##      SCS          .        .        .        .  0.02041        .        .
    ##      UST          .        .        .        .  0.04651        .        .
    ##      KON          .  0.00368  0.00368        .  0.01103        .  0.00735
    ##      FRI          .        .        .        .  0.00800        .        .
    ##      BAR          .        .        .        .        .        .        .
    ##      LIM          .        .        .        .        .        .        .
    ##      BRU          .        .        .        .        .        .        .
    ##      HAN          .        .        .        .        .        .        .
    ##      GEN          .        .        .        .        .        .  0.09091
    ##      SCHIF        .        .        .        .        .        .        .
    ##      Total  0.00030  0.00030  0.00820  0.00243  0.07475  0.00061  0.00213
    ##           Locus
    ## Population aest01_1 aeste10_1 aest18_1 aest24_1 aest25_1 aest29_1 aest31_1
    ##      FRE    0.13433   0.11194        .  0.01493        .  0.00373  0.11194
    ##      ALD    0.07042         .        .  0.04225        .        .  0.07042
    ##      RIE    0.04478   0.01493        .        .        .        .  0.02985
    ##      TRO    0.08108         .        .  0.02703        .        .  0.02703
    ##      SCG    0.25000   0.04167        .  0.04167        .        .  0.12500
    ##      BOH    0.22222   0.03704  0.01852  0.05556        .  0.01852  0.07407
    ##      BOB    0.06765   0.00634  0.00211  0.00634  0.00423        .  0.03171
    ##      FRB    0.10000   0.20000        .        .        .  0.05000  0.10000
    ##      UEB    0.20098   0.02451        .  0.02451  0.00490  0.01471  0.15686
    ##      SCL    0.09396   0.01007  0.00336  0.00336        .  0.00671  0.06040
    ##      SCD    0.13008   0.02439        .  0.02439        .  0.01626  0.07317
    ##      WSL          .         .        .        .        .        .        .
    ##      BUR    0.19322   0.00339        .  0.02034        .  0.00678  0.11864
    ##      NEU    0.10526   0.02632        .        .  0.02632        .  0.02632
    ##      SCS    0.10204         .        .  0.02041        .  0.02041  0.02041
    ##      UST    0.11628   0.02326        .  0.02326        .  0.02326  0.04651
    ##      KON    0.20221         .        .  0.01103        .  0.00368  0.13235
    ##      FRI    0.03200   0.00800        .  0.00800        .        .  0.00800
    ##      BAR    0.50000         .        .        .        .        .        .
    ##      LIM    0.06494         .        .  0.03030        .  0.00866  0.05628
    ##      BRU    0.19214         .        .  0.07860        .  0.02620  0.18341
    ##      HAN    0.10959         .        .  0.02740        .  0.01370  0.01370
    ##      GEN    0.09091         .        .        .        .        .  0.09091
    ##      SCHIF        .         .        .        .        .        .        .
    ##      Total  0.11638   0.01702  0.00091  0.01914  0.00122  0.00729  0.07718
    ##           Locus
    ## Population    Mean
    ##      FRE   0.02825
    ##      ALD   0.01308
    ##      RIE   0.00640
    ##      TRO   0.00965
    ##      SCG   0.04762
    ##      BOH   0.03571
    ##      BOB   0.01012
    ##      FRB   0.04643
    ##      UEB   0.04342
    ##      SCL   0.05849
    ##      SCD   0.02323
    ##      WSL         .
    ##      BUR   0.02785
    ##      NEU   0.01692
    ##      SCS   0.01312
    ##      UST   0.01993
    ##      KON   0.02679
    ##      FRI   0.00457
    ##      BAR   0.03571
    ##      LIM   0.01144
    ##      BRU   0.03431
    ##      HAN   0.01174
    ##      GEN   0.01948
    ##      SCHIF       .
    ##      Total 0.02342

``` r
#for more stuff (removal of individuals or loci etc.): grunwaldlab.github.io/population_genetics_in_r/locus_stats.html
```

    ## `summarise()` has grouped output by 'Pop', 'Month', 'SamplingYear'. You can
    ## override using the `.groups` argument.

![](Truffles-First-Steps_files/figure-gfm/treemap-1.png)<!-- -->

![](Truffles-First-Steps_files/figure-gfm/examine%20genetic%20data%20as%20genind-1.png)<!-- -->

    ## 
    ## // Number of individuals: 3291
    ## // Group sizes: 268 71 67 37 24 54 473 20 204 298 123 283 295 38 49 43 272 125 2 231 229 73 11 1
    ## // Number of alleles per locus: 7 7 6 12 23 7 8 8 8 5 9 5 9 8
    ## // Number of alleles per group: 42 35 44 23 31 37 55 39 42 52 37 25 49 26 34 38 42 43 19 61 57 28 20 14
    ## // Percentage of missing data: 2.34 %
    ## // Observed heterozygosity: 0

![](Truffles-First-Steps_files/figure-gfm/examine%20genetic%20data%20as%20genind-2.png)<!-- -->

    ## NULL

\##Clone correction from
<https://grunwaldlab.github.io/Population_Genetics_in_R/Population_Strata.html>
When dealing with clonal populations, analyses are typically conducted
with and without clone correction. Clone correction is a method of
censoring a data set such that only one individual per MLG is
represented per population (Milgroom, 1996; Grünwald et al., 2003;
Grünwald & Hoheisel, 2006). This technique is commonly used with the
index of association and genotypic diversity measures since clone
corrected populations approximate behavior of sexual populations. Since
we want to only observe unique genotypes per population, clone
correction requires specification of the stratifications at which clones
should be censored. This section will show how to clone correct at a
specific stratification and also compare the results with uncorrected
data.

Question: Will allelic diversity increase or decrease with
clone-censored data?

The graph shows a decrease of diversity for most markers when
clone-correcting the data (increase of Simpson index means a decrease of
genotypic diversity).

\#PCA

``` r
x.pops <- tab(myData_genind, freq=TRUE, NA.method="mean")
pca.pops <- dudi.pca(df = x.pops, center = TRUE, scale = FALSE, scannf = FALSE, nf = 2)
s.class(pca.pops$li, fac=pop(myData_genind), col=funky(15))
```

![](Truffles-First-Steps_files/figure-gfm/first%20pca-1.png)<!-- -->

``` r
#clonecorrected PCA
x.pops_cc <- tab(cc_myData_genind, freq=TRUE, NA.method="mean")
pca.pops_cc <- dudi.pca(df = x.pops_cc, center = TRUE, scale = FALSE, scannf = FALSE, nf = 2)
s.class(pca.pops_cc$li, fac=pop(cc_myData_genind), col=funky(15))
```

![](Truffles-First-Steps_files/figure-gfm/first%20pca-2.png)<!-- -->

## Genotype accumulation curve

A genotype accumulation curve is a tool that allows you to assess how
much power you have to discriminate between unique individuals given a
random sample of n loci. This analysis is particularly important for
clonal organisms to confirm that a plateau has been reached in the
number of loci necessary to discriminate individuals. We specified
sample = 1000 in our function call. This means that for each boxplot, n
loci were randomly sampled 1000 times in order to create the
distribution. Since this data has been curated, we can see that we have
reached the plateau with 13 loci.

``` r
gac <- genotype_curve(myData, sample = 1000, quiet = TRUE)
```

![](Truffles-First-Steps_files/figure-gfm/genotype%20accumulation%20curve-1.png)<!-- -->

\##The index of association The index of association (IA) was originally
proposed by Brown et al. (Brown, Feldman & Nevo, 1980) and implemented
in the poppr R package (Kamvar, Tabima & Grünwald, 2014) using a
permutation approach to assess if loci are linked as described
previously by Agapow and Burt \[@\]. Agapow and Burt also described the
index r¯d that accounts for the number of loci sampled that is less
biased and will be used here. The data we will use in this chapter are
populations of Phytophthora infestans from North and South America (Goss
et al., 2014). We will use the index of association to test the
hypothesis that Mexico is the putative origin of P. infestans where
populations are expected to be sexual while populations in South America
are expected to be clonal. Next, we will analyze the North American
population with the index of association and use 999 permutations of the
data in order to give us a p-value. Note that the p-value is calculated
with the original observation included.

``` r
setPop(myData) <-~Pop
LD_BUR <- popsub(myData, "BUR")
ia(LD_BUR, sample = 999)
```

![](Truffles-First-Steps_files/figure-gfm/index%20of%20association%20and%20LD-1.png)<!-- -->

    ##        Ia      p.Ia     rbarD      p.rD 
    ## 5.2302076 0.0010000 0.4408573 0.0010000

``` r
setPop(cc_myData) <- ~Pop
LD_BUR_cc <- popsub(cc_myData,"BUR")
ia(LD_BUR_cc,sample=999)
```

![](Truffles-First-Steps_files/figure-gfm/index%20of%20association%20and%20LD-2.png)<!-- -->

    ##        Ia      p.Ia     rbarD      p.rD 
    ## 4.3625975 0.0010000 0.3413024 0.0010000

``` r
setPop(myData) <- ~Pop/SamplingYear
popdata <- poppr(myData)
popdata
```

    ##            Pop    N MLG  eMLG       SE      H     G lambda   E.5    Hexp
    ## 1     FRE_2011   27   5  3.37 8.04e-01 1.0415  2.08 0.5185 0.587 0.20512
    ## 2     ALD_2011   25   7  3.90 1.03e+00 1.1814  2.08 0.5184 0.477 0.15381
    ## 3     RIE_2011   20  17  9.16 7.56e-01 2.7616 14.29 0.9300 0.896 0.40576
    ## 4     TRO_2011   17   6  4.60 8.30e-01 1.4469  3.32 0.6990 0.714 0.15847
    ## 5     SCG_2011    9   5  5.00 0.00e+00 1.5230  4.26 0.7654 0.910 0.29960
    ## 6     BOH_2011   12  10  8.50 5.84e-01 2.2103  8.00 0.8750 0.862 0.27641
    ## 7     BOB_2011    6   5  5.00 0.00e+00 1.5607  4.50 0.7778 0.930 0.30000
    ## 8     FRB_2011    4   2  2.00 0.00e+00 0.5623  1.60 0.3750 0.795 0.17857
    ## 9     UEB_2011   27   7  4.81 9.06e-01 1.6201  4.03 0.7517 0.747 0.11905
    ## 10    ALD_2012   13   5  4.08 7.30e-01 1.0438  1.99 0.4970 0.537 0.12179
    ## 11    TRO_2012   10   3  3.00 0.00e+00 0.8979  2.17 0.5400 0.807 0.16667
    ## 12    FRE_2012   21   4  2.43 8.21e-01 0.5671  1.35 0.2585 0.457 0.20272
    ## 13    RIE_2012   11  10  9.18 3.86e-01 2.2719  9.31 0.8926 0.955 0.37792
    ## 14    SCG_2012    7   5  5.00 0.00e+00 1.4751  3.77 0.7347 0.821 0.02041
    ## 15    BOB_2012   41  22  7.02 1.30e+00 2.5485  6.81 0.8531 0.492 0.23468
    ## 16    BOH_2012   19  17  9.47 5.98e-01 2.7985 15.70 0.9363 0.953 0.32479
    ## 17    SCL_2012  106  30  5.50 1.32e+00 2.2257  4.39 0.7723 0.411 0.12831
    ## 18    SCD_2012   28  25  9.64 5.40e-01 3.1837 23.06 0.9566 0.953 0.36200
    ## 19    UEB_2012    7   7  7.00 0.00e+00 1.9459  7.00 0.8571 1.000 0.29048
    ## 20    WSL_2012    1   1  1.00 0.00e+00 0.0000  1.00 0.0000   NaN     NaN
    ## 21    FRB_2012    5   4  4.00 0.00e+00 1.3322  3.57 0.7200 0.922 0.43810
    ## 22    ALD_2013   18   6  4.20 8.98e-01 1.3031  2.70 0.6296 0.634 0.08357
    ## 23    RIE_2013   17  14  8.85 7.95e-01 2.5578 11.56 0.9135 0.887 0.39828
    ## 24    TRO_2013   10   4  4.00 0.00e+00 0.9404  1.92 0.4800 0.591 0.13651
    ## 25    BOB_2013   56  33  8.52 1.06e+00 3.2006 17.42 0.9426 0.697 0.32954
    ## 26    SCG_2013    2   2  2.00 0.00e+00 0.6931  2.00 0.5000 1.000 0.07692
    ## 27    SCL_2013   13  13 10.00 7.30e-08 2.5649 13.00 0.9231 1.000 0.29238
    ## 28    FRE_2013   12   5  4.48 5.84e-01 1.2343  2.57 0.6111 0.645 0.28254
    ## 29    UEB_2013   56  15  5.91 1.16e+00 2.1602  6.40 0.8438 0.704 0.14545
    ## 30    FRB_2013    5   5  5.00 0.00e+00 1.6094  5.00 0.8000 1.000 0.44762
    ## 31    BUR_2013   74  17  4.21 1.32e+00 1.5514  2.24 0.5533 0.333 0.10676
    ## 32    SCD_2013    9   9  9.00 0.00e+00 2.1972  9.00 0.8889 1.000 0.38135
    ## 33    NEU_2013    7   5  5.00 0.00e+00 1.4751  3.77 0.7347 0.821 0.04082
    ## 34    SCS_2013    8   5  5.00 0.00e+00 1.3863  3.20 0.6875 0.733 0.12500
    ## 35    BOH_2013    9   9  9.00 0.00e+00 2.1972  9.00 0.8889 1.000 0.21117
    ## 36    UST_2013   32  22  8.52 1.02e+00 2.8968 13.84 0.9277 0.750 0.37897
    ## 37    KON_2013    5   5  5.00 0.00e+00 1.6094  5.00 0.8000 1.000 0.10714
    ## 38    FRI_2013   54  10  3.76 1.10e+00 1.2763  2.15 0.5350 0.445 0.16651
    ## 39    RIE_2014   13  12  9.42 4.94e-01 2.4583 11.27 0.9112 0.961 0.49833
    ## 40    SCL_2014   60  19  5.05 1.31e+00 1.9374  3.69 0.7289 0.453 0.07239
    ## 41    FRB_2014    5   5  5.00 0.00e+00 1.6094  5.00 0.8000 1.000 0.26905
    ## 42    ALD_2014   10   8  8.00 0.00e+00 2.0253  7.14 0.8600 0.934 0.22500
    ## 43    BOB_2014   39  33  9.46 7.02e-01 3.4013 25.78 0.9612 0.854 0.37496
    ## 44    UEB_2014   47  25  7.22 1.30e+00 2.6871  8.27 0.8791 0.531 0.15447
    ## 45    UST_2014   11   9  8.36 4.81e-01 2.1458  8.07 0.8760 0.936 0.40519
    ## 46    SCS_2014   18   6  3.78 9.72e-01 1.0379  1.86 0.4630 0.473 0.12401
    ## 47    SCG_2014    6   6  6.00 0.00e+00 1.7918  6.00 0.8333 1.000 0.39048
    ## 48    SCD_2014   11   7  6.55 4.98e-01 1.7678  4.84 0.7934 0.790 0.29596
    ## 49    FRE_2014   15   8  5.90 9.09e-01 1.7075  3.81 0.7378 0.623 0.07493
    ## 50    KON_2014   14   9  6.93 8.49e-01 1.9652  5.44 0.8163 0.724 0.22135
    ## 51    FRI_2014    4   3  3.00 0.00e+00 1.0397  2.67 0.6250 0.912 0.34524
    ## 52    BAR_2014    2   2  2.00 0.00e+00 0.6931  2.00 0.5000 1.000 0.38462
    ## 53    BOH_2014   11   7  6.55 4.98e-01 1.7678  4.84 0.7934 0.790 0.25079
    ## 54    BUR_2014   10   8  8.00 0.00e+00 2.0253  7.14 0.8600 0.934 0.46236
    ## 55    NEU_2014    4   3  3.00 0.00e+00 1.0397  2.67 0.6250 0.912 0.14286
    ## 56    ALD_2015    5   4  4.00 0.00e+00 1.3322  3.57 0.7200 0.922 0.07857
    ## 57    BOB_2015   23  22  9.82 3.82e-01 3.0752 21.16 0.9527 0.976 0.41007
    ## 58    KON_2015    2   2  2.00 0.00e+00 0.6931  2.00 0.5000 1.000 0.14286
    ## 59    RIE_2015    6   6  6.00 0.00e+00 1.7918  6.00 0.8333 1.000 0.41190
    ## 60    FRI_2015    2   2  2.00 0.00e+00 0.6931  2.00 0.5000 1.000 0.00000
    ## 61    UEB_2015   39  20  6.92 1.29e+00 2.4896  7.28 0.8626 0.568 0.11977
    ## 62    FRE_2015    8   3  3.00 0.00e+00 0.7356  1.68 0.4062 0.630 0.05867
    ## 63    FRB_2015    1   1  1.00 0.00e+00 0.0000  1.00 0.0000   NaN     NaN
    ## 64    SCL_2015   11   5  4.64 4.81e-01 1.1596  2.28 0.5620 0.586 0.07229
    ## 65    SCD_2015   25  20  8.57 1.02e+00 2.7889 11.36 0.9120 0.679 0.27263
    ## 66    BOH_2015    3   3  3.00 0.00e+00 1.0986  3.00 0.6667 1.000 0.45455
    ## 67    BUR_2015    7   5  5.00 0.00e+00 1.4751  3.77 0.7347 0.821 0.21156
    ## 68    NEU_2015    1   1  1.00 0.00e+00 0.0000  1.00 0.0000   NaN     NaN
    ## 69    BOB_2016   60  29  8.04 1.15e+00 3.0057 13.85 0.9278 0.669 0.34123
    ## 70    SCS_2016   20   3  2.00 6.88e-01 0.3944  1.23 0.1850 0.469 0.05000
    ## 71    SCD_2016   47   8  3.21 1.04e+00 1.0144  1.76 0.4319 0.432 0.16128
    ## 72    KON_2016   54  20  6.40 1.27e+00 2.3764  6.69 0.8505 0.582 0.13086
    ## 73    LIM_2016   71  18  5.34 1.25e+00 2.0069  4.12 0.7570 0.484 0.43362
    ## 74    BRU_2016  118  12  3.85 1.18e+00 1.3599  2.16 0.5381 0.402 0.21440
    ## 75    BUR_2016   44   8  4.51 9.47e-01 1.6104  3.98 0.7490 0.745 0.11641
    ## 76    UEB_2016   17   8  6.15 8.70e-01 1.8746  5.25 0.8097 0.771 0.28075
    ## 77    FRE_2016   70   8  3.23 1.03e+00 1.0507  1.83 0.4531 0.445 0.00255
    ## 78    SCL_2016   58  10  4.31 1.08e+00 1.5019  2.80 0.6427 0.515 0.04849
    ## 79    WSL_2016   58   2  1.32 4.66e-01 0.1500  1.07 0.0666 0.441 0.00968
    ## 80    NEU_2016   17   4  3.28 6.61e-01 0.9161  1.89 0.4706 0.593 0.10568
    ## 81    FRI_2016   10   3  3.00 0.00e+00 0.6390  1.52 0.3400 0.576 0.16825
    ## 82    HAN_2016    1   1  1.00 0.00e+00 0.0000  1.00 0.0000   NaN     NaN
    ## 83    HAN_2017   45   9  4.99 9.82e-01 1.7687  4.70 0.7872 0.760 0.29375
    ## 84    LIM_2017   84   8  3.64 8.38e-01 1.3544  3.10 0.6774 0.731 0.38490
    ## 85    BRU_2017   19   8  4.95 1.04e+00 1.4832  2.76 0.6371 0.515 0.35426
    ## 86    SCL_2017   37   7  2.97 1.01e+00 0.8826  1.60 0.3755 0.424 0.08082
    ## 87    KON_2017   30  18  7.12 1.24e+00 2.4295  6.25 0.8400 0.507 0.18889
    ## 88    NEU_2017    2   1  1.00 0.00e+00 0.0000  1.00 0.0000   NaN 0.00000
    ## 89    WSL_2017   61   2  1.16 3.70e-01 0.0836  1.03 0.0322 0.382 0.00468
    ## 90    BOB_2017  146  42  7.33 1.26e+00 2.9727 10.94 0.9086 0.536 0.31121
    ## 91    UEB_2017    6   4  4.00 0.00e+00 1.3297  3.60 0.7222 0.935 0.20000
    ## 92    FRE_2017  115   9  2.56 9.92e-01 0.7834  1.45 0.3118 0.381 0.03058
    ## 93    FRI_2017   52   5  3.32 5.89e-01 1.2211  3.07 0.6746 0.867 0.28383
    ## 94    SCS_2017    3   2  2.00 0.00e+00 0.6365  1.80 0.4444 0.899 0.00000
    ## 95    SCD_2017    3   2  2.00 0.00e+00 0.6365  1.80 0.4444 0.899 0.33333
    ## 96    FRI_2018    3   2  2.00 0.00e+00 0.6365  1.80 0.4444 0.899 0.33333
    ## 97    HAN_2018   24   7  4.53 9.30e-01 1.5127  3.47 0.7118 0.698 0.26341
    ## 98    UEB_2018    2   2  2.00 0.00e+00 0.6931  2.00 0.5000 1.000 0.27273
    ## 99    WSL_2018   14   1  1.00 0.00e+00 0.0000  1.00 0.0000   NaN 0.00000
    ## 100   BOB_2018   74  40  8.81 9.65e-01 3.4244 23.01 0.9565 0.741 0.36130
    ## 101   KON_2018   45  13  5.03 1.23e+00 1.7672  3.22 0.6894 0.457 0.12869
    ## 102   LIM_2018   76  11  3.95 9.22e-01 1.5132  3.53 0.7164 0.713 0.42304
    ## 103   BRU_2018   92  19  4.93 1.25e+00 1.9227  4.02 0.7509 0.516 0.56134
    ## 104   BUR_2018   41   8  3.64 1.00e+00 1.3173  2.78 0.6401 0.651 0.27119
    ## 105   SCL_2018   13   2  1.77 4.21e-01 0.2712  1.17 0.1420 0.531 0.00000
    ## 106   GEN_2018    1   1  1.00 0.00e+00 0.0000  1.00 0.0000   NaN     NaN
    ## 107   HAN_2019    3   1  1.00 0.00e+00 0.0000  1.00 0.0000   NaN 0.00000
    ## 108   WSL_2019   72   2  1.26 4.39e-01 0.1269  1.06 0.0540 0.422 0.00782
    ## 109   KON_2019   31  12  6.51 1.09e+00 2.1906  7.12 0.8595 0.771 0.18021
    ## 110   BUR_2019   21   6  4.30 8.49e-01 1.4468  3.53 0.7166 0.778 0.04592
    ## 111   BOB_2019   28  16  7.38 1.17e+00 2.4653  8.52 0.8827 0.699 0.38813
    ## 112   UEB_2019    3   3  3.00 0.00e+00 1.0986  3.00 0.6667 1.000 0.14286
    ## 113   GEN_2019    1   1  1.00 0.00e+00 0.0000  1.00 0.0000   NaN     NaN
    ## 114   GEN_2020    3   2  2.00 0.00e+00 0.6365  1.80 0.4444 0.899 0.28571
    ## 115   BUR_2020   61   4  2.26 5.69e-01 0.6771  1.63 0.3854 0.648 0.17000
    ## 116   KON_2020   25   7  4.96 8.97e-01 1.6651  4.37 0.7712 0.786 0.24643
    ## 117   WSL_2020   73   3  1.27 4.83e-01 0.1446  1.06 0.0537 0.365 0.01957
    ## 118   NEU_2020    5   3  3.00 0.00e+00 0.9503  2.27 0.5600 0.802 0.14286
    ## 119   WSL_2021    4   2  2.00 0.00e+00 0.5623  1.60 0.3750 0.795 0.28571
    ## 120   BUR_2021   25   3  2.34 5.31e-01 0.6592  1.61 0.3808 0.659 0.31286
    ## 121   KON_2021   28   3  2.31 5.17e-01 0.6649  1.64 0.3903 0.678 0.16239
    ## 122   NEU_2021    2   2  2.00 0.00e+00 0.6931  2.00 0.5000 1.000 0.21429
    ## 123   GEN_2021    1   1  1.00 0.00e+00 0.0000  1.00 0.0000   NaN     NaN
    ## 124   GEN_2022    5   1  1.00 0.00e+00 0.0000  1.00 0.0000   NaN 0.00000
    ## 125   BUR_2022   12   4  3.50 5.84e-01 0.8370  1.71 0.4167 0.546 0.14740
    ## 126   KON_2022   13   6  5.04 7.26e-01 1.4105  2.96 0.6627 0.634 0.12088
    ## 127   KON_2023   25   7  3.65 1.05e+00 1.0824  1.88 0.4672 0.449 0.15795
    ## 128 SCHIF_2016    1   1  1.00 0.00e+00 0.0000  1.00 0.0000   NaN     NaN
    ## 129      Total 3291 725  9.18 8.70e-01 5.0266 47.33 0.9789 0.306 0.65262
    ##          Ia    rbarD   File
    ## 1    9.5017  0.95125 myData
    ## 2    5.0311  0.58947 myData
    ## 3    0.4744  0.04817 myData
    ## 4    2.1882  0.44171 myData
    ## 5    1.8907  0.27233 myData
    ## 6    0.9894  0.12531 myData
    ## 7    3.1449  0.45153 myData
    ## 8    4.0000  1.00000 myData
    ## 9    0.5010  0.16752 myData
    ## 10   3.8976  0.77955 myData
    ## 11   2.7692  0.69822 myData
    ## 12   9.9552  0.99553 myData
    ## 13   0.8005  0.07362 myData
    ## 14   0.0000      NaN myData
    ## 15   2.8460  0.27082 myData
    ## 16   2.2123  0.20462 myData
    ## 17   2.7289  0.32902 myData
    ## 18   0.9847  0.10993 myData
    ## 19   0.0989  0.01239 myData
    ## 20       NA       NA myData
    ## 21   6.3816  0.53228 myData
    ## 22   1.2687  0.26752 myData
    ## 23   1.3089  0.11040 myData
    ## 24   3.5236  0.59443 myData
    ## 25   2.3433  0.20330 myData
    ## 26       NA       NA myData
    ## 27   0.4092  0.04197 myData
    ## 28   9.8665  0.89797 myData
    ## 29   0.4801  0.07186 myData
    ## 30   5.3399  0.44526 myData
    ## 31   5.9404  0.52561 myData
    ## 32   0.5075  0.05691 myData
    ## 33  -0.1667 -0.16667 myData
    ## 34   2.4306  0.48839 myData
    ## 35   0.6302  0.09065 myData
    ## 36   1.3771  0.11007 myData
    ## 37  -0.0580 -0.02903 myData
    ## 38   4.3458  0.65445 myData
    ## 39   1.1958  0.11004 myData
    ## 40   1.3656  0.16330 myData
    ## 41   0.0714  0.01466 myData
    ## 42   2.1623  0.27368 myData
    ## 43   1.6149  0.13661 myData
    ## 44   1.3403  0.14752 myData
    ## 45   0.0348  0.00348 myData
    ## 46   4.0745  0.52053 myData
    ## 47   1.0553  0.10682 myData
    ## 48   2.4394  0.30525 myData
    ## 49   0.4057  0.10608 myData
    ## 50   1.8325  0.23444 myData
    ## 51   4.9873  0.62382 myData
    ## 52       NA       NA myData
    ## 53   0.6079  0.06809 myData
    ## 54   5.0556  0.39003 myData
    ## 55   1.0000  0.33333 myData
    ## 56  -0.3556 -0.35635 myData
    ## 57   1.5026  0.13672 myData
    ## 58       NA       NA myData
    ## 59   0.8642  0.09692 myData
    ## 60       NA       NA myData
    ## 61   0.6324  0.07722 myData
    ## 62   0.2336  0.11690 myData
    ## 63       NA       NA myData
    ## 64   1.8169  0.61465 myData
    ## 65   1.7211  0.19609 myData
    ## 66   0.5000  0.10000 myData
    ## 67   3.3425  0.66859 myData
    ## 68       NA       NA myData
    ## 69   2.3366  0.21474 myData
    ## 70   6.0000  1.00000 myData
    ## 71   4.6628  0.46168 myData
    ## 72   1.4359  0.19982 myData
    ## 73   4.8434  0.37816 myData
    ## 74  11.5334  0.89034 myData
    ## 75   5.4432  0.53563 myData
    ## 76   0.7825  0.11483 myData
    ## 77   0.0000      NaN myData
    ## 78  -0.0194 -0.02076 myData
    ## 79   1.0000  1.00000 myData
    ## 80   3.8181  0.56171 myData
    ## 81   5.4234  0.90945 myData
    ## 82       NA       NA myData
    ## 83   2.2323  0.23210 myData
    ## 84   5.0601  0.42999 myData
    ## 85  10.9848  0.84780 myData
    ## 86   3.8673  0.73405 myData
    ## 87   0.6870  0.06794 myData
    ## 88       NA       NA myData
    ## 89   1.0000  1.00000 myData
    ## 90   2.4507  0.23452 myData
    ## 91   0.2424  0.08162 myData
    ## 92   9.0541  0.91241 myData
    ## 93   3.6298  0.33092 myData
    ## 94      NaN      NaN myData
    ## 95   6.0000  1.00000 myData
    ## 96   6.0000  1.00000 myData
    ## 97   2.3365  0.26187 myData
    ## 98       NA       NA myData
    ## 99      NaN      NaN myData
    ## 100  2.0371  0.19023 myData
    ## 101  1.9848  0.25862 myData
    ## 102  5.2567  0.42524 myData
    ## 103 10.4067  0.80055 myData
    ## 104  6.0446  0.54306 myData
    ## 105     NaN      NaN myData
    ## 106      NA       NA myData
    ## 107     NaN      NaN myData
    ## 108  1.0000  1.00000 myData
    ## 109  0.5557  0.07358 myData
    ## 110  4.5731  0.92425 myData
    ## 111  2.1155  0.19442 myData
    ## 112 -1.0000 -0.50000 myData
    ## 113      NA       NA myData
    ## 114  5.0000  1.00000 myData
    ## 115  5.4932  0.57088 myData
    ## 116  1.5215  0.17689 myData
    ## 117  5.7556  0.63951 myData
    ## 118  1.0000  0.25000 myData
    ## 119  7.0000  1.00000 myData
    ## 120 10.0188  0.84808 myData
    ## 121  4.6991  0.94010 myData
    ## 122      NA       NA myData
    ## 123      NA       NA myData
    ## 124     NaN      NaN myData
    ## 125 10.6349  0.96746 myData
    ## 126  1.4505  0.24883 myData
    ## 127  4.7206  0.55177 myData
    ## 128      NA       NA myData
    ## 129  2.5713  0.19832 myData

``` r
#N = Number of individuals, MLG = Number of multilocus genotypes, eMLG = number of expected MLG at the smallest sample size >= 10 based on rarefaction
#SE = Standard error based on eMLG, H = Shannon-Wiener Index of MLG diversity
#G = Stoddart & Taylor's Index of MLG diversity
# lambda = Simpsons index, E.5 = Evenness, Hexp = Neis Expected Heterozygosity
#Ia = Index of association, rbarD = stand. Index of association
```

We have anything between 5 and 24 alleles per locus. aest28_01 has the
highest Simpson diversity (0.80) and aest10_1 as well as aest36_1 have
the most evenly distributed alleles (0.88).

## Genotypic evenness

Evenness is a measure of the distribution of genotype abundances,
wherein a population with equally abundant genotypes yields a value
equal to 1 and a population dominated by a single genotype is closer to
zero.

``` r
M.tab <- mlg.table(myData)
```

![](Truffles-First-Steps_files/figure-gfm/genotypic%20evenness-1.png)<!-- -->

``` r
#mll(myData)
```

``` r
setPop(myData) <- ~Pop
popdata_pop <- poppr(myData)
popdata_pop
```

    ##      Pop    N MLG eMLG    SE     H     G lambda   E.5   Hexp     Ia  rbarD
    ## 1    FRE  268  26 3.28 1.217 1.252  1.74 0.4267 0.298 0.0738  8.693 0.7373
    ## 2    ALD   71  21 5.02 1.349 1.908  3.01 0.6681 0.351 0.1414  3.223 0.3235
    ## 3    RIE   67  44 9.29 0.771 3.630 31.84 0.9686 0.840 0.4284  0.870 0.0753
    ## 4    TRO   37   9 4.21 1.087 1.440  2.72 0.6326 0.535 0.1523  2.474 0.3468
    ## 5    SCG   24  14 7.62 1.064 2.438  9.29 0.8924 0.793 0.2724  1.856 0.1613
    ## 6    BOH   54  38 8.90 0.986 3.398 19.97 0.9499 0.656 0.2856  1.370 0.1266
    ## 7    BOB  473 162 8.35 1.172 4.020 19.53 0.9488 0.339 0.3370  2.133 0.1885
    ## 8    FRB   20  14 7.91 1.041 2.441  9.09 0.8900 0.772 0.4538  5.403 0.4216
    ## 9    UEB  204  62 7.51 1.293 3.221 11.70 0.9145 0.445 0.1684  0.571 0.0620
    ## 10   SCL  298  65 5.09 1.380 2.264  3.40 0.7057 0.278 0.1073  2.529 0.2410
    ## 11   SCD  123  58 6.92 1.447 3.054  6.20 0.8388 0.258 0.2766  2.332 0.2256
    ## 12   WSL  283   4 1.30 0.502 0.168  1.07 0.0621 0.362 0.0145  4.281 0.5055
    ## 13   BUR  295  41 5.62 1.242 2.326  5.82 0.8280 0.521 0.2541  5.230 0.4409
    ## 14   NEU   38  10 4.26 1.129 1.441  2.53 0.6053 0.475 0.1031  2.179 0.2714
    ## 15   SCS   49  11 3.34 1.162 1.098  1.73 0.4223 0.366 0.0863  4.069 0.3946
    ## 16   UST   43  30 8.87 0.951 3.209 19.06 0.9475 0.760 0.3908  0.947 0.0762
    ## 17   KON  272  54 6.72 1.315 2.818  8.26 0.8789 0.461 0.2304  1.144 0.0998
    ## 18   FRI  125  16 4.14 1.064 1.612  3.40 0.7060 0.598 0.3014  4.365 0.3804
    ## 19   BAR    2   2 2.00 0.000 0.693  2.00 0.5000 1.000 0.3846     NA     NA
    ## 20   LIM  231  24 4.31 1.094 1.747  3.63 0.7243 0.555 0.4141  5.031 0.4012
    ## 21   BRU  229  24 4.56 1.244 1.773  2.96 0.6625 0.402 0.4112 11.277 0.8681
    ## 22   HAN   73  11 4.83 1.001 1.757  4.31 0.7679 0.690 0.2799  2.203 0.2308
    ## 23   GEN   11   3 2.82 0.386 0.600  1.46 0.3140 0.557 0.0792  4.899 0.9800
    ## 24 SCHIF    1   1 1.00 0.000 0.000  1.00 0.0000   NaN    NaN     NA     NA
    ## 25 Total 3291 725 9.18 0.870 5.027 47.33 0.9789 0.306 0.6526  2.571 0.1983
    ##      File
    ## 1  myData
    ## 2  myData
    ## 3  myData
    ## 4  myData
    ## 5  myData
    ## 6  myData
    ## 7  myData
    ## 8  myData
    ## 9  myData
    ## 10 myData
    ## 11 myData
    ## 12 myData
    ## 13 myData
    ## 14 myData
    ## 15 myData
    ## 16 myData
    ## 17 myData
    ## 18 myData
    ## 19 myData
    ## 20 myData
    ## 21 myData
    ## 22 myData
    ## 23 myData
    ## 24 myData
    ## 25 myData

``` r
#N = Number of individuals, MLG = Number of multilocus genotypes, eMLG = number of expected MLG at the smallest sample size >= 10 based on rarefaction
#SE = Standard error based on eMLG, H = Shannon-Wiener Index of MLG diversity
#G = Stoddart & Taylor's Index of MLG diversity
# lambda = Simpsons index, E.5 = Evenness, Hexp = Neis Expected Heterozygosity
#Ia = Index of association, rbarD = stand. Index of association
```

``` r
M.tab <- mlg.table(myData)
```

![](Truffles-First-Steps_files/figure-gfm/genotypic%20evenness%202-1.png)<!-- -->

## Visualize diversity

Diversity measures incorporate both genotypic richness and abundance.
There are three measures of genotypic diversity employed by poppr, the
Shannon-Wiener index (H), Stoddart and Taylor’s index (G), and Simpson’s
index (lambda). In our example, comparing the diversity of BB to FR
shows that H is greater for FR (4.58 vs. 4.4), but G is lower (53.4
vs. 61.7). Thus, our expectation that diversity is lower for FR than BB
is rejected in the case of H, which is likely due to the sensitivity of
the Shannon-Wiener index to genotypic richness in the uneven sample
sizes, and accepted in the case of G. To be fair, the sample size used
to calculate these diversity measures is different and is therefore not
an appropriate comparison.

For an easier statistic to grasp, we have included the Simpson index,
which is simply one minus the sum of squared genotype frequencies. This
measure provides an estimation of the probability that two randomly
selected genotypes are different and scales from 0 (no genotypes are
different) to 1 (all genotypes are different). In the data above, we can
see that lambda is just barely higher in BB than FR (0.984 vs. 0.981).
Since this might be an artifact of sample size, we can explore a
correction of Simpson’s index for sample size by multiplying lambda by
$N/(N - 1)$. Since R is vectorized, we can do this for all of our
populations at once:

``` r
popdata_pop_year <- separate(popdata,Pop,c("Pop","SamplingYear"))
```

    ## Warning: Expected 2 pieces. Missing pieces filled with `NA` in 1 rows [129].

``` r
ggplot(popdata_pop_year,aes(SamplingYear,lambda)) +
  geom_point() +
  facet_wrap(vars(Pop)) +
  theme_light() +
  labs(y="Simpson's index", title ="Simpson's index over the sampling years")
```

![](Truffles-First-Steps_files/figure-gfm/plot%20diversity%20Simpson-1.png)<!-- -->

``` r
#rarefied:
N      <- popdata_pop_year$N      # number of samples
lambda <- popdata_pop_year$lambda # Simpson's index
Simpson_rarefied <- (N/(N - 1)) * lambda              # Corrected Simpson's index

ggplot(popdata_pop_year,aes(SamplingYear,Simpson_rarefied)) +
  geom_point() +
  facet_wrap(vars(Pop)) +
  theme_light() +
  labs(y="Rarefied Simpson's index", title ="Rarefied Simpson's index over the sampling years")
```

    ## Warning: Removed 8 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](Truffles-First-Steps_files/figure-gfm/plot%20diversity%20Simpson-2.png)<!-- -->

``` r
# TruffleYears
setPop(myData) <- ~Pop/TruffleYear

popdata_pop_year <- separate(popdata,Pop,c("Pop","SamplingYear"))
```

    ## Warning: Expected 2 pieces. Missing pieces filled with `NA` in 1 rows [129].

``` r
ggplot(popdata_pop_year,aes(SamplingYear,lambda)) +
  geom_point() +
  facet_wrap(vars(Pop)) +
  theme_light() +
  labs(y="Simpson's index", title ="Simpson's index over the sampling years")
```

![](Truffles-First-Steps_files/figure-gfm/plot%20diversity%20Simpson-3.png)<!-- -->

``` r
#rarefied:
N      <- popdata_pop_year$N      # number of samples
lambda <- popdata_pop_year$lambda # Simpson's index
Simpson_rarefied <- (N/(N - 1)) * lambda              # Corrected Simpson's index

ggplot(popdata_pop_year,aes(SamplingYear,Simpson_rarefied)) +
  geom_point() +
  facet_wrap(vars(Pop)) +
  theme_light() +
  labs(y="Rarefied Simpson's index", title ="Rarefied Simpson's index over the sampling years")
```

    ## Warning: Removed 8 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](Truffles-First-Steps_files/figure-gfm/plot%20diversity%20Simpson-4.png)<!-- -->

Assess GST

Assessing genetic diversity almost always starts with an analysis of a
parameter such as GST . There are lengthy debates as to what measure of
differentiation is better (Meirmans & Hedrick, 2011). Instead of going
into that lengthy debate, it would be more worthwhile to point you into
the direction of a package dedicated to Modern Methods of
Differentiation called mmod. We will use the data set nancycats
containing 17 colonies of cats collected from Nancy, France. As cats
tend to stay within small groups, we expect to see some population
differentiation. In terms of these diversity measures, an index of GST=0
indicates no differentiation, whereas GST=1 indicates that populations
are segregating for differing alleles. Now we will use Hendrick’s
standardized GST to assess population structure among these populations
(Hedrick, 2005).

``` r
GST <- Gst_Hedrick(myData_genind)
barplot(GST$per.locus)
```

![](Truffles-First-Steps_files/figure-gfm/hedrick%20Gst-1.png)<!-- -->
Very high differentiation!

``` r
setPop(myData) <- ~Pop/SamplingYear
mlg.crosspop <- mlg.crosspop(myData,df=T)
```

    ## MLG.1: (90 inds) BUR_2016 BUR_2018 BUR_2019 BUR_2020 BUR_2021 BUR_2022
    ## MLG.2: (16 inds) BUR_2016 BUR_2019
    ## MLG.3: (20 inds) BUR_2016 BUR_2018 BUR_2019 BUR_2020
    ## MLG.4: (3 inds) BUR_2016 BUR_2022
    ## MLG.8: (5 inds) LIM_2016 BRU_2016 LIM_2018
    ## MLG.9: (2 inds) BRU_2016 BRU_2018
    ## MLG.15: (5 inds) BOB_2016 BOB_2018
    ## MLG.19: (56 inds) BUR_2014 BUR_2016 BUR_2018 BUR_2020 BUR_2021
    ## MLG.21: (2 inds) BUR_2018 BUR_2019
    ## MLG.25: (2 inds) ALD_2011 ALD_2014
    ## MLG.29: (4 inds) FRB_2011 FRB_2012
    ## MLG.30: (2 inds) RIE_2013 RIE_2014
    ## MLG.31: (4 inds) RIE_2011 RIE_2013 RIE_2014
    ## MLG.34: (29 inds) HAN_2017 HAN_2018 HAN_2019
    ## MLG.35: (5 inds) HAN_2017 HAN_2018
    ## MLG.37: (6 inds) WSL_2016 WSL_2017 WSL_2019 WSL_2020
    ## MLG.38: (275 inds) SCD_2016 WSL_2016 WSL_2017 WSL_2018 WSL_2019 WSL_2020 
    ## WSL_2021
    ## MLG.40: (12 inds) HAN_2017 HAN_2018
    ## MLG.47: (3 inds) KON_2014 KON_2015
    ## MLG.58: (129 inds) LIM_2016 BRU_2016 BRU_2017 BRU_2018
    ## MLG.59: (3 inds) BRU_2016 BRU_2017 BRU_2018
    ## MLG.60: (9 inds) BRU_2016 BRU_2017 BRU_2018
    ## MLG.61: (14 inds) BRU_2016 BRU_2017 BRU_2018
    ## MLG.62: (11 inds) BRU_2016 BRU_2017 BRU_2018
    ## MLG.69: (54 inds) FRI_2013 FRI_2016 FRI_2017 FRI_2018
    ## MLG.76: (2 inds) TRO_2011 TRO_2012
    ## MLG.77: (8 inds) ALD_2011 TRO_2011 TRO_2012
    ## MLG.79: (4 inds) BOB_2017 BOB_2019
    ## MLG.87: (5 inds) UEB_2012 UEB_2014 UEB_2015
    ## MLG.88: (3 inds) UEB_2013 UEB_2015
    ## MLG.90: (3 inds) UEB_2012 UEB_2013 UEB_2014
    ## MLG.91: (38 inds) UEB_2011 BOH_2012 UEB_2012 UEB_2013 UEB_2014 UEB_2015
    ## MLG.93: (3 inds) UEB_2014 UEB_2015
    ## MLG.94: (39 inds) UEB_2011 UEB_2013 UEB_2014 UEB_2015
    ## MLG.97: (3 inds) UEB_2014 UEB_2015
    ## MLG.101: (16 inds) UEB_2011 BOH_2012 UEB_2013 UEB_2015
    ## MLG.104: (2 inds) RIE_2011 RIE_2012
    ## MLG.121: (3 inds) BOB_2012 BOB_2015 BOB_2018
    ## MLG.123: (7 inds) BOB_2013 BOB_2016 BOB_2017 BOB_2018
    ## MLG.125: (8 inds) BOB_2012 BOB_2013 BOB_2017 BOB_2018
    ## MLG.128: (5 inds) BOB_2012 BOB_2016 BOB_2017 BOB_2018
    ## MLG.133: (4 inds) UEB_2016 UEB_2017 UEB_2018
    ## MLG.135: (3 inds) UEB_2017 UEB_2018
    ## MLG.145: (2 inds) RIE_2011 RIE_2014
    ## MLG.152: (4 inds) UST_2013 UST_2014
    ## MLG.156: (2 inds) UEB_2013 UEB_2019
    ## MLG.157: (3 inds) UEB_2011 UEB_2012
    ## MLG.158: (6 inds) UEB_2011 BOH_2012 UEB_2013 UEB_2014 UEB_2019
    ## MLG.163: (2 inds) BOB_2017 BOB_2019
    ## MLG.164: (14 inds) UEB_2011 UEB_2013 BOB_2014 UEB_2019
    ## MLG.166: (6 inds) UEB_2011 UEB_2014
    ## MLG.171: (3 inds) BOB_2016 BOB_2018
    ## MLG.172: (40 inds) LIM_2016 LIM_2017 LIM_2018
    ## MLG.175: (9 inds) GEN_2018 GEN_2020 GEN_2021 GEN_2022
    ## MLG.177: (3 inds) ALD_2011 RIE_2011
    ## MLG.179: (2 inds) BOB_2016 BOB_2017
    ## MLG.184: (4 inds) FRE_2011 FRE_2012
    ## MLG.190: (2 inds) BOB_2013 BOB_2015
    ## MLG.193: (5 inds) BOB_2016 BOB_2018
    ## MLG.197: (4 inds) UEB_2016 UEB_2017
    ## MLG.202: (54 inds) BUR_2013 BUR_2014 BUR_2015
    ## MLG.206: (3 inds) BUR_2013 BUR_2015
    ## MLG.208: (3 inds) BUR_2013 BUR_2015
    ## MLG.216: (3 inds) BUR_2013 BUR_2014
    ## MLG.235: (2 inds) LIM_2017 BRU_2018
    ## MLG.251: (21 inds) TRO_2011 TRO_2012 TRO_2013
    ## MLG.260: (38 inds) SCS_2013 SCS_2014 SCS_2016 SCS_2017 SCHIF_2016
    ## MLG.263: (3 inds) SCS_2013 SCS_2016 SCS_2017
    ## MLG.272: (3 inds) LIM_2016 LIM_2018
    ## MLG.276: (27 inds) BRU_2016 BRU_2018
    ## MLG.277: (3 inds) BRU_2016 BRU_2018
    ## MLG.293: (2 inds) LIM_2016 LIM_2018
    ## MLG.298: (2 inds) SCL_2012 SCL_2014
    ## MLG.304: (2 inds) BOB_2017 BOB_2019
    ## MLG.307: (2 inds) BOB_2017 BOB_2019
    ## MLG.308: (27 inds) BOB_2011 BOB_2012 BOB_2013 BOB_2014 BOB_2017 BOB_2018 
    ## BOB_2019
    ## MLG.309: (23 inds) BOB_2011 BOB_2012 BOB_2013 BOB_2014 BOB_2015 BOB_2016 
    ## BOB_2017
    ## MLG.311: (3 inds) BOB_2016 BOB_2017
    ## MLG.319: (2 inds) BOB_2017 BOB_2018
    ## MLG.321: (51 inds) BOB_2013 BOB_2014 BOB_2016 BOB_2017 BOB_2018 BOB_2019
    ## MLG.323: (15 inds) BOB_2012 BOB_2013 BOB_2015 BOB_2016 BOB_2017 BOB_2018 
    ## BOB_2019
    ## MLG.324: (4 inds) BOB_2013 BOB_2016 BOB_2019
    ## MLG.327: (19 inds) BOB_2012 BOB_2013 BOB_2014 BOB_2015 BOB_2017 BOB_2018
    ## MLG.328: (79 inds) BOH_2011 BOB_2011 BOB_2012 BOB_2013 BOB_2014 BOB_2015 
    ## BOB_2016 BOB_2017 BOB_2018
    ## MLG.330: (3 inds) SCD_2012 SCD_2014
    ## MLG.331: (2 inds) SCD_2012 SCD_2015
    ## MLG.333: (2 inds) BOB_2013 BOB_2018
    ## MLG.334: (4 inds) BOB_2012 BOB_2014 BOB_2016
    ## MLG.335: (3 inds) BOB_2016 BOB_2018
    ## MLG.344: (7 inds) SCG_2011 SCL_2012 SCG_2014 SCL_2017
    ## MLG.348: (47 inds) KON_2017 KON_2020 KON_2021 KON_2022 KON_2023
    ## MLG.361: (2 inds) WSL_2020 WSL_2021
    ## MLG.367: (4 inds) BOB_2013 BOB_2016
    ## MLG.373: (3 inds) BOB_2017 BOB_2019
    ## MLG.393: (5 inds) SCG_2011 SCG_2012
    ## MLG.394: (2 inds) SCG_2012 SCG_2013
    ## MLG.397: (2 inds) SCG_2011 SCG_2013
    ## MLG.398: (3 inds) FRB_2011 FRE_2013
    ## MLG.408: (8 inds) BUR_2018 BUR_2021 BUR_2022
    ## MLG.425: (11 inds) HAN_2017 HAN_2018
    ## MLG.432: (73 inds) KON_2013 KON_2014 KON_2016 KON_2017 KON_2018 KON_2019 
    ## KON_2020 KON_2022 KON_2023
    ## MLG.433: (2 inds) KON_2017 KON_2019
    ## MLG.434: (3 inds) ALD_2011 ALD_2012
    ## MLG.436: (12 inds) KON_2014 KON_2016 KON_2017 KON_2019 KON_2022
    ## MLG.437: (11 inds) KON_2016 KON_2018
    ## MLG.439: (101 inds) LIM_2016 LIM_2017 LIM_2018 BRU_2018
    ## MLG.441: (7 inds) FRI_2013 FRI_2014 FRI_2017
    ## MLG.443: (2 inds) LIM_2016 LIM_2018
    ## MLG.444: (4 inds) LIM_2016 LIM_2017 LIM_2018
    ## MLG.460: (3 inds) UEB_2012 UEB_2013
    ## MLG.469: (3 inds) RIE_2011 RIE_2012 RIE_2014
    ## MLG.476: (3 inds) KON_2016 KON_2017 KON_2019
    ## MLG.479: (23 inds) NEU_2013 NEU_2014 NEU_2016 NEU_2017 NEU_2020 NEU_2021
    ## MLG.483: (5 inds) NEU_2015 NEU_2016 NEU_2020 NEU_2021
    ## MLG.486: (3 inds) NEU_2013 NEU_2016
    ## MLG.493: (19 inds) SCL_2012 SCL_2014 SCL_2016 SCL_2017
    ## MLG.496: (5 inds) SCD_2016 SCD_2017
    ## MLG.498: (2 inds) SCD_2014 SCD_2015
    ## MLG.505: (34 inds) SCL_2012 SCL_2014 SCL_2016 SCL_2017
    ## MLG.506: (2 inds) SCL_2012 SCL_2013
    ## MLG.511: (10 inds) HAN_2016 HAN_2017 HAN_2018
    ## MLG.513: (2 inds) RIE_2011 RIE_2014
    ## MLG.514: (4 inds) RIE_2011 RIE_2013
    ## MLG.518: (6 inds) SCL_2012 SCL_2014
    ## MLG.525: (2 inds) BOH_2014 BOH_2015
    ## MLG.529: (157 inds) SCL_2012 SCL_2013 SCL_2014 SCL_2015 SCD_2016 SCL_2016 
    ## SCL_2017 SCL_2018
    ## MLG.530: (2 inds) SCL_2014 SCL_2016
    ## MLG.535: (11 inds) SCL_2014 SCL_2016 SCL_2017
    ## MLG.537: (7 inds) SCL_2015 SCL_2016 SCL_2018
    ## MLG.550: (2 inds) RIE_2012 RIE_2014
    ## MLG.552: (4 inds) RIE_2012 RIE_2013 RIE_2014
    ## MLG.554: (5 inds) BOB_2013 BOB_2016 BOB_2017 BOB_2018
    ## MLG.555: (5 inds) BOB_2013 BOB_2015 BOB_2017
    ## MLG.556: (5 inds) BOB_2013 BOB_2015 BOB_2017 BOB_2018
    ## MLG.560: (2 inds) KON_2014 KON_2018
    ## MLG.561: (20 inds) KON_2017 KON_2018 KON_2019 KON_2020 KON_2023
    ## MLG.562: (8 inds) KON_2015 KON_2016 KON_2018 KON_2020 KON_2022
    ## MLG.563: (6 inds) KON_2018 KON_2019
    ## MLG.564: (19 inds) KON_2014 KON_2016 KON_2017 KON_2018 KON_2019 KON_2020
    ## MLG.565: (4 inds) KON_2016 KON_2018 KON_2022
    ## MLG.568: (2 inds) KON_2017 KON_2018
    ## MLG.573: (2 inds) KON_2016 KON_2017
    ## MLG.574: (2 inds) KON_2016 KON_2018
    ## MLG.575: (13 inds) KON_2016 KON_2019 KON_2020 KON_2021
    ## MLG.577: (2 inds) KON_2016 KON_2017
    ## MLG.578: (6 inds) KON_2016 KON_2017 KON_2019 KON_2021
    ## MLG.595: (2 inds) BOB_2012 BOB_2016
    ## MLG.596: (7 inds) BOB_2017 BOB_2019
    ## MLG.598: (17 inds) BOB_2013 BOB_2016 BOB_2017 BOB_2018 BOB_2019
    ## MLG.599: (3 inds) BOB_2017 BOB_2018
    ## MLG.616: (2 inds) SCD_2014 SCD_2015
    ## MLG.617: (49 inds) SCD_2012 SCL_2014 SCD_2014 SCD_2015 SCD_2016 SCD_2017
    ## MLG.620: (3 inds) SCD_2013 SCD_2015 SCD_2016
    ## MLG.621: (3 inds) SCD_2012 SCD_2015 SCD_2016
    ## MLG.622: (9 inds) BOH_2011 BOH_2012 BOH_2013 BOH_2014
    ## MLG.623: (2 inds) BOH_2011 BOH_2012
    ## MLG.626: (3 inds) BOH_2012 BOH_2013 BOH_2014
    ## MLG.628: (3 inds) BOH_2012 BOH_2013
    ## MLG.634: (57 inds) LIM_2016 LIM_2017 LIM_2018 BRU_2018
    ## MLG.635: (3 inds) LIM_2016 LIM_2017
    ## MLG.641: (38 inds) FRI_2013 FRI_2014 FRI_2016 FRI_2017 FRI_2018
    ## MLG.643: (2 inds) SCL_2014 SCG_2014
    ## MLG.649: (4 inds) RIE_2012 RIE_2013
    ## MLG.656: (3 inds) RIE_2011 RIE_2012 RIE_2013
    ## MLG.660: (11 inds) BRU_2016 BRU_2017 BRU_2018
    ## MLG.676: (4 inds) FRB_2012 FRB_2013 FRB_2014
    ## MLG.685: (4 inds) ALD_2011 ALD_2013 ALD_2015
    ## MLG.689: (5 inds) ALD_2012 ALD_2013
    ## MLG.694: (40 inds) ALD_2011 ALD_2012 ALD_2013 ALD_2014 ALD_2015
    ## MLG.709: (3 inds) FRE_2013 FRE_2014
    ## MLG.710: (202 inds) FRE_2011 FRE_2012 FRE_2013 FRE_2014 FRE_2015 FRE_2016 
    ## FRE_2017
    ## MLG.712: (5 inds) FRE_2011 FRE_2013
    ## MLG.714: (9 inds) FRE_2011 FRE_2014 FRE_2016 FRE_2017
    ## MLG.715: (4 inds) FRE_2016 FRE_2017
    ## MLG.716: (3 inds) FRE_2016 FRE_2017
    ## MLG.717: (13 inds) FRE_2016 FRE_2017
    ## MLG.718: (3 inds) FRE_2016 FRE_2017

``` r
crosspop <- mlg.crosspop %>%
  separate_wider_delim(Population,delim="_",names=c("Pop","SamplingYear"), too_few="align_end")

ggplot(subset(crosspop, Pop %in% "BOB"),aes(SamplingYear,MLG)) + geom_point() + geom_line(aes(group=MLG)) +
  labs(title ="Distribution of MLGs in BOB (Bohlingen Buche) over the years")
```

![](Truffles-First-Steps_files/figure-gfm/plot%20persistence-1.png)<!-- -->

``` r
ggplot(subset(crosspop, Pop %in% "SCL"),aes(SamplingYear,MLG)) + geom_point()  + geom_line(aes(group=MLG)) +
  labs(title ="Distribution of MLGs in SCL (Schlieren Grillplatz) over the years")
```

![](Truffles-First-Steps_files/figure-gfm/plot%20persistence-2.png)<!-- -->

``` r
ggplot(subset(crosspop, Pop %in% "ALD"),aes(SamplingYear,MLG)) + geom_point()  + geom_line(aes(group=MLG)) +
  labs(title ="Distribution of MLGs in ALD (Aldingen) over the years")
```

![](Truffles-First-Steps_files/figure-gfm/plot%20persistence-3.png)<!-- -->

``` r
ggplot(subset(crosspop, Pop %in% "BOH"),aes(SamplingYear,MLG)) + geom_point()  + geom_line(aes(group=MLG)) +
  labs(title ="Distribution of MLGs in BOH (Bohlingen Hasel) over the years")
```

![](Truffles-First-Steps_files/figure-gfm/plot%20persistence-4.png)<!-- -->

``` r
ggplot(subset(crosspop, Pop %in% "WSL"),aes(SamplingYear,MLG)) + geom_point()  + geom_line(aes(group=MLG)) +
  labs(title ="Distribution of MLGs in WSL over the years")
```

![](Truffles-First-Steps_files/figure-gfm/plot%20persistence-5.png)<!-- -->

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-kamvarPopprPackageGenetic2014" class="csl-entry">

Kamvar, Zhian N., Javier F. Tabima, and Niklaus J. Grünwald. 2014.
“Poppr: an R package for genetic analysis of populations with clonal,
partially clonal, and/or sexual reproduction.” *PeerJ* 2: e281.
<https://doi.org/10.7717/peerj.281>.

</div>

<div id="ref-legendreComparisonMantelTest2010" class="csl-entry">

Legendre, Pierre, and Marie-Josée Fortin. 2010. “Comparison of the
Mantel Test and Alternative Approaches for Detecting Complex
Multivariate Relationships in the Spatial Analysis of Genetic Data.”
*Molecular Ecology Resources* 10 (5): 831–44.
<https://doi.org/10.1111/j.1755-0998.2010.02866.x>.

</div>

<div id="ref-molinierSSRbasedIdentificationGenetic2016"
class="csl-entry">

Molinier, Virgine, Claude Murat, Martina Peter, Armelle Gollotte,
Herminia De la Varga, Barbara Meier, Simon Egli, Beatrice Belfiori,
Francesco Paolocci, and Daniel Wipf. 2016. “SSR-Based Identification of
Genetic Groups Within European Populations of &Lt;em&gt;Tuber
Aestivum&lt;/Em&gt; Vittad.” *Mycorrhiza*, 99–110.
<https://doi.org/10.1007/s00572-015-0649-0>.

</div>

<div id="ref-molinierFinescaleGeneticStructure2016" class="csl-entry">

Molinier, Virginie, Claude Murat, Andri Baltensweiler, Ulf Büntgen,
Francis Martin, Barbara Meier, Barbara Moser, et al. 2016. “Fine-scale
genetic structure of natural Tuber aestivum sites in southern Germany.”
*Mycorrhiza* 26 (8): 895–907.
<https://doi.org/10.1007/s00572-016-0719-y>.

</div>

<div id="ref-molinierFirstIdentificationPolymorphic2013"
class="csl-entry">

Molinier, Virginie, Claude Murat, Emmanuelle Morin, Armelle Gollotte,
Daniel Wipf, and Francis Martin. 2013. “First Identification of
Polymorphic Microsatellite Markers in the Burgundy Truffle, Tuber
Aestivum (Tuberaceae)1.” *Applications in Plant Sciences* 1 (2):
apps.1200220. <https://doi.org/10.3732/apps.1200220>.

</div>

<div id="ref-staubliHiddenFairyRings2022" class="csl-entry">

Staubli, Florian, Lea Imola, Benjamin Dauphin, Virginie Molinier,
Stephanie Pfister, Yasmine Piñuela, Laura Schürz, et al. 2022. “Hidden
Fairy Rings and Males—Genetic Patterns of Natural Burgundy Truffle (
*Tuber Aestivum* Vittad.) Populations Reveal New Insights into Its Life
Cycle.” *Environmental Microbiology* 24 (12): 6376–91.
<https://doi.org/10.1111/1462-2920.16131>.

</div>

<div id="ref-steidingerFallSummerTruffle2022" class="csl-entry">

Steidinger, Brian S., Ulf Büntgen, Uli Stobbe, Willy Tegel, Ludger
Sproll, Matthias Haeni, Barbara Moser, et al. 2022. “The Fall of the
Summer Truffle: Recurring Hot, Dry Summers Result in Declining Fruitbody
Production of Tuber Aestivum in Central Europe.” *Global Change Biology*
28 (24): 7376–90. <https://doi.org/10.1111/gcb.16424>.

</div>

</div>
