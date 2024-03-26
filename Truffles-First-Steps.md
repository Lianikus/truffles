Truffle Monitoring
================
Lia Baumann
2024-03-05

## Basis

The basis for the research about Truffle Monitoring data is described in
detail in these publications: Virginie Molinier et al. (2013) Virgine
Molinier et al. (2016) Virginie Molinier et al. (2016) Staubli et al.
(2022) Steidinger et al. (2022) Legendre and Fortin (2010)
(**kamvarPopprPackageGenetic2014?**)

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
| BOB          | Bohlingen Buche       | 477 |
| BOH          | Bohlingen Hasel       |  54 |
| BUR          | Bursins               | 284 |
| FRB          | Freiburg-Brumaleweg   |  20 |
| FRE          | Freiburg-Wittnau      | 268 |
| FRI          | Frick                 | 125 |
| GEN          | Genolier              |   5 |
| KON          | Köniz                 | 234 |
| NEU          | Neuchatel             |  38 |
| RIE          | Rietheim              |  67 |
| SCD          | Schlieren Dreieck     | 123 |
| SCG          | Schlieren Graben      |  24 |
| SCL          | Schlieren Grillplatz  | 299 |
| SCS          | Schiffenensee         |  50 |
| TRO          | Trossingen            |  37 |
| UEB          | Bohlingen Ueberlingen | 204 |
| UST          | Uster                 |  43 |
| WSL          | WSL                   | 283 |

Number of observations per site

![](Truffles-First-Steps_files/figure-gfm/Data%20structure%20per%20site-1.png)<!-- -->

Evtl. später: Standorte mit weniger als 10 Samples werden entfernt. Dies
lässt danach noch folgende Standorte für die Auswertung zu:

\#Idee: alle Observationen auf einer Zeitachse und nach Standort
unterteilt anzeigen

![](Truffles-First-Steps_files/figure-gfm/all%20observations-1.png)<!-- -->

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

    ## 
    ## This is a genclone object
    ## -------------------------
    ## Genotype information:
    ## 
    ##     676 original multilocus genotypes 
    ##    2708 haploid individuals
    ##      14 codominant loci
    ## 
    ## Population information:
    ## 
    ##       1 stratum - Pop
    ##      20 populations defined - ALD, BAR, BOB, ..., UEB, UST, WSL

Show missing data

![](Truffles-First-Steps_files/figure-gfm/missing%20data-1.png)<!-- -->

    ##           Locus
    ## Population aest01_1 aest06_1 aest07_1 aest10_1 aest15_1 aest18_1 aest24_1
    ##      ALD    0.07042        .        .        .        .        .  0.04225
    ##      BAR    0.50000        .        .        .        .        .        .
    ##      BOB    0.06709  0.00210        .  0.00210        .  0.00210  0.00629
    ##      BOH    0.22222        .        .  0.01852        .  0.01852  0.05556
    ##      BUR    0.19366        .        .  0.00352        .        .  0.01761
    ##      FRB    0.10000        .        .  0.20000        .        .        .
    ##      FRE    0.13433        .        .  0.11194        .        .  0.01493
    ##      FRI    0.03200        .        .  0.00800        .        .  0.00800
    ##      GEN          .        .        .        .        .        .        .
    ##      KON    0.21795        .  0.00427        .  0.00427        .  0.00855
    ##      NEU    0.07895        .        .  0.02632        .        .        .
    ##      RIE    0.04478        .        .  0.01493        .        .        .
    ##      SCD    0.13008        .        .  0.02439        .        .  0.02439
    ##      SCG    0.25000        .        .  0.04167        .        .  0.04167
    ##      SCL    0.09365        .        .  0.01003  0.00334  0.00334  0.00334
    ##      SCS    0.10000        .        .        .        .        .  0.02000
    ##      TRO    0.08108        .        .        .        .        .  0.02703
    ##      UEB    0.20098        .        .  0.02451  0.12255        .  0.02451
    ##      UST    0.11628        .        .  0.02326        .        .  0.02326
    ##      WSL          .        .        .        .        .        .        .
    ##      Total  0.11374  0.00037  0.00037  0.01957  0.00997  0.00111  0.01256
    ##           Locus
    ## Population aest25_1 aest26_1 aest28_1 aest29_1 aest31_1 aest35_1 aest36_1
    ##      ALD          .        .        .        .  0.07042        .        .
    ##      BAR          .        .        .        .        .        .        .
    ##      BOB    0.00419        .  0.01887        .  0.03145        .  0.00210
    ##      BOH          .        .  0.07407  0.01852  0.07407        .        .
    ##      BUR          .  0.02465  0.00704  0.00704  0.11620  0.00704  0.00352
    ##      FRB          .        .  0.20000  0.05000  0.10000        .        .
    ##      FRE          .        .  0.01866  0.00373  0.11194        .        .
    ##      FRI          .        .  0.00800        .  0.00800        .        .
    ##      GEN          .        .        .        .        .        .        .
    ##      KON          .        .  0.01282  0.00427  0.14530        .  0.00855
    ##      NEU    0.02632        .  0.05263        .  0.02632        .        .
    ##      RIE          .        .        .        .  0.02985        .        .
    ##      SCD          .        .  0.05691  0.01626  0.07317        .        .
    ##      SCG          .        .  0.20833        .  0.12500        .        .
    ##      SCL          .        .  0.63545  0.00669  0.05686        .        .
    ##      SCS          .        .  0.02000  0.02000  0.02000        .        .
    ##      TRO          .        .        .        .  0.02703        .        .
    ##      UEB    0.00490        .  0.05882  0.01471  0.15686        .        .
    ##      UST          .        .  0.04651  0.02326  0.04651        .        .
    ##      WSL          .        .        .        .        .        .        .
    ##      Total  0.00148  0.00258  0.09121  0.00554  0.07090  0.00074  0.00148
    ##           Locus
    ## Population    Mean
    ##      ALD   0.01308
    ##      BAR   0.03571
    ##      BOB   0.00973
    ##      BOH   0.03439
    ##      BUR   0.02716
    ##      FRB   0.04643
    ##      FRE   0.02825
    ##      FRI   0.00457
    ##      GEN         .
    ##      KON   0.02900
    ##      NEU   0.01504
    ##      RIE   0.00640
    ##      SCD   0.02323
    ##      SCG   0.04762
    ##      SCL   0.05805
    ##      SCS   0.01286
    ##      TRO   0.00965
    ##      UEB   0.04342
    ##      UST   0.01993
    ##      WSL         .
    ##      Total 0.02369

Clone correction

    ## 
    ## This is a genclone object
    ## -------------------------
    ## Genotype information:
    ## 
    ##    676 original multilocus genotypes 
    ##    690 haploid individuals
    ##     14 codominant loci
    ## 
    ## Population information:
    ## 
    ##      1 stratum - Pop
    ##     20 populations defined - ALD, BAR, BOB, ..., UEB, UST, WSL

    ##           summary
    ## locus      allele      1-D     Hexp Evenness
    ##   aest01_1      .  0.03161  0.03058  0.08373
    ##   aest06_1      .  0.14621  0.14571  0.05797
    ##   aest07_1      .  0.06025  0.05966  0.00091
    ##   aest10_1      .  0.04087  0.04005  0.05927
    ##   aest15_1      .  0.14840  0.14824  0.12356
    ##   aest18_1      .  0.01211  0.01135  0.03571
    ##   aest24_1      .  0.06898  0.06849  0.06678
    ##   aest25_1      .  0.02541  0.02480 -0.00906
    ##   aest26_1      . -0.02265 -0.02347  0.01175
    ##   aest28_1      . -0.03892 -0.03995 -0.00317
    ##   aest29_1      .  0.05603  0.05537  0.04754
    ##   aest31_1      . -0.02566 -0.02656  0.04152
    ##   aest35_1      .  0.05544  0.05496  0.11352
    ##   aest36_1      .  0.01563  0.01507  0.07482
    ##   mean          .  0.04098  0.04031  0.05035

![](Truffles-First-Steps_files/figure-gfm/compare%20diversity%20between%20corrected%20and%20uncorrected-1.png)<!-- -->

The graph shows a decrease of diversity for most markers when
clone-correcting the data (increase of Simpson index means a decrease of
genotypic diversity).

``` r
popdata <- poppr(myData)
popdata
```

    ##      Pop    N MLG eMLG    SE     H     G lambda   E.5   Hexp    Ia  rbarD
    ## 1    ALD   71  21 5.02 1.349 1.908  3.01 0.6681 0.351 0.1414 3.223 0.3235
    ## 2    BAR    2   2 2.00 0.000 0.693  2.00 0.5000 1.000 0.3846    NA     NA
    ## 3    BOB  477 172 8.50 1.133 4.134 21.95 0.9544 0.341 0.3425 2.157 0.1854
    ## 4    BOH   54  38 8.90 0.986 3.398 19.97 0.9499 0.656 0.2871 1.341 0.1240
    ## 5    BUR  284  41 5.83 1.251 2.404  6.37 0.8431 0.534 0.2610 4.980 0.4192
    ## 6    FRB   20  14 7.91 1.041 2.441  9.09 0.8900 0.772 0.4538 5.403 0.4216
    ## 7    FRE  268  26 3.28 1.217 1.252  1.74 0.4267 0.298 0.0738 8.693 0.7373
    ## 8    FRI  125  16 4.14 1.064 1.612  3.40 0.7060 0.598 0.3014 4.365 0.3804
    ## 9    GEN    5   2 2.00 0.000 0.500  1.47 0.3200 0.725 0.1714 5.000 1.0000
    ## 10   KON  234  49 6.90 1.297 2.852  8.70 0.8851 0.472 0.2168 1.070 0.0922
    ## 11   NEU   38  11 4.37 1.179 1.492  2.55 0.6080 0.450 0.1089 2.088 0.2622
    ## 12   RIE   67  44 9.29 0.771 3.630 31.84 0.9686 0.840 0.4284 0.870 0.0753
    ## 13   SCD  123  58 6.92 1.447 3.054  6.20 0.8388 0.258 0.2766 2.332 0.2256
    ## 14   SCG   24  14 7.62 1.064 2.438  9.29 0.8924 0.793 0.2724 1.856 0.1613
    ## 15   SCL  299  65 5.08 1.380 2.258  3.38 0.7041 0.278 0.1080 2.629 0.2494
    ## 16   SCS   50  11 3.30 1.157 1.082  1.71 0.4152 0.364 0.0846 4.079 0.3958
    ## 17   TRO   37   9 4.21 1.087 1.440  2.72 0.6326 0.535 0.1523 2.474 0.3468
    ## 18   UEB  204  62 7.51 1.293 3.221 11.70 0.9145 0.445 0.1684 0.571 0.0621
    ## 19   UST   43  30 8.87 0.951 3.209 19.06 0.9475 0.760 0.3908 0.947 0.0762
    ## 20   WSL  283   5 1.30 0.505 0.173  1.07 0.0621 0.351 0.0145 4.281 0.5054
    ## 21 Total 2708 676 9.04 0.945 4.971 38.79 0.9742 0.264 0.6329 2.459 0.1895
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

``` r
#N = Number of individuals, MLG = Number of multilocus genotypes, eMLG = number of expected MLG at the smallest sample size >= 10 based on rarefaction
#SE = Standard error based on eMLG, H = Shannon-Wiener Index of MLG diversity
#G = Stoddart & Taylor's Index of MLG diversity
# lambda = Simpsons index, E.5 = Evenness, Hexp = Neis Expected Heterozygosity
#Ia = Index of association, rbarD = stand. Index of association
```

<div id="refs" class="references csl-bib-body hanging-indent">

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
