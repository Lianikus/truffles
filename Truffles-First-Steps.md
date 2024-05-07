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

    ## /// GENIND OBJECT /////////
    ## 
    ##  // 2,708 individuals; 14 loci; 127 alleles; size: 1.5 Mb
    ## 
    ##  // Basic content
    ##    @tab:  2708 x 127 matrix of allele counts
    ##    @loc.n.all: number of alleles per locus (range: 5-24)
    ##    @loc.fac: locus factor for the 127 columns of @tab
    ##    @all.names: list of allele names for each locus
    ##    @ploidy: ploidy of each individual  (range: 1-1)
    ##    @type:  codom
    ##    @call: read.genalex(genalex = "Daten_Genalex.csv", ploidy = 1, genclone = FALSE, 
    ##     sep = ";")
    ## 
    ##  // Optional content
    ##    @pop: population of each individual (group size range: 2-477)
    ##    @strata: a data frame with 1 columns ( Pop )

![](Truffles-First-Steps_files/figure-gfm/load%20genalex%20data%20and%20first%20examinations-1.png)<!-- -->

    ## 
    ## // Number of individuals: 2708
    ## // Group sizes: 71 2 477 54 284 20 268 125 5 234 38 67 123 24 299 50 37 204 43 283
    ## // Number of alleles per locus: 9 7 6 9 6 5 10 5 12 24 9 10 7 8
    ## // Number of alleles per group: 35 19 68 38 51 39 42 43 20 40 28 44 37 31 54 34 23 43 38 26
    ## // Percentage of missing data: 2.37 %
    ## // Observed heterozygosity: 0

    ## 
    ## // Number of individuals: 2708
    ## // Group sizes: 71 2 477 54 284 20 268 125 5 234 38 67 123 24 299 50 37 204 43 283
    ## // Number of alleles per locus: 9 7 6 9 6 5 10 5 12 24 9 10 7 8
    ## // Number of alleles per group: 35 19 68 38 51 39 42 43 20 40 28 44 37 31 54 34 23 43 38 26
    ## // Percentage of missing data: 2.37 %
    ## // Observed heterozygosity: 0

    ## NULL

![](Truffles-First-Steps_files/figure-gfm/load%20genalex%20data%20and%20first%20examinations-2.png)<!-- -->

Show missing data

Clone correction

The graph shows a decrease of diversity for most markers when
clone-correcting the data (increase of Simpson index means a decrease of
genotypic diversity).

Next step: I want to compare genetic diversity across years and show it
for populations. So I need to combine the sample and assign the sampling
date. I have several ideas how to do it. 1. –\> Join Genalex file with
monitoring file on the sampling column. –\> then read it in as a
genclone and genind file again. 2. –\> Filter monitoring file with
genalex file to exclude the rows not matching –\> choose columns that I
want to process as a new file and read in as genclone and genind file
again.

``` r
genalex_dates <- read.csv("Daten_Genalex.csv",sep=";", header=FALSE)
T_all_dates <- T_all %>%
  mutate(Sampling_month = month(Sampling_date)) %>%
  select(Code_Analyses_2024,Sampling_month, Sampling_year)
genalex_dates <- left_join(genalex_dates,T_all_dates,by=c("V1"="Code_Analyses_2024"))

genalex_dates <- genalex_dates %>%
  unite(col="pop_date",c("V1","Sampling_month","Sampling_year"),sep="_")
#correct column names
genalex_dates$pop_date[genalex_dates$pop_date =="2708_NA_NA"] <- "2708"
genalex_dates$pop_date[genalex_dates$pop_date =="_NA_NA"] <- ""
genalex_dates$pop_date[genalex_dates$pop_date =="pop_NA_NA"] <- "pop"

#jetzt als csv Datei neu reinladen
#write.table(genalex_dates, col.names=FALSE, sep=",", "C:/Users/liaba/OneDrive - Eidg. Forschungsanstalt WSL/R/truffles/genalex_dates.csv")

microsats_dates <- read.genalex("genalex_dates.csv",ploidy=1)
splitStrata(microsats_dates) <- ~Pop/Month/Year
microsats_dates
```

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
    ##       3 strata - Pop, Month, Year
    ##     506 populations defined - 
    ## ALD_3_2011, ALD_6_2011, ALD_7_2011, ..., WSL_11_2020, WSL_6_2021, WSL_7_2021

``` r
#treemap:
monstrata <- strata(microsats_dates) %>%     
  group_by(Pop,Month,Year) %>%
  summarize(Count = n())
```

    ## `summarise()` has grouped output by 'Pop', 'Month'. You can override using the
    ## `.groups` argument.

``` r
monstrata
```

    ## # A tibble: 506 × 4
    ## # Groups:   Pop, Month [174]
    ##    Pop   Month Year  Count
    ##    <fct> <fct> <fct> <int>
    ##  1 ALD   3     2011      4
    ##  2 ALD   3     2012      2
    ##  3 ALD   3     2013      2
    ##  4 ALD   6     2011      2
    ##  5 ALD   7     2011      2
    ##  6 ALD   7     2014      1
    ##  7 ALD   8     2011      3
    ##  8 ALD   8     2014      2
    ##  9 ALD   10    2011      2
    ## 10 ALD   10    2012      1
    ## # ℹ 496 more rows

``` r
library(treemap)
treemap(dtf = monstrata, index = nameStrata(microsats_dates), vSize = "Count",
        type = "categorical", vColor = "Pop", title = "Truffles")
```

![](Truffles-First-Steps_files/figure-gfm/add%20dates%20into%20microsat%20file-1.png)<!-- -->

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
gac <- genotype_curve(microsats_dates, sample = 1000, quiet = TRUE)
```

![](Truffles-First-Steps_files/figure-gfm/genotype%20accumulation%20curve-1.png)<!-- -->

We have anything between 5 and 24 alleles per locus. aest28_01 has the
highest Simpson diversity (0.80) and aest10_1 as well as aest36_1 have
the most evenly distributed alleles (0.88).

``` r
(microsats_lt <- locus_table(microsats_dates))
```

    ## 
    ## allele = Number of observed alleles
    ## 1-D = Simpson index
    ## Hexp = Nei's 1978 gene diversity
    ## ------------------------------------------

    ##           summary
    ## locus      allele   1-D  Hexp Evenness
    ##   aest01_1   9.00  0.74  0.74     0.87
    ##   aest06_1   7.00  0.66  0.66     0.62
    ##   aest07_1   6.00  0.63  0.63     0.80
    ##   aest10_1   9.00  0.77  0.77     0.88
    ##   aest15_1   6.00  0.35  0.35     0.57
    ##   aest18_1   5.00  0.71  0.71     0.86
    ##   aest24_1  10.00  0.52  0.52     0.51
    ##   aest25_1   5.00  0.59  0.60     0.77
    ##   aest26_1  12.00  0.73  0.73     0.85
    ##   aest28_1  24.00  0.80  0.80     0.76
    ##   aest29_1   9.00  0.67  0.67     0.77
    ##   aest31_1  10.00  0.65  0.65     0.86
    ##   aest35_1   7.00  0.52  0.52     0.70
    ##   aest36_1   8.00  0.53  0.53     0.88
    ##   mean       9.07  0.63  0.63     0.76

``` r
info_table(microsats_dates, type = "missing", plot = TRUE)
```

![](Truffles-First-Steps_files/figure-gfm/allele%20frequencies,%20missing%20data,%20ploidy-1.png)<!-- -->

    ##              Locus
    ## Population    aest01_1 aest06_1 aest07_1 aest10_1 aest15_1 aest18_1 aest24_1
    ##   ALD_3_2011         .        .        .        .        .        .        .
    ##   ALD_6_2011         .        .        .        .        .        .        .
    ##   ALD_7_2011         .        .        .        .        .        .  0.50000
    ##   ALD_8_2011         .        .        .        .        .        .        .
    ##   ALD_10_2011        .        .        .        .        .        .        .
    ##   ALD_11_2011        .        .        .        .        .        .        .
    ##   ALD_12_2011        .        .        .        .        .        .        .
    ##   ALD_1_2012         .        .        .        .        .        .        .
    ##   ALD_2_2012         .        .        .        .        .        .        .
    ##   ALD_3_2012   0.50000        .        .        .        .        .        .
    ##   ALD_10_2012        .        .        .        .        .        .        .
    ##   ALD_11_2012        .        .        .        .        .        .        .
    ##   ALD_12_2012        .        .        .        .        .        .        .
    ##   ALD_1_2013         .        .        .        .        .        .        .
    ##   ALD_2_2013         .        .        .        .        .        .        .
    ##   ALD_3_2013   0.50000        .        .        .        .        .        .
    ##   ALD_10_2013        .        .        .        .        .        .        .
    ##   ALD_4_2011         .        .        .        .        .        .        .
    ##   ALD_11_2013        .        .        .        .        .        .        .
    ##   ALD_12_2013  0.25000        .        .        .        .        .        .
    ##   ALD_1_2014   0.50000        .        .        .        .        .  1.00000
    ##   ALD_2_2014         .        .        .        .        .        .        .
    ##   ALD_5_2011         .        .        .        .        .        .        .
    ##   ALD_7_2014         .        .        .        .        .        .        .
    ##   ALD_8_2014         .        .        .        .        .        .        .
    ##   ALD_9_2014         .        .        .        .        .        .        .
    ##   ALD_10_2014        .        .        .        .        .        .        .
    ##   ALD_1_2015         .        .        .        .        .        .        .
    ##   ALD_12_2015  0.25000        .        .        .        .        .        .
    ##   BAR_5_2014   1.00000        .        .        .        .        .        .
    ##   BAR_6_2014         .        .        .        .        .        .        .
    ##   BOB_8_2011         .        .        .        .        .        .        .
    ##   BOB_2_2012         .        .        .        .        .        .        .
    ##   BOB_1_2015         .        .        .        .        .        .        .
    ##   BOB_3_2015   0.25000        .        .        .        .        .        .
    ##   BOB_6_2015         .        .        .        .        .        .        .
    ##   BOB_8_2015   0.23077        .        .        .        .        .        .
    ##   BOB_3_2012         .        .        .        .        .        .        .
    ##   BOB_10_2015  0.50000        .        .        .        .        .        .
    ##   BOB_1_2016         .        .        .        .        .        .        .
    ##   BOB_2_2016   1.00000        .        .        .        .        .        .
    ##   BOB_11_2012        .        .        .        .        .        .        .
    ##   BOB_9_2014         .        .        .        .        .        .  0.12500
    ##   BOB_10_2013  0.02564        .        .        .        .  0.02564        .
    ##   BOB_5_2012         .        .        .        .        .        .        .
    ##   BOB_7_2012   0.14286        .        .        .        .        .        .
    ##   BOB_6_2016   1.00000        .        .        .        .        .        .
    ##   BOB_7_2016         .        .        .        .        .        .        .
    ##   BOB_8_2016   0.28571        .        .        .        .        .        .
    ##   BOB_9_2016         .        .        .        .        .        .        .
    ##   BOB_10_2016        .        .        .        .        .        .        .
    ##   BOB_9_2011         .        .        .        .        .        .        .
    ##   BOB_11_2016  0.07692  0.07692        .        .        .        .        .
    ##   BOB_8_2012         .        .        .        .        .        .        .
    ##   BOB_12_2016        .        .        .        .        .        .        .
    ##   BOB_9_2012   0.18182        .        .        .        .        .        .
    ##   BOB_6_2017         .        .        .        .        .        .        .
    ##   BOB_7_2017         .        .        .        .        .        .        .
    ##   BOB_8_2017         .        .        .        .        .        .        .
    ##   BOB_9_2017         .        .        .        .        .        .        .
    ##   BOB_10_2011        .        .        .        .        .        .        .
    ##   BOB_10_2017  0.05556        .        .        .        .        .        .
    ##   BOB_11_2017        .        .        .        .        .        .        .
    ##   BOB_10_2012        .        .        .        .        .        .  0.14286
    ##   BOB_11_2011        .        .        .        .        .        .        .
    ##   BOB_12_2017        .        .        .        .        .        .        .
    ##   BOB_1_2018         .        .        .        .        .        .        .
    ##   BOB_3_2018         .        .        .        .        .        .        .
    ##   BOB_4_2018         .        .        .        .        .        .        .
    ##   BOB_6_2018         .        .        .        .        .        .        .
    ##   BOB_8_2018         .        .        .        .        .        .        .
    ##   BOB_1_2013         .        .        .        .        .        .        .
    ##   BOB_9_2018         .        .        .        .        .        .        .
    ##   BOB_2_2013   1.00000        .        .  1.00000        .        .        .
    ##   BOB_5_2013   0.50000        .        .        .        .        .        .
    ##   BOB_10_2018  0.38095        .        .        .        .        .        .
    ##   BOB_12_2011        .        .        .        .        .        .        .
    ##   BOB_6_2013   1.00000        .        .        .        .        .        .
    ##   BOB_8_2013         .        .        .        .        .        .        .
    ##   BOB_9_2019         .        .        .        .        .        .        .
    ##   BOB_9_2013         .        .        .        .        .        .        .
    ##   BOB_12_2013        .        .        .        .        .        .        .
    ##   BOB_1_2014         .        .        .        .        .        .        .
    ##   BOB_2_2014         .        .        .        .        .        .        .
    ##   BOB_3_2014         .        .        .        .        .        .        .
    ##   BOB_7_2014   0.25000        .        .        .        .        .  0.12500
    ##   BOB_8_2014         .        .        .        .        .        .        .
    ##   BOB_10_2014  0.33333        .        .        .        .        .        .
    ##   BOB_11_2014  0.27273        .        .        .        .        .        .
    ##   BOB_12_2014        .        .        .        .        .        .        .
    ##   BOH_5_2011         .        .        .        .        .        .        .
    ##   BOH_10_2011        .        .        .        .        .        .        .
    ##   BOH_12_2011        .        .        .        .        .        .        .
    ##   BOH_3_2012   1.00000        .        .        .        .        .        .
    ##   BOH_8_2012         .        .        .        .        .        .        .
    ##   BOH_9_2012   0.33333        .        .        .        .        .        .
    ##   BOH_6_2011         .        .        .        .        .        .        .
    ##   BOH_10_2012        .        .        .        .        .        .        .
    ##   BOH_11_2012  0.33333        .        .  0.16667        .  0.16667  0.16667
    ##   BOH_8_2011         .        .        .        .        .        .        .
    ##   BOH_9_2013   0.40000        .        .        .        .        .  0.20000
    ##   BOH_10_2013  0.50000        .        .        .        .        .        .
    ##   BOH_9_2011         .        .        .        .        .        .        .
    ##   BOH_11_2013        .        .        .        .        .        .        .
    ##   BOH_7_2014   0.50000        .        .        .        .        .        .
    ##   BOH_8_2014         .        .        .        .        .        .        .
    ##   BOH_9_2014         .        .        .        .        .        .        .
    ##   BOH_11_2014        .        .        .        .        .        .        .
    ##   BOH_9_2015   1.00000        .        .        .        .        .        .
    ##   BOH_10_2015  1.00000        .        .        .        .        .  1.00000
    ##   BUR_7_2013         .        .        .        .        .        .        .
    ##   BUR_9_2013         .        .        .        .        .        .        .
    ##   BUR_8_2016   0.83333        .        .        .        .        .  0.16667
    ##   BUR_9_2016   0.21053        .        .        .        .        .        .
    ##   BUR_10_2016  0.61538        .        .        .        .        .  0.07692
    ##   BUR_12_2016  1.00000        .        .        .        .        .        .
    ##   BUR_6_2018         .        .        .        .        .        .        .
    ##   BUR_8_2018   0.12500        .        .        .        .        .        .
    ##   BUR_9_2018   0.09091        .        .        .        .        .        .
    ##   BUR_10_2013  0.11364        .        .        .        .        .  0.04545
    ##   BUR_10_2018  0.05263        .        .        .        .        .        .
    ##   BUR_11_2018  0.50000        .        .        .        .        .        .
    ##   BUR_7_2019   1.00000        .        .        .        .        .        .
    ##   BUR_8_2019         .        .        .        .        .        .        .
    ##   BUR_9_2019   0.85714        .        .  0.14286        .        .  0.14286
    ##   BUR_10_2019  0.75000        .        .        .        .        .        .
    ##   BUR_11_2019  1.00000        .        .        .        .        .        .
    ##   BUR_12_2019  0.80000        .        .        .        .        .        .
    ##   BUR_1_2020   1.00000        .        .        .        .        .        .
    ##   BUR_5_2020         .        .        .        .        .        .        .
    ##   BUR_6_2020         .        .        .        .        .        .        .
    ##   BUR_7_2020         .        .        .        .        .        .        .
    ##   BUR_8_2020         .        .        .        .        .        .        .
    ##   BUR_10_2020        .        .        .        .        .        .        .
    ##   BUR_11_2020        .        .        .        .        .        .        .
    ##   BUR_12_2020        .        .        .        .        .        .        .
    ##   BUR_8_2021         .        .        .        .        .        .        .
    ##   BUR_9_2021         .        .        .        .        .        .        .
    ##   BUR_10_2021        .        .        .        .        .        .        .
    ##   BUR_8_2013         .        .        .        .        .        .        .
    ##   BUR_11_2021        .        .        .        .        .        .        .
    ##   BUR_12_2021        .        .        .        .        .        .        .
    ##   BUR_11_2013  0.11111        .        .        .        .        .        .
    ##   BUR_12_2013        .        .        .        .        .        .        .
    ##   BUR_8_2014         .        .        .        .        .        .        .
    ##   BUR_9_2014         .        .        .        .        .        .        .
    ##   BUR_10_2014  0.33333        .        .        .        .        .        .
    ##   BUR_11_2014        .        .        .        .        .        .        .
    ##   BUR_12_2014  0.50000        .        .        .        .        .        .
    ##   BUR_9_2015   0.20000        .        .        .        .        .        .
    ##   BUR_11_2015        .        .        .        .        .        .        .
    ##   BUR_7_2016   1.00000        .        .        .        .        .        .
    ##   FRB_10_2011        .        .        .        .        .        .        .
    ##   FRB_2_2013         .        .        .  0.50000        .        .        .
    ##   FRB_10_2013        .        .        .        .        .        .        .
    ##   FRB_12_2013        .        .        .        .        .        .        .
    ##   FRB_1_2014         .        .        .  1.00000        .        .        .
    ##   FRB_2_2014         .        .        .        .        .        .        .
    ##   FRB_5_2014   1.00000        .        .        .        .        .        .
    ##   FRB_8_2014         .        .        .        .        .        .        .
    ##   FRB_12_2014  1.00000        .        .        .        .        .        .
    ##   FRB_5_2015         .        .        .  1.00000        .        .        .
    ##   FRB_12_2011        .        .        .        .        .        .        .
    ##   FRB_10_2012        .        .        .        .        .        .        .
    ##   FRB_11_2012        .        .        .  0.33333        .        .        .
    ##   FRE_1_2011         .        .        .        .        .        .        .
    ##   FRE_6_2011         .        .        .  0.50000        .        .        .
    ##   FRE_7_2016   0.12903        .        .  0.12903        .        .  0.03226
    ##   FRE_7_2011         .        .        .        .        .        .        .
    ##   FRE_10_2016  0.32258        .        .  0.22581        .        .  0.03226
    ##   FRE_11_2016        .        .        .        .        .        .        .
    ##   FRE_7_2017   0.04082        .        .  0.06122        .        .  0.02041
    ##   FRE_8_2011         .        .        .        .        .        .        .
    ##   FRE_10_2011        .        .        .        .        .        .        .
    ##   FRE_3_2011         .        .        .        .        .        .        .
    ##   FRE_9_2017   0.13953        .        .  0.06977        .        .        .
    ##   FRE_11_2011  0.33333        .        .        .        .        .        .
    ##   FRE_10_2017  0.34783        .        .  0.21739        .        .        .
    ##   FRE_12_2011        .        .        .  1.00000        .        .        .
    ##   FRE_1_2012         .        .        .        .        .        .        .
    ##   FRE_3_2012         .        .        .        .        .        .        .
    ##   FRE_6_2012         .        .        .        .        .        .        .
    ##   FRE_7_2012         .        .        .        .        .        .        .
    ##   FRE_8_2012   0.25000        .        .        .        .        .        .
    ##   FRE_9_2012         .        .        .        .        .        .        .
    ##   FRE_10_2012        .        .        .        .        .        .        .
    ##   FRE_11_2012        .        .        .        .        .        .        .
    ##   FRE_12_2012        .        .        .        .        .        .        .
    ##   FRE_1_2013         .        .        .        .        .        .        .
    ##   FRE_4_2011         .        .        .  0.50000        .        .        .
    ##   FRE_6_2013         .        .        .  1.00000        .        .        .
    ##   FRE_7_2013         .        .        .        .        .        .        .
    ##   FRE_8_2013         .        .        .        .        .        .        .
    ##   FRE_9_2013         .        .        .        .        .        .        .
    ##   FRE_10_2013        .        .        .  0.33333        .        .        .
    ##   FRE_11_2013        .        .        .        .        .        .        .
    ##   FRE_12_2013        .        .        .        .        .        .        .
    ##   FRE_2_2014         .        .        .  0.50000        .        .        .
    ##   FRE_4_2014   0.50000        .        .        .        .        .        .
    ##   FRE_5_2014   0.33333        .        .        .        .        .        .
    ##   FRE_7_2014         .        .        .        .        .        .        .
    ##   FRE_8_2014   0.50000        .        .  0.50000        .        .  0.50000
    ##   FRE_5_2011         .        .        .        .        .        .        .
    ##   FRE_10_2014        .        .        .        .        .        .        .
    ##   FRE_11_2014        .        .        .        .        .        .        .
    ##   FRE_12_2014        .        .        .        .        .        .        .
    ##   FRE_3_2015         .        .        .        .        .        .        .
    ##   FRE_8_2015         .        .        .        .        .        .        .
    ##   FRE_10_2015  1.00000        .        .        .        .        .        .
    ##   FRE_11_2015        .        .        .        .        .        .        .
    ##   FRE_12_2015        .        .        .        .        .        .        .
    ##   FRE_9_2015         .        .        .        .        .        .        .
    ##   FRI_10_2013  0.08333        .        .        .        .        .        .
    ##   FRI_10_2017        .        .        .        .        .        .        .
    ##   FRI_11_2017        .        .        .        .        .        .        .
    ##   FRI_12_2017        .        .        .        .        .        .        .
    ##   FRI_1_2018         .        .        .        .        .        .        .
    ##   FRI_11_2013        .        .        .        .        .        .        .
    ##   FRI_12_2013        .        .        .        .        .        .        .
    ##   FRI_3_2014         .        .        .        .        .        .        .
    ##   FRI_9_2014         .        .        .        .        .        .        .
    ##   FRI_1_2015         .        .        .  0.50000        .        .        .
    ##   FRI_9_2016         .        .        .        .        .        .        .
    ##   FRI_11_2016  0.11111        .        .        .        .        .  0.11111
    ##   FRI_8_2017         .        .        .        .        .        .        .
    ##   FRI_9_2017         .        .        .        .        .        .        .
    ##   GEN_NA_NA          .        .        .        .        .        .        .
    ##   GEN_8_2020         .        .        .        .        .        .        .
    ##   GEN_12_2021        .        .        .        .        .        .        .
    ##   GEN_1_2022         .        .        .        .        .        .        .
    ##   KON_9_2013         .        .        .        .        .        .        .
    ##   KON_10_2014        .        .        .        .        .        .        .
    ##   KON_12_2016  0.46154        .        .        .        .        .        .
    ##   KON_2_2017   0.80000        .        .        .        .        .        .
    ##   KON_7_2017   1.00000        .  0.25000        .        .        .        .
    ##   KON_9_2017   0.20000        .        .        .        .        .        .
    ##   KON_10_2017  0.40000        .        .        .        .        .        .
    ##   KON_11_2017        .        .        .        .        .        .        .
    ##   KON_12_2017        .        .        .        .        .        .        .
    ##   KON_11_2014        .        .        .        .        .        .        .
    ##   KON_1_2018         .        .        .        .        .        .        .
    ##   KON_7_2018   1.00000        .        .        .        .        .        .
    ##   KON_8_2018         .        .        .        .        .        .        .
    ##   KON_9_2018   0.07692        .        .        .        .        .        .
    ##   KON_10_2018  0.10526        .        .        .        .        .        .
    ##   KON_12_2014        .        .        .        .        .        .        .
    ##   KON_12_2018        .        .        .        .        .        .        .
    ##   KON_7_2019   0.50000        .        .        .        .        .        .
    ##   KON_8_2019   0.66667        .        .        .        .        .        .
    ##   KON_9_2019         .        .        .        .        .        .        .
    ##   KON_10_2019        .        .        .        .        .        .        .
    ##   KON_11_2013        .        .        .        .        .        .        .
    ##   KON_1_2015         .        .        .        .        .        .        .
    ##   KON_11_2019  0.11111        .        .        .        .        .        .
    ##   KON_3_2016   1.00000        .        .        .        .        .        .
    ##   KON_1_2020         .        .        .        .        .        .        .
    ##   KON_7_2020         .        .        .        .        .        .        .
    ##   KON_8_2020         .        .        .        .        .        .        .
    ##   KON_9_2020         .        .        .        .        .        .        .
    ##   KON_10_2020        .        .        .        .        .        .        .
    ##   KON_8_2016   0.50000        .        .        .        .        .        .
    ##   KON_11_2020        .        .        .        .        .        .        .
    ##   KON_12_2020        .        .        .        .        .        .        .
    ##   KON_8_2021         .        .        .        .        .        .        .
    ##   KON_10_2021  0.12500        .        .        .        .        .        .
    ##   KON_11_2021        .        .        .        .        .        .        .
    ##   KON_12_2021        .        .        .        .        .        .        .
    ##   KON_9_2016   0.83333        .        .        .        .        .        .
    ##   KON_10_2016  0.47826        .        .        .  0.04348        .  0.04348
    ##   KON_12_2013        .        .        .        .        .        .        .
    ##   KON_2_2014   1.00000        .        .        .        .        .        .
    ##   KON_9_2014         .        .        .        .        .        .        .
    ##   KON_11_2016  0.55556        .        .        .        .        .  0.11111
    ##   NEU_8_2013         .        .        .        .        .        .        .
    ##   NEU_10_2014        .        .        .        .        .        .        .
    ##   NEU_11_2014        .        .        .        .        .        .        .
    ##   NEU_10_2015        .        .        .        .        .        .        .
    ##   NEU_8_2016   0.50000        .        .        .        .        .        .
    ##   NEU_9_2016         .        .        .        .        .        .        .
    ##   NEU_9_2013         .        .        .        .        .        .        .
    ##   NEU_10_2016  0.20000        .        .        .        .        .        .
    ##   NEU_10_2013        .        .        .        .        .        .        .
    ##   NEU_2_2017         .        .        .        .        .        .        .
    ##   NEU_11_2013        .        .        .        .        .        .        .
    ##   NEU_7_2020         .        .        .        .        .        .        .
    ##   NEU_8_2020         .        .        .        .        .        .        .
    ##   NEU_9_2020         .        .        .        .        .        .        .
    ##   NEU_10_2021        .        .        .        .        .        .        .
    ##   NEU_11_2021        .        .        .        .        .        .        .
    ##   NEU_12_2013  0.50000        .        .  0.50000        .        .        .
    ##   RIE_3_2011         .        .        .        .        .        .        .
    ##   RIE_8_2011         .        .        .        .        .        .        .
    ##   RIE_10_2011        .        .        .        .        .        .        .
    ##   RIE_11_2011        .        .        .        .        .        .        .
    ##   RIE_12_2011        .        .        .        .        .        .        .
    ##   RIE_2_2012         .        .        .        .        .        .        .
    ##   RIE_3_2012         .        .        .        .        .        .        .
    ##   RIE_10_2012        .        .        .        .        .        .        .
    ##   RIE_11_2012        .        .        .        .        .        .        .
    ##   RIE_12_2012        .        .        .        .        .        .        .
    ##   RIE_1_2013         .        .        .        .        .        .        .
    ##   RIE_2_2013         .        .        .        .        .        .        .
    ##   RIE_3_2013         .        .        .        .        .        .        .
    ##   RIE_10_2013        .        .        .        .        .        .        .
    ##   RIE_11_2013        .        .        .        .        .        .        .
    ##   RIE_12_2013  0.25000        .        .        .        .        .        .
    ##   RIE_1_2014   0.33333        .        .        .        .        .        .
    ##   RIE_4_2011         .        .        .        .        .        .        .
    ##   RIE_2_2014         .        .        .        .        .        .        .
    ##   RIE_9_2014         .        .        .        .        .        .        .
    ##   RIE_10_2014        .        .        .        .        .        .        .
    ##   RIE_11_2014        .        .        .        .        .        .        .
    ##   RIE_12_2014        .        .        .        .        .        .        .
    ##   RIE_1_2015         .        .        .        .        .        .        .
    ##   RIE_2_2015         .        .        .        .        .        .        .
    ##   RIE_10_2015  1.00000        .        .        .        .        .        .
    ##   RIE_12_2015        .        .        .        .        .        .        .
    ##   RIE_6_2011         .        .        .        .        .        .        .
    ##   RIE_7_2011         .        .        .  1.00000        .        .        .
    ##   SCD_7_2012         .        .        .        .        .        .        .
    ##   SCD_9_2012   0.06667        .        .        .        .        .  0.06667
    ##   SCD_10_2016        .        .        .        .        .        .        .
    ##   SCD_11_2016        .        .        .        .        .        .        .
    ##   SCD_12_2016        .        .        .        .        .        .        .
    ##   SCD_10_2017        .        .        .        .        .        .        .
    ##   SCD_10_2012        .        .        .        .        .        .        .
    ##   SCD_11_2012  0.14286        .        .        .        .        .        .
    ##   SCD_12_2012        .        .        .        .        .        .        .
    ##   SCD_8_2013   0.66667        .        .        .        .        .        .
    ##   SCD_9_2013   0.50000        .        .        .        .        .        .
    ##   SCD_10_2013  0.33333        .        .        .        .        .        .
    ##   SCD_11_2013        .        .        .        .        .        .        .
    ##   SCD_1_2014         .        .        .        .        .        .        .
    ##   SCD_3_2014         .        .        .        .        .        .        .
    ##   SCD_7_2014         .        .        .        .        .        .        .
    ##   SCD_9_2014         .        .        .        .        .        .        .
    ##   SCD_11_2014  1.00000        .        .  0.50000        .        .  0.50000
    ##   SCD_7_2015         .        .        .        .        .        .        .
    ##   SCD_8_2015   0.25000        .        .        .        .        .        .
    ##   SCD_9_2015   0.16667        .        .  0.33333        .        .        .
    ##   SCD_10_2015  0.75000        .        .        .        .        .  0.25000
    ##   SCD_11_2015        .        .        .        .        .        .        .
    ##   SCD_2_2016   1.00000        .        .        .        .        .        .
    ##   SCD_8_2016   0.66667        .        .        .        .        .        .
    ##   SCD_9_2016         .        .        .        .        .        .        .
    ##   SCG_3_2011         .        .        .        .        .        .        .
    ##   SCG_2_2012         .        .        .        .        .        .        .
    ##   SCG_3_2012   0.50000        .        .  0.50000        .        .        .
    ##   SCG_7_2011         .        .        .        .        .        .        .
    ##   SCG_11_2011        .        .        .        .        .        .        .
    ##   SCG_10_2012        .        .        .        .        .        .        .
    ##   SCG_11_2012  0.66667        .        .        .        .        .  0.33333
    ##   SCG_1_2013   0.50000        .        .        .        .        .        .
    ##   SCG_1_2014         .        .        .        .        .        .        .
    ##   SCG_6_2014         .        .        .        .        .        .        .
    ##   SCG_11_2014  1.00000        .        .        .        .        .        .
    ##   SCL_7_2012   0.04762        .        .  0.01190        .        .        .
    ##   SCL_9_2014         .        .        .        .        .        .        .
    ##   SCL_10_2014  0.50000        .        .        .        .        .        .
    ##   SCL_6_2015   0.50000        .        .        .        .        .        .
    ##   SCL_7_2015         .        .        .        .        .        .        .
    ##   SCL_9_2012         .        .        .        .        .        .        .
    ##   SCL_10_2012  0.14286        .        .        .        .        .        .
    ##   SCL_8_2016   0.53333        .        .        .        .        .        .
    ##   SCL_9_2016   0.20000        .        .        .        .        .        .
    ##   SCL_10_2016  0.16129        .        .  0.03226  0.03226        .        .
    ##   SCL_11_2012        .        .        .        .        .        .        .
    ##   SCL_4_2012         .        .        .        .        .        .        .
    ##   SCL_11_2016        .        .        .        .        .        .        .
    ##   SCL_12_2016        .        .        .        .        .        .        .
    ##   SCL_2_2017         .        .        .        .        .        .        .
    ##   SCL_6_2017   0.33333        .        .        .        .        .        .
    ##   SCL_12_2012        .        .        .        .        .        .        .
    ##   SCL_7_2017         .        .        .        .        .        .        .
    ##   SCL_1_2013         .        .        .        .        .        .        .
    ##   SCL_8_2017         .        .        .        .        .        .        .
    ##   SCL_7_2013   0.12500        .        .  0.12500        .        .        .
    ##   SCL_9_2017         .        .        .        .        .        .        .
    ##   SCL_10_2017        .        .        .        .        .        .        .
    ##   SCL_11_2017        .        .        .        .        .        .        .
    ##   SCL_12_2017        .        .        .        .        .        .        .
    ##   SCL_6_2018         .        .        .        .        .        .        .
    ##   SCL_7_2018         .        .        .        .        .        .        .
    ##   SCL_8_2018   0.14286        .        .        .        .        .        .
    ##   SCL_8_2013         .        .        .        .        .        .  0.50000
    ##   SCL_10_2013        .        .        .        .        .        .        .
    ##   SCL_1_2014         .        .        .        .        .  1.00000        .
    ##   SCL_4_2014         .        .        .        .        .        .        .
    ##   SCL_6_2014         .        .        .        .        .        .        .
    ##   SCL_7_2014   0.04167        .        .        .        .        .        .
    ##   SCL_8_2014         .        .        .        .        .        .        .
    ##   SCS_8_2013         .        .        .        .        .        .        .
    ##   SCS_1_2014   0.50000        .        .        .        .        .        .
    ##   SCS_5_2014         .        .        .        .        .        .        .
    ##   SCS_7_2014         .        .        .        .        .        .        .
    ##   SCS_8_2014         .        .        .        .        .        .        .
    ##   SCS_9_2014         .        .        .        .        .        .        .
    ##   SCS_11_2014        .        .        .        .        .        .        .
    ##   SCS_12_2014        .        .        .        .        .        .  0.33333
    ##   SCS_1_2016         .        .        .        .        .        .        .
    ##   SCS_8_2016   0.16667        .        .        .        .        .        .
    ##   SCS_9_2016         .        .        .        .        .        .        .
    ##   SCS_10_2013  0.25000        .        .        .        .        .        .
    ##   SCS_10_2016  0.16667        .        .        .        .        .        .
    ##   SCS_9_2017   0.33333        .        .        .        .        .        .
    ##   SCS_NA_2016        .        .        .        .        .        .        .
    ##   SCS_11_2013        .        .        .        .        .        .        .
    ##   TRO_3_2011         .        .        .        .        .        .        .
    ##   TRO_11_2011        .        .        .        .        .        .        .
    ##   TRO_12_2011        .        .        .        .        .        .        .
    ##   TRO_1_2012         .        .        .        .        .        .        .
    ##   TRO_2_2012         .        .        .        .        .        .        .
    ##   TRO_3_2012         .        .        .        .        .        .        .
    ##   TRO_11_2012        .        .        .        .        .        .        .
    ##   TRO_12_2012        .        .        .        .        .        .        .
    ##   TRO_1_2013   0.25000        .        .        .        .        .        .
    ##   TRO_2_2013   0.50000        .        .        .        .        .  0.50000
    ##   TRO_3_2013         .        .        .        .        .        .        .
    ##   TRO_4_2011   0.50000        .        .        .        .        .        .
    ##   TRO_8_2011         .        .        .        .        .        .        .
    ##   TRO_10_2011        .        .        .        .        .        .        .
    ##   UEB_7_2012         .        .        .        .        .        .        .
    ##   UEB_8_2015   0.42857        .        .        .        .        .        .
    ##   UEB_9_2015   0.15385        .        .        .        .        .  0.15385
    ##   UEB_3_2013         .        .        .        .        .        .        .
    ##   UEB_10_2015        .        .        .        .        .        .        .
    ##   UEB_10_2011        .        .        .        .        .        .        .
    ##   UEB_4_2013         .        .        .        .        .        .        .
    ##   UEB_5_2013         .        .        .        .        .        .        .
    ##   UEB_10_2013        .        .        .        .        .        .        .
    ##   UEB_8_2013         .        .        .        .        .        .        .
    ##   UEB_9_2013         .        .        .        .        .        .        .
    ##   UEB_7_2016   1.00000        .        .        .  1.00000        .        .
    ##   UEB_8_2016   1.00000        .        .        .  1.00000        .        .
    ##   UEB_9_2016   1.00000        .        .        .  1.00000        .        .
    ##   UEB_10_2016  1.00000        .        .        .  1.00000        .        .
    ##   UEB_11_2013  0.16667        .        .        .        .        .        .
    ##   UEB_11_2016  1.00000        .        .        .  1.00000        .        .
    ##   UEB_12_2016  1.00000        .        .        .  1.00000        .        .
    ##   UEB_6_2017   1.00000        .        .        .  1.00000        .        .
    ##   UEB_8_2012         .        .        .  0.50000        .        .        .
    ##   UEB_7_2017   1.00000        .        .        .  1.00000        .        .
    ##   UEB_1_2018   1.00000        .        .        .  1.00000        .        .
    ##   UEB_9_2019         .        .        .        .        .        .        .
    ##   UEB_12_2013  0.50000        .        .        .        .        .        .
    ##   UEB_1_2014         .        .        .        .        .        .        .
    ##   UEB_2_2014         .        .        .        .        .        .        .
    ##   UEB_3_2014         .        .        .        .        .        .        .
    ##   UEB_5_2014         .        .        .  0.33333        .        .  0.33333
    ##   UEB_6_2014   0.50000        .        .        .        .        .        .
    ##   UEB_9_2012         .        .        .        .        .        .        .
    ##   UEB_7_2014   0.18182        .        .  0.18182        .        .  0.18182
    ##   UEB_8_2014         .        .        .        .        .        .        .
    ##   UEB_9_2014   0.14286        .        .        .        .        .        .
    ##   UEB_10_2014  0.33333        .        .        .        .        .        .
    ##   UEB_11_2014  0.16667        .        .        .        .        .        .
    ##   UEB_12_2014        .        .        .        .        .        .        .
    ##   UEB_1_2015   0.20000        .        .        .        .        .        .
    ##   UEB_2_2013         .        .        .        .        .        .        .
    ##   UEB_2_2015         .        .        .        .        .        .        .
    ##   UEB_3_2015   0.12500        .        .        .        .        .        .
    ##   UEB_5_2015         .        .        .        .        .        .        .
    ##   UEB_6_2015         .        .        .  1.00000        .        .        .
    ##   UEB_7_2015   1.00000        .        .        .        .        .        .
    ##   UST_9_2013   0.12500        .        .        .        .        .        .
    ##   UST_10_2013  0.20000        .        .        .        .        .        .
    ##   UST_11_2013  0.16667        .        .  0.08333        .        .        .
    ##   UST_12_2013        .        .        .        .        .        .        .
    ##   UST_1_2014         .        .        .        .        .        .  0.33333
    ##   UST_3_2014         .        .        .        .        .        .        .
    ##   UST_8_2014         .        .        .        .        .        .        .
    ##   UST_9_2014         .        .        .        .        .        .        .
    ##   WSL_9_2012         .        .        .        .        .        .        .
    ##   WSL_8_2016         .        .        .        .        .        .        .
    ##   WSL_9_2016         .        .        .        .        .        .        .
    ##   WSL_10_2016        .        .        .        .        .        .        .
    ##   WSL_11_2016        .        .        .        .        .        .        .
    ##   WSL_12_2016        .        .        .        .        .        .        .
    ##   WSL_3_2017         .        .        .        .        .        .        .
    ##   WSL_7_2017         .        .        .        .        .        .        .
    ##   WSL_8_2017         .        .        .        .        .        .        .
    ##   WSL_9_2017         .        .        .        .        .        .        .
    ##   WSL_10_2017        .        .        .        .        .        .        .
    ##   WSL_11_2017        .        .        .        .        .        .        .
    ##   WSL_1_2018         .        .        .        .        .        .        .
    ##   WSL_3_2018         .        .        .        .        .        .        .
    ##   WSL_5_2018         .        .        .        .        .        .        .
    ##   WSL_6_2018         .        .        .        .        .        .        .
    ##   WSL_7_2018         .        .        .        .        .        .        .
    ##   WSL_8_2018         .        .        .        .        .        .        .
    ##   WSL_9_2018         .        .        .        .        .        .        .
    ##   WSL_6_2019         .        .        .        .        .        .        .
    ##   WSL_7_2019         .        .        .        .        .        .        .
    ##   WSL_8_2019         .        .        .        .        .        .        .
    ##   WSL_9_2019         .        .        .        .        .        .        .
    ##   WSL_10_2019        .        .        .        .        .        .        .
    ##   WSL_4_2020         .        .        .        .        .        .        .
    ##   WSL_5_2020         .        .        .        .        .        .        .
    ##   WSL_6_2020         .        .        .        .        .        .        .
    ##   WSL_7_2020         .        .        .        .        .        .        .
    ##   WSL_8_2020         .        .        .        .        .        .        .
    ##   WSL_9_2020         .        .        .        .        .        .        .
    ##   WSL_10_2020        .        .        .        .        .        .        .
    ##   WSL_11_2020        .        .        .        .        .        .        .
    ##   WSL_6_2021         .        .        .        .        .        .        .
    ##   WSL_7_2021         .        .        .        .        .        .        .
    ##   Total        0.11374  0.00037  0.00037  0.01957  0.00997  0.00111  0.01256
    ##              Locus
    ## Population    aest25_1 aest26_1 aest28_1 aest29_1 aest31_1 aest35_1 aest36_1
    ##   ALD_3_2011         .        .        .        .        .        .        .
    ##   ALD_6_2011         .        .        .        .        .        .        .
    ##   ALD_7_2011         .        .        .        .        .        .        .
    ##   ALD_8_2011         .        .        .        .        .        .        .
    ##   ALD_10_2011        .        .        .        .        .        .        .
    ##   ALD_11_2011        .        .        .        .        .        .        .
    ##   ALD_12_2011        .        .        .        .        .        .        .
    ##   ALD_1_2012         .        .        .        .        .        .        .
    ##   ALD_2_2012         .        .        .        .        .        .        .
    ##   ALD_3_2012         .        .        .        .        .        .        .
    ##   ALD_10_2012        .        .        .        .        .        .        .
    ##   ALD_11_2012        .        .        .        .        .        .        .
    ##   ALD_12_2012        .        .        .        .        .        .        .
    ##   ALD_1_2013         .        .        .        .        .        .        .
    ##   ALD_2_2013         .        .        .        .        .        .        .
    ##   ALD_3_2013         .        .        .        .        .        .        .
    ##   ALD_10_2013        .        .        .        .        .        .        .
    ##   ALD_4_2011         .        .        .        .        .        .        .
    ##   ALD_11_2013        .        .        .        .        .        .        .
    ##   ALD_12_2013        .        .        .        .  0.25000        .        .
    ##   ALD_1_2014         .        .        .        .  0.50000        .        .
    ##   ALD_2_2014         .        .        .        .        .        .        .
    ##   ALD_5_2011         .        .        .        .        .        .        .
    ##   ALD_7_2014         .        .        .        .  1.00000        .        .
    ##   ALD_8_2014         .        .        .        .        .        .        .
    ##   ALD_9_2014         .        .        .        .        .        .        .
    ##   ALD_10_2014        .        .        .        .        .        .        .
    ##   ALD_1_2015         .        .        .        .        .        .        .
    ##   ALD_12_2015        .        .        .        .  0.50000        .        .
    ##   BAR_5_2014         .        .        .        .        .        .        .
    ##   BAR_6_2014         .        .        .        .        .        .        .
    ##   BOB_8_2011         .        .        .        .        .        .        .
    ##   BOB_2_2012         .        .        .        .        .        .        .
    ##   BOB_1_2015         .        .  0.33333        .  0.33333        .        .
    ##   BOB_3_2015         .        .        .        .        .        .        .
    ##   BOB_6_2015         .        .        .        .        .        .        .
    ##   BOB_8_2015         .        .  0.07692        .  0.15385        .        .
    ##   BOB_3_2012         .        .        .        .        .        .        .
    ##   BOB_10_2015        .        .  0.50000        .  0.50000        .        .
    ##   BOB_1_2016         .        .        .        .        .        .        .
    ##   BOB_2_2016         .        .  1.00000        .        .        .        .
    ##   BOB_11_2012        .        .        .        .        .        .        .
    ##   BOB_9_2014         .        .  0.12500        .  0.12500        .        .
    ##   BOB_10_2013        .        .  0.02564        .        .        .        .
    ##   BOB_5_2012         .        .        .        .        .        .        .
    ##   BOB_7_2012         .        .        .        .        .        .        .
    ##   BOB_6_2016         .        .        .        .  1.00000        .        .
    ##   BOB_7_2016         .        .        .        .        .        .        .
    ##   BOB_8_2016         .        .        .        .  0.14286        .        .
    ##   BOB_9_2016         .        .        .        .        .        .        .
    ##   BOB_10_2016        .        .        .        .        .        .        .
    ##   BOB_9_2011         .        .        .        .        .        .        .
    ##   BOB_11_2016        .        .        .        .        .        .  0.07692
    ##   BOB_8_2012         .        .        .        .        .        .        .
    ##   BOB_12_2016        .        .        .        .        .        .        .
    ##   BOB_9_2012   0.09091        .        .        .        .        .        .
    ##   BOB_6_2017         .        .        .        .        .        .        .
    ##   BOB_7_2017         .        .        .        .        .        .        .
    ##   BOB_8_2017         .        .        .        .        .        .        .
    ##   BOB_9_2017         .        .        .        .        .        .        .
    ##   BOB_10_2011        .        .        .        .        .        .        .
    ##   BOB_10_2017        .        .        .        .        .        .        .
    ##   BOB_11_2017        .        .        .        .        .        .        .
    ##   BOB_10_2012        .        .        .        .        .        .        .
    ##   BOB_11_2011        .        .        .        .        .        .        .
    ##   BOB_12_2017        .        .        .        .        .        .        .
    ##   BOB_1_2018         .        .        .        .        .        .        .
    ##   BOB_3_2018         .        .        .        .        .        .        .
    ##   BOB_4_2018         .        .        .        .        .        .        .
    ##   BOB_6_2018         .        .        .        .        .        .        .
    ##   BOB_8_2018         .        .        .        .  0.05263        .        .
    ##   BOB_1_2013         .        .        .        .        .        .        .
    ##   BOB_9_2018         .        .        .        .        .        .        .
    ##   BOB_2_2013         .        .        .        .        .        .        .
    ##   BOB_5_2013         .        .  0.50000        .  0.50000        .        .
    ##   BOB_10_2018        .        .        .        .  0.23810        .        .
    ##   BOB_12_2011        .        .        .        .        .        .        .
    ##   BOB_6_2013         .        .        .        .        .        .        .
    ##   BOB_8_2013         .        .        .        .        .        .        .
    ##   BOB_9_2019         .        .        .        .        .        .        .
    ##   BOB_9_2013         .        .        .        .        .        .        .
    ##   BOB_12_2013        .        .  0.20000        .        .        .        .
    ##   BOB_1_2014         .        .        .        .        .        .        .
    ##   BOB_2_2014         .        .        .        .        .        .        .
    ##   BOB_3_2014         .        .        .        .        .        .        .
    ##   BOB_7_2014   0.12500        .        .        .  0.12500        .        .
    ##   BOB_8_2014         .        .        .        .        .        .        .
    ##   BOB_10_2014        .        .        .        .        .        .        .
    ##   BOB_11_2014        .        .  0.09091        .        .        .        .
    ##   BOB_12_2014        .        .        .        .        .        .        .
    ##   BOH_5_2011         .        .        .        .        .        .        .
    ##   BOH_10_2011        .        .        .        .        .        .        .
    ##   BOH_12_2011        .        .        .        .        .        .        .
    ##   BOH_3_2012         .        .        .        .        .        .        .
    ##   BOH_8_2012         .        .        .        .        .        .        .
    ##   BOH_9_2012         .        .        .  0.16667  0.33333        .        .
    ##   BOH_6_2011         .        .        .        .        .        .        .
    ##   BOH_10_2012        .        .        .        .        .        .        .
    ##   BOH_11_2012        .        .        .        .        .        .        .
    ##   BOH_8_2011         .        .        .        .        .        .        .
    ##   BOH_9_2013         .        .        .        .        .        .        .
    ##   BOH_10_2013        .        .        .        .        .        .        .
    ##   BOH_9_2011         .        .  0.33333        .        .        .        .
    ##   BOH_11_2013        .        .        .        .        .        .        .
    ##   BOH_7_2014         .        .  0.50000        .        .        .        .
    ##   BOH_8_2014         .        .        .        .        .        .        .
    ##   BOH_9_2014         .        .        .        .        .        .        .
    ##   BOH_11_2014        .        .        .        .        .        .        .
    ##   BOH_9_2015         .        .  1.00000        .  0.50000        .        .
    ##   BOH_10_2015        .        .        .        .  1.00000        .        .
    ##   BUR_7_2013         .        .        .        .        .        .        .
    ##   BUR_9_2013         .        .        .        .        .        .        .
    ##   BUR_8_2016         .        .        .        .  0.83333        .        .
    ##   BUR_9_2016         .        .        .        .        .        .        .
    ##   BUR_10_2016        .        .        .        .  0.23077        .        .
    ##   BUR_12_2016        .        .        .        .  1.00000        .        .
    ##   BUR_6_2018         .        .        .        .        .        .        .
    ##   BUR_8_2018         .        .        .        .  0.25000        .        .
    ##   BUR_9_2018         .        .        .        .        .  0.09091  0.09091
    ##   BUR_10_2013        .        .        .  0.02273  0.09091        .        .
    ##   BUR_10_2018        .  0.10526        .        .        .        .        .
    ##   BUR_11_2018        .        .        .        .  0.50000        .        .
    ##   BUR_7_2019         .        .        .        .        .        .        .
    ##   BUR_8_2019         .        .        .        .        .        .        .
    ##   BUR_9_2019         .        .        .        .  0.42857        .        .
    ##   BUR_10_2019        .        .        .        .  0.50000        .        .
    ##   BUR_11_2019        .        .        .        .  0.50000  0.25000        .
    ##   BUR_12_2019        .        .        .        .  0.40000        .        .
    ##   BUR_1_2020         .        .        .        .  1.00000        .        .
    ##   BUR_5_2020         .        .        .        .        .        .        .
    ##   BUR_6_2020         .        .        .        .        .        .        .
    ##   BUR_7_2020         .        .        .        .        .        .        .
    ##   BUR_8_2020         .        .  0.04348        .        .        .        .
    ##   BUR_10_2020        .        .        .        .        .        .        .
    ##   BUR_11_2020        .        .        .        .        .        .        .
    ##   BUR_12_2020        .        .        .        .        .        .        .
    ##   BUR_8_2021         .        .        .        .        .        .        .
    ##   BUR_9_2021         .  0.33333        .        .        .        .        .
    ##   BUR_10_2021        .  0.23077        .        .        .        .        .
    ##   BUR_8_2013         .        .        .        .        .        .        .
    ##   BUR_11_2021        .  0.33333        .        .        .        .        .
    ##   BUR_12_2021        .        .        .        .        .        .        .
    ##   BUR_11_2013        .        .        .  0.11111  0.11111        .        .
    ##   BUR_12_2013        .        .        .        .        .        .        .
    ##   BUR_8_2014         .        .        .        .        .        .        .
    ##   BUR_9_2014         .        .        .        .        .        .        .
    ##   BUR_10_2014        .        .        .        .        .        .        .
    ##   BUR_11_2014        .        .        .        .        .        .        .
    ##   BUR_12_2014        .        .  0.50000        .  0.50000        .        .
    ##   BUR_9_2015         .        .        .        .        .        .        .
    ##   BUR_11_2015        .        .        .        .        .        .        .
    ##   BUR_7_2016         .        .        .        .  1.00000        .        .
    ##   FRB_10_2011        .        .        .        .        .        .        .
    ##   FRB_2_2013         .        .        .        .        .        .        .
    ##   FRB_10_2013        .        .        .        .        .        .        .
    ##   FRB_12_2013        .        .        .        .        .        .        .
    ##   FRB_1_2014         .        .        .  1.00000        .        .        .
    ##   FRB_2_2014         .        .        .        .        .        .        .
    ##   FRB_5_2014         .        .  1.00000        .  1.00000        .        .
    ##   FRB_8_2014         .        .  1.00000        .        .        .        .
    ##   FRB_12_2014        .        .  1.00000        .        .        .        .
    ##   FRB_5_2015         .        .        .        .        .        .        .
    ##   FRB_12_2011        .        .        .        .        .        .        .
    ##   FRB_10_2012        .        .        .        .        .        .        .
    ##   FRB_11_2012        .        .  0.33333        .  0.33333        .        .
    ##   FRE_1_2011         .        .        .        .        .        .        .
    ##   FRE_6_2011         .        .        .        .        .        .        .
    ##   FRE_7_2016         .        .        .        .  0.16129        .        .
    ##   FRE_7_2011         .        .        .        .        .        .        .
    ##   FRE_10_2016        .        .        .        .  0.35484        .        .
    ##   FRE_11_2016        .        .        .        .        .        .        .
    ##   FRE_7_2017         .        .        .        .  0.04082        .        .
    ##   FRE_8_2011         .        .        .        .        .        .        .
    ##   FRE_10_2011        .        .        .        .        .        .        .
    ##   FRE_3_2011         .        .        .        .  0.33333        .        .
    ##   FRE_9_2017         .        .        .        .  0.09302        .        .
    ##   FRE_11_2011        .        .        .        .        .        .        .
    ##   FRE_10_2017        .        .        .        .  0.21739        .        .
    ##   FRE_12_2011        .        .        .        .        .        .        .
    ##   FRE_1_2012         .        .        .        .        .        .        .
    ##   FRE_3_2012         .        .        .        .        .        .        .
    ##   FRE_6_2012         .        .        .        .        .        .        .
    ##   FRE_7_2012         .        .        .        .        .        .        .
    ##   FRE_8_2012         .        .        .        .        .        .        .
    ##   FRE_9_2012         .        .        .        .        .        .        .
    ##   FRE_10_2012        .        .        .        .        .        .        .
    ##   FRE_11_2012        .        .        .        .        .        .        .
    ##   FRE_12_2012        .        .        .        .        .        .        .
    ##   FRE_1_2013         .        .        .        .        .        .        .
    ##   FRE_4_2011         .        .        .        .        .        .        .
    ##   FRE_6_2013         .        .        .  1.00000        .        .        .
    ##   FRE_7_2013         .        .        .        .        .        .        .
    ##   FRE_8_2013         .        .        .        .        .        .        .
    ##   FRE_9_2013         .        .        .        .        .        .        .
    ##   FRE_10_2013        .        .        .        .        .        .        .
    ##   FRE_11_2013        .        .        .        .        .        .        .
    ##   FRE_12_2013        .        .        .        .        .        .        .
    ##   FRE_2_2014         .        .        .        .        .        .        .
    ##   FRE_4_2014         .        .        .        .        .        .        .
    ##   FRE_5_2014         .        .  0.33333        .        .        .        .
    ##   FRE_7_2014         .        .        .        .        .        .        .
    ##   FRE_8_2014         .        .  1.00000        .        .        .        .
    ##   FRE_5_2011         .        .        .        .        .        .        .
    ##   FRE_10_2014        .        .        .        .        .        .        .
    ##   FRE_11_2014        .        .        .        .  1.00000        .        .
    ##   FRE_12_2014        .        .        .        .        .        .        .
    ##   FRE_3_2015         .        .  0.50000        .  0.50000        .        .
    ##   FRE_8_2015         .        .        .        .        .        .        .
    ##   FRE_10_2015        .        .  1.00000        .        .        .        .
    ##   FRE_11_2015        .        .        .        .        .        .        .
    ##   FRE_12_2015        .        .        .        .        .        .        .
    ##   FRE_9_2015         .        .        .        .        .        .        .
    ##   FRI_10_2013        .        .        .        .  0.02778        .        .
    ##   FRI_10_2017        .        .        .        .        .        .        .
    ##   FRI_11_2017        .        .        .        .        .        .        .
    ##   FRI_12_2017        .        .        .        .        .        .        .
    ##   FRI_1_2018         .        .        .        .        .        .        .
    ##   FRI_11_2013        .        .        .        .        .        .        .
    ##   FRI_12_2013        .        .        .        .        .        .        .
    ##   FRI_3_2014         .        .        .        .        .        .        .
    ##   FRI_9_2014         .        .        .        .        .        .        .
    ##   FRI_1_2015         .        .        .        .        .        .        .
    ##   FRI_9_2016         .        .        .        .        .        .        .
    ##   FRI_11_2016        .        .  0.11111        .        .        .        .
    ##   FRI_8_2017         .        .        .        .        .        .        .
    ##   FRI_9_2017         .        .        .        .        .        .        .
    ##   GEN_NA_NA          .        .        .        .        .        .        .
    ##   GEN_8_2020         .        .        .        .        .        .        .
    ##   GEN_12_2021        .        .        .        .        .        .        .
    ##   GEN_1_2022         .        .        .        .        .        .        .
    ##   KON_9_2013         .        .        .        .        .        .        .
    ##   KON_10_2014        .        .        .        .        .        .        .
    ##   KON_12_2016        .        .        .        .  0.30769        .        .
    ##   KON_2_2017         .        .        .        .  0.40000        .        .
    ##   KON_7_2017         .        .  0.25000        .  0.75000        .  0.25000
    ##   KON_9_2017         .        .        .        .  0.20000        .        .
    ##   KON_10_2017        .        .        .        .        .        .        .
    ##   KON_11_2017        .        .        .        .        .        .        .
    ##   KON_12_2017        .        .        .        .        .        .        .
    ##   KON_11_2014        .        .        .        .        .        .        .
    ##   KON_1_2018         .        .        .        .        .        .        .
    ##   KON_7_2018         .        .        .        .  1.00000        .        .
    ##   KON_8_2018         .        .        .        .        .        .        .
    ##   KON_9_2018         .        .        .        .        .        .        .
    ##   KON_10_2018        .        .        .        .        .        .        .
    ##   KON_12_2014        .        .        .        .        .        .        .
    ##   KON_12_2018        .        .        .        .        .        .        .
    ##   KON_7_2019         .        .        .        .  0.50000        .        .
    ##   KON_8_2019         .        .        .        .  0.66667        .        .
    ##   KON_9_2019         .        .        .        .  0.33333        .        .
    ##   KON_10_2019        .        .        .        .        .        .        .
    ##   KON_11_2013        .        .        .        .  0.33333        .        .
    ##   KON_1_2015         .        .        .        .        .        .        .
    ##   KON_11_2019        .        .        .        .        .        .        .
    ##   KON_3_2016         .        .        .        .  1.00000        .        .
    ##   KON_1_2020         .        .        .        .        .        .        .
    ##   KON_7_2020         .        .        .        .        .        .        .
    ##   KON_8_2020         .        .        .        .        .        .        .
    ##   KON_9_2020         .        .        .        .        .        .        .
    ##   KON_10_2020        .        .        .        .        .        .        .
    ##   KON_8_2016         .        .        .        .  0.50000        .        .
    ##   KON_11_2020        .        .        .        .        .        .        .
    ##   KON_12_2020        .        .        .        .        .        .        .
    ##   KON_8_2021         .        .        .        .        .        .        .
    ##   KON_10_2021        .        .        .        .  0.12500        .        .
    ##   KON_11_2021        .        .        .        .        .        .        .
    ##   KON_12_2021        .        .        .        .        .        .        .
    ##   KON_9_2016         .        .        .        .  0.50000        .        .
    ##   KON_10_2016        .        .        .        .  0.30435        .  0.04348
    ##   KON_12_2013        .        .  1.00000        .  1.00000        .        .
    ##   KON_2_2014         .        .        .        .        .        .        .
    ##   KON_9_2014         .        .  0.50000  0.50000        .        .        .
    ##   KON_11_2016        .        .        .        .  0.33333        .        .
    ##   NEU_8_2013         .        .        .        .        .        .        .
    ##   NEU_10_2014        .        .  0.33333        .        .        .        .
    ##   NEU_11_2014        .        .        .        .        .        .        .
    ##   NEU_10_2015        .        .        .        .        .        .        .
    ##   NEU_8_2016         .        .        .        .        .        .        .
    ##   NEU_9_2016         .        .        .        .        .        .        .
    ##   NEU_9_2013         .        .        .        .        .        .        .
    ##   NEU_10_2016        .        .        .        .        .        .        .
    ##   NEU_10_2013  1.00000        .  1.00000        .  1.00000        .        .
    ##   NEU_2_2017         .        .        .        .        .        .        .
    ##   NEU_11_2013        .        .        .        .        .        .        .
    ##   NEU_7_2020         .        .        .        .        .        .        .
    ##   NEU_8_2020         .        .        .        .        .        .        .
    ##   NEU_9_2020         .        .        .        .        .        .        .
    ##   NEU_10_2021        .        .        .        .        .        .        .
    ##   NEU_11_2021        .        .        .        .        .        .        .
    ##   NEU_12_2013        .        .        .        .        .        .        .
    ##   RIE_3_2011         .        .        .        .        .        .        .
    ##   RIE_8_2011         .        .        .        .        .        .        .
    ##   RIE_10_2011        .        .        .        .        .        .        .
    ##   RIE_11_2011        .        .        .        .        .        .        .
    ##   RIE_12_2011        .        .        .        .        .        .        .
    ##   RIE_2_2012         .        .        .        .        .        .        .
    ##   RIE_3_2012         .        .        .        .        .        .        .
    ##   RIE_10_2012        .        .        .        .        .        .        .
    ##   RIE_11_2012        .        .        .        .        .        .        .
    ##   RIE_12_2012        .        .        .        .        .        .        .
    ##   RIE_1_2013         .        .        .        .        .        .        .
    ##   RIE_2_2013         .        .        .        .        .        .        .
    ##   RIE_3_2013         .        .        .        .        .        .        .
    ##   RIE_10_2013        .        .        .        .        .        .        .
    ##   RIE_11_2013        .        .        .        .        .        .        .
    ##   RIE_12_2013        .        .        .        .        .        .        .
    ##   RIE_1_2014         .        .        .        .        .        .        .
    ##   RIE_4_2011         .        .        .        .        .        .        .
    ##   RIE_2_2014         .        .        .        .        .        .        .
    ##   RIE_9_2014         .        .        .        .        .        .        .
    ##   RIE_10_2014        .        .        .        .        .        .        .
    ##   RIE_11_2014        .        .        .        .        .        .        .
    ##   RIE_12_2014        .        .        .        .        .        .        .
    ##   RIE_1_2015         .        .        .        .  0.33333        .        .
    ##   RIE_2_2015         .        .        .        .        .        .        .
    ##   RIE_10_2015        .        .        .        .  1.00000        .        .
    ##   RIE_12_2015        .        .        .        .        .        .        .
    ##   RIE_6_2011         .        .        .        .        .        .        .
    ##   RIE_7_2011         .        .        .        .        .        .        .
    ##   SCD_7_2012         .        .        .        .        .        .        .
    ##   SCD_9_2012         .        .        .  0.06667  0.06667        .        .
    ##   SCD_10_2016        .        .  0.03125        .        .        .        .
    ##   SCD_11_2016        .        .        .        .        .        .        .
    ##   SCD_12_2016        .        .        .        .        .        .        .
    ##   SCD_10_2017        .        .        .        .        .        .        .
    ##   SCD_10_2012        .        .        .        .        .        .        .
    ##   SCD_11_2012        .        .        .        .        .        .        .
    ##   SCD_12_2012        .        .        .        .        .        .        .
    ##   SCD_8_2013         .        .        .        .        .        .        .
    ##   SCD_9_2013         .        .        .        .        .        .        .
    ##   SCD_10_2013        .        .        .        .        .        .        .
    ##   SCD_11_2013        .        .        .        .        .        .        .
    ##   SCD_1_2014         .        .        .        .        .        .        .
    ##   SCD_3_2014         .        .        .        .        .        .        .
    ##   SCD_7_2014         .        .        .        .        .        .        .
    ##   SCD_9_2014         .        .  0.25000        .        .        .        .
    ##   SCD_11_2014        .        .        .        .  0.50000        .        .
    ##   SCD_7_2015         .        .        .        .        .        .        .
    ##   SCD_8_2015         .        .  0.25000        .  0.50000        .        .
    ##   SCD_9_2015         .        .  0.33333        .  0.33333        .        .
    ##   SCD_10_2015        .        .  0.25000  0.25000  0.50000        .        .
    ##   SCD_11_2015        .        .        .        .        .        .        .
    ##   SCD_2_2016         .        .  1.00000        .        .        .        .
    ##   SCD_8_2016         .        .        .        .  0.33333        .        .
    ##   SCD_9_2016         .        .        .        .        .        .        .
    ##   SCG_3_2011         .        .        .        .        .        .        .
    ##   SCG_2_2012         .        .        .        .        .        .        .
    ##   SCG_3_2012         .        .  1.00000        .        .        .        .
    ##   SCG_7_2011         .        .        .        .        .        .        .
    ##   SCG_11_2011        .        .        .        .        .        .        .
    ##   SCG_10_2012        .        .        .        .        .        .        .
    ##   SCG_11_2012        .        .        .        .  0.33333        .        .
    ##   SCG_1_2013         .        .        .        .        .        .        .
    ##   SCG_1_2014         .        .  1.00000        .  0.50000        .        .
    ##   SCG_6_2014         .        .        .        .        .        .        .
    ##   SCG_11_2014        .        .  0.50000        .  0.50000        .        .
    ##   SCL_7_2012         .        .  0.54762  0.01190  0.02381        .        .
    ##   SCL_9_2014         .        .  0.66667        .        .        .        .
    ##   SCL_10_2014        .        .  0.50000        .        .        .        .
    ##   SCL_6_2015         .        .  1.00000        .  0.50000        .        .
    ##   SCL_7_2015         .        .  0.77778        .        .        .        .
    ##   SCL_9_2012         .        .  0.20000        .        .        .        .
    ##   SCL_10_2012        .        .        .        .        .        .        .
    ##   SCL_8_2016         .        .  0.93333        .  0.06667        .        .
    ##   SCL_9_2016         .        .  0.90000        .  0.10000        .        .
    ##   SCL_10_2016        .        .  0.77419        .  0.19355        .        .
    ##   SCL_11_2012        .        .        .        .        .        .        .
    ##   SCL_4_2012         .        .  0.40000        .        .        .        .
    ##   SCL_11_2016        .        .  1.00000        .        .        .        .
    ##   SCL_12_2016        .        .  1.00000        .        .        .        .
    ##   SCL_2_2017         .        .  1.00000        .        .        .        .
    ##   SCL_6_2017         .        .  1.00000        .        .        .        .
    ##   SCL_12_2012        .        .        .        .        .        .        .
    ##   SCL_7_2017         .        .  0.80000        .        .        .        .
    ##   SCL_1_2013         .        .        .        .        .        .        .
    ##   SCL_8_2017         .        .  1.00000        .        .        .        .
    ##   SCL_7_2013         .        .  0.12500        .  0.12500        .        .
    ##   SCL_9_2017         .        .  0.40000        .        .        .        .
    ##   SCL_10_2017        .        .  1.00000        .        .        .        .
    ##   SCL_11_2017        .        .  1.00000        .        .        .        .
    ##   SCL_12_2017        .        .  1.00000        .        .        .        .
    ##   SCL_6_2018         .        .  1.00000        .        .        .        .
    ##   SCL_7_2018         .        .  1.00000        .        .        .        .
    ##   SCL_8_2018         .        .  1.00000        .  0.14286        .        .
    ##   SCL_8_2013         .        .  0.50000        .  0.50000        .        .
    ##   SCL_10_2013        .        .  0.50000  0.50000        .        .        .
    ##   SCL_1_2014         .        .  1.00000        .        .        .        .
    ##   SCL_4_2014         .        .  0.71429        .  0.14286        .        .
    ##   SCL_6_2014         .        .  0.60000        .        .        .        .
    ##   SCL_7_2014         .        .  0.45833        .  0.08333        .        .
    ##   SCL_8_2014         .        .  0.60000        .        .        .        .
    ##   SCS_8_2013         .        .        .        .        .        .        .
    ##   SCS_1_2014         .        .  0.50000  0.50000  0.50000        .        .
    ##   SCS_5_2014         .        .        .        .        .        .        .
    ##   SCS_7_2014         .        .        .        .        .        .        .
    ##   SCS_8_2014         .        .        .        .        .        .        .
    ##   SCS_9_2014         .        .        .        .        .        .        .
    ##   SCS_11_2014        .        .        .        .        .        .        .
    ##   SCS_12_2014        .        .        .        .        .        .        .
    ##   SCS_1_2016         .        .        .        .        .        .        .
    ##   SCS_8_2016         .        .        .        .        .        .        .
    ##   SCS_9_2016         .        .        .        .        .        .        .
    ##   SCS_10_2013        .        .        .        .        .        .        .
    ##   SCS_10_2016        .        .        .        .        .        .        .
    ##   SCS_9_2017         .        .        .        .        .        .        .
    ##   SCS_NA_2016        .        .        .        .        .        .        .
    ##   SCS_11_2013        .        .        .        .        .        .        .
    ##   TRO_3_2011         .        .        .        .        .        .        .
    ##   TRO_11_2011        .        .        .        .        .        .        .
    ##   TRO_12_2011        .        .        .        .        .        .        .
    ##   TRO_1_2012         .        .        .        .        .        .        .
    ##   TRO_2_2012         .        .        .        .        .        .        .
    ##   TRO_3_2012         .        .        .        .        .        .        .
    ##   TRO_11_2012        .        .        .        .        .        .        .
    ##   TRO_12_2012        .        .        .        .        .        .        .
    ##   TRO_1_2013         .        .        .        .        .        .        .
    ##   TRO_2_2013         .        .        .        .  0.50000        .        .
    ##   TRO_3_2013         .        .        .        .        .        .        .
    ##   TRO_4_2011         .        .        .        .        .        .        .
    ##   TRO_8_2011         .        .        .        .        .        .        .
    ##   TRO_10_2011        .        .        .        .        .        .        .
    ##   UEB_7_2012         .        .        .        .        .        .        .
    ##   UEB_8_2015         .        .  0.42857        .        .        .        .
    ##   UEB_9_2015         .        .  0.07692        .  0.07692        .        .
    ##   UEB_3_2013         .        .        .        .        .        .        .
    ##   UEB_10_2015        .        .        .  1.00000        .        .        .
    ##   UEB_10_2011        .        .        .        .        .        .        .
    ##   UEB_4_2013         .        .        .        .        .        .        .
    ##   UEB_5_2013         .        .        .        .        .        .        .
    ##   UEB_10_2013        .        .        .        .        .        .        .
    ##   UEB_8_2013         .        .        .        .        .        .        .
    ##   UEB_9_2013         .        .        .        .        .        .        .
    ##   UEB_7_2016         .        .        .        .  1.00000        .        .
    ##   UEB_8_2016         .        .        .        .  1.00000        .        .
    ##   UEB_9_2016         .        .        .        .  1.00000        .        .
    ##   UEB_10_2016        .        .        .        .  1.00000        .        .
    ##   UEB_11_2013        .        .        .        .        .        .        .
    ##   UEB_11_2016        .        .        .        .  1.00000        .        .
    ##   UEB_12_2016        .        .        .        .  1.00000        .        .
    ##   UEB_6_2017         .        .        .        .  1.00000        .        .
    ##   UEB_8_2012         .        .        .        .        .        .        .
    ##   UEB_7_2017         .        .        .        .  1.00000        .        .
    ##   UEB_1_2018         .        .        .        .  1.00000        .        .
    ##   UEB_9_2019         .        .        .        .        .        .        .
    ##   UEB_12_2013        .        .        .        .        .        .        .
    ##   UEB_1_2014         .        .        .        .        .        .        .
    ##   UEB_2_2014         .        .        .        .        .        .        .
    ##   UEB_3_2014         .        .        .        .        .        .        .
    ##   UEB_5_2014         .        .  0.33333        .  0.33333        .        .
    ##   UEB_6_2014         .        .        .        .        .        .        .
    ##   UEB_9_2012         .        .        .        .        .        .        .
    ##   UEB_7_2014         .        .  0.27273  0.18182  0.09091        .        .
    ##   UEB_8_2014         .        .        .        .        .        .        .
    ##   UEB_9_2014         .        .        .        .        .        .        .
    ##   UEB_10_2014        .        .        .        .  0.33333        .        .
    ##   UEB_11_2014        .        .        .        .  0.16667        .        .
    ##   UEB_12_2014        .        .        .        .        .        .        .
    ##   UEB_1_2015         .        .  0.20000        .        .        .        .
    ##   UEB_2_2013         .        .        .        .        .        .        .
    ##   UEB_2_2015         .        .        .        .        .        .        .
    ##   UEB_3_2015         .        .  0.25000        .  0.12500        .        .
    ##   UEB_5_2015         .        .        .        .        .        .        .
    ##   UEB_6_2015   1.00000        .  1.00000        .  1.00000        .        .
    ##   UEB_7_2015         .        .        .        .        .        .        .
    ##   UST_9_2013         .        .        .        .        .        .        .
    ##   UST_10_2013        .        .        .        .        .        .        .
    ##   UST_11_2013        .        .        .  0.08333  0.08333        .        .
    ##   UST_12_2013        .        .  0.50000        .        .        .        .
    ##   UST_1_2014         .        .  0.33333        .  0.33333        .        .
    ##   UST_3_2014         .        .        .        .        .        .        .
    ##   UST_8_2014         .        .        .        .        .        .        .
    ##   UST_9_2014         .        .        .        .        .        .        .
    ##   WSL_9_2012         .        .        .        .        .        .        .
    ##   WSL_8_2016         .        .        .        .        .        .        .
    ##   WSL_9_2016         .        .        .        .        .        .        .
    ##   WSL_10_2016        .        .        .        .        .        .        .
    ##   WSL_11_2016        .        .        .        .        .        .        .
    ##   WSL_12_2016        .        .        .        .        .        .        .
    ##   WSL_3_2017         .        .        .        .        .        .        .
    ##   WSL_7_2017         .        .        .        .        .        .        .
    ##   WSL_8_2017         .        .        .        .        .        .        .
    ##   WSL_9_2017         .        .        .        .        .        .        .
    ##   WSL_10_2017        .        .        .        .        .        .        .
    ##   WSL_11_2017        .        .        .        .        .        .        .
    ##   WSL_1_2018         .        .        .        .        .        .        .
    ##   WSL_3_2018         .        .        .        .        .        .        .
    ##   WSL_5_2018         .        .        .        .        .        .        .
    ##   WSL_6_2018         .        .        .        .        .        .        .
    ##   WSL_7_2018         .        .        .        .        .        .        .
    ##   WSL_8_2018         .        .        .        .        .        .        .
    ##   WSL_9_2018         .        .        .        .        .        .        .
    ##   WSL_6_2019         .        .        .        .        .        .        .
    ##   WSL_7_2019         .        .        .        .        .        .        .
    ##   WSL_8_2019         .        .        .        .        .        .        .
    ##   WSL_9_2019         .        .        .        .        .        .        .
    ##   WSL_10_2019        .        .        .        .        .        .        .
    ##   WSL_4_2020         .        .        .        .        .        .        .
    ##   WSL_5_2020         .        .        .        .        .        .        .
    ##   WSL_6_2020         .        .        .        .        .        .        .
    ##   WSL_7_2020         .        .        .        .        .        .        .
    ##   WSL_8_2020         .        .        .        .        .        .        .
    ##   WSL_9_2020         .        .        .        .        .        .        .
    ##   WSL_10_2020        .        .        .        .        .        .        .
    ##   WSL_11_2020        .        .        .        .        .        .        .
    ##   WSL_6_2021         .        .        .        .        .        .        .
    ##   WSL_7_2021         .        .        .        .        .        .        .
    ##   Total        0.00148  0.00258  0.09121  0.00554  0.07090  0.00074  0.00148
    ##              Locus
    ## Population       Mean
    ##   ALD_3_2011        .
    ##   ALD_6_2011        .
    ##   ALD_7_2011  0.03571
    ##   ALD_8_2011        .
    ##   ALD_10_2011       .
    ##   ALD_11_2011       .
    ##   ALD_12_2011       .
    ##   ALD_1_2012        .
    ##   ALD_2_2012        .
    ##   ALD_3_2012  0.03571
    ##   ALD_10_2012       .
    ##   ALD_11_2012       .
    ##   ALD_12_2012       .
    ##   ALD_1_2013        .
    ##   ALD_2_2013        .
    ##   ALD_3_2013  0.03571
    ##   ALD_10_2013       .
    ##   ALD_4_2011        .
    ##   ALD_11_2013       .
    ##   ALD_12_2013 0.03571
    ##   ALD_1_2014  0.14286
    ##   ALD_2_2014        .
    ##   ALD_5_2011        .
    ##   ALD_7_2014  0.07143
    ##   ALD_8_2014        .
    ##   ALD_9_2014        .
    ##   ALD_10_2014       .
    ##   ALD_1_2015        .
    ##   ALD_12_2015 0.05357
    ##   BAR_5_2014  0.07143
    ##   BAR_6_2014        .
    ##   BOB_8_2011        .
    ##   BOB_2_2012        .
    ##   BOB_1_2015  0.04762
    ##   BOB_3_2015  0.01786
    ##   BOB_6_2015        .
    ##   BOB_8_2015  0.03297
    ##   BOB_3_2012        .
    ##   BOB_10_2015 0.10714
    ##   BOB_1_2016        .
    ##   BOB_2_2016  0.14286
    ##   BOB_11_2012       .
    ##   BOB_9_2014  0.02679
    ##   BOB_10_2013 0.00549
    ##   BOB_5_2012        .
    ##   BOB_7_2012  0.01020
    ##   BOB_6_2016  0.14286
    ##   BOB_7_2016        .
    ##   BOB_8_2016  0.03061
    ##   BOB_9_2016        .
    ##   BOB_10_2016       .
    ##   BOB_9_2011        .
    ##   BOB_11_2016 0.01648
    ##   BOB_8_2012        .
    ##   BOB_12_2016       .
    ##   BOB_9_2012  0.01948
    ##   BOB_6_2017        .
    ##   BOB_7_2017        .
    ##   BOB_8_2017        .
    ##   BOB_9_2017        .
    ##   BOB_10_2011       .
    ##   BOB_10_2017 0.00397
    ##   BOB_11_2017       .
    ##   BOB_10_2012 0.01020
    ##   BOB_11_2011       .
    ##   BOB_12_2017       .
    ##   BOB_1_2018        .
    ##   BOB_3_2018        .
    ##   BOB_4_2018        .
    ##   BOB_6_2018        .
    ##   BOB_8_2018  0.00376
    ##   BOB_1_2013        .
    ##   BOB_9_2018        .
    ##   BOB_2_2013  0.14286
    ##   BOB_5_2013  0.10714
    ##   BOB_10_2018 0.04422
    ##   BOB_12_2011       .
    ##   BOB_6_2013  0.07143
    ##   BOB_8_2013        .
    ##   BOB_9_2019        .
    ##   BOB_9_2013        .
    ##   BOB_12_2013 0.01429
    ##   BOB_1_2014        .
    ##   BOB_2_2014        .
    ##   BOB_3_2014        .
    ##   BOB_7_2014  0.04464
    ##   BOB_8_2014        .
    ##   BOB_10_2014 0.02381
    ##   BOB_11_2014 0.02597
    ##   BOB_12_2014       .
    ##   BOH_5_2011        .
    ##   BOH_10_2011       .
    ##   BOH_12_2011       .
    ##   BOH_3_2012  0.07143
    ##   BOH_8_2012        .
    ##   BOH_9_2012  0.05952
    ##   BOH_6_2011        .
    ##   BOH_10_2012       .
    ##   BOH_11_2012 0.05952
    ##   BOH_8_2011        .
    ##   BOH_9_2013  0.04286
    ##   BOH_10_2013 0.03571
    ##   BOH_9_2011  0.02381
    ##   BOH_11_2013       .
    ##   BOH_7_2014  0.07143
    ##   BOH_8_2014        .
    ##   BOH_9_2014        .
    ##   BOH_11_2014       .
    ##   BOH_9_2015  0.17857
    ##   BOH_10_2015 0.21429
    ##   BUR_7_2013        .
    ##   BUR_9_2013        .
    ##   BUR_8_2016  0.13095
    ##   BUR_9_2016  0.01504
    ##   BUR_10_2016 0.06593
    ##   BUR_12_2016 0.14286
    ##   BUR_6_2018        .
    ##   BUR_8_2018  0.02679
    ##   BUR_9_2018  0.01948
    ##   BUR_10_2013 0.01948
    ##   BUR_10_2018 0.01128
    ##   BUR_11_2018 0.07143
    ##   BUR_7_2019  0.07143
    ##   BUR_8_2019        .
    ##   BUR_9_2019  0.11224
    ##   BUR_10_2019 0.08929
    ##   BUR_11_2019 0.12500
    ##   BUR_12_2019 0.08571
    ##   BUR_1_2020  0.14286
    ##   BUR_5_2020        .
    ##   BUR_6_2020        .
    ##   BUR_7_2020        .
    ##   BUR_8_2020  0.00311
    ##   BUR_10_2020       .
    ##   BUR_11_2020       .
    ##   BUR_12_2020       .
    ##   BUR_8_2021        .
    ##   BUR_9_2021  0.02381
    ##   BUR_10_2021 0.01648
    ##   BUR_8_2013        .
    ##   BUR_11_2021 0.02381
    ##   BUR_12_2021       .
    ##   BUR_11_2013 0.02381
    ##   BUR_12_2013       .
    ##   BUR_8_2014        .
    ##   BUR_9_2014        .
    ##   BUR_10_2014 0.02381
    ##   BUR_11_2014       .
    ##   BUR_12_2014 0.10714
    ##   BUR_9_2015  0.01429
    ##   BUR_11_2015       .
    ##   BUR_7_2016  0.14286
    ##   FRB_10_2011       .
    ##   FRB_2_2013  0.03571
    ##   FRB_10_2013       .
    ##   FRB_12_2013       .
    ##   FRB_1_2014  0.14286
    ##   FRB_2_2014        .
    ##   FRB_5_2014  0.21429
    ##   FRB_8_2014  0.07143
    ##   FRB_12_2014 0.14286
    ##   FRB_5_2015  0.07143
    ##   FRB_12_2011       .
    ##   FRB_10_2012       .
    ##   FRB_11_2012 0.07143
    ##   FRE_1_2011        .
    ##   FRE_6_2011  0.03571
    ##   FRE_7_2016  0.03226
    ##   FRE_7_2011        .
    ##   FRE_10_2016 0.06682
    ##   FRE_11_2016       .
    ##   FRE_7_2017  0.01166
    ##   FRE_8_2011        .
    ##   FRE_10_2011       .
    ##   FRE_3_2011  0.02381
    ##   FRE_9_2017  0.02159
    ##   FRE_11_2011 0.02381
    ##   FRE_10_2017 0.05590
    ##   FRE_12_2011 0.07143
    ##   FRE_1_2012        .
    ##   FRE_3_2012        .
    ##   FRE_6_2012        .
    ##   FRE_7_2012        .
    ##   FRE_8_2012  0.01786
    ##   FRE_9_2012        .
    ##   FRE_10_2012       .
    ##   FRE_11_2012       .
    ##   FRE_12_2012       .
    ##   FRE_1_2013        .
    ##   FRE_4_2011  0.03571
    ##   FRE_6_2013  0.14286
    ##   FRE_7_2013        .
    ##   FRE_8_2013        .
    ##   FRE_9_2013        .
    ##   FRE_10_2013 0.02381
    ##   FRE_11_2013       .
    ##   FRE_12_2013       .
    ##   FRE_2_2014  0.03571
    ##   FRE_4_2014  0.03571
    ##   FRE_5_2014  0.04762
    ##   FRE_7_2014        .
    ##   FRE_8_2014  0.17857
    ##   FRE_5_2011        .
    ##   FRE_10_2014       .
    ##   FRE_11_2014 0.07143
    ##   FRE_12_2014       .
    ##   FRE_3_2015  0.07143
    ##   FRE_8_2015        .
    ##   FRE_10_2015 0.14286
    ##   FRE_11_2015       .
    ##   FRE_12_2015       .
    ##   FRE_9_2015        .
    ##   FRI_10_2013 0.00794
    ##   FRI_10_2017       .
    ##   FRI_11_2017       .
    ##   FRI_12_2017       .
    ##   FRI_1_2018        .
    ##   FRI_11_2013       .
    ##   FRI_12_2013       .
    ##   FRI_3_2014        .
    ##   FRI_9_2014        .
    ##   FRI_1_2015  0.03571
    ##   FRI_9_2016        .
    ##   FRI_11_2016 0.02381
    ##   FRI_8_2017        .
    ##   FRI_9_2017        .
    ##   GEN_NA_NA         .
    ##   GEN_8_2020        .
    ##   GEN_12_2021       .
    ##   GEN_1_2022        .
    ##   KON_9_2013        .
    ##   KON_10_2014       .
    ##   KON_12_2016 0.05495
    ##   KON_2_2017  0.08571
    ##   KON_7_2017  0.17857
    ##   KON_9_2017  0.02857
    ##   KON_10_2017 0.02857
    ##   KON_11_2017       .
    ##   KON_12_2017       .
    ##   KON_11_2014       .
    ##   KON_1_2018        .
    ##   KON_7_2018  0.14286
    ##   KON_8_2018        .
    ##   KON_9_2018  0.00549
    ##   KON_10_2018 0.00752
    ##   KON_12_2014       .
    ##   KON_12_2018       .
    ##   KON_7_2019  0.07143
    ##   KON_8_2019  0.09524
    ##   KON_9_2019  0.02381
    ##   KON_10_2019       .
    ##   KON_11_2013 0.02381
    ##   KON_1_2015        .
    ##   KON_11_2019 0.00794
    ##   KON_3_2016  0.14286
    ##   KON_1_2020        .
    ##   KON_7_2020        .
    ##   KON_8_2020        .
    ##   KON_9_2020        .
    ##   KON_10_2020       .
    ##   KON_8_2016  0.07143
    ##   KON_11_2020       .
    ##   KON_12_2020       .
    ##   KON_8_2021        .
    ##   KON_10_2021 0.01786
    ##   KON_11_2021       .
    ##   KON_12_2021       .
    ##   KON_9_2016  0.09524
    ##   KON_10_2016 0.06522
    ##   KON_12_2013 0.14286
    ##   KON_2_2014  0.07143
    ##   KON_9_2014  0.07143
    ##   KON_11_2016 0.07143
    ##   NEU_8_2013        .
    ##   NEU_10_2014 0.02381
    ##   NEU_11_2014       .
    ##   NEU_10_2015       .
    ##   NEU_8_2016  0.03571
    ##   NEU_9_2016        .
    ##   NEU_9_2013        .
    ##   NEU_10_2016 0.01429
    ##   NEU_10_2013 0.21429
    ##   NEU_2_2017        .
    ##   NEU_11_2013       .
    ##   NEU_7_2020        .
    ##   NEU_8_2020        .
    ##   NEU_9_2020        .
    ##   NEU_10_2021       .
    ##   NEU_11_2021       .
    ##   NEU_12_2013 0.07143
    ##   RIE_3_2011        .
    ##   RIE_8_2011        .
    ##   RIE_10_2011       .
    ##   RIE_11_2011       .
    ##   RIE_12_2011       .
    ##   RIE_2_2012        .
    ##   RIE_3_2012        .
    ##   RIE_10_2012       .
    ##   RIE_11_2012       .
    ##   RIE_12_2012       .
    ##   RIE_1_2013        .
    ##   RIE_2_2013        .
    ##   RIE_3_2013        .
    ##   RIE_10_2013       .
    ##   RIE_11_2013       .
    ##   RIE_12_2013 0.01786
    ##   RIE_1_2014  0.02381
    ##   RIE_4_2011        .
    ##   RIE_2_2014        .
    ##   RIE_9_2014        .
    ##   RIE_10_2014       .
    ##   RIE_11_2014       .
    ##   RIE_12_2014       .
    ##   RIE_1_2015  0.02381
    ##   RIE_2_2015        .
    ##   RIE_10_2015 0.14286
    ##   RIE_12_2015       .
    ##   RIE_6_2011        .
    ##   RIE_7_2011  0.07143
    ##   SCD_7_2012        .
    ##   SCD_9_2012  0.01905
    ##   SCD_10_2016 0.00223
    ##   SCD_11_2016       .
    ##   SCD_12_2016       .
    ##   SCD_10_2017       .
    ##   SCD_10_2012       .
    ##   SCD_11_2012 0.01020
    ##   SCD_12_2012       .
    ##   SCD_8_2013  0.04762
    ##   SCD_9_2013  0.03571
    ##   SCD_10_2013 0.02381
    ##   SCD_11_2013       .
    ##   SCD_1_2014        .
    ##   SCD_3_2014        .
    ##   SCD_7_2014        .
    ##   SCD_9_2014  0.01786
    ##   SCD_11_2014 0.17857
    ##   SCD_7_2015        .
    ##   SCD_8_2015  0.07143
    ##   SCD_9_2015  0.08333
    ##   SCD_10_2015 0.14286
    ##   SCD_11_2015       .
    ##   SCD_2_2016  0.14286
    ##   SCD_8_2016  0.07143
    ##   SCD_9_2016        .
    ##   SCG_3_2011        .
    ##   SCG_2_2012        .
    ##   SCG_3_2012  0.14286
    ##   SCG_7_2011        .
    ##   SCG_11_2011       .
    ##   SCG_10_2012       .
    ##   SCG_11_2012 0.09524
    ##   SCG_1_2013  0.03571
    ##   SCG_1_2014  0.10714
    ##   SCG_6_2014        .
    ##   SCG_11_2014 0.14286
    ##   SCL_7_2012  0.04592
    ##   SCL_9_2014  0.04762
    ##   SCL_10_2014 0.07143
    ##   SCL_6_2015  0.14286
    ##   SCL_7_2015  0.05556
    ##   SCL_9_2012  0.01429
    ##   SCL_10_2012 0.01020
    ##   SCL_8_2016  0.10952
    ##   SCL_9_2016  0.08571
    ##   SCL_10_2016 0.08525
    ##   SCL_11_2012       .
    ##   SCL_4_2012  0.02857
    ##   SCL_11_2016 0.07143
    ##   SCL_12_2016 0.07143
    ##   SCL_2_2017  0.07143
    ##   SCL_6_2017  0.09524
    ##   SCL_12_2012       .
    ##   SCL_7_2017  0.05714
    ##   SCL_1_2013        .
    ##   SCL_8_2017  0.07143
    ##   SCL_7_2013  0.03571
    ##   SCL_9_2017  0.02857
    ##   SCL_10_2017 0.07143
    ##   SCL_11_2017 0.07143
    ##   SCL_12_2017 0.07143
    ##   SCL_6_2018  0.07143
    ##   SCL_7_2018  0.07143
    ##   SCL_8_2018  0.09184
    ##   SCL_8_2013  0.10714
    ##   SCL_10_2013 0.07143
    ##   SCL_1_2014  0.14286
    ##   SCL_4_2014  0.06122
    ##   SCL_6_2014  0.04286
    ##   SCL_7_2014  0.04167
    ##   SCL_8_2014  0.04286
    ##   SCS_8_2013        .
    ##   SCS_1_2014  0.14286
    ##   SCS_5_2014        .
    ##   SCS_7_2014        .
    ##   SCS_8_2014        .
    ##   SCS_9_2014        .
    ##   SCS_11_2014       .
    ##   SCS_12_2014 0.02381
    ##   SCS_1_2016        .
    ##   SCS_8_2016  0.01190
    ##   SCS_9_2016        .
    ##   SCS_10_2013 0.01786
    ##   SCS_10_2016 0.01190
    ##   SCS_9_2017  0.02381
    ##   SCS_NA_2016       .
    ##   SCS_11_2013       .
    ##   TRO_3_2011        .
    ##   TRO_11_2011       .
    ##   TRO_12_2011       .
    ##   TRO_1_2012        .
    ##   TRO_2_2012        .
    ##   TRO_3_2012        .
    ##   TRO_11_2012       .
    ##   TRO_12_2012       .
    ##   TRO_1_2013  0.01786
    ##   TRO_2_2013  0.10714
    ##   TRO_3_2013        .
    ##   TRO_4_2011  0.03571
    ##   TRO_8_2011        .
    ##   TRO_10_2011       .
    ##   UEB_7_2012        .
    ##   UEB_8_2015  0.06122
    ##   UEB_9_2015  0.03297
    ##   UEB_3_2013        .
    ##   UEB_10_2015 0.07143
    ##   UEB_10_2011       .
    ##   UEB_4_2013        .
    ##   UEB_5_2013        .
    ##   UEB_10_2013       .
    ##   UEB_8_2013        .
    ##   UEB_9_2013        .
    ##   UEB_7_2016  0.21429
    ##   UEB_8_2016  0.21429
    ##   UEB_9_2016  0.21429
    ##   UEB_10_2016 0.21429
    ##   UEB_11_2013 0.01190
    ##   UEB_11_2016 0.21429
    ##   UEB_12_2016 0.21429
    ##   UEB_6_2017  0.21429
    ##   UEB_8_2012  0.03571
    ##   UEB_7_2017  0.21429
    ##   UEB_1_2018  0.21429
    ##   UEB_9_2019        .
    ##   UEB_12_2013 0.03571
    ##   UEB_1_2014        .
    ##   UEB_2_2014        .
    ##   UEB_3_2014        .
    ##   UEB_5_2014  0.09524
    ##   UEB_6_2014  0.03571
    ##   UEB_9_2012        .
    ##   UEB_7_2014  0.07792
    ##   UEB_8_2014        .
    ##   UEB_9_2014  0.01020
    ##   UEB_10_2014 0.04762
    ##   UEB_11_2014 0.02381
    ##   UEB_12_2014       .
    ##   UEB_1_2015  0.02857
    ##   UEB_2_2013        .
    ##   UEB_2_2015        .
    ##   UEB_3_2015  0.03571
    ##   UEB_5_2015        .
    ##   UEB_6_2015  0.28571
    ##   UEB_7_2015  0.07143
    ##   UST_9_2013  0.00893
    ##   UST_10_2013 0.01429
    ##   UST_11_2013 0.02976
    ##   UST_12_2013 0.03571
    ##   UST_1_2014  0.07143
    ##   UST_3_2014        .
    ##   UST_8_2014        .
    ##   UST_9_2014        .
    ##   WSL_9_2012        .
    ##   WSL_8_2016        .
    ##   WSL_9_2016        .
    ##   WSL_10_2016       .
    ##   WSL_11_2016       .
    ##   WSL_12_2016       .
    ##   WSL_3_2017        .
    ##   WSL_7_2017        .
    ##   WSL_8_2017        .
    ##   WSL_9_2017        .
    ##   WSL_10_2017       .
    ##   WSL_11_2017       .
    ##   WSL_1_2018        .
    ##   WSL_3_2018        .
    ##   WSL_5_2018        .
    ##   WSL_6_2018        .
    ##   WSL_7_2018        .
    ##   WSL_8_2018        .
    ##   WSL_9_2018        .
    ##   WSL_6_2019        .
    ##   WSL_7_2019        .
    ##   WSL_8_2019        .
    ##   WSL_9_2019        .
    ##   WSL_10_2019       .
    ##   WSL_4_2020        .
    ##   WSL_5_2020        .
    ##   WSL_6_2020        .
    ##   WSL_7_2020        .
    ##   WSL_8_2020        .
    ##   WSL_9_2020        .
    ##   WSL_10_2020       .
    ##   WSL_11_2020       .
    ##   WSL_6_2021        .
    ##   WSL_7_2021        .
    ##   Total       0.02369

``` r
setPop(microsats_dates) <- ~Pop/Year
popdata <- poppr(microsats_dates)
#N = Number of individuals, MLG = Number of multilocus genotypes, eMLG = number of expected MLG at the smallest sample size >= 10 based on rarefaction
#SE = Standard error based on eMLG, H = Shannon-Wiener Index of MLG diversity
#G = Stoddart & Taylor's Index of MLG diversity
# lambda = Simpsons index, E.5 = Evenness, Hexp = Neis Expected Heterozygosity
#Ia = Index of association, rbarD = stand. Index of association
```

## Genotypic evenness

Evenness is a measure of the distribution of genotype abundances,
wherein a population with equally abundant genotypes yields a value
equal to 1 and a population dominated by a single genotype is closer to
zero.

``` r
M.tab <- mlg.table(microsats_dates)
```

![](Truffles-First-Steps_files/figure-gfm/genotypic%20evenness-1.png)<!-- -->

## Visualize diversity

``` r
popdata_pop_year <- separate(popdata,Pop,c("Pop","Year"))
```

    ## Warning: Expected 2 pieces. Missing pieces filled with `NA` in 1 rows [114].

``` r
ggplot(popdata_pop_year,aes(Year,lambda)) +
  geom_point() +
  facet_wrap(vars(Pop)) +
  theme_light() +
  labs(y="Simpson's index", title ="Simpson's index over the years")
```

![](Truffles-First-Steps_files/figure-gfm/plot%20diversity-1.png)<!-- -->

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
