Persistence of MLGs with all markers
================
Lia Baumann
2024-05-29

## Genetic information

Erklärung der Datenherkunft:

1.  myData = Daten_Genalex.csv

Aus myData werden danach alle missing (NA) gefiltert:

Tuaest_allMarkersOnly \<- myData %\>% missingno(“geno”, cutoff = 0)

Info zur Tabelle: N = Number of individuals, MLG = Number of multilocus
genotypes, eMLG = number of expected MLG at the smallest sample size \>=
10 based on rarefaction SE = Standard error based on eMLG, H =
Shannon-Wiener Index of MLG diversity G = Stoddart & Taylor’s Index of
MLG diversity lambda = Simpsons index, E.5 = Evenness, Hexp = Neis
Expected Heterozygosity Ia = Index of association, rbarD = stand. Index
of association

– Wegen der wenigen Samples habe ich SCHIF und BAR danach ausgelassen. –

Hier nochmals die gleiche Tabelle, jedoch ohne Populationen mit kleiner
Sample Size (n\<20):

## Persistence over the years (MLG per population)

1.  Import der klonkorrigierten Datei “cc_TUAEST_ALLMARKERSONLY”, welche
    mit MLGSim generiert wurde (enthält Sample, MLG, PSex sowie die
    Marker):

cc_FullTable_Tuaest_allMarkersOnly \<-
read.csv(“C:/Users/liaba/OneDrive - Eidg. Forschungsanstalt
WSL/R/truffles/cc_TUAEST_ALLMARKERSONLY_HWE.csv”)

2.  MLG Duplikate werden entfernt

cc_distinctMLGs_allMarkersOnly \<-
cc_FullTable_Tuaest_allMarkersOnly\[!duplicated(cc_FullTable_Tuaest_allMarkersOnly\$MLG),\]

3.  T_all kommt aus der Masterliste
    Tuaest_Master_genMonit_haploid_MP_Mai2024.xlsx

Hier wird zuerst die Spalte “Code_Analyses2” vom Jahr separiert, und
danach verjoint, damit die MLGs zusammen mit den Sample-Informationen
verfügbar sind:

T_all \<- separate_wider_delim(T_all,Code_Analyses2, delim=“\_”, names =
c(“Code”,“Yr”))

T_all_withMLGs_allMarkersOnly \<-
inner_join(T_all,cc_FullTable_Tuaest_allMarkersOnly,by=c(“Code”=“Sample”))
%\>% select(Code, MLG, Site_1_abrev, truffle_year,
Sampling_year,Sampling_date)

4.  Nun können wir damit die Grafiken generieren:

## Looking at persistence over years

–\> bei LIM und BRU kommen die MLGs doppelt vor bzw. an beiden Orten…
sollte man diese vielleicht zusammen nehmen?

``` r
# Count the number of unique years each individual was found
year_count_k6 <- T_all_withMLGs_withClusters_allMarkersOnly %>%
  group_by(MLG) %>%
  summarize(years_found = n_distinct(Sampling_year))

proportion <- year_count_k6 %>%
  mutate(found_once = ifelse(years_found == 1, "Once", "Multiple")) %>%
  count(found_once) %>%
  mutate(proportion = n / sum(n))
proportion
```

    ## # A tibble: 2 × 3
    ##   found_once     n proportion
    ##   <chr>      <int>      <dbl>
    ## 1 Multiple     126      0.281
    ## 2 Once         323      0.719

``` r
# Join back to the original data to classify each sample as "Once" or "Multiple"
MLGcount_with_classification <- T_all_withMLGs_withClusters_allMarkersOnly %>%
  left_join(year_count_k6, by = "MLG") %>%
  mutate(found_once = ifelse(years_found == 1, "Once", "Multiple")) %>%
  count(found_once) %>%
  mutate(proportion = n / sum(n))
MLGcount_with_classification #over all samples
```

    ## # A tibble: 2 × 3
    ##   found_once     n proportion
    ##   <chr>      <int>      <dbl>
    ## 1 Multiple    2258      0.856
    ## 2 Once         381      0.144

``` r
# Count the number of samples for each individual
sample_count <- T_all_withMLGs_withClusters_allMarkersOnly %>%
  group_by(MLG) %>%
  summarize(total_samples = n())
sample_count
```

    ## # A tibble: 449 × 2
    ##    MLG    total_samples
    ##    <chr>          <int>
    ##  1 MLG1              56
    ##  2 MLG10              1
    ##  3 MLG100             2
    ##  4 MLG101             1
    ##  5 MLG102             4
    ##  6 MLG103             1
    ##  7 MLG104             1
    ##  8 MLG105             3
    ##  9 MLG106             4
    ## 10 MLG107             1
    ## # ℹ 439 more rows

``` r
# Count how many individuals produced each sample count
individuals_per_sample_count <- sample_count %>%
  count(total_samples, name = "number_of_individuals") %>%
  mutate(total_samples_contributed = total_samples * number_of_individuals,
         proportion_of_total_samples = total_samples_contributed / sum(total_samples_contributed)*100)
individuals_per_sample_count
```

    ## # A tibble: 40 × 4
    ##    total_samples number_of_individuals total_samples_contributed
    ##            <int>                 <int>                     <int>
    ##  1             1                   291                       291
    ##  2             2                    47                        94
    ##  3             3                    27                        81
    ##  4             4                    16                        64
    ##  5             5                    12                        60
    ##  6             6                     6                        36
    ##  7             7                     4                        28
    ##  8             8                     3                        24
    ##  9             9                     2                        18
    ## 10            10                     1                        10
    ## # ℹ 30 more rows
    ## # ℹ 1 more variable: proportion_of_total_samples <dbl>

``` r
persistence_bob <- ggplot(subset(T_all_withMLGs_allMarkersOnly, Site_1_abrev %in% "BOB"),aes(Sampling_year,MLG)) + geom_point() + geom_line(aes(group=MLG)) +
  labs(title ="distribution of MLGs in BOB (Bohlingen Buche) over the years") + scale_x_continuous(breaks=2010:2023)
persistence_bob
```

![](Persistence-of-MLGs-with-all-markers_files/figure-gfm/MLGs%20per%20year%20per%20surface-1.png)<!-- -->

![](Persistence-of-MLGs-with-all-markers_files/figure-gfm/more%20MLGs%20without%20code%20display-1.png)<!-- -->

    ## png 
    ##   2

![](Persistence-of-MLGs-with-all-markers_files/figure-gfm/more%20MLGs%20without%20code%20display-2.png)<!-- -->

    ## png 
    ##   2

![](Persistence-of-MLGs-with-all-markers_files/figure-gfm/more%20MLGs%20without%20code%20display-3.png)<!-- -->

    ## png 
    ##   2

![](Persistence-of-MLGs-with-all-markers_files/figure-gfm/more%20MLGs%20without%20code%20display-4.png)<!-- -->

    ## png 
    ##   2

![](Persistence-of-MLGs-with-all-markers_files/figure-gfm/more%20MLGs%20without%20code%20display-5.png)<!-- -->

    ## png 
    ##   2

![](Persistence-of-MLGs-with-all-markers_files/figure-gfm/more%20MLGs%20without%20code%20display-6.png)<!-- -->

    ## png 
    ##   2

![](Persistence-of-MLGs-with-all-markers_files/figure-gfm/more%20MLGs%20without%20code%20display-7.png)<!-- -->![](Persistence-of-MLGs-with-all-markers_files/figure-gfm/more%20MLGs%20without%20code%20display-8.png)<!-- -->![](Persistence-of-MLGs-with-all-markers_files/figure-gfm/more%20MLGs%20without%20code%20display-9.png)<!-- -->![](Persistence-of-MLGs-with-all-markers_files/figure-gfm/more%20MLGs%20without%20code%20display-10.png)<!-- -->

    ## png 
    ##   2

![](Persistence-of-MLGs-with-all-markers_files/figure-gfm/more%20MLGs%20without%20code%20display-11.png)<!-- -->![](Persistence-of-MLGs-with-all-markers_files/figure-gfm/more%20MLGs%20without%20code%20display-12.png)<!-- -->![](Persistence-of-MLGs-with-all-markers_files/figure-gfm/more%20MLGs%20without%20code%20display-13.png)<!-- -->![](Persistence-of-MLGs-with-all-markers_files/figure-gfm/more%20MLGs%20without%20code%20display-14.png)<!-- -->![](Persistence-of-MLGs-with-all-markers_files/figure-gfm/more%20MLGs%20without%20code%20display-15.png)<!-- -->![](Persistence-of-MLGs-with-all-markers_files/figure-gfm/more%20MLGs%20without%20code%20display-16.png)<!-- -->![](Persistence-of-MLGs-with-all-markers_files/figure-gfm/more%20MLGs%20without%20code%20display-17.png)<!-- -->![](Persistence-of-MLGs-with-all-markers_files/figure-gfm/more%20MLGs%20without%20code%20display-18.png)<!-- -->![](Persistence-of-MLGs-with-all-markers_files/figure-gfm/more%20MLGs%20without%20code%20display-19.png)<!-- -->![](Persistence-of-MLGs-with-all-markers_files/figure-gfm/more%20MLGs%20without%20code%20display-20.png)<!-- -->
