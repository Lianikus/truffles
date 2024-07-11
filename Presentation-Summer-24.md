Presentation Summer 2024
================
Lia Baumann
2024-07-10

## Nei’s Genetic diversity

# myData_genind_allMarkersOnly locus table mean over all data

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

# myData_genind_allMarkersOnly locus table per Population (Strata: Population / SamplingYear)

``` r
setPop(myData_genind_allMarkersOnly) <- ~Pop/SamplingYear

# create a data frame and transpose the result;
# use sapply to iterate over all populations calculated by seppop.
# then use the locus_table function on the level Genotype

locus_table_SY <- data.frame(t(
  sapply(seppop(myData_genind_allMarkersOnly),
    function(ls) poppr::locus_table(ls,lev="genotype"))))

# choose only the columns corresponding to Hexp values and rename them with their loci

locus_table_SY_hexp <- 
  select(locus_table_SY, X31:X45) %>%
  rename(aest06_1=X31,aest07_1=X32, aest15_1=X33, aest26_1=X34,aest28_1=X35, aest35_1=X36, aest36_1=X37, aest01_1=X38, aest10_1=X39, aest18_1=X40, aest24_1=X41, aest25_1=X42, aest29_1=X43, aest31_1=X44, mean=X45) %>%
  tibble::rownames_to_column("Pop") %>%
  separate(.,Pop,sep="_", into=c("Pop","SamplingYear"))
                                      
#check with one population:
#locus_table(myData_genind_AllMarkersOnly, lev="genotype",population="UEB_2015")
#looks good

as.tibble(locus_table_SY_hexp)
```

    ## # A tibble: 121 × 17
    ##    Pop   SamplingYear aest06_1 aest07_1 aest15_1 aest26_1 aest28_1 aest35_1
    ##    <chr> <chr>           <dbl>    <dbl>    <dbl>    <dbl>    <dbl>    <dbl>
    ##  1 FRE   2011            0.257    0.257    0.257    0        0.257    0    
    ##  2 ALD   2011            0.373    0.159    0        0.391    0.301    0    
    ##  3 RIE   2011            0.608    0.105    0        0.661    0.649    0    
    ##  4 TRO   2011            0.458    0        0        0        0        0.233
    ##  5 SCG   2011            0        0        0        0.806    0.556    0    
    ##  6 BOH   2011            0.182    0.545    0        0.6      0.491    0.182
    ##  7 BOB   2011            0        0.333    0        0.533    0.8      0    
    ##  8 FRB   2011            0.5      0        0        0.5      0        0.5  
    ##  9 UEB   2011            0        0.462    0        0        0        0    
    ## 10 ALD   2012            0.303    0        0        0.303    0.303    0.318
    ## # ℹ 111 more rows
    ## # ℹ 9 more variables: aest36_1 <dbl>, aest01_1 <dbl>, aest10_1 <dbl>,
    ## #   aest18_1 <dbl>, aest24_1 <dbl>, aest25_1 <dbl>, aest29_1 <dbl>,
    ## #   aest31_1 <dbl>, mean <dbl>

``` r
ggplot(locus_table_SY_hexp,aes(x=SamplingYear, y=mean)) +
  geom_point() +
  facet_wrap(vars(Pop))
```

![](Presentation-Summer-24_files/figure-gfm/Hexp%20table%20Pop/SamplingYear-1.png)<!-- -->

# myData_genind_allMarkersOnly locus table per Population (Strata: Population / TruffleYear)

``` r
setPop(myData_genind_allMarkersOnly) <- ~Pop/TruffleYear

# create a data frame and transpose the result;
# use sapply to iterate over all populations calculated by seppop.
# then use the locus_table function on the level Genotype

locus_table_TY <- data.frame(t(
  sapply(seppop(myData_genind_allMarkersOnly),
    function(ls) poppr::locus_table(ls,lev="genotype"))))

# choose only the columns corresponding to Hexp values and rename them with their loci

locus_table_TY_hexp <- 
  select(locus_table_TY, X31:X45) %>%
  rename(aest06_1=X31,aest07_1=X32, aest15_1=X33, aest26_1=X34,aest28_1=X35, aest35_1=X36, aest36_1=X37, aest01_1=X38, aest10_1=X39, aest18_1=X40, aest24_1=X41, aest25_1=X42, aest29_1=X43, aest31_1=X44, mean=X45) %>%
    tibble::rownames_to_column("Pop") %>%
  separate(.,Pop,sep="_", into=c("Pop","TruffleYear"))
                                      
#check with one population:
#locus_table(myData_genind_AllMarkersOnly, lev="genotype",population="UEB_2015")
#looks good

as.tibble(locus_table_TY_hexp)
```

    ## # A tibble: 122 × 17
    ##    Pop   TruffleYear aest06_1 aest07_1 aest15_1 aest26_1 aest28_1 aest35_1
    ##    <chr> <chr>          <dbl>    <dbl>    <dbl>    <dbl>    <dbl>    <dbl>
    ##  1 FRE   2010           0        0        0        0        0        0    
    ##  2 ALD   2010           0        0        0        0        0        0    
    ##  3 RIE   2010           0.333    0        0        0.733    0.6      0    
    ##  4 TRO   2010           0.6      0        0        0        0        0.6  
    ##  5 SCG   2010           0        0        0        0.667    0.667    0    
    ##  6 ALD   2011           0.467    0.153    0        0.453    0.397    0.08 
    ##  7 BOH   2011           0.182    0.545    0        0.6      0.491    0.182
    ##  8 FRE   2011           0.395    0.395    0.395    0        0.416    0    
    ##  9 RIE   2011           0.647    0.221    0        0.676    0.684    0    
    ## 10 SCG   2011           0        0        0        0.571    0.571    0    
    ## # ℹ 112 more rows
    ## # ℹ 9 more variables: aest36_1 <dbl>, aest01_1 <dbl>, aest10_1 <dbl>,
    ## #   aest18_1 <dbl>, aest24_1 <dbl>, aest25_1 <dbl>, aest29_1 <dbl>,
    ## #   aest31_1 <dbl>, mean <dbl>

``` r
ggplot(locus_table_TY_hexp,aes(x=TruffleYear, y=mean)) +
  geom_point() +
  facet_wrap(vars(Pop))
```

![](Presentation-Summer-24_files/figure-gfm/Hexp%20table%20Pop/TruffleYear-1.png)<!-- -->
