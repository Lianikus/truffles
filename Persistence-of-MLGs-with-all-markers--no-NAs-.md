Persistence of MLGs with all markers (no NAs)
================
Lia Baumann
2024-05-29

## Genetic information

Erklärung der Datenherkunft:

1.  myData = Daten_Genalex.csv

Aus myData werden danach alle missing (NA) gefiltert:

Tuaest_allMarkersOnly \<- myData %\>% missingno(“geno”, cutoff = 0)

``` r
setPop(Tuaest_allMarkersOnly) <- ~Pop
poppr_withoutNA <- poppr(Tuaest_allMarkersOnly)
kable(select(poppr_withoutNA,Pop,N,MLG,eMLG,SE,H,G,lambda,E.5,Hexp,Ia,rbarD))
```

| Pop | N | MLG | eMLG | SE | H | G | lambda | E.5 | Hexp | Ia | rbarD |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| FRE | 214 | 6 | 1.539399 | 0.6722074 | 0.3010990 | 1.121406 | 0.1082627 | 0.3455497 | 0.0674671 | 9.8312168 | 0.8435767 |
| ALD | 62 | 12 | 4.094020 | 1.1922532 | 1.4301026 | 2.310096 | 0.5671176 | 0.4120929 | 0.1400242 | 3.6838855 | 0.3652128 |
| RIE | 62 | 39 | 9.175120 | 0.8155850 | 3.5065466 | 28.264706 | 0.9646202 | 0.8432482 | 0.4263051 | 0.8823601 | 0.0778223 |
| TRO | 34 | 6 | 3.539940 | 0.9044355 | 1.1637449 | 2.312000 | 0.5674740 | 0.5958486 | 0.1529157 | 2.4640775 | 0.3408221 |
| SCG | 15 | 6 | 5.138861 | 0.7032333 | 1.6170532 | 4.411765 | 0.7733333 | 0.8448681 | 0.2891156 | 2.2825279 | 0.2888403 |
| BOH | 40 | 27 | 8.363669 | 1.1159864 | 3.0081359 | 12.903226 | 0.9225000 | 0.6183617 | 0.2779304 | 1.3529540 | 0.1363888 |
| BOB | 429 | 125 | 8.053505 | 1.2226228 | 3.7299217 | 16.156703 | 0.9381062 | 0.3726217 | 0.3379152 | 2.1478566 | 0.1899626 |
| FRB | 13 | 7 | 5.846154 | 0.7692308 | 1.7118451 | 4.567568 | 0.7810651 | 0.7859511 | 0.5064103 | 6.0738307 | 0.4700290 |
| UEB | 153 | 29 | 6.165380 | 1.2801106 | 2.4477034 | 6.814847 | 0.8532616 | 0.5505564 | 0.1454126 | 0.7315988 | 0.0909901 |
| SCL | 99 | 39 | 6.359011 | 1.3989524 | 2.6490635 | 6.114161 | 0.8364453 | 0.3891822 | 0.1762376 | 2.2112657 | 0.2492060 |
| SCD | 100 | 40 | 6.000809 | 1.4565731 | 2.5226599 | 4.163197 | 0.7598000 | 0.2759798 | 0.2686003 | 2.5489256 | 0.2569022 |
| WSL | 283 | 4 | 1.300608 | 0.5021267 | 0.1679418 | 1.066190 | 0.0620809 | 0.3619560 | 0.0144994 | 4.2811571 | 0.5055186 |
| BUR | 223 | 19 | 3.898280 | 0.9326109 | 1.5637569 | 3.503276 | 0.7145529 | 0.6628151 | 0.2315650 | 4.7673727 | 0.4420210 |
| SCS | 42 | 6 | 2.190476 | 0.9047313 | 0.5566229 | 1.283843 | 0.2210884 | 0.3811146 | 0.0505226 | 3.9530106 | 0.4511154 |
| NEU | 32 | 6 | 3.119230 | 0.9106796 | 0.9606248 | 1.835125 | 0.4550781 | 0.5176412 | 0.1031106 | 2.5052006 | 0.3431153 |
| UST | 35 | 22 | 8.356705 | 1.0497898 | 2.8767637 | 13.764045 | 0.9273469 | 0.7617272 | 0.3733493 | 1.1158778 | 0.0960655 |
| KON | 212 | 25 | 5.421408 | 1.2057463 | 2.1288833 | 5.213921 | 0.8082058 | 0.5690278 | 0.2292319 | 1.2911268 | 0.1186808 |
| FRI | 120 | 11 | 3.802453 | 0.9470008 | 1.4376132 | 3.138622 | 0.6813889 | 0.6661060 | 0.2952481 | 4.2495522 | 0.3842205 |
| BAR | 1 | 1 | 1.000000 | 0.0000000 | 0.0000000 | 1.000000 | 0.0000000 | NaN | NaN | NA | NA |
| BRU | 178 | 12 | 2.962104 | 0.9294881 | 1.0212038 | 1.836541 | 0.4554980 | 0.4708833 | 0.4326432 | 12.1705207 | 0.9362151 |
| LIM | 216 | 15 | 3.771861 | 0.9168608 | 1.4707184 | 3.179501 | 0.6854853 | 0.6501392 | 0.4155131 | 5.1943548 | 0.4137422 |
| HAN | 65 | 7 | 4.054125 | 0.7804336 | 1.4532618 | 3.494624 | 0.7138462 | 0.7612424 | 0.2884272 | 2.1145302 | 0.2207899 |
| GEN | 10 | 2 | 2.000000 | 0.0000000 | 0.3250830 | 1.219512 | 0.1800000 | 0.5714298 | 0.0857143 | 5.0000000 | 1.0000000 |
| SCHIF | 1 | 1 | 1.000000 | 0.0000000 | 0.0000000 | 1.000000 | 0.0000000 | NaN | NaN | NA | NA |
| Total | 2639 | 449 | 8.910634 | 0.9757313 | 4.5410836 | 34.529314 | 0.9710391 | 0.3613369 | 0.6444092 | 3.0044991 | 0.2324626 |

Info zur Tabelle: N = Number of individuals, MLG = Number of multilocus
genotypes, eMLG = number of expected MLG at the smallest sample size \>=
10 based on rarefaction SE = Standard error based on eMLG, H =
Shannon-Wiener Index of MLG diversity G = Stoddart & Taylor’s Index of
MLG diversity lambda = Simpsons index, E.5 = Evenness, Hexp = Neis
Expected Heterozygosity Ia = Index of association, rbarD = stand. Index
of association

– Wegen der wenigen Samples habe ich SCHIF und BAR danach ausgelassen. –

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
ggplot(subset(T_all_withMLGs_allMarkersOnly, Site_1_abrev %in% "BOB"),aes(Sampling_year,MLG)) + geom_point() + geom_line(aes(group=MLG)) +
  labs(title ="distribution of MLGs in BOB (Bohlingen Buche) over the years") + scale_x_continuous(breaks=2010:2023)
```

![](Persistence-of-MLGs-with-all-markers--no-NAs-_files/figure-gfm/MLGs%20per%20year%20per%20surface-1.png)<!-- -->

![](Persistence-of-MLGs-with-all-markers--no-NAs-_files/figure-gfm/more%20MLGs%20without%20code%20display-1.png)<!-- -->![](Persistence-of-MLGs-with-all-markers--no-NAs-_files/figure-gfm/more%20MLGs%20without%20code%20display-2.png)<!-- -->![](Persistence-of-MLGs-with-all-markers--no-NAs-_files/figure-gfm/more%20MLGs%20without%20code%20display-3.png)<!-- -->![](Persistence-of-MLGs-with-all-markers--no-NAs-_files/figure-gfm/more%20MLGs%20without%20code%20display-4.png)<!-- -->![](Persistence-of-MLGs-with-all-markers--no-NAs-_files/figure-gfm/more%20MLGs%20without%20code%20display-5.png)<!-- -->![](Persistence-of-MLGs-with-all-markers--no-NAs-_files/figure-gfm/more%20MLGs%20without%20code%20display-6.png)<!-- -->![](Persistence-of-MLGs-with-all-markers--no-NAs-_files/figure-gfm/more%20MLGs%20without%20code%20display-7.png)<!-- -->![](Persistence-of-MLGs-with-all-markers--no-NAs-_files/figure-gfm/more%20MLGs%20without%20code%20display-8.png)<!-- -->![](Persistence-of-MLGs-with-all-markers--no-NAs-_files/figure-gfm/more%20MLGs%20without%20code%20display-9.png)<!-- -->![](Persistence-of-MLGs-with-all-markers--no-NAs-_files/figure-gfm/more%20MLGs%20without%20code%20display-10.png)<!-- -->![](Persistence-of-MLGs-with-all-markers--no-NAs-_files/figure-gfm/more%20MLGs%20without%20code%20display-11.png)<!-- -->![](Persistence-of-MLGs-with-all-markers--no-NAs-_files/figure-gfm/more%20MLGs%20without%20code%20display-12.png)<!-- -->![](Persistence-of-MLGs-with-all-markers--no-NAs-_files/figure-gfm/more%20MLGs%20without%20code%20display-13.png)<!-- -->![](Persistence-of-MLGs-with-all-markers--no-NAs-_files/figure-gfm/more%20MLGs%20without%20code%20display-14.png)<!-- -->![](Persistence-of-MLGs-with-all-markers--no-NAs-_files/figure-gfm/more%20MLGs%20without%20code%20display-15.png)<!-- -->![](Persistence-of-MLGs-with-all-markers--no-NAs-_files/figure-gfm/more%20MLGs%20without%20code%20display-16.png)<!-- -->![](Persistence-of-MLGs-with-all-markers--no-NAs-_files/figure-gfm/more%20MLGs%20without%20code%20display-17.png)<!-- -->![](Persistence-of-MLGs-with-all-markers--no-NAs-_files/figure-gfm/more%20MLGs%20without%20code%20display-18.png)<!-- -->![](Persistence-of-MLGs-with-all-markers--no-NAs-_files/figure-gfm/more%20MLGs%20without%20code%20display-19.png)<!-- -->![](Persistence-of-MLGs-with-all-markers--no-NAs-_files/figure-gfm/more%20MLGs%20without%20code%20display-20.png)<!-- -->
