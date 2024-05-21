Code Archiv

#Create genalex file with dates of sampling:
I want to compare genetic diversity across years and show it for populations. So I need to combine the sample and assign the sampling date.
I have several ideas how to do it.
1.
--> Join Genalex file with monitoring file on the sampling column.
--> then read it in as a genclone and genind file again.
2.
--> Filter monitoring file with genalex file to exclude the rows not matching
--> choose columns that I want to process as a new file and read in as genclone and genind file again.

```{r add dates into microsat file, eval=FALSE, include=FALSE}

##Achtung hier wird eine Tabelle generiert, nur bei Datenänderungen laufen lassen!!

#genalex_dates <- read.csv("Daten_Genalex.csv",sep=";", header=FALSE)
#T_all_dates <- T_all %>%
#  mutate(Sampling_month = month(Sampling_date)) %>%
#  select(Code_Analyses2,Sampling_month)
#genalex_dates <- left_join(genalex_dates,T_all_dates,by=c("V1"="Code_Analyses_2024"))

#genalex_dates <- genalex_dates %>%
#  unite(col="pop_date",c("V2","Sampling_month","Sampling_year"),sep="_")
#correct column names
#genalex_dates$pop_date[genalex_dates$pop_date =="2708_NA_NA"] <- "2708"
#genalex_dates$pop_date[genalex_dates$pop_date =="_NA_NA"] <- ""
#genalex_dates$pop_date[genalex_dates$pop_date =="pop_NA_NA"] <- "pop"

#jetzt als csv Datei neu reinladen
#write.table(genalex_dates, col.names=FALSE, sep=",", "C:/Users/liaba/OneDrive - Eidg. Forschungsanstalt WSL/R/truffles/genalex_dates.csv")
```


Evtl. später:
  Standorte mit weniger als 10 Samples werden entfernt.
Dies lässt danach noch folgende Standorte für die Auswertung zu:
  
  ```{r Standorte mit weniger als 10 Samples entfernen, include=FALSE}
# of observations per site
T_all_corr1 <- subset(T_all,Site_1_abrev!="BAR")
T_all_corr2 <- subset(T_all_corr1, Site_1_abrev != "GEN")
```