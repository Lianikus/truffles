Ben's Input

install.packages("dartR.data")
install.packages("dartR")
install.packages("SNPRelate")
library(dartR)
library(dartR.data)
library(SNPRelate)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SNPRelate")
gl.impute(myData_genind, method="frequency")
myData_genlight <- as.genlight(myData_genind)
head(myData_genind$tab,10)
myData_missingMean <- missingno(myData, type ="mean")
rownames(myData_missingMean$tab)
View(myData_genind$tab)
myData_missingMean$tab[15,]

myData_genlight <- gi2gl(myData_genind)

# access mean
mean(myData_genlight$gen[15,"aest31_1.313"])