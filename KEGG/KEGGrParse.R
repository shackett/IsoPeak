setwd("/Users/seanhackett/Desktop/KEGG/")

source(file = "peakLib.R")

options(digits = 13)

elFreq <- data.frame(element = c("H", "H", "O", "O", "O", "C", "C", "N", "N", "Na", "Fe", "Fe", "Fe", "Fe", "K", "K", "K", "S", "S", "S", "S", "P", "Cl", "Cl", "Br", "Br", "F"), MW = c( 1.00782503207, 2.0141017778, 15.99491461956, 16.99913170, 17.9991610, 12.0000000, 13.0033548378, 14.0030740048, 15.0001088982, 22.9897692809, 53.9396105, 55.9349375, 56.9353940, 57.9332756, 38.96370668, 39.96399848, 40.96182576, 31.97207100, 32.97145876, 33.96786690, 35.96708076, 30.973762, 34.968853, 36.965903, 78.918336, 80.916290, 18.998403), abund = c(0.999885, 0.000115, 0.99757, 0.00038, 0.00205, 0.9893, 0.0107, 0.99636, 0.00364, 1.000, 0.05845, 0.91754, 0.02119, 0.00282, 0.932581, 0.000117, 0.067302, 0.9499, 0.0075, 0.0425, 0.0001, 1.000, 0.7577, 0.2423, 0.5069, 0.4931, 1))

#isotopic labeling expected abundance
### 99% N-15

NlabFreq <- data.frame(element = c("H", "H", "O", "O", "O", "C", "C", "N", "N", "Na", "Fe", "Fe", "Fe", "Fe", "K", "K", "K", "S", "S", "S", "S", "P", "Cl", "Cl", "Br", "Br", "F"), MW = c( 1.00782503207, 2.0141017778, 15.99491461956, 16.99913170, 17.9991610, 12.0000000, 13.0033548378, 14.0030740048, 15.0001088982, 22.9897692809, 53.9396105, 55.9349375, 56.9353940, 57.9332756, 38.96370668, 39.96399848, 40.96182576, 31.97207100, 32.97145876, 33.96786690, 35.96708076, 30.973762, 34.968853, 36.965903, 78.918336, 80.916290, 18.998403), abund = c(0.999885, 0.000115, 0.99757, 0.00038, 0.00205, 0.9893, 0.0107, 0.01, 0.99, 1.000, 0.05845, 0.91754, 0.02119, 0.00282, 0.932581, 0.000117, 0.067302, 0.9499, 0.0075, 0.0425, 0.0001, 1.000, 0.7577, 0.2423, 0.5069, 0.4931, 1))

KEGGcomp <- read.table("KEGGyTheKEGG.fq", sep = "\t", colClasses = c("character", "character"))

#remove compounds without a provided formula
KEGGm <- KEGGcomp[!(KEGGcomp[,2] %in% ""),]

#remove polymers
KEGGm <- KEGGm[grep(")", KEGGm[,2], invert = TRUE),]

#remove chemical motifs
KEGGm <-  KEGGm[grep("R", KEGGm[,2], invert = TRUE),]

#remove ion-complexed/adduct chemicals and chemicals with elements which shouldn't be present or present besides as adducts

KEGGm <- KEGGm[grep(". ", KEGGm[,2], invert = TRUE),]

oddEls <- union(union(union(union(union(union(grep("Na", KEGGm[,2]),grep("Ni", KEGGm[,2])), union(union(grep("Cd", KEGGm[,2]), grep("As", KEGGm[,2])), grep("Cl", KEGGm[,2]))), union(union(union(grep("Hg", KEGGm[,2]), grep("X", KEGGm[,2])), union(grep("K", KEGGm[,2]), grep("Cu", KEGGm[,2]))), union(union(grep("Zn", KEGGm[,2]), grep("Mo", KEGGm[,2])), union(grep("W", KEGGm[,2]), grep("F", KEGGm[,2]))))), union(grep("Mn", KEGGm[,2]), grep("Ca", KEGGm[,2]))), union(union(grep("I", KEGGm[,2]), grep("Mg", KEGGm[,2])), union(union(grep("Co", KEGGm[,2]), grep("Se", KEGGm[,2])), grep("Br", KEGGm[,2])))), union(union(union(grep("Sr", KEGGm[,2]), grep("Ag", KEGGm[,2])), union(grep("Al", KEGGm[,2]), grep("Be", KEGGm[,2]))), union(grep("Ba", KEGGm[,2]), grep("Bi", KEGGm[,2]))))


KEGGm2 <- KEGGm[c(1:length(KEGGm[,2]))[!(c(1:length(KEGGm[,2])) %in% oddEls)],]

### determine unique molecular weights

uniqueCF <- unique(KEGGm2[,2])







KEGGunLab <- NULL
KEGGnLab <- NULL

for (i in 1:length(knownRT[,1])){

chem.formula <- knownRT[i,2]
chem.weight <- chemweight(elFreq, chem.formula, -1, 1)
chem.weightlab <- chemweight(NlabFreq, chem.formula, -1, 1)

MzRT <- rbind(MzRT, data.frame(compound = knownRT$compound[i], chem.weight, RT = knownRT$rt[i]))
NlabMzRT <- rbind(NlabMzRT, data.frame(compound = knownRT$compound[i], chem.weightlab, RT = knownRT$rt[i]))
}


