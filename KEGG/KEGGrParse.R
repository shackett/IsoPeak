setwd("/Users/seanhackett/Desktop/RabinowitzLab/IsoPeak/")

source(file = "peakLib.R")
knowns <- read.table(file = "/Users/seanhackett/Desktop/RabinowitzLab/IsoPeak/knownRTchemForm.csv", sep = ",", header = TRUE, colClasses = c("character", "character", "character", "numeric"))

options(digits = 13)

elFreq <- data.frame(element = c("H", "H", "O", "O", "O", "C", "C", "N", "N", "Na", "Fe", "Fe", "Fe", "Fe", "K", "K", "K", "S", "S", "S", "S", "P", "Cl", "Cl", "Br", "Br", "F"), MW = c( 1.00782503207, 2.0141017778, 15.99491461956, 16.99913170, 17.9991610, 12.0000000, 13.0033548378, 14.0030740048, 15.0001088982, 22.9897692809, 53.9396105, 55.9349375, 56.9353940, 57.9332756, 38.96370668, 39.96399848, 40.96182576, 31.97207100, 32.97145876, 33.96786690, 35.96708076, 30.973762, 34.968853, 36.965903, 78.918336, 80.916290, 18.998403), abund = c(0.999885, 0.000115, 0.99757, 0.00038, 0.00205, 0.9893, 0.0107, 0.99636, 0.00364, 1.000, 0.05845, 0.91754, 0.02119, 0.00282, 0.932581, 0.000117, 0.067302, 0.9499, 0.0075, 0.0425, 0.0001, 1.000, 0.7577, 0.2423, 0.5069, 0.4931, 1))

#isotopic labeling expected abundance
### 99% N-15

NlabFreq <- data.frame(element = c("H", "H", "O", "O", "O", "C", "C", "N", "N", "Na", "Fe", "Fe", "Fe", "Fe", "K", "K", "K", "S", "S", "S", "S", "P", "Cl", "Cl", "Br", "Br", "F"), MW = c( 1.00782503207, 2.0141017778, 15.99491461956, 16.99913170, 17.9991610, 12.0000000, 13.0033548378, 14.0030740048, 15.0001088982, 22.9897692809, 53.9396105, 55.9349375, 56.9353940, 57.9332756, 38.96370668, 39.96399848, 40.96182576, 31.97207100, 32.97145876, 33.96786690, 35.96708076, 30.973762, 34.968853, 36.965903, 78.918336, 80.916290, 18.998403), abund = c(0.999885, 0.000115, 0.99757, 0.00038, 0.00205, 0.9893, 0.0107, 0.01, 0.99, 1.000, 0.05845, 0.91754, 0.02119, 0.00282, 0.932581, 0.000117, 0.067302, 0.9499, 0.0075, 0.0425, 0.0001, 1.000, 0.7577, 0.2423, 0.5069, 0.4931, 1))

KEGGcomp <- read.table("KEGG/KEGGyTheKEGG.fq", sep = "\t", colClasses = c("character", "character"))

#remove compounds without a provided formula
KEGGm <- KEGGcomp[!(KEGGcomp[,2] %in% ""),]

#remove polymers
KEGGm <- KEGGm[grep(")", KEGGm[,2], invert = TRUE),]

#remove chemical motifs
KEGGm <-  KEGGm[grep("R", KEGGm[,2], invert = TRUE),]

#remove ion-complexed/adduct chemicals and chemicals with elements which shouldn't be present or present besides as adducts

KEGGm <- KEGGm[grep(". ", KEGGm[,2], invert = TRUE),]

unexpectedEL <- c("Na", "Ni", "Cd", "As", "Cl", "Hg", "X", "K", "Cu", "Zn", "Mo", "W", "F", "Mn", "Ca", "I", "Mg", "Co", "Se", "Sr", "Ag", "Al", "B", "Ge", "Si", "Te", "Au", "Sn", "Pt")

unEX <- NULL
for(i in 1:length(unexpectedEL)){
	unEX <- union(unEX, grep(unexpectedEL[i], KEGGm[,2]))}
	

KEGGm2 <- KEGGm[c(1:length(KEGGm[,2]))[!(c(1:length(KEGGm[,2])) %in% unEX)],]

### determine unique molecular weights

uniqueCF <- unique(KEGGm2[,2])

#remove molecular formulas that are the same as knowns or that have no hydrogen (causes an error with molecular abundance and don't contain metabolites)

uniqueCF <- uniqueCF[grep("H", uniqueCF, invert = FALSE)]

knowns$formula[knowns$formula %in% "P2H4O7"] <- "H4P2O7"

#knowns not in KEGG - mostly off by PO3 - knowns[!(knowns$formula %in% uniqueCF),]

uniqueCF <- uniqueCF[!(uniqueCF %in% knowns$formula)]

uniqueCF <- uniqueCF[!(uniqueCF %in% c("H"))] 


clen = function(x){length(unlist(strsplit(x, split = "")))}

CFcom <- rep(NA, times = length(uniqueCF))

for (i in 1:length(uniqueCF)){
vcoms <- KEGGm2[,1][KEGGm2[,2] %in% uniqueCF[i]]
namelen <- unlist(lapply(vcoms, clen))
CFcom[i] <- paste(c(vcoms[namelen == min(namelen)][1], length(vcoms)), collapse = "-")
}

KEGGunLab <- NULL
KEGGnLab <- NULL

workingUnLab <- NULL
workingnLab <- NULL

length(CFcom)

for (i in 1:length(CFcom)){

chem.formula <- uniqueCF[i]
chem.weight <- chemweight(elFreq, chem.formula, -1, 1)
chem.weightlab <- chemweight(NlabFreq, chem.formula, -1, 1)

workingUnLab <- rbind(workingUnLab, data.frame(compound = CFcom[i], chem.weight))
workingnLab <- rbind(workingnLab, data.frame(compound = CFcom[i], chem.weightlab))
}

KEGGunLab <- rbind(KEGGunLab, workingUnLab)
KEGGnLab <- rbind(KEGGnLab, workingnLab)


#save(KEGGunLab, KEGGnLab, file = "KEGGcompounds.R")
#load("KEGGcompounds.R")

### Determine exact mass of each compound under the two labeling conditions

exactUn <- NULL 
exactN <- NULL
for(i in 1:length(CFcom)){
exactUn <- rbind(exactUn, KEGGunLab[KEGGunLab$compound == CFcom[i],][KEGGunLab[KEGGunLab$compound == CFcom[i],3]== max(KEGGunLab[KEGGunLab$compound == CFcom[i],]$prob),])
exactN <- rbind(exactN, KEGGnLab[KEGGnLab$compound == CFcom[i],][KEGGnLab[KEGGnLab$compound == CFcom[i],3]== max(KEGGnLab[KEGGnLab$compound == CFcom[i],]$prob),])
}

#save(KEGGunLab, KEGGnLab, exactUn, exactN, file = "KEGG/KEGGcompounds.R")

load("KEGG/KEGGcompounds.R")

keggcom <- unique(exactUn$compound)
	
KEGGcombo <- NULL

for(i in 1:length(keggcom)){
sublist <- KEGGunLab[KEGGunLab$compound %in% keggcom[i],c(1,3,4)]
nsublist <- KEGGnLab[KEGGnLab$compound %in% keggcom[i],c(1,3,4)]

jointmass <- union(sublist$mass, nsublist$mass)
if(length(jointmass[!(jointmass %in% sublist$mass)]) != 0){
sublist  <- rbind(sublist, data.frame(compound = sublist$compound[1], mass = jointmass[!(jointmass %in% sublist$mass)], prob = 0)) 
sublist <- sublist[order(sublist$mass),]
	}
	
if(length(jointmass[!(jointmass %in% nsublist$mass)]) != 0){
nsublist  <- rbind(nsublist, data.frame(compound = nsublist$compound[1], mass = jointmass[!(jointmass %in% nsublist$mass)], prob = 0)) 
nsublist <- nsublist[order(nsublist$mass),]
	}	

KEGGcombo <- rbind(KEGGcombo, data.frame(compound = sublist$compound, mass = sublist$mass, uLabp = sublist$prob, nLabp = nsublist$prob))}

save(KEGGcombo, exactUn, exactN, file = "KEGG/KEGGcombo.R")


