options(digits = 13)

CETUSUSED <- FALSE

if(CETUSUSED == TRUE){setwd("/Genomics/grid/users/shackett/ISOpeakFinder/")}else{
setwd("/Users/seanhackett/Desktop/RabinowitzLab/IsoPeak/")}

source("pfLibrary.R")
load("Saved_Filez/knownMzRTre.R")
load("Saved_Filez/knownProbs.R")

#Constants
isoCov <- 0.7
Threshold <- 300
Signal.Thresh <- 2*Threshold

################ Prepare Peak sizes ##################################################

allPeaks <- read.table("Saved_Filez/4.23.11aligned.csv", sep = ",", header = TRUE)

samples <- c("N.1", "N.2", "P.1", "P.2")
sampleclass <- c("N", "N", "P", "P")
replicates <- c(1,1,2,2)
nsamples <- length(samples)
nuniqueSamp <- length(unique(replicates))

blanks <- c("blank.1", "blank.2", "blank.3")
compounds <- unique(MzRTrefine$compound)

### Compare Harmonic means of blanks and samples - baseline so approach is valid where samples have no value

allPeaks[,names(allPeaks) %in% c(samples, blanks)][allPeaks[,names(allPeaks) %in% c(samples, blanks)] < Threshold] <- Threshold

SampM = apply(allPeaks[,names(allPeaks) %in% samples], 1, harmM)
BLM = apply(allPeaks[,names(allPeaks) %in% blanks], 1, harmM)
allPeaks <- cbind(allPeaks, BLM)

valS <- allPeaks[SampM > 2*BLM,]

peaksizeMat <- valS[,colnames(allPeaks) %in% samples] - valS$BLM
peaksizeMat[peaksizeMat < Threshold] <- Threshold

probMat <- matrix(data =NA, ncol = length(sampleclass), nrow = length(combinedProbs[,1]))
for(i in 1:length(sampleclass)){
if(sampleclass[i] == "P"){probMat[,i] <- combinedProbs$uLabp}else{probMat[,i] <- combinedProbs$nLabp}
}	

coToiso <- matrix(data = NA, nrow = length(combinedProbs[,1]), ncol = length(compounds))
colnames(coToiso) <- compounds

for(i in 1:length(compounds)){
	coToiso[,i] <- combinedProbs[,1] %in% compounds[i]
	}

##################laod MZ offset, RT scaling and SD scaling from analysis of compounds ###################
load("SpectrumScale.R")

nanneal <- j

RTvals = valS$medRt
MZvals = valS$medMz

RTn = 20
RTpos <- range(RTvals)[1] + c(0:(RTn-1))*(diff(range(RTvals))/(RTn-1))
RTpoints <- RTpos*RTtrack[nanneal,]
RTcoefs <- summary(lm(RTpoints ~ RTpos + I(RTpos^2) + I(RTpos^3)))$coef[,1]
RTeval <- RTcoefs[1] + RTcoefs[2]*RTvals + RTcoefs[3]*RTvals^2 + RTcoefs[4]*RTvals^3
	

HETbase <- 2.5
hetR <- c(floor(range(log(valS[,colnames(valS) %in% samples], base = 2))[1]), ceiling(range(log(valS[,colnames(valS) %in% samples], base = HETbase))[2]))

nhet <- length(hetR[1]:hetR[2])

SDcoefeval <- log(HETbase^(hetR[1]:hetR[2])*SDtrack[nanneal,], base = HETbase)

SDpoints <- c(hetR[1]:hetR[2])
SDlmMat <- matrix(summary(lm(SDcoefeval ~ SDpoints + I(SDpoints^2) + I(SDpoints^3)))$coef[,1])


#import adduct list from ISOsetup.R
# - adducts of a compound identified by ISOpeakFinder, with a (-) charge.  Adducts are expected to be proportionally related to their derivatized peak (X-H) by a scaling factor to be determined

load("Saved_Filez/negAdducts.R")

combinedAdds <- combinedAdds[-1,]

#import KEGG compounds from ISOsetup.R
# - compounds without a known RT (& aren't known to exist in the sample) but with known isotopic distribution

#load("KEGG/KEGGcompounds.R")
load("KEGG/KEGGcombo.R")

#


##############

#peak vals
pMZ = MZvals
pRT = RTeval

npeaks  <- length(pMZ)
nstd <- length(combinedProbs[,1])

#compare against all standards - combinedProbs$mass

MZe <- log(dnorm(sapply(combinedProbs$mass, masserror, standard = pMZ) + MZoffsetrack[nanneal], mean = 0, sd = 1), base = 2)
RTe <- log(dnorm(sapply(combinedProbs$RT, RTdiff, standard = pRT), mean = 0, sd = pRT*RT.SDtrack[nanneal]), base = 2)
SIZe <- t(log(probMat %*% t(peaksizeMat), base = 2))
posL <- MZe + RTe + SIZe

#LIKeval <- cbind(matrix(MZe + RTe + SIZe, ncol = length(combinedProbs[,1]), nrow = length(peakpos), byrow = FALSE), -50)

#peakLIK <- lapply(c(1:length(compounds)), GaussianLik, coToiso, posL, peaksizeMat, probMat, npeaks, nhet+1, SDlmMat, HETbase)

LIKform <- NULL
for (znum in 1:length(peakLIK)){
LIKform <- rbind(LIKform, data.frame(peakLIK[znum[1]]))	}

LIKmat <- matrix(0, ncol = nstd, nrow = npeaks)

for(z in 1:length(LIKform[,1])){
LIKmat[LIKform[z,1],LIKform[z,3]] <- LIKform[z,4]}

MZlikE <- sum(apply(LIKmat, 1, max))/npeaks


com <- 169

GauS.w.Adduct(com, coToiso, posL, peaksizeMat, probMat, npeaks, h, SDlmMat, HETbase, pMZ, pRT, combinedProbs, combinedAdds, ADDUCT.USE = TRUE, ADDUCT.OUT = FALSE, RT.UNKNOWN = FALSE, isoCov)

peakLIK <- lapply(c(1:length(compounds)), GauS.w.Adduct, coToiso, posL, peaksizeMat, probMat, npeaks, h, SDlmMat, HETbase, pMZ, pRT, combinedProbs, combinedAdds, ADDUCT.USE = TRUE, ADDUCT.OUT = TRUE, RT.UNKNOWN = FALSE, isoCov)





