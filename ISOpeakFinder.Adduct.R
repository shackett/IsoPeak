options(digits = 13)

CETUSUSED <- FALSE

if(CETUSUSED == TRUE){setwd("/Genomics/grid/users/shackett/ISOpeakFinder/")}else{
setwd("/Users/seanhackett/Desktop/RabinowitzLab/IsoPeak/")}

source("pfLibrary.R")
source("Adduct.Lib.R")
load("Saved_Filez/knownMzRTre.R")
load("Saved_Filez/knownProbs.R")



############### Constants ################
### From prepare peak sizes

peakQualThresh <- 0.8

####From adduct.lib.R ###

isoCov <- 0.7
Threshold <- 300
Signal.Thresh <- 2*Threshold
Supthresh <- -20 
RT.coel.SD <- 0.1

### From Paramter evaluation by SA ###

#nanneal - Number of rounds of Metropolis-Hastings sampling 
nanneal <- 20
#RTn - the number of locally altered polynomials used to assess the RT scaling of the 
RTn = 20
#Degree of the fitted polynomial for RT alignment
RTpolyD = 4
#HETbase - the base used for forming the locally altered polynomials used to assess the relationship between peak abundance and SD 
HETbase <- 1.8
#Degree of the fitted polynomial for assessing heteroscedasticity
HETpolyD = 4
#ppm error (systematic) for mass determination
MZ.syst.error = 2

################ Prepare Peak sizes ##################################################

allPeaks <- read.table("Saved_Filez/4.23.11aligned.csv", sep = ",", header = TRUE)

#filter by peak quality
allPeaks <- allPeaks[allPeaks$maxQuality > peakQualThresh,]

###Define sample names and labeling status: N = N15, P = natural.  Sample class points to expected isotope abundances in combinedProbs

samples <- c("N.1", "N.2", "P.1", "P.2")
sampleclass <- c("N", "N", "P", "P")
replicates <- c(1,1,2,2)
nsamples <- length(samples)
nuniqueSamp <- length(unique(replicates))

blanks <- c("blank.1", "blank.2", "blank.3")
compounds <- unique(MzRTrefine$compound)

### Compare means of blanks and samples - baseline so approach is valid where samples have no value

allPeaks[,names(allPeaks) %in% c(samples, blanks)][allPeaks[,names(allPeaks) %in% c(samples, blanks)] < Threshold] <- Threshold

#switched harmonic to arithmetic mean
SampM = apply(allPeaks[,names(allPeaks) %in% samples], 1, mean)
BLM = apply(allPeaks[,names(allPeaks) %in% blanks], 1, mean)
allPeaks <- cbind(allPeaks, BLM)

valS <- allPeaks[SampM > 2*BLM,]

peaksizeMat <- valS[,colnames(allPeaks) %in% samples] - valS$BLM
peaksizeMat[peaksizeMat < Threshold] <- Threshold

#Make an nsamples x mSTD matrix for taking the dot product abundance

probMat <- matrix(data =NA, ncol = length(sampleclass), nrow = length(combinedProbs[,1]))
for(i in 1:length(sampleclass)){
if(sampleclass[i] == "P"){probMat[,i] <- combinedProbs$uLabp}else{probMat[,i] <- combinedProbs$nLabp}
}	

#Isotopic variants corresponding to each compound

coToiso <- matrix(data = NA, nrow = length(combinedProbs[,1]), ncol = length(compounds))
colnames(coToiso) <- compounds

for(i in 1:length(compounds)){
	coToiso[,i] <- combinedProbs[,1] %in% compounds[i]
	}

#import adduct list from ISOsetup.R
# - adducts of a compound identified by ISOpeakFinder, with a (-) charge.  Adducts are expected to be proportionally related to their derivatized peak (X-H) by a scaling factor to be determined

load("Saved_Filez/negAdducts.R")

#remove the -H form as this is taken as the default ion
combinedAdds <- combinedAdds[-1,]

################### Parameter Evaluation by Simulated Annealing ############################

# Establish annealing schedule # must be a multiple of 20

anneals <- rep(NA, times = nanneal)
for(j in 1:nanneal){
	anneals[j] = 1*j^-0.05
	}
anneals <- (anneals - min(anneals))*(1/(max(anneals) - min(anneals)))

#### Fluxuating variance #####

annealvar <- rep(c(rep(1, times = 5), rep(0.2, times = 5), rep(0.05, times = 5), rep(0.01, times = 5)), times = nanneal/20)

compoundRT <- NULL
for(i in 1:length(compounds)){
compoundRT <- c(compoundRT, combinedProbs[combinedProbs$compound %in% compounds[i],]$RT[1])
}

# Determine the retention times and peak sizes that will be permuted to get the respective alignments and sd(peaksize) fxns.

RTvals = valS$medRt
MZvals = valS$medMz

RTpos <- range(RTvals)[1] + c(0:(RTn-1))*(diff(range(RTvals))/(RTn-1))

hetR <- c(floor(range(log(valS[,colnames(valS) %in% samples], base = HETbase))[1]), ceiling(range(log(valS[,colnames(valS) %in% samples], base = HETbase))[2]))
nhet <- length(hetR[1]:hetR[2])

#track the parameter estimates across iterations

RTtrack <- matrix(NA, ncol = RTn, nrow = nanneal)
MZoffsetrack  <- rep(NA, times = nanneal)
SDtrack <- matrix(NA, ncol = nhet, nrow = nanneal)
RT.SDtrack <- rep(NA, times = nanneal)
RT.coel.SDtrack <- rep(NA, times = nanneal)
MZ.SDtrack <- rep(NA, times = nanneal)

#track the likelihood (more accurately collective support) across iterations

RTliktrack <- matrix(NA, ncol = RTn, nrow = nanneal)
MZliktrack  <- rep(NA, times = nanneal)
SDliktrack <- matrix(NA, ncol = nhet, nrow = nanneal)
RT.SDliktrack <- rep(NA, times = nanneal)
RT.coel.SDliktrack <- rep(NA, times = nanneal)
MZ.SDliktrack <- rep(NA, times = nanneal)

#parameters from previous iteration

RTcoefsO <- rep(1, times = RTn)
MZcoefO <- 0
SDcoefO <- rep(1, times = nhet)
RT.SDcoefO <- 0.1
RT.coel.SDO <- RT.coel.SD
MZ.SDO <- 1

#sampled parameters to be evalutated

RTcoefsE <- rep(NA, times = RTn)
MZcoefE <- NA
SDcoefE <- rep(NA, times = nhet)
RT.SDcoefE <- NA
RT.coel.SDE <- NA
MZ.SDE <- NA

#the likelihood from previous iteration 

RTlik  <- rep(NA, times = RTn)
MZlik <- NA
peakSDlik <- rep(NA, times = nhet)
RT.peakSDlik <- NA
RT.coel.SDlik <- NA
MZ.SDlik <- NA

#the likelihood upon changing the evaluation parameter

RTlikE <- rep(NA, times = RTn)
MZlikE <- NA
peakSDlikE <- rep(NA, times = nhet) 
RT.peakSDlikE <- NA
RT.coel.SDlikE <- NA
MZ.SDlikE <- NA


for(j in 1:nanneal){

#set evaluation parameters equal to the initial parameter estimates for cycle 1
if(j == 1){MZcoefE <- MZcoefO; RTcoefsE <- RTcoefsO; SDcoefE <- SDcoefO; RT.SDcoefE <- RT.SDcoefO; RT.coel.SDE <- RT.coel.SDO; MZ.SDE <- MZ.SDO
	
	}else{

#for cycle j>1, sample evaluation parameters as a fxn of the current paramter values, for a j vs j-1 comparision of parameter fits

MZcoefE <- rnorm(1, MZcoefO, 2*annealvar[j])
	
RTcoefsE <- (rnorm(n = length(RTcoefsO), mean = RTcoefsO, sd = 0.1*annealvar[j]) + rnorm(n = length(RTcoefsO), mean = 1, sd = 0.1*annealvar[j])*annealvar[j])/(annealvar[j]+1)

SDcoefE <- SDcoefO*runif(nhet, 1-annealvar[j]*(2/5), 1+annealvar[j]*(2/5))

RT.SDcoefE <- RT.SDcoefO*runif(1, 1-annealvar[j]*(2/5), 1+annealvar[j]*(2/5))

RT.coel.SDE <- RT.coel.SDO*runif(1, 1-annealvar[j]*(2/5), 1+annealvar[j]*(2/5))

MZ.SDE <- MZ.SDO*runif(1, 1-annealvar[j]*(2/5), 1+annealvar[j]*(2/5))
}

### Determine alignment between centroided samples and standard RTs ###

RTcoefsMat <- matrix(data = RTcoefsO, ncol = RTn+1, nrow = RTn)
diag(RTcoefsMat) <- RTcoefsE

RTeval <- matrix(data = NA, ncol = RTn+1, nrow = length(RTvals))

for(i in 1:(RTn+1)){

RTpoints <- RTpos*RTcoefsMat[,i]
RTcoefs <- lm(RTpoints ~ polym(RTpos, degree = RTpolyD, raw = TRUE))$coef
RTeval[,i] <- poly.fit.vec(RTvals, RTcoefs, RTpolyD)

}	
		
### Determine SD of a peak given its intensity - Heteroscedasticity ###

SDcoefMat <- matrix(data = log((HETbase^(hetR[1]:hetR[2]))*SDcoefO, base = HETbase), ncol = nhet+1, nrow = nhet)
diag(SDcoefMat) <- log(HETbase^(hetR[1]:hetR[2])*SDcoefE, base = HETbase)

SDlmMat <- matrix(data = NA, ncol = nhet+1, nrow = HETpolyD+1)

for(i in 1:(nhet+1)){
	
SDpoints <- c(hetR[1]:hetR[2])
SDlmMat[,i] <- lm(SDcoefMat[,i] ~ polym(SDpoints, degree = HETpolyD, raw = TRUE))$coef
	
	}

	
###### M/Z eval #######	Sample MZ systematic error from current estimate with an sd of the instruments precision

pMZ = MZvals
pRT = RTeval[,RTn+1]

npeaks  <- length(pMZ)
nstd <- length(combinedProbs[,1])

MZe <- log(dnorm(sapply(combinedProbs$mass, masserror, standard = pMZ) + MZcoefE, mean = 0, sd = MZ.SDO), base = 2)
RTe <- log(dnorm(sapply(combinedProbs$RT, RTdiff, standard = pRT), mean = 0, sd = pRT*RT.SDcoefO), base = 2)
SIZe <- t(log(probMat %*% t(peaksizeMat), base = 2))
posL <- MZe + RTe + SIZe


peakLIK <- lapply(c(1:length(compounds)), GauS.w.Adduct, coToiso, posL, peaksizeMat, probMat, npeaks, nhet + 1, HETbase, pMZ, pRT, combinedProbs, combinedAdds, ADDUCT.USE = TRUE, ADDUCT.OUT = FALSE, RT.UNKNOWN, isoCov, RT.coel.SD, Supthresh, SDlmMat, MZcoefE, MZ.SDO, HETpolyD)

#GauS.w.Adduct(203, coToiso, posL, peaksizeMat, probMat, npeaks, nhet + 1, HETbase, pMZ, pRT, combinedProbs, combinedAdds, ADDUCT.USE = TRUE, ADDUCT.OUT = FALSE, RT.UNKNOWN, isoCov, RT.coel.SD, Supthresh, SDlmMat, MZcoefE, MZ.SDO, HETpolyD)


LIKform <- NULL
for (znum in 1:length(peakLIK)){
LIKform <- rbind(LIKform, data.frame(peakLIK[znum[1]]))	}

LIKmat <- matrix(0, ncol = nstd, nrow = npeaks)

for(z in 1:length(LIKform[,1])){
LIKmat[LIKform[z,1],LIKform[z,3]] <- LIKform[z,4]}

MZlikE <- sum(apply(LIKmat, 1, max))/npeaks

##### sd(RT) ####### determine the coefficient that scales the sd of the RT variation (b/w standards and samples) as a linear fxn of RT.

MZe <- log(dnorm(sapply(combinedProbs$mass, masserror, standard = pMZ) + MZcoefO, mean = 0, sd = MZ.SDO), base = 2)
RTe <- log(dnorm(sapply(combinedProbs$RT, RTdiff, standard = pRT), mean = 0, sd = pRT*RT.SDcoefE), base = 2)
SIZe <- t(log(probMat %*% t(peaksizeMat), base = 2))
posL <- MZe + RTe + SIZe

peakLIK <- lapply(c(1:length(compounds)), GauS.w.Adduct, coToiso, posL, peaksizeMat, probMat, npeaks, nhet + 1, HETbase, pMZ, pRT, combinedProbs, combinedAdds, ADDUCT.USE = TRUE, ADDUCT.OUT = FALSE, RT.UNKNOWN, isoCov, RT.coel.SD, Supthresh, SDlmMat, MZcoefO, MZ.SDO, HETpolyD)


LIKform <- NULL
for (znum in 1:length(peakLIK)){
LIKform <- rbind(LIKform, data.frame(peakLIK[znum[1]]))	}

LIKmat <- matrix(0, ncol = nstd, nrow = npeaks)

for(z in 1:length(LIKform[,1])){
LIKmat[LIKform[z,1],LIKform[z,3]] <- LIKform[z,4]}

RT.peakSDlikE <- sum(apply(LIKmat, 1, max))/npeaks


##### RT scaling ######## Alter points defining the relationship between the standards RT and the RT of the centroided samples, then create a polynomial of degree RTpolyD fitting these points.  This polynomial is then used to transform the list of standards for comparison with the sample peaks.

for(i in 1:RTn){
	
MZe <- log(dnorm(sapply(combinedProbs$mass, masserror, standard = MZvals) + MZcoefO, mean = 0, sd = MZ.SDO), base = 2)
RTe <- log(dnorm(sapply(combinedProbs$RT, RTdiff, standard = RTeval[,i]), mean = 0, sd = pRT*RT.SDcoefO), base = 2)
SIZe <- t(log(probMat %*% t(peaksizeMat), base = 2))
posL <- MZe + RTe + SIZe

peakLIK <- lapply(c(1:length(compounds)), GauS.w.Adduct, coToiso, posL, peaksizeMat, probMat, npeaks, nhet + 1, HETbase, pMZ, pRT, combinedProbs, combinedAdds, ADDUCT.USE = TRUE, ADDUCT.OUT = FALSE, RT.UNKNOWN, isoCov, RT.coel.SD, Supthresh, SDlmMat, MZcoefO, MZ.SDO, HETpolyD)

LIKform <- NULL
for (znum in 1:length(peakLIK)){
LIKform <- rbind(LIKform, data.frame(peakLIK[znum[1]]))	}

if(dim(LIKform)[1] == 0){RTlikE[i] <- 0}else{

LIKmat <- matrix(0, ncol = nstd, nrow = npeaks)

for(z in 1:length(LIKform[,1])){
LIKmat[LIKform[z,1],LIKform[z,3]] <- LIKform[z,4]}

RTlikE[i] <- sum(apply(LIKmat, 1, max))/npeaks
}}

######### Heteroscedasticity determination ############# Alter points defining the relationship between a peak's size and its sd(variation in sd(peak) == heteroscedasticity).  Fit a polynomial transforming the sparse peak sizes to their sd. 


for(i in 1:nhet){

### all x all
npeaks  <- length(pMZ)
nstd <- length(combinedProbs[,1])

pMZ = MZvals
pRT = RTeval[,RTn+1]

#compare against all standards - combinedProbs$mass

MZe <- log(dnorm(sapply(combinedProbs$mass, masserror, standard = pMZ) + MZcoefO, mean = 0, sd = 2), base = 2)
RTe <- log(dnorm(sapply(combinedProbs$RT, RTdiff, standard = pRT), mean = 0, sd = pRT*RT.SDcoefO), base = 2)
SIZe <- t(log(probMat %*% t(peaksizeMat), base = 2))
posL <- MZe + RTe + SIZe

#LIKeval <- cbind(matrix(MZe + RTe + SIZe, ncol = length(combinedProbs[,1]), nrow = length(peakpos), byrow = FALSE), -50)

peakLIK <- lapply(c(1:length(compounds)), GaussianLik, coToiso, posL, peaksizeMat, probMat, npeaks, i, SDlmMat, HETbase)

LIKform <- NULL
for (znum in 1:length(peakLIK)){
LIKform <- rbind(LIKform, data.frame(peakLIK[znum[1]]))	}

LIKmat <- matrix(0, ncol = nstd, nrow = npeaks)

for(z in 1:length(LIKform[,1])){
LIKmat[LIKform[z,1],LIKform[z,3]] <- LIKform[z,4]}

peakSDlikE[i] <- sum(apply(LIKmat, 1, max))/npeaks
}


######## Determine acceptance prob - alpha #####

if(j == 1){
MZcoefO <- MZcoefE
MZlik <- MZlikE	
RTcoefsO <- RTcoefsE
RTlik <- RTlikE
SDcoefO <- SDcoefE
peakSDlik <- peakSDlikE
RT.SDcoefO <- RT.SDcoefE
RT.peakSDlik <- RT.peakSDlikE

RTtrack[j,] <- RTcoefsO
MZoffsetrack[j]  <- MZcoefO
SDtrack[j,] <- SDcoefO
RTliktrack[j,] <- RTlikE
MZliktrack[j]  <- MZlikE
SDliktrack[j,] <- peakSDlikE
RT.SDtrack[j] <- RT.SDcoefE
RT.SDliktrack[j] <- RT.peakSDlikE

next
}

MZalpha = 10^(MZlikE - MZlik)
if(MZalpha > 1){
MZaccept <- TRUE}else{if(abs(1-MZalpha) < runif(1, 0, 1)*anneals[j]){MZaccept <- TRUE} else {MZaccept <- FALSE}}

if(MZaccept == TRUE){
MZcoefO <- MZcoefE
MZlik <- MZlikE	
MZoffsetrack[j]  <- MZcoefE
MZliktrack[j]  <- MZlikE
}else{
MZcoefO <- MZcoefO
MZlik <- MZlik	
MZoffsetrack[j]  <- MZcoefO
MZliktrack[j]  <- MZlik}


RT.SDalpha = 10^(RT.peakSDlikE - RT.peakSDlik)
if(RT.SDalpha > 1){
RT.SDaccept <- TRUE}else{if(abs(1-RT.SDalpha) < runif(1, 0, 1)*anneals[j]){RT.SDaccept <- TRUE} else {RT.SDaccept <- FALSE}}

if(RT.SDaccept == TRUE){
RT.SDcoefO <- RT.SDcoefE
RT.peakSDlik <- RT.peakSDlikE
RT.SDtrack[j] <- RT.SDcoefE
RT.SDliktrack[j] <- RT.peakSDlikE
}else{
RT.SDcoefO <- RT.SDcoefO
RT.peakSDlik <- RT.peakSDlik
RT.SDtrack[j] <- RT.SDcoefO
RT.SDliktrack[j] <- RT.peakSDlik}


RTalpha = 10^(RTlikE - RTlik)	
RTaccept = rep(NA, times = RTn)

for(i in 1: RTn){
if(RTalpha[i] > 1){
RTaccept[i] <- TRUE}else{if(abs(1-RTalpha[i]) < runif(1, 0, 1)*anneals[j]){RTaccept[i] <- TRUE} else {RTaccept[i] <- FALSE}}

if(RTaccept[i] == TRUE){
RTcoefsO[i] <- RTcoefsE[i]
RTlik[i] <- RTlikE[i]	
RTtrack[j,i]  <- RTcoefsE[i]
RTliktrack[j,i]  <- RTlikE[i]
}else{
RTcoefsO[i] <- RTcoefsO[i]
RTlik[i] <- RTlik[i]	
RTtrack[j,i]  <- RTcoefsO[i]
RTliktrack[j,i]  <- RTlik[i]}}	

SDalpha = 10^(peakSDlikE - peakSDlik)	
SDaccept = rep(NA, times = nhet)

for(i in 1: nhet){
if(SDalpha[i] > 1){
SDaccept[i] <- TRUE}else{if(abs(1-SDalpha[i]) < runif(1, 0, 1)*anneals[j]){SDaccept[i] <- TRUE} else {SDaccept[i] <- FALSE}}

if(SDaccept[i] == TRUE){
SDcoefO[i] <- SDcoefE[i]
peakSDlik[i] <- peakSDlikE[i]	
SDtrack[j,i]  <- SDcoefE[i]
SDliktrack[j,i]  <- peakSDlikE[i]
}else{
SDcoefO[i] <- SDcoefO[i]
peakSDlik[i] <- peakSDlik[i]	
SDtrack[j,i]  <- SDcoefO[i]
SDliktrack[j,i]  <- peakSDlik[i]}}	
}

save(j, RTtrack, RTliktrack, MZoffsetrack, MZliktrack, SDtrack, SDliktrack, RT.SDtrack, RT.SDliktrack, file = "SpectrumScale.R")











