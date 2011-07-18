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

###### constants #######
	
combined <- combinedProbs


#### M/Z offset - ~3ppm
#### Scaling b/w 0.7 and 1.5

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

#track the likelihood (more accurately collective support) across iterations

RTliktrack <- matrix(NA, ncol = RTn, nrow = nanneal)
MZliktrack  <- rep(NA, times = nanneal)
SDliktrack <- matrix(NA, ncol = nhet, nrow = nanneal)
RT.SDliktrack <- rep(NA, times = nanneal)
RT.coel.SDliktrack <- rep(NA, times = nanneal)

#parameters from previous iteration

RTcoefsO <- rep(1, times = RTn)
MZcoefO <- 0
SDcoefO <- rep(1, times = nhet)
RT.SDcoefO <- 0.1
RT.coel.SDO <- RT.coel.SD

#sampled parameters to be evalutated

RTcoefsE <- rep(NA, times = RTn)
MZcoefE <- NA
SDcoefE <- rep(NA, times = nhet)
RT.SDcoefE <- NA
RT.coel.SDE <- NA

#the likelihood from previous iteration 

RTlik  <- rep(NA, times = RTn)
MZlik <- NA
peakSDlik <- rep(NA, times = nhet)
RT.peakSDlik <- NA
RT.coel.SDlik <- NA

#the likelihood upon changing the evaluation parameter

RTlikE <- rep(NA, times = RTn)
MZlikE <- NA
peakSDlikE <- rep(NA, times = nhet) 
RT.peakSDlikE <- NA
RT.coel.SDlikE <- NA



for(j in 1:nanneal){

#ppm offset - systematic difference between peaks and standards
if(j == 1){MZcoefE <- MZcoefO; RTcoefsE <- RTcoefsO; SDcoefE <- SDcoefO; RT.SDcoefE <- RT.SDcoefO; RT.coel.SDE <- RT.coel.SDO
	
	}else{

MZcoefE <- rnorm(1, MZcoefO, 2*annealvar[j])
	
RTcoefsE <- (rnorm(n = length(RTcoefsO), mean = RTcoefsO, sd = 0.1*annealvar[j]) + rnorm(n = length(RTcoefsO), mean = 1, sd = 0.1*annealvar[j])*annealvar[j])/(annealvar[j]+1)

SDcoefE <- SDcoefO*runif(nhet, 1-annealvar[j]*(2/5), 1+annealvar[j]*(2/5))

RT.SDcoefE <- RT.SDcoefO*runif(1, 1-annealvar[j]*(2/5), 1+annealvar[j]*(2/5))

RT.coel.SDE <- RT.coel.SDO*runif(1, 1-annealvar[j]*(2/5), 1+annealvar[j]*(2/5))

}

RTcoefsMat <- matrix(data = RTcoefsO, ncol = RTn+1, nrow = RTn)
diag(RTcoefsMat) <- c(1:20)
diag(RTcoefsMat) <- RTcoefsE

RTeval <- matrix(data = NA, ncol = RTn+1, nrow = length(RTvals))

for(i in 1:(RTn+1)){

RTpoints <- RTpos*RTcoefsMat[,i] + rnorm(length(RTpos*RTcoefsMat[,i]), mean = 0, sd = 5)
RTeval[,i] <- predict(smooth.spline(RTpos*RTcoefsMat[,i], RTpoints, df = RTpolyD, cv = TRUE, all.knots=TRUE), RTvals)$y

	}	
		
### Determine SD of a peak given its intensity - Heteroscedasticity ###

SDcoefMat <- matrix(data = log((HETbase^(hetR[1]:hetR[2]))*SDcoefO, base = HETbase), ncol = nhet+1, nrow = nhet)
diag(SDcoefMat) <- log(HETbase^(hetR[1]:hetR[2])*SDcoefE, base = HETbase)

sd.knot <- NULL
sd.nk <- NULL
sd.min <- NULL
sd.range <- NULL
sd.coef	<- NULL

for(i in 1:(nhet+1)){

SDpoints <- c(hetR[1]:hetR[2])
#SDlmMat[,i] <- summary(lm(SDcoefMat[,i] ~ SDpoints + I(SDpoints^2) + I(SDpoints^3)))$coef[,1]
SDspline <- smooth.spline(SDpoints, SDcoefMat[,i]+rnorm(length(SDcoefMat[,1]),0,5), df = HETpolyD, cv = TRUE, all.knots=TRUE)$fit

sd.knot <- cbind(sd.knot, SDspline$knot)
sd.nk <- c(sd.nk, SDspline$nk)
sd.min <- c(sd.min, SDspline$min)
sd.range <- c(sd.range, SDspline$range)
sd.coef	<- cbind(sd.coef, SDspline$coef)
}






	
###### M/Z eval #######	

pMZ = MZvals
pRT = RTeval[,RTn+1]

npeaks  <- length(pMZ)
nstd <- length(combinedProbs[,1])

#compare against all standards - combinedProbs$mass

MZe <- log(dnorm(sapply(combinedProbs$mass, masserror, standard = pMZ) + MZcoefE, mean = 0, sd = 2), base = 2)
RTe <- log(dnorm(sapply(combinedProbs$RT, RTdiff, standard = pRT), mean = 0, sd = pRT*RT.SDcoefO), base = 2)
SIZe <- t(log(probMat %*% t(peaksizeMat), base = 2))
posL <- MZe + RTe + SIZe

#LIKeval <- cbind(matrix(MZe + RTe + SIZe, ncol = length(combinedProbs[,1]), nrow = length(peakpos), byrow = FALSE), -50)

peakLIK <- lapply(c(1:length(compounds)), GaussianLik, coToiso, posL, peaksizeMat, probMat, npeaks, nhet+1, SDlmMat, HETbase)

LIKform <- NULL
for (znum in 1:length(peakLIK)){
LIKform <- rbind(LIKform, data.frame(peakLIK[znum[1]]))	}

LIKmat <- matrix(0, ncol = nstd, nrow = npeaks)

for(z in 1:length(LIKform[,1])){
LIKmat[LIKform[z,1],LIKform[z,3]] <- LIKform[z,4]}

MZlikE <- sum(apply(LIKmat, 1, max))/npeaks

##### sd(RT) #######

MZe <- log(dnorm(sapply(combinedProbs$mass, masserror, standard = pMZ) + MZcoefO, mean = 0, sd = 2), base = 2)
RTe <- log(dnorm(sapply(combinedProbs$RT, RTdiff, standard = pRT), mean = 0, sd = pRT*RT.SDcoefE), base = 2)
SIZe <- t(log(probMat %*% t(peaksizeMat), base = 2))
posL <- MZe + RTe + SIZe

#LIKeval <- cbind(matrix(MZe + RTe + SIZe, ncol = length(combinedProbs[,1]), nrow = length(peakpos), byrow = FALSE), -50)

peakLIK <- lapply(c(1:length(compounds)), GaussianLik, coToiso, posL, peaksizeMat, probMat, npeaks, nhet+1, SDlmMat, HETbase)

LIKform <- NULL
for (znum in 1:length(peakLIK)){
LIKform <- rbind(LIKform, data.frame(peakLIK[znum[1]]))	}

LIKmat <- matrix(0, ncol = nstd, nrow = npeaks)

for(z in 1:length(LIKform[,1])){
LIKmat[LIKform[z,1],LIKform[z,3]] <- LIKform[z,4]}

RT.peakSDlikE <- sum(apply(LIKmat, 1, max))/npeaks


##### RT scaling ########

for(i in 1:RTn){
	
MZe <- log(dnorm(sapply(combinedProbs$mass, masserror, standard = MZvals) + MZcoefO, mean = 0, sd = 2), base = 2)
RTe <- log(dnorm(sapply(combinedProbs$RT, RTdiff, standard = RTeval[,i]), mean = 0, sd = pRT*RT.SDcoefO), base = 2)
SIZe <- t(log(probMat %*% t(peaksizeMat), base = 2))
posL <- MZe + RTe + SIZe

#LIKeval <- cbind(matrix(MZe + RTe + SIZe, ncol = length(combinedProbs[,1]), nrow = length(peakpos), byrow = FALSE), -50)

peakLIK <- lapply(c(1:length(compounds)), GaussianLik, coToiso, posL, peaksizeMat, probMat, npeaks, nhet+1, SDlmMat, HETbase)

LIKform <- NULL
for (znum in 1:length(peakLIK)){
LIKform <- rbind(LIKform, data.frame(peakLIK[znum[1]]))	}

if(dim(LIKform)[1] == 0){RTlikE[i] <- 0}else{

LIKmat <- matrix(0, ncol = nstd, nrow = npeaks)

for(z in 1:length(LIKform[,1])){
LIKmat[LIKform[z,1],LIKform[z,3]] <- LIKform[z,4]}

RTlikE[i] <- sum(apply(LIKmat, 1, max))/npeaks
}}

######### Heteroscedasticity determination #############



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











