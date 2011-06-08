#qsub -l 1day -cwd -sync n Rscript //Genomics/grid/users/shackett/ISOpeakFinder/ISOpeakFinder.R


options(digits = 13)

CETUSUSED <- FALSE

if(CETUSUSED == TRUE){setwd("/Genomics/grid/users/shackett/ISOpeakFinder/")}else{
setwd("/Users/seanhackett/Desktop/RabinowitzLab/IsoPeak/")}

load("knownMzRTre.R")
source("pfLibrary.R")


allPeaks <- read.table("4.23.11aligned.csv", sep = ",", header = TRUE)

#Constants
isoCov <- 0.7

samples <- c("N.1", "N.2", "P.1", "P.2")
sampleclass <- c("N", "N", "P", "P")
replicates <- c(1,1,2,2)
nsamples <- length(samples)
nuniqueSamp <- length(unique(replicates))

blanks <- c("blank.1", "blank.2", "blank.3")
compounds <- unique(MzRTrefine$compound)

### Compare Harmonic means of blanks and samples - baseline so approach is valid where samples have no value

allPeaks[,names(allPeaks) %in% c(samples, blanks)][allPeaks[,names(allPeaks) %in% c(samples, blanks)] < 300] <- 300

SampM = apply(allPeaks[,names(allPeaks) %in% samples], 1, harmM)
BLM = apply(allPeaks[,names(allPeaks) %in% blanks], 1, harmM)
allPeaks <- cbind(allPeaks, BLM)

valS <- allPeaks[SampM > 2*BLM,]

peaksizeMat <- valS[,colnames(allPeaks) %in% samples] - valS$BLM
peaksizeMat[peaksizeMat < 300] <- 300

#make n-isotopic labeled and non-labeled indexes equivalent

combinedProbs <- NULL

for(i in 1:length(compounds)){
sublist <- MzRTrefine[MzRTrefine$compound %in% compounds[i],]
nsublist <- nisoMzRTrefine[nisoMzRTrefine$compound %in% compounds[i],]

jointmass <- union(sublist$mass, nsublist$mass)
if(length(jointmass[!(jointmass %in% sublist$mass)]) != 0){
sublist  <- rbind(sublist, data.frame(compound = sublist$compound[1], mass = jointmass[!(jointmass %in% sublist$mass)], prob = 0, RT = sublist$RT[1])) 
sublist <- sublist[order(sublist$mass),]
	}
	
if(length(jointmass[!(jointmass %in% nsublist$mass)]) != 0){
nsublist  <- rbind(nsublist, data.frame(compound = nsublist$compound[1], mass = jointmass[!(jointmass %in% nsublist$mass)], prob = 0, RT = nsublist$RT[1])) 
nsublist <- nsublist[order(nsublist$mass),]
	}	

combinedProbs <- rbind(combinedProbs, data.frame(compound = sublist$compound, mass = sublist$mass, RT = sublist$RT, uLabp = sublist$prob, nLabp = nsublist$prob))}

#save(combinedProbs, file = "knownProbs.R")

probMat <- matrix(data =NA, ncol = length(sampleclass), nrow = length(combinedProbs[,1]))
for(i in 1:length(sampleclass)){
if(sampleclass[i] == "P"){probMat[,i] <- combinedProbs$uLabp}else{probMat[,i] <- combinedProbs$nLabp}
}	

coToiso <- matrix(data = NA, nrow = length(combinedProbs[,1]), ncol = length(compounds))
colnames(coToiso) <- compounds

for(i in 1:length(compounds)){
	coToiso[,i] <- combinedProbs[,1] %in% compounds[i]
	}
	




#### Simulated Annealing ####

######### Setup ####################

# Establish annealing schedule # must be a multiple of 20
nanneal <- 20
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

RTn = 20
RTpos <- range(RTvals)[1] + c(0:(RTn-1))*(diff(range(RTvals))/(RTn-1))
RTinc <- (diff(range(RTvals))/(RTn-1))/2

RTtrack <- matrix(NA, ncol = RTn, nrow = nanneal)
MZoffsetrack  <- rep(NA, times = nanneal)
RTliktrack <- matrix(NA, ncol = RTn, nrow = nanneal)
MZliktrack  <- rep(NA, times = nanneal)

RTcoefsO <- rep(1, times = RTn)
RTcoefsE <- rep(NA, times = RTn)
MZcoefO <- 0
MZcoefE <- NA

RTlik  <- rep(NA, times = RTn)
MZlik <- NA

RTlikE <- rep(NA, times = RTn)
MZlikE <- NA

HETbase <- 2.5
hetR <- c(floor(range(log(valS[,colnames(valS) %in% samples], base = 2))[1]), ceiling(range(log(valS[,colnames(valS) %in% samples], base = HETbase))[2]))

nhet <- length(hetR[1]:hetR[2])

SDliktrack <- matrix(NA, ncol = nhet, nrow = nanneal)
SDtrack <- matrix(NA, ncol = nhet, nrow = nanneal)
peakSDlik <- rep(NA, times = nhet)
peakSDlikE <- rep(NA, times = nhet) 
SDcoefO <- rep(1, times = nhet)
SDcoefE <- rep(NA, times = nhet)

save(valS, HETbase, RTpos, samples, nanneal, file = "PeakFparams.R")

for(j in 1:nanneal){

#ppm offset - systematic difference between peaks and standards
if(j == 1){MZcoefE <- 0; RTcoefsE <- rep(1, times = RTn); SDcoefE <- rep(1, times = nhet)
	}else{
MZcoefE <- rnorm(1, MZcoefO, 2*annealvar[j])
	
RTcoefsE <- (rnorm(n = length(RTcoefsO), mean = RTcoefsO, sd = 0.1*annealvar[j]) + rnorm(n = length(RTcoefsO), mean = 1, sd = 0.1*annealvar[j]))/2

SDcoefE <- SDcoefO*runif(nhet, 1-annealvar[j]*(2/5), 1+annealvar[j]*(2/5))
}
#C1-C20 are alternative H, C21 is the null for LRT
RTcoefsMat <- matrix(data = RTcoefsO, ncol = RTn+1, nrow = RTn)
diag(RTcoefsMat) <- RTcoefsE

RTeval <- matrix(data = NA, ncol = RTn+1, nrow = length(RTvals))

for(i in 1:(RTn+1)){

RTpoints <- RTpos*RTcoefsMat[,i]
RTcoefs <- summary(lm(RTpoints ~ RTpos + I(RTpos^2) + I(RTpos^3)))$coef[,1]
RTeval[,i] <- RTcoefs[1] + RTcoefs[2]*RTvals + RTcoefs[3]*RTvals^2 + RTcoefs[4]*RTvals^3
	
	}	

### NULL looks at points around RTpos for each comparison ###
### ALT looks at points around RTpos*diag(RTcoefsMat)

RTalt <- matrix(data = NA, ncol = RTn, nrow = length(RTvals))

for(i in 1:length(RTpos)){
	RTalt[,i] <- ifelse(abs(RTvals - RTpos[i]*diag(RTcoefsMat)[i]) < RTinc + 0.2, TRUE, FALSE)
	}
colnames(RTalt) <- RTpos			

### Determine SD of a peak given its intensity - Heteroscedasticity ###

SDcoefMat <- matrix(data = (hetR[1]:hetR[2])*SDcoefO, ncol = nhet+1, nrow = nhet)
diag(SDcoefMat) <- log(HETbase^(hetR[1]:hetR[2])*SDcoefE, base = HETbase)

SDlmMat <- matrix(data = NA, ncol = nhet+1, nrow = 4)

for(i in 1:(nhet+1)){

SDpoints <- c(hetR[1]:hetR[2])
SDlmMat[,i] <- summary(lm(SDcoefMat[,i] ~ SDpoints + I(SDpoints^2) + I(SDpoints^3)))$coef[,1]
	
	}	

#evalINT <- seq(hetR[1], hetR[2], by = 0.1)
#plot(VarlmMat[1,2] + VarlmMat[2,2]*evalINT + VarlmMat[3,2]*evalINT^2 + VarlmMat[4,2]*evalINT^3 ~ evalINT)
#points(evalINT+1 ~ evalINT, type = "l")




	
###### M/Z eval #######	

### all x all

pMZ = MZvals
pRT = RTeval[,RTn+1]

npeaks  <- length(pMZ)
nstd <- length(combinedProbs[,1])

#compare against all standards - combinedProbs$mass

MZe <- log(dnorm(sapply(combinedProbs$mass, masserror, standard = pMZ) + MZcoefE, mean = 0, sd = 1), base = 2)
RTe <- log(dnorm(sapply(combinedProbs$RT, RTdiff, standard = pRT), mean = 0, sd = 0.5), base = 2)
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

##### RT eval ########

### comparing fits around RT region x all peaks ###

for(i in 1:RTn){

peakpos <- c(1:length(RTalt[,1]))[RTalt[,i]]

if(length(peakpos) == 0){RTlikE[i] <- -50; next}


peaks = MZvals[peakpos]	
pRT = RTeval[peakpos,i]
pSize = peaksizeMat[peakpos,]

npeaks <- length(pRT)

		
#compare against all standards - combinedProbs$mass
	
MZe <- log(dnorm(sapply(combinedProbs$mass, masserror, standard = peaks) + MZcoefO, mean = 0, sd = 3), base = 2)
RTe <- log(dnorm(sapply(combinedProbs$RT, RTdiff, standard = pRT), mean = 0, sd = 0.5), base = 2)
SIZe <- t(log(probMat %*% t(pSize), base = 2))
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

MZe <- log(dnorm(sapply(combinedProbs$mass, masserror, standard = pMZ) + MZcoefO, mean = 0, sd = 3), base = 2)
RTe <- log(dnorm(sapply(combinedProbs$RT, RTdiff, standard = pRT), mean = 0, sd = 0.5), base = 2)
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
		
RTtrack[j,] <- RTcoefsO
MZoffsetrack[j]  <- MZcoefO
SDtrack[j,] <- SDcoefO
RTliktrack[j,] <- RTlikE
MZliktrack[j]  <- MZlikE
SDliktrack[j,] <- peakSDlikE

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

save(j, RTtrack, RTliktrack, MZoffsetrack, MZliktrack, SDtrack, SDliktrack, file = "SpectrumScale.R")





