options(digits = 13)

CETUSUSED <- FALSE

if(CETUSUSED == TRUE){setwd("/Genomics/grid/users/shackett/ISOpeakFinder/")}else{
setwd("/Users/seanhackett/Desktop/RabinowitzLab/IsoPeak/")}

source("pfLibrary.R")
load("Saved_Filez/knownMzRTre.R")
load("Saved_Filez/knownProbs.R")

#Constants
isoCov <- 0.7

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

allPeaks[,names(allPeaks) %in% c(samples, blanks)][allPeaks[,names(allPeaks) %in% c(samples, blanks)] < 300] <- 300

SampM = apply(allPeaks[,names(allPeaks) %in% samples], 1, harmM)
BLM = apply(allPeaks[,names(allPeaks) %in% blanks], 1, harmM)
allPeaks <- cbind(allPeaks, BLM)

valS <- allPeaks[SampM > 2*BLM,]

peaksizeMat <- valS[,colnames(allPeaks) %in% samples] - valS$BLM
peaksizeMat[peaksizeMat < 300] <- 300

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

peakLIK <- lapply(c(1:length(compounds)), GaussianLik, coToiso, posL, peaksizeMat, probMat, npeaks, nhet+1, SDlmMat, HETbase)

LIKform <- NULL
for (znum in 1:length(peakLIK)){
LIKform <- rbind(LIKform, data.frame(peakLIK[znum[1]]))	}

LIKmat <- matrix(0, ncol = nstd, nrow = npeaks)

for(z in 1:length(LIKform[,1])){
LIKmat[LIKform[z,1],LIKform[z,3]] <- LIKform[z,4]}

MZlikE <- sum(apply(LIKmat, 1, max))/npeaks



com <- 169

############# Determine likelihood of a set of peaks corresponding to a compouds isotopes given the current evidence from peak size and position (posL) and the observed ratios of isotopes ###############

GauS.w.Adduct <- function(com, coToiso, posL, peaksizeMat, probMat, npeaks, h, SDlmMat, HETbase, pMZ, pRT, combinedAdds, ADDUCT.OUT, RT.UNKNOWN){
#single peaks have the right ratio
#print(com)
subprob <- probMat[coToiso[,com],]
posLsub <- posL[,coToiso[,com]]

nfacs <- apply(posLsub, 2, sumthresh, thresh = -20, nvec = npeaks)
if(length(unlist(nfacs)) != 0){

stacker <- NULL
for(i in 1:length(nfacs)){
	#if(length(unlist(nfacs[i]) > 0)){stacker <- rbind(stacker, data.frame(peaks = unlist(nfacs[i]), iso = i))}
	if(length(unlist(nfacs[i]) > 0)){stacker <- rbind(stacker, data.frame(peaks = unlist(nfacs[i]), iso = i), data.frame(peaks = NA, iso = i))}else{stacker <- rbind(stacker, data.frame(peaks = NA, iso = i))		}}

factlevels = unstack(stacker)
levgrid <- expand.grid(factlevels)
definedP <- ifelse(is.na(levgrid), 0, 1)

Nf <- combinedProbs$uLabp[coToiso[,com]]
Pf <- combinedProbs$nLabp[coToiso[,com]]

#######

validP <- levgrid[apply(definedP, 1, validPerms, Nf, Pf, isoCov),]

if(length(validP[,1]) != 0){

colZ <- stacker[!is.na(stacker)[,1],]

probMatsub <- probMat[coToiso[,com],]



###### require that peaks in the same permutation have matching RT w/ sd = 0.2 min.

RT.weights <- apply(cbind(Nf, Pf), 1, mean)

colZ.RT <- pRT[colZ[,1]]

RTgrid <- validP

for(i in 1:length(colZ[,1])){
	RTgrid[RTgrid == colZ[i,1]] <- colZ.RT[i]
	}
RTgrid[is.na(RTgrid)] <- 0	
	

permRTs <- apply(RTgrid*matrix(RT.weights, ncol = length(validP[1,]), nrow = length(validP[,1]), byrow = TRUE), 1, sum, na.rm = TRUE)/(ifelse(is.na(validP), 0, 1)%*%RT.weights)

logL.RT <- log(gausD(peaks=RTgrid, expected = matrix(permRTs, ncol = length(validP[1,]), nrow = length(validP[,1])), sdP = matrix(0.1, ncol = length(validP[1,]), nrow = length(validP[,1]))), base = 2) * ifelse(is.na(validP), 0, 1)
logL.RT[is.nan(logL.RT)] <- NA
logL.RT.perm <- apply(logL.RT, 1, sum, na.rm = TRUE)

######

validRT <- validP[logL.RT.perm > 0,]

if(length(validRT[,1]) != 0){

RT.perm <- permRTs[logL.RT.perm > 0]

nuMiso <- length(levgrid[1,])
npeakISO <- length(colZ[,1]) 
nperms <- length(validRT[,1])

indies <- c(1:nsamples)

PMAT <- matrix(data = unlist(t(peaksizeMat[colZ[,1],])), ncol = nsamples*npeakISO, nrow = nperms, byrow = TRUE)
#colnames(PMAT) <- rep(colZ[,1], each = nsamples)
colnames(PMAT) <- rep(indies, times = npeakISO)

PPROB <- matrix(data = unlist(t(as.data.frame(probMatsub[colZ[,2],]))), ncol = nsamples*npeakISO, nrow = nperms, byrow = TRUE)
colnames(PPROB) <- rep(indies, times = npeakISO)

sdPMAT <- matrix(data = peakSD(unlist(t(peaksizeMat[colZ[,1],])), h = 1, SDlmMat, HETbase), ncol = nsamples*npeakISO, nrow = nperms, byrow = TRUE)
colnames(sdPMAT) <- rep(indies, times = npeakISO)

aMAT <- matrix(data = unlist(peaksizeMat[colZ[,1],])/peakSD(unlist(peaksizeMat[colZ[,1],]), 1, SDlmMat, HETbase), ncol = nsamples*npeakISO, nrow = nperms, byrow = TRUE)
colnames(sdPMAT) <- rep(indies, times = npeakISO)

posLinfo <- NULL
for(i in 1:length(colZ[,1])){
	posLinfo <- c(posLinfo, posLsub[colZ[i,1], colZ[i,2]])
}

POSLMAT <- matrix(rep(posLinfo, each = nsamples), ncol = nsamples*npeakISO, nrow = nperms, byrow = TRUE) 


#convert validRT to peak-sample x nperm format using colZ

isoRep <- table(colZ[,2])*nsamples
gridEXP <- validRT[,as.numeric(rep(names(isoRep), times = isoRep))]
gridCOM <- matrix(data = (-1*rep(colZ[,1], each = nsamples)), ncol = nsamples*npeakISO, nrow = nperms, byrow = TRUE)

USED <- ifelse(gridEXP+gridCOM == 0 & !is.na(gridEXP*gridCOM), 1, 0)

# Take the sum of P for each samples x perm

Psum <- USED*PPROB
PAsum <- USED*PPROB*aMAT
colnames(Psum) <- rep(indies, times = npeakISO)
colnames(PAsum) <- rep(indies, times = npeakISO)

MUcolz <- NULL
if(nperms == 1){

for(i in 1:nsamples){ 
DEN <- Psum[,colnames(Psum) == i]*PAsum[,colnames(PAsum) == i]
NUM <- PMAT[,colnames(PMAT) == i]*USED[,colnames(PMAT) == i]*aMAT[,colnames(PMAT) == i]*PPROB[,colnames(PMAT) == i]
MUcolz <- cbind(MUcolz, NUM/DEN)
}
MUcolz <- ifelse(is.nan(MUcolz), 0, MUcolz)
for(i in 1:nuniqueSamp){	
MUcolz[,c(1:nsamples)[replicates == i]] <- mean(MUcolz[,c(1:nsamples)[replicates == i]])
}}else{

for(i in 1:nsamples){ 
DEN <- apply(Psum[,colnames(Psum) == i], 1, sum)*apply(PAsum[,colnames(PAsum) == i], 1, sum)
NUM <- apply(PMAT[,colnames(PMAT) == i]*USED[,colnames(PMAT) == i]*aMAT[,colnames(PMAT) == i]*PPROB[,colnames(PMAT) == i], 1, sum)
MUcolz <- cbind(MUcolz, NUM/DEN)
}
MUcolz <- ifelse(is.nan(MUcolz), 0, MUcolz)
for(i in 1:nuniqueSamp){	
MUcolz[,c(1:nsamples)[replicates == i]] <- apply(MUcolz[,c(1:nsamples)[replicates == i]], 1, mean)
}}

##################

REPMOD <- USED*5
EVAL <- (log(gausD(PMAT, matrix(MUcolz, ncol = nsamples*npeakISO, nrow = nperms, byrow = FALSE)*USED*PPROB, sdPMAT), base = 10) + POSLMAT)*USED+ REPMOD


################# identify peaks that fit with both isotopes ################

sharedPeaks <- names((apply(table(stacker[!is.na(stacker[,1]),]), 1, sum) != 1)[(apply(table(stacker[!is.na(stacker[,1]),]), 1, sum) != 1) == TRUE])
if(length(sharedPeaks) > 0){
sharedPeakPairs <- NULL
sharedPnumVec <- NULL
sharedPeakPairsFULL <- NULL
for(i in 1:length(sharedPeaks)){
for(j in 2:sum((stacker$peaks[!is.na(stacker$peaks)] %in% sharedPeaks[i]) == TRUE)){
	
table(stacker)[row.names(table(stacker)) == sharedPeaks[i],]
	
combos <- combn(c(1:nuMiso)[table(stacker)[row.names(table(stacker)) == sharedPeaks[i],] == 1], j)
tMat <- matrix(data = 0, ncol = nuMiso, nrow = length(combos[1,]))
combosFULL <- combn(c(1:npeakISO)[stacker$peaks[!is.na(stacker$peaks)] %in% sharedPeaks[i]], j)
tMatFULL <- matrix(data = 0, ncol = npeakISO, nrow = length(combos[1,]))		
		
	
for(k in 1:length(combos[1,])){	
tMat[k,combos[,1]] <- (1/j)	
tMatFULL[k, combosFULL[,1]] <- 1}
}
sharedPeakPairs <- rbind(sharedPeakPairs, tMat)	
sharedPnumVec <- c(sharedPnumVec, rep(sharedPeaks[i], times = length(tMat[,1])))
sharedPeakPairsFULL <- rbind(sharedPeakPairsFULL, tMatFULL)
}

rTemp <- validRT
rTemp[is.na(rTemp)] <- 0

Tcombos <- (as.matrix(rTemp) %*% t(sharedPeakPairs)) == matrix(as.numeric(sharedPnumVec), nrow = length(rTemp[,1]), ncol = length(sharedPeakPairs[,1]), byrow = TRUE)

Tcombos <- ifelse(Tcombos, 1, 0)

for(k in 1:nsamples){
	
LIKadjust <- log(exp((Tcombos %*% sharedPeakPairsFULL)*EVAL[,colnames(EVAL) == k]) / (Tcombos*(exp(EVAL[,colnames(EVAL) == k]) %*% t(sharedPeakPairsFULL)))%*% sharedPeakPairsFULL)

EVAL[,colnames(EVAL) == k] <- EVAL[,colnames(EVAL) == k] + ifelse(abs(LIKadjust) == Inf, 1, LIKadjust)

	}}

#######################################

#EVALsum[order(EVALsum, decreasing = TRUE)  == 1]
 
EVALsum <- apply(EVAL, 1, sum)
peakeval <- matrix(EVAL[order(EVALsum, decreasing = TRUE)[1:min(10, nperms)],], nrow = nsamples*npeakISO, byrow = TRUE)
rownames(peakeval) <- rep(c(1:npeakISO), each = length(indies))

RT.perm.eval <- RT.perm[order(EVALsum, decreasing = TRUE)[1:min(10, nperms)]]
MU.perm.eval <- MUcolz[order(EVALsum, decreasing = TRUE)[1:min(10, nperms)],]


outz <- matrix(sapply(c(1:npeakISO), factcond, peakeval), ncol = npeakISO)
outz <- rbind(outz, rep(0, times = length(outz[1,])))


output <- data.frame(colZ, standard = c(1:length(coToiso[,1]))[coToiso[,com]][colZ[,2]], value = apply(outz, 2, max))
output <- output[output$value != 0,]

if(length(output[,1]) != 0){
	
if(ADDUCT.OUT == FALSE){output}else{


############ Look for adducts of each peak in stacker

parentpeaks <- unique(stacker$peaks[!is.na(stacker$peaks)][(apply(outz, 2, sum) != 0)])
STD <- combinedProbs[combinedProbs$compound %in% compounds[com],][unique(stacker$iso[!is.na(stacker$peaks)][(apply(outz, 2, sum) != 0)]),]


transM <- MZtransform(combinedAdds, STD)

addMZmat <- (matrix(STD$mass, ncol = length(STD$mass), nrow = length(transM[,1]), byrow = TRUE) + matrix(transM$add, ncol = length(STD$mass), nrow = length(transM[,1])))*matrix(transM$scale, ncol = length(STD$mass), nrow = length(transM[,1]))

adductL <- length(unlist(t(addMZmat)))

#points adduct type to corresponding columns of adduct isotopic variants * parent peak

adduct.name <- unique(transM[,1])
adduct.to.pbya <- matrix(data = NA, ncol = length(adduct.name), nrow = adductL)

for(i in 1:length(adduct.name)){

adduct.to.pbya[,i] <- rep(transM[,1] == adduct.name[i], times = length(STD[,1]))

}

addMZmat <- matrix(unlist(t(matrix(addMZmat, ncol = length(STD$mass), nrow = length(transM[,1]), byrow = TRUE))), ncol = 4, nrow = adductL)

addPmat <- matrix(sampleclass, ncol = nsamples, nrow = adductL, byrow = TRUE)

addPmat <- ifelse(addPmat == "N", 1, 0)*matrix(unlist((matrix(STD$nLabp, ncol = length(STD$mass), nrow = length(transM[,1]), byrow = TRUE)*matrix(transM$nLab, ncol = length(STD$mass), nrow = length(transM[,1])))), ncol = 4, nrow = adductL) + ifelse(addPmat == "N", 0, 1)*matrix(unlist((matrix(STD$uLabp, ncol = length(STD$mass), nrow = length(transM[,1]), byrow = TRUE)*matrix(transM$unLab, ncol = length(STD$mass), nrow = length(transM[,1])))), ncol = 4, nrow = adductL) 

#matrix(unlist((matrix(STD$nLabp, ncol = length(STD$mass), nrow = length(transM[,1]), byrow = TRUE) + matrix(transM$nLab, ncol = length(STD$mass), nrow = length(transM[,1]))))

#matrix(paste(t(matrix(STD$mass, ncol = length(STD$mass), nrow = length(transM[,1]), byrow = TRUE)), t(matrix(transM$adduct, ncol = length(STD$mass), nrow = length(transM[,1])))), ncol = 4, nrow = adductL)




#PeakAddMZmat <- matrix(data = rep(pMZ, each = nsamples), ncol = nsamples*length(pMZ), nrow = length(addPmat[,1]), byrow = TRUE)

#STDaddMZmat <- matrix(data = unlist(t(matrix(addMZmat, ncol = length(STD$mass), nrow = length(transM[,1]), byrow = TRUE))), ncol = nsamples*length(pMZ), nrow = length(transM[,1]))



addMZeval <- t(log(dnorm(sapply(addMZmat[,1], masserror, standard = pMZ) + MZoffsetrack[nanneal], mean = 0, sd = 1), base = 2))

addRTeval <- log(dnorm(sapply(RT.perm.eval, RTdiff, standard = pRT), mean = 0, sd = 0.1), base = 2)

addSIZEeval <- log(addPmat %*% t(peaksizeMat), base = 2)

overall.add.output <- NULL

for(k in 1:length(MU.perm.eval[,1])){
	
addLik <- addMZeval + addSIZEeval + matrix(addRTeval[,k], ncol = length(peaksizeMat[,1]), nrow = length(addPmat[,1]), byrow = TRUE)
add.output <- NULL

for(add in 1:length(adduct.name)){
	
ind.addL <- addLik[adduct.to.pbya[,add],]
add.nfacs <- apply(ind.addL, 1, sumthresh, thresh = -20, nvec = npeaks)
	
if(length(unlist(add.nfacs)) != 0){
	
#print(add)}}	
	
stacker <- NULL
for(i in 1:length(add.nfacs)){
	if(length(unlist(add.nfacs[i]) > 0)){stacker <- rbind(stacker, data.frame(peaks = unlist(add.nfacs[i]), iso = i), data.frame(peaks = NA, iso = i))}else{stacker <- rbind(stacker, data.frame(peaks = NA, iso = i))	}}

factlevels = unstack(stacker)
levgrid <- expand.grid(factlevels)
definedP <- ifelse(is.na(levgrid), 0, 1)

validP <- levgrid[definedP%*%apply(addPmat[adduct.to.pbya[,add],], 1, max) > isoCov,]





if(length(validP[,1]) != 0){

colZ <- stacker[!is.na(stacker)[,1],]

add.probMatsub <- addPmat[adduct.to.pbya[,add],]
nuMiso <- length(levgrid[1,])
npeakISO <- length(colZ[,1]) 
nperms <- length(validP[,1])
posLsub <- t(addLik[adduct.to.pbya[,add],])
			
PMAT <- matrix(data = unlist(t(peaksizeMat[colZ[,1],])), ncol = nsamples*npeakISO, nrow = nperms, byrow = TRUE)
#colnames(PMAT) <- rep(colZ[,1], each = nsamples)
colnames(PMAT) <- rep(indies, times = npeakISO)

PPROB <- matrix(data = unlist(t(as.data.frame(add.probMatsub[colZ[,2],]))), ncol = nsamples*npeakISO, nrow = nperms, byrow = TRUE)
colnames(PPROB) <- rep(indies, times = npeakISO)

sdPMAT <- matrix(data = peakSD(unlist(t(peaksizeMat[colZ[,1],])), h = 1, SDlmMat, HETbase), ncol = nsamples*npeakISO, nrow = nperms, byrow = TRUE)
colnames(sdPMAT) <- rep(indies, times = npeakISO)

aMAT <- matrix(data = unlist(peaksizeMat[colZ[,1],])/peakSD(unlist(peaksizeMat[colZ[,1],]), 1, SDlmMat, HETbase), ncol = nsamples*npeakISO, nrow = nperms, byrow = TRUE)
colnames(sdPMAT) <- rep(indies, times = npeakISO)

posLinfo <- NULL
for(i in 1:length(colZ[,1])){
	posLinfo <- c(posLinfo, posLsub[colZ[i,1], colZ[i,2]])
}

POSLMAT <- matrix(rep(posLinfo, each = nsamples), ncol = nsamples*npeakISO, nrow = nperms, byrow = TRUE) 


#convert validRT to peak-sample x nperm format using colZ

isoRep <- table(colZ[,2])*nsamples
gridEXP <- validP[,as.numeric(rep(names(isoRep), times = isoRep))]
gridCOM <- matrix(data = (-1*rep(colZ[,1], each = nsamples)), ncol = nsamples*npeakISO, nrow = nperms, byrow = TRUE)

USED <- ifelse(gridEXP+gridCOM == 0 & !is.na(gridEXP*gridCOM), 1, 0)

# Take the sum of P for each samples x perm

Psum <- USED*PPROB
PAsum <- USED*PPROB*aMAT
colnames(Psum) <- rep(indies, times = npeakISO)
colnames(PAsum) <- rep(indies, times = npeakISO)

add.MUcolz <- NULL
if(nperms == 1){

for(i in 1:nsamples){ 
DEN <- Psum[,colnames(Psum) == i]*PAsum[,colnames(PAsum) == i]
NUM <- PMAT[,colnames(PMAT) == i]*USED[,colnames(PMAT) == i]*aMAT[,colnames(PMAT) == i]*PPROB[,colnames(PMAT) == i]
add.MUcolz <- cbind(add.MUcolz, NUM/DEN)
}
add.MUcolz <- ifelse(is.nan(add.MUcolz), 0, add.MUcolz)
for(i in 1:nuniqueSamp){	
add.MUcolz[,c(1:nsamples)[replicates == i]] <- mean(add.MUcolz[,c(1:nsamples)[replicates == i]])
}}else{

for(i in 1:nsamples){ 
DEN <- apply(Psum[,colnames(Psum) == i], 1, sum)*apply(PAsum[,colnames(PAsum) == i], 1, sum)
NUM <- apply(PMAT[,colnames(PMAT) == i]*USED[,colnames(PMAT) == i]*aMAT[,colnames(PMAT) == i]*PPROB[,colnames(PMAT) == i], 1, sum)
add.MUcolz <- cbind(add.MUcolz, NUM/DEN)
}
add.MUcolz <- ifelse(is.nan(add.MUcolz), 0, add.MUcolz)
for(i in 1:nuniqueSamp){	
add.MUcolz[,c(1:nsamples)[replicates == i]] <- apply(add.MUcolz[,c(1:nsamples)[replicates == i]], 1, mean)
}}

add.fract = apply(add.MUcolz / matrix(MU.perm.eval[k,], nrow = nperms, ncol = nsamples, byrow = TRUE), 1, mean)


REPMOD <- USED*5
EVAL <- (log(gausD(PMAT, matrix(rep(MU.perm.eval[k,], times = npeakISO*nperms)*rep(add.fract, each = nsamples*npeakISO), ncol = nsamples*npeakISO, nrow = nperms, byrow = TRUE)*USED*PPROB, sdPMAT), base = 10) + POSLMAT)*USED+ REPMOD


EVALsum <- apply(EVAL, 1, sum)
peakeval <- matrix(EVAL[order(EVALsum, decreasing = TRUE)[1:min(10, nperms)],], nrow = nsamples*npeakISO, byrow = TRUE)
rownames(peakeval) <- rep(c(1:npeakISO), each = length(indies))

outz <- matrix(sapply(c(1:npeakISO), factcond, peakeval), ncol = npeakISO)
outz <- rbind(outz, rep(0, times = length(outz[1,])))


output <- data.frame(parentpk = k, addnum = add, adduct = adduct.name[add], colZ, standard = c(1:length(coToiso[,1]))[coToiso[,com]][colZ[,2]], value = apply(outz, 2, max))
output <- output[output$value != 0,]
add.output <- rbind(add.output, output)

}}}

overall.add.output <- rbind(overall.add.output, add.output)
}





}}}}}}
	
#for expected i isotopes w/ known MZ and l adducts, generate a m*l matrix of expected M/Z including the global offset 
		
MZtransform <- function(combinedAdds, STD){
STDmz <- STD$mass
abundSTDmz <- STD[(STD$uLabp + STD$nLabp) > 0.1,]

transM <- NULL
for(i in 1:length(combinedAdds[,1])){
	if(combinedAdds$nmol[i] == 1){
		transM <- rbind(transM, data.frame(adduct = combinedAdds$adduct[i], add = combinedAdds$weightch[i], scale = abs(1/combinedAdds$charge[i]), unLab = combinedAdds$uLabp[i], nLab = combinedAdds$nLabp[i]))
	}else{
		transM <- rbind(transM, data.frame(adduct = combinedAdds$adduct[i], add = combinedAdds$weightch[i] + combinedAdds$nmol[i]*abundSTDmz$mass, scale = rep(abs(1/combinedAdds$charge[i])), unLab = combinedAdds$uLabp[i]*abundSTDmz$uLabp^(combinedAdds$nmol[i]), nLab = combinedAdds$nLabp[i]*abundSTDmz$nLabp^(combinedAdds$nmol[i])))
		}}
transM}





