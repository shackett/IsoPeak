############# Determine likelihood of a set of peaks corresponding to a compouds isotopes given the current evidence from peak size and position (posL) and the observed ratios of isotopes ###############

GauS.w.Adduct <- function(com, coToiso, posL, peaksizeMat, probMat, npeaks, h, HETbase, pMZ, pRT, combinedProbs, combinedAdds, ADDUCT.USE, ADDUCT.OUT, RT.UNKNOWN, isoCov, RT.coel.SD, Supthresh, SDlmMat, MZoffset, MZ.SD, HETpolyD){

#Description of constants

#Supthresh - The support threshold for LIKe above which a possible association between a peak and standard will be evaluated.
#default = -20 

#isoCov - Only those peak-STD permutations which have peaks corresponding to sum(p) >isoCov of the expected compound for any labeling state will be accepted as validP permutations.  For adducts the sum of expected isotopic variant abundances across the labeling conditions must be greater than isoCov.
#default = 0.70

#RT.coel.SD - The support for peak-STD permutations - all peaks within a permutation should coelute and have a similar RT.  LIKe is modified by this factor and compared to Supthresh
#default = 0.1

#Signal.Thresh - Minimum detectable signal - Taken as a function of Threshold - This variable only acts to determine whether combinations of standards should be penalized for the absence of an expected signal above the Signal.Thresh.
#default = Threshold*2


print(com)

#subset of attributes for a single compound's isotopic variants
subprob <- probMat[coToiso[,com],]
posLsub <- posL[,coToiso[,com]]

#peaks which are associated with 1+ standards
nfacs <- apply(posLsub, 2, sumthresh, thresh = Supthresh, nvec = npeaks)
if(length(unlist(nfacs)) != 0){

stacker <- NULL
for(i in 1:length(nfacs)){
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

###### require that peaks in the same permutation have matching RT w/ sd = 0.1 min.

RT.weights <- apply(cbind(Nf, Pf), 1, mean)

colZ.RT <- pRT[colZ[,1]]

RTgrid <- validP

for(i in 1:length(colZ[,1])){
	RTgrid[RTgrid == colZ[i,1]] <- colZ.RT[i]
	}
RTgrid[is.na(RTgrid)] <- 0	
	

permRTs <- apply(RTgrid*matrix(RT.weights, ncol = length(validP[1,]), nrow = length(validP[,1]), byrow = TRUE), 1, sum, na.rm = TRUE)/(ifelse(is.na(validP), 0, 1)%*%RT.weights)

logL.RT <- log(gausD(peaks=RTgrid, expected = matrix(permRTs, ncol = length(validP[1,]), nrow = length(validP[,1])), sdP = matrix(RT.coel.SD, ncol = length(validP[1,]), nrow = length(validP[,1]))), base = 2) * ifelse(is.na(validP), 0, 1)
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

sdPMAT = matrix(unlist(t(peakSD(peaksizeMat[colZ[,1],], h, SDlmMat, HETbase, HETpolyD))), ncol = nsamples*npeakISO, nrow = nperms, byrow = TRUE)
colnames(sdPMAT) <- rep(indies, times = npeakISO)

aMAT <- PMAT/sdPMAT



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
DEN <- sum(Psum[,colnames(Psum) == i])*sum(PAsum[,colnames(PAsum) == i])
NUM <- sum(PMAT[,colnames(PMAT) == i]*USED[,colnames(PMAT) == i]*aMAT[,colnames(PMAT) == i]*PPROB[,colnames(PMAT) == i])
MUcolz <- c(MUcolz, NUM/DEN)
}
MUcolz <- ifelse(is.nan(MUcolz), 0, MUcolz)
for(i in 1:nuniqueSamp){	
MUcolz[c(1:nsamples)[replicates == i]] <- mean(MUcolz[c(1:nsamples)[replicates == i]])
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

absent.used.setup <- ifelse(is.na(validRT), 1, 0)

absent.used <- matrix(data = 0, ncol = nuMiso*nsamples, nrow = nperms)
colnames(absent.used) <- rep(c(1:nuMiso), each = nsamples)

for(i in 1:nuMiso){absent.used[,colnames(absent.used) == i] <- absent.used.setup[,i]}

absent.peaksize = matrix(data = Threshold, ncol = nuMiso*nsamples, nrow = nperms)

absent.sd = sdPMAT = matrix(peakSD(matrix(Threshold), h, SDlmMat, HETbase, HETpolyD), ncol = nuMiso*nsamples, nrow = nperms)

absent.prob = matrix(unlist(t(subprob)), ncol = nuMiso*nsamples, nrow = nperms, byrow = TRUE)

#evaluate how unlikely absent peaks are to be missing, and then round peaks with a low penalty to 0 (so as to not penalize absent peaks where the expected abundance is minute), set threshold where penalty applies to MUcolz*absent.prob > 2*threshold

absent.EVAL = (log(gausD(absent.peaksize, matrix(MUcolz, ncol = nuMiso*nsamples, nrow = nperms, byrow = FALSE)*absent.used*absent.prob, absent.sd), base = 10))*absent.used*ifelse(matrix(MUcolz, ncol = nuMiso*nsamples, nrow = nperms, byrow = FALSE)*absent.prob < Signal.Thresh, 0, 1)

EVAL <- cbind(EVAL, absent.EVAL)

################# identify peaks that fit with both isotopes ################

#sharedPeaks <- names((apply(table(stacker[!is.na(stacker[,1]),]), 1, sum) != 1)[(apply(table(stacker[!is.na(stacker[,1]),]), 1, sum) != 1) == TRUE])
#if(length(sharedPeaks) > 0){
#sharedPeakPairs <- NULL
#sharedPnumVec <- NULL
#sharedPeakPairsFULL <- NULL
#for(i in 1:length(sharedPeaks)){
#for(j in 2:sum((stacker$peaks[!is.na(stacker$peaks)] %in% sharedPeaks[i]) == TRUE)){
	
#table(stacker)[row.names(table(stacker)) == sharedPeaks[i],]
	
#combos <- combn(c(1:nuMiso)[table(stacker)[row.names(table(stacker)) == sharedPeaks[i],] == 1], j)
#tMat <- matrix(data = 0, ncol = nuMiso, nrow = length(combos[1,]))
#combosFULL <- combn(c(1:npeakISO)[stacker$peaks[!is.na(stacker$peaks)] %in% sharedPeaks[i]], j)
#tMatFULL <- matrix(data = 0, ncol = npeakISO, nrow = length(combos[1,]))		
		
	
#for(k in 1:length(combos[1,])){	
#tMat[k,combos[,1]] <- (1/j)	
#tMatFULL[k, combosFULL[,1]] <- 1}
#}
#sharedPeakPairs <- rbind(sharedPeakPairs, tMat)	
#sharedPnumVec <- c(sharedPnumVec, rep(sharedPeaks[i], times = length(tMat[,1])))
#sharedPeakPairsFULL <- rbind(sharedPeakPairsFULL, tMatFULL)
#}

#rTemp <- validRT
#rTemp[is.na(rTemp)] <- 0

#Tcombos <- (as.matrix(rTemp) %*% t(sharedPeakPairs)) == matrix(as.numeric(sharedPnumVec), nrow = length(rTemp[,1]), ncol = length(sharedPeakPairs[,1]), byrow = TRUE)

#Tcombos <- ifelse(Tcombos, 1, 0)

#for(k in 1:nsamples){
	
#LIKadjust <- log(exp((Tcombos %*% sharedPeakPairsFULL)*EVAL[,colnames(EVAL) == k]) / (Tcombos*(exp(EVAL[,colnames(EVAL) == k]) %*% t(sharedPeakPairsFULL)))%*% sharedPeakPairsFULL)

#EVAL[,colnames(EVAL) == k] <- EVAL[,colnames(EVAL) == k] + ifelse(abs(LIKadjust) == Inf, 1, LIKadjust)

#	}}

#######################################

#EVALsum[order(EVALsum, decreasing = TRUE)  == 1]
 
EVALsum <- apply(EVAL, 1, sum)
posEVAL <- length(EVALsum[EVALsum > 0])

peakeval <- matrix(EVAL[order(EVALsum, decreasing = TRUE)[1:min(10, posEVAL, nperms)],], nrow = nsamples*npeakISO + nuMiso*nsamples, byrow = TRUE)
#remove peaks with negative support
rownames(peakeval) <- c(rep(c(1:npeakISO), each = length(indies)), rep("absent", times = nsamples*nuMiso))

#penalize peak-std matches fractionaly by p/sum(p)

if(posEVAL !=0){
if(min(10, posEVAL, nperms) == 1){
	
eval.correction <- matrix(unlist(t(probMatsub[colZ[,2],])), ncol = min(10, posEVAL, nperms), nrow = nsamples*npeakISO)*(USED[order(EVALsum, decreasing = TRUE)[1:min(10, posEVAL, nperms)],])
eval.correction <- eval.correction*rep(sum(peakeval[rownames(peakeval) == "absent",]), each = nsamples*npeakISO)/rep(apply(eval.correction, 2, sum), each = nsamples*npeakISO)
	
	}else{

eval.correction <- matrix(unlist(t(probMatsub[colZ[,2],])), ncol = min(10, posEVAL, nperms), nrow = nsamples*npeakISO)*t(USED[order(EVALsum, decreasing = TRUE)[1:min(10, posEVAL, nperms)],])
eval.correction <- eval.correction*rep(apply(peakeval[rownames(peakeval) == "absent",], 2, sum), each = nsamples*npeakISO)/rep(apply(eval.correction, 2, sum), each = nsamples*npeakISO)
	}
	
peakeval <- peakeval[rownames(peakeval) != "absent",] + eval.correction
rownames(peakeval) <- c(rep(c(1:npeakISO), each = length(indies)))

RT.perm.eval <- RT.perm[order(EVALsum, decreasing = TRUE)[1:min(10, nperms, posEVAL)]]

if(is.matrix(MUcolz) == FALSE){MU.perm.eval <- t(matrix(MUcolz))}else{
MU.perm.eval <- MUcolz[order(EVALsum, decreasing = TRUE)[1:min(10, nperms, posEVAL)],]}

par.outz <- matrix(sapply(c(1:npeakISO), factcond, peakeval), ncol = npeakISO)
par.outz <- rbind(par.outz, rep(0, times = length(par.outz[1,])))

par.output <- data.frame(colZ, standard = c(1:length(coToiso[,1]))[coToiso[,com]][colZ[,2]], value = apply(par.outz, 2, max))
par.output <- par.output[par.output$value != 0,]

}else{par.output <- NULL}

if(length(par.output[,1]) != 0 & !is.null(par.output)){
	
if(ADDUCT.USE == FALSE){output}else{


############ Look for adducts of each peak in stacker ########

STD <- combinedProbs[unique(par.output$standard),]
#STD <- combinedProbs[combinedProbs$compound %in% compounds[com],][unique(stacker$iso[!is.na(stacker$peaks)][(apply(par.outz, 2, sum) != 0)]),]

transM <- MZtransform(combinedAdds, STD, sub.combinedProbs = combinedProbs[coToiso[,com],])

addMZmat <- (matrix(STD$mass, ncol = length(STD$mass), nrow = length(transM[,1]), byrow = TRUE) + matrix(transM$add, ncol = length(STD$mass), nrow = length(transM[,1])))*matrix(transM$scale, ncol = length(STD$mass), nrow = length(transM[,1]))

peakAssoc <- rep(rownames(STD), each = length(transM[,1]))

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


# MZ, RT diff for the adducts of the peaks observed in the initial survey

addMZeval <- t(log(dnorm(sapply(addMZmat[,1], masserror, standard = pMZ) + MZoffset, mean = 0, sd = MZ.SD), base = 2))

addRTeval <- log(dnorm(sapply(RT.perm.eval, RTdiff, standard = pRT), mean = 0, sd = RT.coel.SD), base = 2)

addSIZEeval <- log(addPmat %*% t(peaksizeMat), base = 2)

overall.add.output <- NULL

par.peak.perms <- ifelse(is.matrix(MU.perm.eval) == TRUE, length(MU.perm.eval[,1]), 1)



for(k in 1:par.peak.perms){
	
addLik <- addMZeval + addSIZEeval + matrix(addRTeval[,k], ncol = length(peaksizeMat[,1]), nrow = length(addPmat[,1]), byrow = TRUE)
add.output <- NULL

for(add in 1:length(adduct.name)){

ind.addL <- addLik[adduct.to.pbya[,add],]
if(is.matrix(ind.addL) == TRUE){
add.nfacs <- apply(ind.addL, 1, sumthresh, thresh = Supthresh, nvec = npeaks)}else{add.nfacs <- sumthresh(ind.addL, thresh = Supthresh, nvec = npeaks)}
	
if(length(unlist(add.nfacs)) != 0){
	
#print(add)}}
	
stacker <- NULL
for(i in 1:length(add.nfacs)){
	if(length(unlist(add.nfacs[i]) > 0)){stacker <- rbind(stacker, data.frame(peaks = unlist(add.nfacs[i]), iso = i), data.frame(peaks = NA, iso = i))}else{stacker <- rbind(stacker, data.frame(peaks = NA, iso = i))	}}

factlevels = unstack(stacker)
levgrid <- expand.grid(factlevels)
definedP <- ifelse(is.na(levgrid), 0, 1)

validP <- as.matrix(levgrid[definedP%*%apply(matrix(addPmat[adduct.to.pbya[,add],], ncol = nsamples),1, max) > isoCov,], ncol = length(definedP[1,]))

if(length(validP[,1]) != 0){

add.colZ <- stacker[!is.na(stacker)[,1],]

add.probMatsub <- matrix(addPmat[adduct.to.pbya[,add],], ncol = nsamples)
nuMiso <- length(levgrid[1,])
npeakISO <- length(add.colZ[,1])
nperms <- ifelse(is.matrix(validP), length(validP[,1]), 1)
posLsub <- matrix(t(addLik[adduct.to.pbya[,add],]), ncol = nuMiso)
			
PMAT <- matrix(data = unlist(t(peaksizeMat[add.colZ[,1],])), ncol = nsamples*npeakISO, nrow = nperms, byrow = TRUE)
colnames(PMAT) <- rep(indies, times = npeakISO)

PPROB <- matrix(data = unlist(t(add.probMatsub[add.colZ[,2],])), ncol = nsamples*npeakISO, nrow = nperms, byrow = TRUE)
colnames(PPROB) <- rep(indies, times = npeakISO)

sdPMAT = matrix(unlist(t(peakSD(peaksizeMat[add.colZ[,1],], h, SDlmMat, HETbase, HETpolyD))), ncol = nsamples*npeakISO, nrow = nperms, byrow = TRUE)
colnames(sdPMAT) <- rep(indies, times = npeakISO)

aMAT <- PMAT/sdPMAT


posLinfo <- NULL
for(i in 1:length(add.colZ[,1])){
	posLinfo <- c(posLinfo, posLsub[add.colZ[i,1], add.colZ[i,2]])
}

POSLMAT <- matrix(rep(posLinfo, each = nsamples), ncol = nsamples*npeakISO, nrow = nperms, byrow = TRUE) 


#convert validRT to peak-sample x nperm format using add.colZ

isoRep <- table(add.colZ[,2])*nsamples
gridEXP <- validP[,as.numeric(rep(names(isoRep), times = isoRep))]
gridCOM <- matrix(data = (-1*rep(add.colZ[,1], each = nsamples)), ncol = nsamples*npeakISO, nrow = nperms, byrow = TRUE)

USED <- ifelse(gridEXP+gridCOM == 0 & !is.na(gridEXP*gridCOM), 1, 0)

# Take the sum of P for each samples x perm

Psum <- USED*PPROB
PAsum <- USED*PPROB*aMAT
colnames(Psum) <- rep(indies, times = npeakISO)
colnames(PAsum) <- rep(indies, times = npeakISO)

#############

add.MUcolz <- NULL
if(nperms == 1){

for(i in 1:nsamples){ 
DEN <- sum(Psum[,colnames(Psum) == i])*sum(PAsum[,colnames(PAsum) == i])
NUM <- sum(PMAT[,colnames(PMAT) == i]*USED[,colnames(PMAT) == i]*aMAT[,colnames(PMAT) == i]*PPROB[,colnames(PMAT) == i])
add.MUcolz <- c(add.MUcolz, NUM/DEN)
}
add.MUcolz <- ifelse(is.nan(add.MUcolz), 0, add.MUcolz)
for(i in 1:nuniqueSamp){	
add.MUcolz[c(1:nsamples)[replicates == i]] <- mean(add.MUcolz[c(1:nsamples)[replicates == i]])
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

if(par.peak.perms != 1){MU.perm <- MU.perm.eval[k,]}else{MU.perm <- MU.perm.eval}


add.fract = add.MUcolz / matrix(MU.perm, nrow = nperms, ncol = nsamples, byrow = TRUE)
add.fract[is.nan(add.fract)] <- NA
add.fract = apply(add.fract, 1, mean, na.rm = TRUE)

REPMOD <- USED*5
EVAL <- (log(gausD(PMAT, matrix(rep(MU.perm, times = npeakISO*nperms)*rep(add.fract, each = nsamples*npeakISO), ncol = nsamples*npeakISO, nrow = nperms, byrow = TRUE)*USED*PPROB, sdPMAT), base = 10) + POSLMAT)*USED+ REPMOD

### Evaluate penalties for absent - large signal adducts

absent.used.setup <- ifelse(is.na(validP), 1, 0)

absent.used <- matrix(data = 0, ncol = nuMiso*nsamples, nrow = nperms)
colnames(absent.used) <- rep(c(1:nuMiso), each = nsamples)

for(i in 1:nuMiso){absent.used[,colnames(absent.used) == i] <- absent.used.setup[,i]}

absent.peaksize = matrix(data = Threshold, ncol = nuMiso*nsamples, nrow = nperms)

absent.sd = matrix(data = peakSD(matrix(Threshold), h, SDlmMat, HETbase, HETpolyD), ncol = nuMiso*nsamples, nrow = nperms)

absent.prob = matrix(unlist(t(add.probMatsub)), ncol = nuMiso*nsamples, nrow = nperms, byrow = TRUE)

#evaluate how unlikely absent peaks are to be missing, and then round peaks with a low penalty to 0 (so as to not penalize absent peaks where the expected abundance is minute), set threshold where penalty applies to MUcolz*absent.prob > 2*threshold

absent.EVAL = (log(gausD(absent.peaksize, matrix(rep(MU.perm, times = nuMiso*nperms)*rep(add.fract, each = nsamples*nuMiso), ncol = nsamples*nuMiso, nrow = nperms, byrow = TRUE)*absent.used*absent.prob, absent.sd), base = 10))*absent.used*ifelse(matrix(rep(MU.perm, times = nuMiso*nperms)*rep(add.fract, each = nsamples*nuMiso), ncol = nsamples*nuMiso, nrow = nperms, byrow = TRUE)*absent.prob < Signal.Thresh, 0, 1)

EVAL <- cbind(EVAL, absent.EVAL)

EVALsum <- apply(EVAL, 1, sum)
EVALsum[is.nan(EVALsum)] <- 0
posEVAL <- length(EVALsum[EVALsum > 0])

peakeval <- matrix(EVAL[order(EVALsum, decreasing = TRUE)[1:min(10, nperms, posEVAL)],], nrow = nsamples*npeakISO + nuMiso*nsamples, byrow = TRUE)
rownames(peakeval) <- c(rep(c(1:npeakISO), each = length(indies)), rep("absent", times = nsamples*nuMiso))
#####


EVALsum <- apply(EVAL, 1, sum)
posEVAL <- length(EVALsum[EVALsum > 0])

peakeval <- matrix(EVAL[order(EVALsum, decreasing = TRUE)[1:min(10, posEVAL, nperms)],], nrow = nsamples*npeakISO + nuMiso*nsamples, byrow = TRUE)
#remove peaks with negative support
rownames(peakeval) <- c(rep(c(1:npeakISO), each = length(indies)), rep("absent", times = nsamples*nuMiso))

#penalize peak-std matches fractionaly by p/sum(p)

if(posEVAL !=0){
if(min(10, posEVAL, nperms) == 1){
	
eval.correction <- matrix(unlist(t(add.probMatsub[add.colZ[,2],])), ncol = min(10, posEVAL, nperms), nrow = nsamples*npeakISO)*(USED[order(EVALsum, decreasing = TRUE)[1:min(10, posEVAL, nperms)],])
eval.correction <- eval.correction*rep(sum(peakeval[rownames(peakeval) == "absent",]), each = nsamples*npeakISO)/rep(apply(eval.correction, 2, sum), each = nsamples*npeakISO)
	
	}else{

eval.correction <- matrix(unlist(t(add.probMatsub[add.colZ[,2],])), ncol = min(10, posEVAL, nperms), nrow = nsamples*npeakISO)*t(USED[order(EVALsum, decreasing = TRUE)[1:min(10, posEVAL, nperms)],])
eval.correction <- eval.correction*rep(apply(peakeval[rownames(peakeval) == "absent",], 2, sum), each = nsamples*npeakISO)/rep(apply(eval.correction, 2, sum), each = nsamples*npeakISO)
}

#####

peakeval <- peakeval[rownames(peakeval) != "absent",] + eval.correction
rownames(peakeval) <- c(rep(c(1:npeakISO), each = length(indies)))


outz <- matrix(sapply(c(1:npeakISO), factcond, peakeval), ncol = npeakISO)
frac.out <- add.fract[order(EVALsum, decreasing = TRUE)[1:min(10, nperms, posEVAL)]][apply(outz, 1, sum) == max(apply(outz, 1, sum))]
outz <- rbind(outz, rep(0, times = length(outz[1,])))


output <- data.frame(parentpkcombos = k, addnum = add, adduct = adduct.name[add], add.colZ, standard = c(1:length(adduct.to.pbya[,1]))[adduct.to.pbya[,add]][add.colZ[,2]], value = apply(outz, 2, max), abundance.fraction = round(frac.out, digits = 4))
output <- output[output$value != 0,]
output <- output[-is.na(output$value),]
if(length(output[,1]) == 0){output <- NULL}

add.output <- rbind(add.output, output)

}}}}

overall.add.output <- rbind(overall.add.output, add.output)
}


if(!is.null(overall.add.output)){

#identify adducts by parent-peak permutation (with the associated mu), parental peak origin and type of adduct, adduct 

Adduct.add <- matrix(data = 0, nrow = par.peak.perms, ncol = length(STD[,1]))
rownames(Adduct.add) <- c(1:par.peak.perms)
colnames(Adduct.add) <- rownames(STD)

for(i in 1:length(overall.add.output[,1])){	
	Adduct.add[overall.add.output$parentpkcombos[i],peakAssoc[overall.add.output$standard][i]] <- Adduct.add[overall.add.output$parentpkcombos[i],peakAssoc[overall.add.output$standard][i]] + overall.add.output$value[i]
	}
	
par.outz <- par.outz[-length(par.outz[,1]),]

stdnames <- rownames(combinedProbs[coToiso[,com],])
colznames <- stdnames[colZ$iso]

#make an adduct modification matrix that corresponds to the peak-correspondences in par.outz

add.outz <- matrix(data = 0, ncol = ifelse(is.matrix(par.outz), length(par.outz[1,]), length(par.outz)), nrow = par.peak.perms)

for(i in 1:length(colznames)){
	
if(is.matrix(par.outz) == TRUE){add.match <- ifelse(par.outz[,i] == 0, 0, 1) * Adduct.add[,colnames(Adduct.add) == colznames[i]]}else{
add.match <- ifelse(par.outz[i] == 0, 0, 1)*Adduct.add[,colnames(Adduct.add) == colznames[i]]}
if(length(add.match) == 0){add.outz[,i] <- rep(0, times = par.peak.perms)}else{add.outz[,i] <- add.match}
	
	}

#total confidence in parental peaks is the confidence in the peak and the combined support from adducts that agree with its abundance and located standards.

total.outz <- par.outz + add.outz

if(ADDUCT.OUT == FALSE){
	
output <- data.frame(colZ, standard = c(1:length(coToiso[,1]))[coToiso[,com]][colZ[,2]], value = apply(total.outz, 2, max))
output <- output[output$value != 0,]


	}
if(ADDUCT.OUT == TRUE){
	
peak.add <- paste(overall.add.output$standard, overall.add.output$peaks, sep = "-")
unique.peak.add <- unique(peak.add)

#take the best adduct for each peak-adduct correspondence across all parental peak combinations

top.adducts <- NULL

for(i in 1:length(unique.peak.add)){
	
top.adducts <- rbind(top.adducts, overall.add.output[overall.add.output$value == max(overall.add.output$value[peak.add == unique.peak.add[i]]),][1,])
	
	}

output <- data.frame(colZ, type = "STD", parent = c(1:length(coToiso[,1]))[coToiso[,com]][colZ[,2]], adduct = NA, value = apply(total.outz, 2, max))
	
output <- rbind(output, data.frame(peaks = top.adducts$peaks, iso = top.adducts$iso, type = "ADDUCT", parent = peakAssoc[top.adducts$standard], adduct = top.adducts$standard, value = top.adducts$value))	
}

output <- output[output$value != 0,]
output

}else{
	
if(ADDUCT.OUT == FALSE){
	
output <- data.frame(colZ, standard = c(1:length(coToiso[,1]))[coToiso[,com]][colZ[,2]], value = apply(par.outz, 2, max))
output <- output[output$value != 0,]


	}
if(ADDUCT.OUT == TRUE){	
	
	
output <- data.frame(colZ, type = "STD", standard = c(1:length(coToiso[,1]))[coToiso[,com]][colZ[,2]], adduct = NA, value = apply(par.outz, 2, max))
output <- output[output$value != 0,]
	
}
output

}}}}}}}

	
#for expected i isotopes w/ known MZ and l adducts, generate a m*l matrix of expected M/Z including the global offset 

	
MZtransform <- function(combinedAdds, STD, sub.combinedProbs){

STDmz <- STD$mass
abundSTDmz <- sub.combinedProbs[(sub.combinedProbs$uLabp + sub.combinedProbs$nLabp) > 0.1,]

transM <- NULL
for(i in 1:length(combinedAdds[,1])){
	if(combinedAdds$nmol[i] == 1){
		transM <- rbind(transM, data.frame(adduct = combinedAdds$adduct[i], add = combinedAdds$weightch[i], scale = abs(1/combinedAdds$charge[i]), unLab = combinedAdds$uLabp[i], nLab = combinedAdds$nLabp[i]))
	}else{
		transM <- rbind(transM, data.frame(adduct = combinedAdds$adduct[i], add = combinedAdds$weightch[i] + combinedAdds$nmol[i]*abundSTDmz$mass, scale = rep(abs(1/combinedAdds$charge[i])), unLab = combinedAdds$uLabp[i]*abundSTDmz$uLabp^(combinedAdds$nmol[i]), nLab = combinedAdds$nLabp[i]*abundSTDmz$nLabp^(combinedAdds$nmol[i])))
		}}
transM}

sum.by.k <- function(i, overall.add.output){
sum(overall.add.output$value[overall.add.output$parentpkcombos == i])
}


SDlmfit.spline <- function(vec, knot, nk, min, range, coef, HETbase, h){

SDspline <- smooth.spline(SDpoints, SDcoefMat[,i], df = HETpolyD, cv = TRUE, all.knots=TRUE)$fit
SDspline$knot <- sd.knot[,h]
SDspline$nk <- sd.nk[h]
SDspline$min <- sd.min[h]
SDspline$range <- sd.range[h]
SDspline$coef <- sd.coef[,h]

predict(SDspline, vec)$y

	}




