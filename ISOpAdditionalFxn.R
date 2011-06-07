sumByISO <- function(unUSEDa, colZ){
sumMat <- matrix(data = FALSE, ncol = length(colZ[,1]), nrow = length(unique(colZ$iso)))
for(i in 1:length(colZ[,1])){
	sumMat[unique(colZ$iso) %in% colZ[i,2],i] <- TRUE
	}

output <- matrix(data = 0, ncol = length(unique(colZ$iso)), nrow = nperms)
for(i in 1:length(unique(colZ$iso))){
if(table(colZ$iso)[i] == 1){output[,i] <- unUSEDa[,sumMat[i,]]}else{output[,i] <- apply(unUSEDa[,sumMat[i,]], 1, sum)}}
ifelse(output == 1, 0, 1)
}


GaussianLik2 <- function(com, coToiso, posL, peaksizeMat, probMat, npeaks, h, SDlmMat, HETbase){
#single peaks have the right ratio
print(com)
subprob <- probMat[coToiso[,com],]
posLsub <- posL[,coToiso[,com]]

nfacs <- apply(posLsub, 2, sumthresh, thresh = -60, nvec = npeaks)
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

validP <- levgrid[apply(definedP, 1, validPerms, Nf, Pf, isoCov),]

if(length(validP[,1]) != 0){

colZ <- stacker[!is.na(stacker)[,1],]

probMatsub <- probMat[coToiso[,com],]


nuMiso <- length(levgrid[1,])
npeakISO <- length(colZ[,1]) 
nperms <- length(validP[,1])

indies <- c(1:nsamples)

PMAT <- matrix(data = unlist(t(peaksizeMat[colZ[,1],])), ncol = nsamples*npeakISO, nrow = nperms, byrow = TRUE)
#colnames(PMAT) <- rep(colZ[,1], each = nsamples)
colnames(PMAT) <- rep(indies, times = npeakISO)

PPROB <- matrix(data = unlist(t(as.data.frame(probMatsub[colZ[,2],]))), ncol = nsamples*npeakISO, nrow = nperms, byrow = TRUE)
colnames(PPROB) <- rep(indies, times = npeakISO)

sdPMAT <- matrix(data = peakSD(unlist(t(peaksizeMat[colZ[,1],])), h, SDlmMat, HETbase), ncol = nsamples*npeakISO, nrow = nperms, byrow = TRUE)
colnames(sdPMAT) <- rep(indies, times = npeakISO)

aMAT <- matrix(data = unlist(peaksizeMat[colZ[,1],])/peakSD(unlist(peaksizeMat[colZ[,1],]), h, SDlmMat, HETbase), ncol = nsamples*npeakISO, nrow = nperms, byrow = TRUE)
colnames(sdPMAT) <- rep(indies, times = npeakISO)

posLinfo <- NULL
for(i in 1:length(colZ[,1])){
	posLinfo <- c(posLinfo, posLsub[colZ[i,1], colZ[i,2]])
}

POSLMAT <- matrix(rep(posLinfo, each = nsamples), ncol = nsamples*npeakISO, nrow = nperms, byrow = TRUE) 


#convert validP to peak-sample x nperm format using colZ

isoRep <- table(colZ[,2])*nsamples
gridEXP <- validP[,as.numeric(rep(names(isoRep), times = isoRep))]
gridP <- matrix(data = (rep(colZ[,1], each = nsamples)), ncol = nsamples*npeakISO, nrow = nperms, byrow = TRUE)
USED <- ifelse(gridEXP-gridP == 0 & !is.na(gridEXP-gridP), 1, 0)


#USED <- ifelse(gridEXP*gridINV == 1 & !is.na(gridEXP*gridINV ), 1, 0)


unUSEDa <- USED[,rep(c(1:nsamples), times = npeakISO) == 1]
unUSED <- sumByISO(unUSEDa, colZ)

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

unUSEDf <- NULL
for(i in 1:length(unique(colZ$iso))){
for(k in 1:nsamples){
	unUSEDf <- cbind(unUSEDf, unUSED[,i])
	}}

unUSEDf <- ifelse(unUSEDf == 0, NA, 1)

unUSEDlik <- log(dnorm(matrix(MUcolz, ncol = nsamples*length(unique(colZ$iso)), nrow = nperms, byrow = FALSE)*matrix(unlist(t(probMatsub[unique(colZ$iso),])), ncol = nsamples*length(unique(colZ$iso)), nrow = nperms, byrow = TRUE), mean = 0, peakSD(300, h, SDlmMat, HETbase))*unUSEDf, base = 10)

##

EVAL <- (log(gausD(PMAT, matrix(MUcolz, ncol = nsamples*npeakISO, nrow = nperms, byrow = FALSE)*USED*PPROB, sdPMAT), base = 10) + POSLMAT)*USED
EVALsum <- apply(EVAL, 1, sum) + apply(unUSEDlik, 1, sum, na.rm = TRUE)
peakeval <- matrix(EVAL[EVALsum >= max(EVALsum) - 10,], nrow = nsamples*npeakISO, byrow = TRUE)
rownames(peakeval) <- rep(c(1:npeakISO), each = length(indies))

outz <- matrix(sapply(c(1:npeakISO), factcond, peakeval), ncol = npeakISO)
#outz <- rbind(outz, rep(0, times = length(outz[1,])))

output <- data.frame(colZ, standard = c(1:length(coToiso[,1]))[coToiso[,com]][colZ[,2]], value = apply(outz, 2, max))
output <- output[output$value != 0,]

if(length(output[,1]) != 0){
output}

}}}