gausD <- function(peaks, expected, sdP){
	1/(sqrt(2*pi)*sdP) * exp(-((peaks - expected)^2/(2*(sdP^2))))
	}

validPerms <- function(vec, Nf, Pf, isoCov){
	if(sum(vec*Nf) > isoCov | sum(vec*Pf) > isoCov){TRUE}else{FALSE}
	}
	
factcond <- function(i, peakeval){
	apply(matrix(peakeval[rownames(peakeval) == i,], nrow = nsamples), 2, sum)
	}
	
peakSD <- function(peaksizeMat, h, SDlmMat, HETbase){
	sdCoef <- SDlmMat[,h]
	logSize <- log(peaksizeMat, base = HETbase)
	HETbase^(sdCoef[1] + sdCoef[2]*logSize + sdCoef[3]*logSize^2 + sdCoef[4]*logSize^3)}
	
harmM <- function(x){length(x)/(sum(1/x))}

masserror <- function(peakvec, standard){
	((peakvec - standard)/standard)*10^6}
	
RTdiff <- function(peakvec, standard){
	peakvec - standard}
	
sumthresh <- function(vec, thresh, nvec){
	c(1:nvec)[vec > thresh]
	}



############# Determine likelihood of a set of peaks corresponding to a compouds isotopes given the current evidence from peak size and position (posL) and the observed ratios of isotopes ###############

GaussianLik <- function(com, coToiso, posL, peaksizeMat, probMat, npeaks, h, SDlmMat, HETbase){
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
USED <- ifelse(is.na(gridEXP), 0, 1)



#gridINV <- matrix(data = (1/rep(colZ[,1], each = nsamples)), ncol = nsamples*npeakISO, nrow = nperms, byrow = TRUE)
#USED <- ifelse(gridEXP*gridINV == 1 & !is.na(gridEXP*gridINV ), 1, 0)

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






REPMOD <- USED*10
EVAL <- (log(gausD(PMAT, matrix(MUcolz, ncol = nsamples*npeakISO, nrow = nperms, byrow = FALSE)*USED*PPROB, sdPMAT), base = 10) + POSLMAT)*USED + REPMOD
EVALsum <- apply(EVAL, 1, sum)
peakeval <- matrix(EVAL[EVALsum >= max(EVALsum) - 10,], nrow = nsamples*npeakISO, byrow = TRUE)
rownames(peakeval) <- rep(c(1:npeakISO), each = length(indies))

outz <- matrix(sapply(c(1:npeakISO), factcond, peakeval), ncol = npeakISO)
outz <- rbind(outz, rep(0, times = length(outz[1,])))


output <- data.frame(colZ, standard = c(1:length(coToiso[,1]))[coToiso[,com]][colZ[,2]], value = apply(outz, 2, max))
output <- output[output$value != 0,]

if(length(output[,1]) != 0){
output}

}}}
		

	
			