setwd("/Users/seanhackett/Desktop/RabinowitzLab/IsoPeak/")

source(file = "peakLib.R")

options(digits = 13)

elFreq <- data.frame(element = c("H", "H", "O", "O", "O", "C", "C", "N", "N", "Na", "Fe", "Fe", "Fe", "Fe", "K", "K", "K", "S", "S", "S", "S", "P", "Cl", "Cl", "Br", "Br", "F"), MW = c( 1.00782503207, 2.0141017778, 15.99491461956, 16.99913170, 17.9991610, 12.0000000, 13.0033548378, 14.0030740048, 15.0001088982, 22.9897692809, 53.9396105, 55.9349375, 56.9353940, 57.9332756, 38.96370668, 39.96399848, 40.96182576, 31.97207100, 32.97145876, 33.96786690, 35.96708076, 30.973762, 34.968853, 36.965903, 78.918336, 80.916290, 18.998403), abund = c(0.999885, 0.000115, 0.99757, 0.00038, 0.00205, 0.9893, 0.0107, 0.99636, 0.00364, 1.000, 0.05845, 0.91754, 0.02119, 0.00282, 0.932581, 0.000117, 0.067302, 0.9499, 0.0075, 0.0425, 0.0001, 1.000, 0.7577, 0.2423, 0.5069, 0.4931, 1))

#isotopic labeling expected abundance
### 99% N-15

NlabFreq <- data.frame(element = c("H", "H", "O", "O", "O", "C", "C", "N", "N", "Na", "Fe", "Fe", "Fe", "Fe", "K", "K", "K", "S", "S", "S", "S", "P", "Cl", "Cl", "Br", "Br", "F"), MW = c( 1.00782503207, 2.0141017778, 15.99491461956, 16.99913170, 17.9991610, 12.0000000, 13.0033548378, 14.0030740048, 15.0001088982, 22.9897692809, 53.9396105, 55.9349375, 56.9353940, 57.9332756, 38.96370668, 39.96399848, 40.96182576, 31.97207100, 32.97145876, 33.96786690, 35.96708076, 30.973762, 34.968853, 36.965903, 78.918336, 80.916290, 18.998403), abund = c(0.999885, 0.000115, 0.99757, 0.00038, 0.00205, 0.9893, 0.0107, 0.01, 0.99, 1.000, 0.05845, 0.91754, 0.02119, 0.00282, 0.932581, 0.000117, 0.067302, 0.9499, 0.0075, 0.0425, 0.0001, 1.000, 0.7577, 0.2423, 0.5069, 0.4931, 1))




####### Known compounds with retention times and chemical formula ######

knowns <- read.table("KNOWNS.csv", sep = ",", header = TRUE, fill = TRUE, colClasses = c("character", "character", "character", "numeric"))

#Fix messed up values
knowns[,3][as.character(knowns[,3]) %in% c("", "?")] <- NA

numerics <- as.numeric(as.character(knowns[,3]))

numericpos <- c(1:length(knowns[,1]))[!is.na(numerics)]

#check when running
knowns[numericpos,]	
knowns[numericpos,4] <- as.numeric(as.character(knowns[numericpos,3]))
knowns[numericpos,3] <- NA


knownRT <- knowns[!is.na(knowns[,4]),]

### Redundant compounds are removed ###

redund <- table(knownRT[,1])[table(knownRT[,1]) > 1]

knownRT[knownRT$compound %in% names(redund),]

# Remove 
removerow <- c(10, 12, 25)

knownRT <- knownRT[-removerow,]

#write.table(knownRT, file = "knownRTchemForm.csv", sep = ",", col.names = TRUE, row.names = FALSE)

MzRT <- NULL
NlabMzRT <- NULL

for (i in 1:length(knownRT[,1])){

chem.formula <- knownRT[i,2]
chem.weight <- chemweight(elFreq, chem.formula, -1, 1)
chem.weightlab <- chemweight(NlabFreq, chem.formula, -1, 1)

MzRT <- rbind(MzRT, data.frame(compound = knownRT$compound[i], chem.weight, RT = knownRT$rt[i]))
NlabMzRT <- rbind(NlabMzRT, data.frame(compound = knownRT$compound[i], chem.weightlab, RT = knownRT$rt[i]))
}

######### Determine MW distribution of adducts ##########

adducts <- read.table("Saved_Filez/ADDUCTS.csv", sep =",", header = TRUE)

negions <- adducts[adducts$charge < 1,]
Unadd <- NULL

for(i in 1:length(negions[,1])){
	
	addform <- as.character(negions$addForm[i])
	remform <- as.character(negions$remove[i])
	
	add.weight <- if(!is.na(addform)){chemweight(elFreq, addform, 1, 1)}else{add.weight <- NA}
	rem.weight <- if(!is.na(remform)){chemweight(elFreq, remform, 0, 1)}else{rem.weight <- NA}
	
	if(is.na(rem.weight)){
		Unadd <- rbind(Unadd, data.frame(adduct = negions[i,1], nmol = negions$nmol[i], charge = negions$charge[i], weightprob = add.weight$prob, weightch = add.weight$mass))
		next}
	if(is.na(add.weight)){
		Unadd <- rbind(Unadd, data.frame(adduct = negions[i,1], nmol = negions$nmol[i], charge = negions$charge[i], weightprob = rem.weight$prob, weightch = ((rem.weight$mass*-1)+elFreq[1,2])))
		next}
	else{
	combos <- expand.grid(1:length(add.weight[,1]), 1:length(rem.weight[,1]))
	if(length(combos[,1]) == 1){Unadd <- rbind(Unadd, data.frame(adduct = negions[i,1], nmol = negions$nmol[i], charge = negions$charge[i], weightprob = add.weight$prob*rem.weight$prob, weightch = add.weight$mass-rem.weight$mass))}else{
		isoCombs <- NULL
		for(k in 1:length(combos[,1])){
		isoCombs <- rbind(isoCombs, data.frame(adduct = negions[i,1], nmol = negions$nmol[i], charge = negions$charge[i], weightprob = add.weight$prob[combos[k,1]]*rem.weight$prob[combos[k,2]], weightch = add.weight$mass[combos[k,1]] - rem.weight$mass[combos[k,2]]))}
		Unadd <- rbind(Unadd, isoCombs[isoCombs$weightprob > 0.001,])
		}}}	

Nadd <- NULL

for(i in 1:length(negions[,1])){
	
	addform <- as.character(negions$addForm[i])
	remform <- as.character(negions$remove[i])
	
	add.weight <- if(!is.na(addform)){chemweight(NlabFreq, addform, 1, 1)}else{add.weight <- NA}
	rem.weight <- if(!is.na(remform)){chemweight(NlabFreq, remform, 0, 1)}else{rem.weight <- NA}
	
	if(is.na(rem.weight)){
		Nadd <- rbind(Nadd, data.frame(adduct = negions[i,1], nmol = negions$nmol[i], charge = negions$charge[i], weightprob = add.weight$prob, weightch = add.weight$mass))
		next}
	if(is.na(add.weight)){
		Nadd <- rbind(Nadd, data.frame(adduct = negions[i,1], nmol = negions$nmol[i], charge = negions$charge[i], weightprob = rem.weight$prob, weightch = ((rem.weight$mass*-1)+elFreq[1,2])))
		next}
	else{
	combos <- expand.grid(1:length(add.weight[,1]), 1:length(rem.weight[,1]))
	if(length(combos[,1]) == 1){Nadd <- rbind(Nadd, data.frame(adduct = negions[i,1], nmol = negions$nmol[i], charge = negions$charge[i], weightprob = add.weight$prob*rem.weight$prob, weightch = add.weight$mass-rem.weight$mass))}else{
		isoCombs <- NULL
		for(k in 1:length(combos[,1])){
		isoCombs <- rbind(isoCombs, data.frame(adduct = negions[i,1], nmol = negions$nmol[i], charge = negions$charge[i], weightprob = add.weight$prob[combos[k,1]]*rem.weight$prob[combos[k,2]], weightch = add.weight$mass[combos[k,1]] - rem.weight$mass[combos[k,2]]))}
		Nadd <- rbind(Nadd, isoCombs[isoCombs$weightprob > 0.001,])
		}}}


###Combine N and unlabeled adduct libraries

adds <- unique(Nadd$adduct)
		
combinedAdds <- NULL

for(i in 1:length(adds)){
sublist <- Unadd[Unadd$adduct %in% adds[i],]
nsublist <- Nadd[Nadd$adduct %in% adds[i],]

jointmass <- union(sublist$weightch, nsublist$weightch)
if(length(jointmass[!(jointmass %in% sublist$weightch)]) != 0){
sublist  <- rbind(sublist, data.frame(adduct = sublist$adduct[1], nmol = sublist$nmol[1], charge = sublist$charge[1], weightch = jointmass[!(jointmass %in% sublist$weightch)], weightprob = 0)) 
sublist <- sublist[order(sublist$weightch),]
	}
	
if(length(jointmass[!(jointmass %in% nsublist$weightch)]) != 0){
nsublist  <- rbind(nsublist, data.frame(adduct = nsublist$adduct[1], nmol = nsublist$nmol[1], charge = nsublist$charge[1], weightch = jointmass[!(jointmass %in% nsublist$weightch)], weightprob = 0))  
nsublist <- nsublist[order(nsublist$weightch),]
	}	

combinedAdds <- rbind(combinedAdds, data.frame(adduct = sublist$adduct, nmol = sublist$nmol, charge = sublist$charge, weightch = sublist$weightch, uLabp = sublist$weightprob, nLabp = nsublist$weightprob))}

		
### Only interested in (-) adducts because exactive generates (-) ions ###

#save(combinedAdds, file = "negAdducts.R")

#save(MzRT, NlabMzRT, file = "knownMzRT.R")

load("knownMzRT.R")

########### Isotopes with similar M/z are pooled ########
##### Rounded if M/z diff < 0.0000001 ####

compounds <- unique(MzRT$compound)
MzRTrefine <- NULL
nisoMzRTrefine <- NULL
for(i in 1:length(compounds)){

isotopes <- MzRT[MzRT$compound %in% compounds[i],]
probround <- round(isotopes$mass, digits = 8)
nisotopes <- NlabMzRT[NlabMzRT$compound %in% compounds[i],]
nprobround <- round(nisotopes$mass, digits = 8)

uniqueprob <- unique(probround)
nuniqueprob <- unique(nprobround)


for(j in 1:length(uniqueprob)){
MzRTrefine <- rbind(MzRTrefine, data.frame(compound = as.character(compounds[i]), mass = uniqueprob[j], prob = sum(isotopes[probround %in% uniqueprob[j],]$prob), RT = isotopes$RT[1]))
	}
for(j in 1:length(nuniqueprob)){
nisoMzRTrefine <- rbind(nisoMzRTrefine, data.frame(compound = as.character(compounds[i]), mass = nuniqueprob[j], prob = sum(nisotopes[nprobround %in% nuniqueprob[j],]$prob), RT = nisotopes$RT[1]))
	}	
	}

#save(MzRTrefine, nisoMzRTrefine, Unadd, file = "knownMzRTre.R")
load("knownMzRTre.R")	
	
############ Plots of Predicted and Observed Peaks #################

allPeaks  <- read.table("allslices.csv", sep = ",", header = TRUE)

library("gplots")

#observed peaks

colors <- redgreen(1000)
minval = min(allPeaks$maxQuality)
incval = (max(allPeaks$maxQuality) - minval)/1000
colbyQual <- colors[round((allPeaks$maxQuality - minval)/incval)]

peakCol <- redgreen(7)
peakcount <- peakCol[allPeaks$goodPeakCount]

peaks <- c("N.1", "N.2", "P.1", "P.2")
blanks <- c("blank.1", "blank.2", "blank.3")

meansamp <- apply(allPeaks[,colnames(allPeaks) %in% peaks], 1, mean)
meanBL <- apply(allPeaks[,colnames(allPeaks) %in% blanks], 1, mean)

plot(log(meansamp) ~ log(meanBL))

peaksize <- meansamp - meanBL

peaksizeMat <- allPeaks[,colnames(allPeaks) %in% peaks] - meanBL
### Baseline
peaksizeMat[peaksizeMat < 300] <- 300
sampleclass <- c("N", "N", "P", "P")

#known metabolites
#MzRT

plot(allPeaks$medMz ~ allPeaks$medRt, col = peakcount, pch = 2, type = "p", cex = peaksize/max(peaksize)*10)
#points(MzRTrefine$mass ~ MzRTrefine$RT, col = "pink", pch = 19, cex = MzRTrefine$prob/max(MzRTrefine$prob))
points(MzRTtop$mass ~ MzRTtop$RT, col = "pink", pch = 19)


####### Determine most abundant expected peak for each compound ########

compounds <- unique(MzRTrefine$compound)

MzRTtop <- NULL

for(i in 1:length(compounds)){
	sublist <- MzRTrefine[MzRTrefine$compound %in% compounds[i],]
	MzRTtop <- rbind(MzRTtop, sublist[sublist[,3] == max(sublist[,3]),])
	}


############ Compared probabilities for Predicted and Observed Peaks ############

#Peakcomp <- matrix(data = NA, ncol = length(allPeaks[,1]), nrow = length(MzRTtop[,1]))
massvec <- allPeaks$medMz
rtvec <- allPeaks$medRt
#peaksize

for(i in 1:length(compounds)){
sublist <- MzRTrefine[MzRTrefine$compound %in% compounds[i],]
nsublist <- nisoMzRTrefine[nisoMzRTrefine$compound %in% compounds[i],]

jointmass <- union(sublist$mass, nsublist$mass)
if(length(jointmass[!(jointmass %in% sublist$mass)]) != 0){
sublist  <- rbind(sublist, data.frame(compound = sublist$compound[1], mass = jointmass[!(jointmass %in% sublist$mass)], prob = 0, RT = sublist$RT)) 
sublist <- sublist[order(sublist$mass),]
	}
	
if(length(jointmass[!(jointmass %in% nsublist$mass)]) != 0){
nsublist  <- rbind(nsublist, data.frame(compound = nsublist$compound[1], mass = jointmass[!(jointmass %in% nsublist$mass)], prob = 0, RT = nsublist$RT[1])) 
nsublist <- nsublist[order(nsublist$mass),]
	}	

nanneal <- 10000
anneals <- rep(NA, times = nanneal)
for(j in 1:nanneal){
	anneals[j] = 1*j^-0.05
	}
anneals <- (anneals - min(anneals))*(1/(max(anneals) - min(anneals)))

annealvar <- rep(c(rep(1, times = 5), rep(0.2, times = 5), rep(0.05, times = 5), rep(0.01, times = 5)), times = nanneal/20)


paroffset <- matrix(data = NA, ncol = 2, nrow = nanneal)
paroffset[1,] <- c(0, 0)
massvec <- allPeaks$medMz
rtvec <- allPeaks$medRt

### sample updated offset from the joint dist N(offset[1], 1) and N(offset[2], offset[2]/13)

previouspt = NA
liktrack = rep(NA, times = nanneal)
for(j in 1:nanneal){
if(j == 1){
	previouspt <- compoundprob(massvec, rtvec, paroffset[1,], peaksizeMat, sampleclass, sublist, nsublist)
	liktrack[j] <- previouspt
	}else{
	offsetdraw <- c(rnorm(1, paroffset[j-1,1], 2*annealvar[j]), rnorm(1, paroffset[j-1,2], (sublist$RT[1]+paroffset[j-1,2])*(annealvar[j]/4)))
		
	currentlik <- compoundprob(massvec, rtvec, offsetdraw, peaksizeMat, sampleclass, sublist, nsublist) + log(dnorm(offsetdraw[1], 0, 1)*dnorm(offsetdraw[2], 0, sublist$RT[1]/13))
	
	alpha = exp(currentlik - previouspt)
	if(alpha > 1){previouspt <- currentlik
		paroffset[j,] <- offsetdraw
		liktrack[j] <- currentlik
		}else{
		if(abs(1-alpha) < runif(1, 0, 1)*anneals[j]){
			previouspt <- currentlik  
			paroffset[j,] <- offsetdraw
			liktrack[j] <- currentlik
		}else{paroffset[j,] <- paroffset[j-1,]
			liktrack[j] <- previouspt
			}}}}


#### Visualization ######

plot(paroffset[,1] ~ paroffset[,2], cex = 0.01, type = "l")
points(jitter(paroffset[,1], amount = 0.01) ~ jitter(paroffset[,2], amount = 0.03))



concatmz <- NULL
concatrt <- NULL
concatp <- NULL
concatcol <- NULL
for(x in 1:length(sublist[,1])){
	concatrt <- c(concatrt, paroffset[,2] + sublist$RT[x])
	concatmz <- c(concatmz, paroffset[,1] + sublist$mass[x])
	concatp <- c(concatp, rep(mean(c(sublist$prob[x], nsublist$prob[x])), times = nanneal))
	concatcol <- c(concatcol, redblue(nanneal))
	}


plot(concatrt ~ concatmz, cex = concatp*5, col = concatcol)
points(allPeaks$medRt ~ allPeaks$medMz, col = rainbow(1, s=1, start = 0.2, end = 1, alpha = 0.4), pch = 16, cex = (log(peaksize/50)))

plot(liktrack)

######################## Fxns #######################

compoundprob <- function(massvec, rtvec, paroffset, peaksizeMat, sampleclass, sublist, nsublist){

### peaks near the offset are defined by mass density N(massvec + offset, 0.01) and rt density sd = cmpdrt/130
floor = log(dnorm(0.4, 0, 0.05)^2*300)

isomerabund <- rep(NA, times = length(sublist[,1]))
for(i in 1:length(sublist[,1])){
parprob <- log(dnorm(massvec, sublist$mass[i] + paroffset[1], 0.15)*dnorm(rtvec, sublist$RT[i] + paroffset[2], (sublist$RT[i] - paroffset[2])/60))

weights <- NULL
for(j in 1:length(sampleclass)){
if(sampleclass[j] == "N"){
weights <- cbind(weights, log(peaksizeMat[,j]*nsublist$prob[i]))
	}
if(sampleclass[j] == "P"){	
weights <- cbind(weights, log(peaksizeMat[,j]*sublist$prob[i]))
	}
}
#floor of the equivalent of N(0.3, 0.05)^2 off of target and baseline signal
maxval <- mapply(max,apply(parprob + weights, 2, max), floor)
isomerabund[i] <- sum(maxval)
}	
sum(isomerabund)
}




# RT probability from N(compoundRT, rt/13)
# MZ probability from N(compoundMz, 1)

#determine appropriate RT probability scaling
retTimes = seq(1, 21, by = 2)
timepts = seq(0, 25, by = 0.01)
retDense <- matrix(data = NA, ncol = length(retTimes), nrow = length(timepts))
for(i in 1:10){
	retDense[,i] <- dnorm(timepts, mean = retTimes[i], sd = retTimes[i]/13)
	}

plot(retDense[,1] ~ timepts, type = "l")
for(i in 2:10){
	lines(retDense[,i] ~ timepts)
	}
	



offsetdraw <- data.frame(x = rnorm(5000, 0, 1), y = rnorm(5000, 0, 5))
plot(offsetdraw[,2] ~ offsetdraw[,1])






