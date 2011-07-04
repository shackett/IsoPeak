setwd("/Users/seanhackett/Desktop/RabinowitzLab/IsoPeak/")

load("SpectrumScale.R")
load("Saved_Filez/PeakFparams.R")

library(gplots)

nanneal <- j

#### Summary Plotz

### RT

RTpoints <- RTpos*RTtrack[nanneal,]
RTcoefs <- summary(lm(RTpoints ~ RTpos + I(RTpos^2) + I(RTpos^3)))$coef[,1]
RTpoly <- RTcoefs[1] + RTcoefs[2]*RTpos + RTcoefs[3]*RTpos^2 + RTcoefs[4]*RTpos^3
	
plot(RTpoints ~ RTpos)
lines(RTpoly ~ RTpos)

library(gplots)
annealCol <- redblue(nanneal)

for(j in 1:nanneal){
RTpoints <- RTpos*RTtrack[j,]
RTcoefs <- summary(lm(RTpoints ~ RTpos + I(RTpos^2) + I(RTpos^3)))$coef[,1]
RTpoly <- RTcoefs[1] + RTcoefs[2]*RTpos + RTcoefs[3]*RTpos^2 + RTcoefs[4]*RTpos^3

if(j == 1){plot(RTpoly ~ RTpos, type = "l", col = annealCol[1])}else{lines(RTpoly ~ RTpos, col = annealCol[j])}}

### MZ offset

plot(MZoffsetrack, col = annealCol)

### SD offset

Hetrange <- c(floor(range(log(valS[,colnames(valS) %in% samples], base = 2))[1]) : ceiling(range(log(valS[,colnames(valS) %in% samples], base = HETbase))[2]))
expHetRange <- HETbase^Hetrange

Hetpoints <- log(HETbase^Hetrange*SDtrack[nanneal,], base = HETbase)
SDcoefs <- summary(lm(Hetpoints ~ Hetrange + I(Hetrange ^2) + I(Hetrange ^3)))$coef[,1]
SDpoly <- SDcoefs[1] + SDcoefs[2]*Hetrange + SDcoefs[3]*Hetrange^2 + SDcoefs[4]*Hetrange^3

plot(Hetpoints ~ Hetrange)
lines(SDpoly ~ Hetrange)

for(j in 1:nanneal){
Hetpoints <- log(HETbase^Hetrange*SDtrack[j,], base = HETbase)
SDcoefs <- summary(lm(Hetpoints ~ Hetrange + I(Hetrange ^2) + I(Hetrange ^3)))$coef[,1]
SDpoly <- SDcoefs[1] + SDcoefs[2]*Hetrange + SDcoefs[3]*Hetrange^2 + SDcoefs[4]*Hetrange^3
expSDpoly <- HETbase^SDpoly	
	
#if(j == 1){plot(SDpoly ~ Hetrange, type = "l", col = annealCol[1])}else{lines(SDpoly ~ Hetrange, col = annealCol[j])}

if(j == 1){plot(expSDpoly ~ expHetRange, type = "l", col = annealCol[1])}else{lines(expSDpoly ~ expHetRange, col = annealCol[j])}
}

#for just the Coefficient of variation

for(j in 1:nanneal){
Hetpoints <- log(HETbase^Hetrange*SDtrack[j,], base = HETbase)
SDcoefs <- summary(lm(Hetpoints ~ Hetrange + I(Hetrange ^2) + I(Hetrange ^3)))$coef[,1]
SDpoly <- SDcoefs[1] + SDcoefs[2]*Hetrange + SDcoefs[3]*Hetrange^2 + SDcoefs[4]*Hetrange^3
expSDpoly <- HETbase^SDpoly	
	
#if(j == 1){plot(SDpoly ~ Hetrange, type = "l", col = annealCol[1])}else{lines(SDpoly ~ Hetrange, col = annealCol[j])}

if(j == 1){plot(expSDpoly/expHetRange ~ expHetRange, type = "l", col = annealCol[1], ylim = c(0,3))}else{lines(expSDpoly/expHetRange ~ expHetRange, col = annealCol[j])}
}

	
#### SD(RT)

plot(RT.SDtrack ~ c(1:nanneal), col = redblue(nanneal))

RT.SDtrack



	
	
		
######## Determine peaks corresponding to the final 

pMZ = MZvals + MZoffsetrack[nanneal]


RTpoints <- RTpos*RTtrack[nanneal,]
RTcoefs <- summary(lm(RTpoints ~ RTpos + I(RTpos^2) + I(RTpos^3)))$coef[,1]
pRT <- RTcoefs[1] + RTcoefs[2]*RTvals + RTcoefs[3]*RTvals^2 + RTcoefs[4]*RTvals^3
pRT <- RTvals
	

#compare against all standards - combinedProbs$mass
	
MZe <- log(dnorm(sapply(combinedProbs$mass, masserror, standard = pMZ), mean = 0, sd = 3))
RTe <- log(dnorm(sapply(combinedProbs$RT, RTdiff, standard = pRT), mean = 0, sd = 0.5))
#SIZe <- t(log(probMat %*% t(peaksizeMat)))
SIZe <- t(probMat %*% t(log(peaksizeMat, base = 2)))
"


#persLik <- matrix(NA, nrow = length(peakpos), ncol = length(combinedProbs[,1]))

### Determine likelihood of best peak-freq match
### applied over compounds