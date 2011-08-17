chemweight <- function(elFreq, chem.formula, Hchange, imerization){
options(warn=-1)
charz <- unlist(strsplit(chem.formula, split = ""))
prevchar <- charz[1]
rem <- rep(0, times = length(charz))
for(i in 2:length(charz)){
	if(prevchar %in% c("F") & charz[i] %in% c("e")){
		charz[i-1] <- "Fe"
		rem[i] <- 1
		next}
	if(prevchar %in% c("N") & charz[i] %in% c("a")){
		charz[i-1] <- "Na"
		rem[i] <- 1
		next}
	if(prevchar %in% c("C") & charz[i] %in% c("l")){
		charz[i-1] <- "Cl"
		rem[i] <- 1
		next}
	if(prevchar %in% c("B") & charz[i] %in% c("r")){
		charz[i-1] <- "Br"
		rem[i] <- 1
		next}			
		prevchar <- charz[i]	
	}
charz <- charz[rem == 0]

if(length(charz) > 1){
ischar <- as.numeric(charz)	
elcount <- as.numeric(unlist(strsplit(paste(as.character(ischar), sep = "", collapse = ""), split = "NA")))
elcount <- elcount[!is.na(elcount)]
if(length(elcount) == 0){elcount <- 0}
characters <- charz[is.na(ischar)]

if(length(elcount != length(characters))){
	pairedel <- data.frame(element = characters, count = rep(NA, times = length(characters)))
	waslet <- NA
	curnum <- 0
	for(i in 1:length(ischar)){
		#letter -> letter
		if(is.na(ischar[i]) & waslet %in% characters){
			pairedel[,2][pairedel[,1] == waslet] <- 1
			waslet <- charz[i]
			}
		#number/NA -> letter
		if(is.na(ischar[i]) & !(waslet %in% characters)){
			curnum <- curnum + 1
			waslet <- charz[i]
			}
		#letter -> number
		if(!is.na(ischar[i]) & waslet %in% characters){
			pairedel[,2][pairedel[,1] == waslet] <- elcount[curnum]
			waslet <- NA
			}
		#number -> number
		if(!is.na(ischar[i]) & !(waslet %in% characters)){			waslet <- NA		
			}}
		if(!is.na(waslet)){
			pairedel[,2][pairedel[,1] == waslet] <- 1
			}
		MW <- pairedel	
		}else{MW <- data.frame(element = characters, count = elcount)}}else{MW <- data.frame(element = charz, count = 1)}

#adjust for ionization and dimer/trimerization
if(sum(MW[,1] == "H") == 0){MW <- rbind(MW, data.frame(element = "H", count = 0))}


MW[MW[,1] == "H",2] <- MW[MW[,1] == "H",2] + Hchange
MW[,2] <- MW[,2]*imerization

if(MW[MW[,1] == "H",2] == 0){MW <- MW[!c(MW[,1] == "H"),]}


isoComp <- NULL
for(i in 1:length(MW[,1])){
isotopes <- elFreq[elFreq[,1] %in% MW[i,1],] 	
n <- MW[i,2]
isoComp <- rbind(isoComp, massrange(n, isotopes))}
	
niso <- table(isoComp[,1])
nameniso <- names(niso)
niso <- niso[niso != 0]

isocombos <- expand.grid(lapply(niso, nisotopes))
isocombos <- isocombos[,order(names(isocombos))]
isoComp <- isoComp[order(isoComp[,1]),]
MW <- MW[order(as.character(MW[,1])),]

if(!is.null(nrow(isocombos))){
possibleisos <- data.frame(formula = rep(NA, times = length(isocombos[,1])), prob = rep(NA, times = length(isocombos[,1])), mass = rep(NA, times = length(isocombos[,1])))

for(j in 1:length(isocombos[,1])){

aggregateForm <- NULL
for(i in 1:length(MW[,1])){
aggregateForm <- rbind(aggregateForm, isoComp[isoComp[,1] %in% MW[i,1],][isocombos[j,i],])	}

isotopeInfo <- combinedFormula(aggregateForm)
possibleisos[j,1] <- as.character(format(isotopeInfo[1]))
possibleisos[j,c(2:3)] <- isotopeInfo[2:3]
}
possibleisos[possibleisos[,2] > 0.001,]
}else{data.frame(formula = isoComp[,2], prob = isoComp[,3], mass = isoComp[,4])}
}

combinedFormula <- function(aggForm){
	formula <- paste(unlist(strsplit(paste(".", aggForm$comp), split = " "))[-1], collapse = "")
	prob <- prod(aggForm$probability)
	mass <- sum(aggForm$mass)
	data.frame(formula, prob, mass)
	}

massrange <- function(n, isotopes){
	#eval where np < 1/20 (min 1)
	if(n == 1){
	data.frame(element = isotopes[1,1], comp = paste(isotopes[,1], round(isotopes[,2]), sep = "_") , probability = isotopes[,3], mass = isotopes[,2])}else{	
	isotopepos = ceiling(n*isotopes[,3]+1)
	isotopepos[isotopepos > n] = n

	#return k counts of each isotope
	
	allperms <- expand.grid(lapply(isotopepos, isotopeapply))
	validperms <- allperms[apply(allperms,1, sum) == n,]
	if (is.null(dim(validperms)) == TRUE){
	isotopiccomp <- paste(paste(isotopes[,1], round(isotopes[2]), sep = "_"), n, sep = "-")
	comp <- data.frame(c(as.character(isotopes[,1]), isotopiccomp, 1, isotopes[2]))
	colnames(comp) <- c("element", "comp", "probability", "mass")
	comp
	}else{		
	probperms <- apply(validperms, 1, dmultinom, prob = isotopes[,3])
	validperms <- validperms[probperms > 0.001,]
	probperms <- probperms[probperms > 0.001]
	
	massperm <- apply(validperms, 1, sumprod, isotopes[,2])
	
	isotopename <- paste(isotopes[,1], round(isotopes[,2]), sep = "_")
	elcomp <- apply(validperms, 1, namenum, name = isotopename)
	
	data.frame(element = isotopes[1,1], comp = elcomp, probability = probperms, mass = massperm)}}}
		
namenum <- function(num, name){
	nonzeroname <- name[num != 0]
	nonzeronum <- num[num != 0]
	paste(unlist(strsplit(paste(".", mapply(paste, nonzeroname, nonzeronum, sep = "-")), split = " "))[-1], collapse = "")
	}
	
	
	
sumprod <- function(p1, p2){
	sum(p1*p2)
	}

nisotopes <- function(n){
	1:n}	
isotopeapply <- function(n){
	0:n}