library(dplyr)
library(evd)
library(tidyverse)
library(multcomp)
library(gridExtra)
library("lmtest")
library("sandwich")
library('extraDistr')
library("permute")
library("xlsx")
library(ggplot2)
library(DescTools)
library(patchwork)

set.seed(218219)

### read input data
metdta <- data.table::fread("metadata.csv", data.table = FALSE)
metdta <- metdta %>% relocate(Subject_ID, Model_ID)
colnames(metdta)[1] = "sampName"

simres <-data.table::fread("inputDiet_net_secretion_fluxes.csv", data.table = FALSE)
simres <- simres[, colSums(is.na(simres)) != nrow(simres)]
colnames(simres)[1] <- "NMPCs"

ag_sec_mets <- simres$NMPCs
faec_mets  <- data.table::fread("metabolite_data.csv", header = TRUE)

# AGORA <-> keggIds
agr <- read.csv(file = "AGORA1_metabolites_abbr.csv")
agr <- agr[-which(agr$keggId %in% agr$keggId[duplicated(agr$keggId)]), ]

a <- gsub(".*EX_","", ag_sec_mets)
b <- gsub("\\[.*", "", a)

for (i in 1:length((b))){
  if(b[i] %in% agr$abbreviation){
    b[i] = agr$keggId[match(b[i],agr$abbreviation)]
  }
}

# AGORA secreted mets with keggIds
keggId_ag <- b[grepl('^_', b)]
keggId_ag <- gsub("_", "", keggId_ag)

# search keggIds in AGORA reconstruction in faecal metabolome
Pattern = paste(keggId_ag, collapse="|")
keepId <- grepl(Pattern, faec_mets$V1)
faec_mets <- faec_mets[keepId, ]

# only consider metabolites >50%
#faec_mets <- faec_mets[ rowSums(faec_mets > 0) > length(faec_mets)/2, ]
faec_mets$V1 <- gsub("_.*", "", faec_mets$V1)

#kgg <- gsub("_.*", "", faec_mets$V1)
#intersect(kgg, gsub("_", "", agr$keggId))
#vmhabb <- agr$abbreviation[which(agr$keggId %in% gsub("^", "_", intersect(kgg, gsub("_", "", agr$keggId))))]

# search keggIds in simulation results
simres$NMPCs <- b
Pattern <- paste(faec_mets$V1, collapse="|")
keepId <- grepl(Pattern, simres$NMPCs)
simres <- simres[keepId, ]

simres$NMPCs <- gsub("_", "EX_", simres$NMPCs)
colnames(simres) <- gsub("sample_", "", colnames(simres))
simres <- simres[, ! names(simres) %in%  colnames(simres)[grepl("_1", colnames(simres))], drop = F]

# sample with only zeros - 10177
#simres <- simres[, ! names(simres) %in% c("10177"), drop = F]

# normcoverage abundances
normcov <-data.table::fread("normCoverage_CRC.csv", data.table = FALSE)
abunspec <- normcov

abunspec <- abunspec[, 1:617]
colnames(abunspec) <- gsub(".*_", "", colnames(abunspec))

# merge
allmerged <- rbind(faec_mets, simres, abunspec, fill=TRUE)
allmerged  <- allmerged  %>%   select_if(~ !any(is.na(.)))
tamgd <- t(allmerged)
tamgd <- as.data.frame(tamgd)
tamgd$sampNames <- rownames(tamgd)

allmerged <- tamgd[ , c(ncol(tamgd), 1:(length(tamgd)-1))]
rownames(allmerged) <- c()
colnames(allmerged) <- c("sampName", faec_mets$V1, simres$NMPCs, abunspec$V1)

# delete metabolites that are <50%
mets <- allmerged[, 2:107]
ts1<- mets[,colSums(mets[,]>0.001) > 173]

exflx <- gsub("C", "EX_C", colnames(ts1))
fluxes <- allmerged[, exflx]
ts2 <- fluxes[,colSums(fluxes[,]>0.001) > 173]

ts1  <- gsub("EX_", "", colnames(ts2))
abuncrct <- allmerged[, 214:575]
ts3 <- abuncrct[,colSums(abuncrct[,]>0) > 34]
keepcols <- c(ts1, colnames(ts2), colnames(ts3))
df <- allmerged[, c("sampName", keepcols)]

# metadata
df <- merge(df, metdta, by = "sampName") # NA's match
df$BMI <- as.numeric(df$BMI)
df$Age <- as.numeric(df$Age)

# filter sign spec
# remove C00695, C00245, keep C00329
df <- df[ , !(names(df) %in% c("C00695", "C00245", "EX_C00695", "EX_C00245"))]
micr <- colnames(df)[grepl("pan", colnames(df))]
n <- length(micr)

#write.csv(df, "C:/Users/faesslerd/Documents/allmerged.csv", row.names = FALSE)
# function that calculates sign: number of sig in vivo species, accuracy, expected accuracy, pval,  size: regression coef, rsq, pval 
vivsil <- function(vivomet, flux, data){
  
  logvi <- log(data[, vivomet])
  logvi[which(is.infinite(logvi))]=NA
  
  
  abun <- micr
  sigvivospec <- length(abun)
  
  flx <- data[, flux]
  
  if(all(as.numeric(data[, micr[1]]) %in% c(0,1))){
    lmV <- lapply(abun, function(x) coeftest(lm(logvi ~ data[, x] + data$Age + data$BMI + data$Stratification + data$Gender), vcov = vcovHC(lm(logvi ~ data[, x] + data$Age + data$BMI + data$Stratification + data$Gender), type="HC1"))[2,1])
    lmS <- lapply(abun, function(x) coeftest(lm(flx  ~ data[, x] + data$Age + data$BMI + data$Stratification + data$Gender), vcov = vcovHC(lm(flx ~ data[, x] + data$Age + data$BMI + data$Stratification + data$Gender), type="HC1"))[2,1]) 
  } else {
     lmV <- lapply(abun, function(x) coeftest(lm(logvi ~ scale(data[, x]) + data$Age + data$BMI + data$Stratification + data$Gender), vcov = vcovHC(lm(logvi ~ scale(data[, x]) + data$Age + data$BMI + data$Stratification + data$Gender), type="HC1"))[2,1])
     lmS <- lapply(abun, function(x) coeftest(lm(flx  ~ scale(data[, x]) + data$Age + data$BMI + data$Stratification + data$Gender), vcov = vcovHC(lm(flx ~ scale(data[, x]) + data$Age + data$BMI + data$Stratification + data$Gender), type="HC1"))[2,1]) 
  }
 
  lmV <- as.numeric(lmV)/1000
  lmS <- as.numeric(lmS)/1000
  
  linmod <- lm(lmV ~ lmS)
  
  coef <-  summary(linmod)$coef[2,1]
  confinter <- confint(linmod)[2,]
  pval <- summary(linmod)$coef[2,4]
  rsq <- summary(linmod)$r.squared
  
  if(sd(unlist(lmS))>0){
    ratertab <- xtabs(~ sign(lmS) + sign(lmV))
    if(length(ratertab) <=2){
      mat <- matrix(rep(0,4), nrow=2)
      mat[1,] <- ratertab[1,]
      ratertab <- mat
    }
    
    acc <- (ratertab[1,1]+ratertab[2,2])/sum(ratertab[,])
    expacc <- 1/sum(ratertab[,])^2*(sum(ratertab[1,])*sum(ratertab[,1])+sum(ratertab[2,])*sum(ratertab[,2]))
    pvalfish <- fisher.test(ratertab)$p.value
  }
  return(c(sigvivospec, acc, expacc, pvalfish, coef, confinter, rsq, pval))
}

### Vspatternplot
vspattern <- function(vivomet, flux, data, micr) {
  
  vi <- df[, vivomet]
  si <- data[, flux]
  
  logvi <- log(vi)
  logvi[which(is.infinite(logvi))]=NA
  
  lmV <- lapply(micr, function(x) coeftest(lm(logvi ~ scale(df[, x]) + df$Age + df$BMI + df$Stratification + df$Gender), vcov = vcovHC(lm(logvi ~ scale(df[, x]) + df$Age + df$BMI + df$Stratification + df$Gender), type="HC1"))[2,1]) 
  lmV <- as.numeric(lmV)
  
  lmS <- lapply(micr, function(x) coeftest(lm(si  ~ scale(df[, x]) + df$Age + df$BMI + df$Stratification + df$Gender), vcov = vcovHC(lm(si ~ scale(df[, x]) + df$Age + df$BMI + df$Stratification + df$Gender), type="HC1"))[2,1]) 
  lmS <- as.numeric(lmS)
  
  plot(lmS, lmV, main=vivomet, xlab="silico", ylab="vivo")
  linmod <- lm(lmV ~ lmS)
  
  abline(lm(lmV ~ lmS))
  print(paste(vivomet, "R^2: "))
  print(summary(linmod)$r.squared)
}

### Total effect abundance
metabolites <- agr$abbreviation[which(gsub("_", "", agr$keggId) %in% setdiff(ts1, c("C00695", "C00245")))]
keggIDs <- gsub("_", "", agr$keggId[which(gsub("_", "", agr$keggId) %in% setdiff(ts1, c("C00695", "C00245")))])

dframe3 <- data.frame(matrix(ncol=298, nrow = 54))
colnames(dframe3) <- c("keggID", "metabolite", gsub("$", "_s", colnames(df)[110:257]), gsub("$", "_v", colnames(df)[110:257]))
dframe3$metabolite <- metabolites
dframe3$keggID <- keggIDs
abun  <- micr
for(i in 1:length(metabolites)){
  print(i)
    vivomet <- dframe3$keggID[i]
    flxname <-gsub("^", "EX_", vivomet)
    vi <- df[, vivomet]
    si <- df[, flxname]
    
    logvi <- log(vi)
    logvi[which(is.infinite(logvi))]=NA
    
    lmS <- lapply(abun, function(x) coeftest(lm(si  ~ scale(df[, x]) + df$Age + df$BMI + df$Stratification + df$Gender), vcov = vcovHC(lm(si ~ scale(df[, x]) + df$Age + df$BMI + df$Stratification + df$Gender), type="HC1"))[2,1]) 
    lmS <- as.numeric(lmS)
    
    lmV <- lapply(abun, function(x) coeftest(lm(logvi ~ scale(df[, x]) + df$Age + df$BMI + df$Stratification + df$Gender), vcov = vcovHC(lm(logvi ~ scale(df[, x]) + df$Age + df$BMI + df$Stratification + df$Gender), type="HC1"))[2,1])
    lmV <- as.numeric(lmV)

    dframe3[i, 3:ncol(dframe3)] <- c(lmS,lmV)
}
# Total effect presence
dframe4 <- data.frame(matrix(ncol=298, nrow = 54))
colnames(dframe4) <- c("keggID", "metabolite", gsub("$", "_s", colnames(dfpres)[110:257]), gsub("$", "_v", colnames(dfpres)[110:257]))
dframe4$metabolite <- metabolites
dframe4$keggID <- keggIDs

for(i in 1:length(metabolites)){
  print(i)
  vivomet <- dframe4$keggID[i]
  flxname <-gsub("^", "EX_", vivomet)
  vi <- dfpres[, vivomet]
  si <- dfpres[, flxname]
  
  logvi <- log(vi)
  logvi[which(is.infinite(logvi))]=NA
  
  lmS <- lapply(abun, function(x) coeftest(lm(si  ~ dfpres[, x] + dfpres$Age + dfpres$BMI + dfpres$Stratification + dfpres$Gender), vcov = vcovHC(lm(si ~ dfpres[, x] + dfpres$Age + dfpres$BMI + dfpres$Stratification + dfpres$Gender), type="HC1"))[2,1]) 
  lmS <- as.numeric(lmS)
  
  lmV <- lapply(abun, function(x) coeftest(lm(logvi ~ dfpres[, x] + dfpres$Age + dfpres$BMI + dfpres$Stratification + dfpres$Gender), vcov = vcovHC(lm(logvi ~ dfpres[, x] + dfpres$Age + dfpres$BMI + dfpres$Stratification + dfpres$Gender), type="HC1"))[2,1])
  lmV <- as.numeric(lmV)
  
  dframe4[i, 3:ncol(dframe4)] <- c(lmS,lmV)
}

######################  all metabolites - species abundance
# check for all 54 metabolites
silmetlist <- vector(mode = "list", length = 55)

for(i in 2:55){
  silmetlist[[i]] <- data.frame(beta=double(),
                                pval=double(), 
                                rsq=double(),
                                pcons=double(),
                                pexp=double(),
                                cohensk=double(),
                                stringsAsFactors=FALSE)
  
  logvi <- log(df[, i])
  logvi[which(is.infinite(logvi))]=NA
  
  metname <- colnames(df)[i]
  abun <- micr
  origflx <- paste("EX_", colnames(df)[i], sep="")
  # permutationtest
  for (u in 1:5000){
    lmV <- lapply(abun, function(x) coeftest(lm(logvi ~ scale(df[, x]) + df$Age + df$BMI + df$Stratification + df$Gender), vcov = vcovHC(lm(logvi ~ scale(df[, x]) + df$Age + df$BMI + df$Stratification + df$Gender), type="HC1"))[2,1]) 
    lmV <- as.numeric(lmV)/1000
    shflx <- df[, origflx][shuffle(length(df[, origflx]))]
    
    lmS <- lapply(abun, function(x) coeftest(lm(shflx  ~ scale(df[, x]) + df$Age + df$BMI + df$Stratification + df$Gender), vcov = vcovHC(lm(shflx ~ scale(df[, x]) + df$Age + df$BMI + df$Stratification + df$Gender), type="HC1"))[2,1]) 
    lmS <- as.numeric(lmS)/1000
    
    linmod <- lm(lmV ~ lmS)
    
    coef <-  summary(linmod)$coef[2,1]
    pval <- summary(linmod)$coef[2,4]
    rsq <- summary(linmod)$r.squared
    
    silmetlist[[i]][u,1:3] <- c(coef, pval, rsq)
    
    ratertab <- xtabs(~ sign(lmS) + sign(lmV))
    if(length(ratertab) <=2){
      mat <- matrix(rep(0,4), nrow=2)
      mat[1,] <- ratertab[1,]
      ratertab <- mat
    }
    
    pcons <- (ratertab[1,1]+ratertab[2,2])/sum(ratertab[,])
    pexp <- 1/sum(ratertab[,])^2*(sum(ratertab[1,])*sum(ratertab[,1])+sum(ratertab[2,])*sum(ratertab[,2]))
    silmetlist[[i]][u,4] <- pcons
    silmetlist[[i]][u,5] <- pexp
    silmetlist[[i]][u,6] <- CohenKappa(ratertab)
  }
  print(i)
}

dfsilshuffl <- data.frame(keggID=character(),
                          NoSiViS=double(),
                          RSQ=double(),
                          RSQG=double(),
                          pcons=double(),
                          pexp=double(),
                          pG = double(),
                          stringsAsFactors=FALSE)

for(i in 2:55){
  logvi <- log(df[, i])
  logvi[which(is.infinite(logvi))]=NA
  origflx <- paste("EX_", colnames(df)[i], sep="")
  
  metname <- colnames(df)[i]
  abun <- micr
  flx <- df[, origflx]
  
  lmV <- lapply(abun, function(x) coeftest(lm(logvi ~ scale(df[, x]) + df$Age + df$BMI + df$Stratification + df$Gender), vcov = vcovHC(lm(logvi ~ scale(df[, x]) + df$Age + df$BMI + df$Stratification + df$Gender), type="HC1"))[2,1]) 
  lmV <- as.numeric(lmV)/1000
  
  lmS <- lapply(abun, function(x) coeftest(lm(flx  ~ scale(df[, x]) + df$Age + df$BMI + df$Stratification + df$Gender), vcov = vcovHC(lm(flx ~ scale(df[, x]) + df$Age + df$BMI + df$Stratification + df$Gender), type="HC1"))[2,1]) 
  lmS <- as.numeric(lmS)/1000
  
  linmod <- lm(lmV ~ lmS)
  
  coef <-  summary(linmod)$coef[2,1]
  pval <- summary(linmod)$coef[2,4]
  rsq <- summary(linmod)$r.squared

  ratertab <- xtabs(~ sign(lmS) + sign(lmV))
  if(length(ratertab) <=2){
    mat <- matrix(rep(0,4), nrow=2)
    mat[1,] <- ratertab[1,]
    ratertab <- mat
  }
  
  pcons <- (ratertab[1,1]+ratertab[2,2])/sum(ratertab[,])
  pexp <- 1/sum(ratertab[,])^2*(sum(ratertab[1,])*sum(ratertab[,1])+sum(ratertab[2,])*sum(ratertab[,2]))
  
  dfsilshuffl[i,] = c(origflx, length(abun) ,rsq, sum(silmetlist[[i]][, 3]>rsq)/5000, pcons, pexp, sum(abs(silmetlist[[i]][, 4]-silmetlist[[i]][, 5])>abs(pcons-pexp))/5000)
  print(i)
}

dfsilshuffl <- na.omit(dfsilshuffl)
rownames(dfsilshuffl) <- c()
sum(as.numeric(dfsilshuffl$RSQG)<100)

dfsilshuffl$keggID[(as.numeric(dfsilshuffl$RSQG)<100)]
unique(sort(as.numeric(dfsilshuffl$RSQG)))

dfsilshuffl$keggID <- gsub("EX_", "", dfsilshuffl$keggID)
keggJoh <- read.csv(file = "keggIDName.csv", sep=";", header=TRUE)
names(keggJoh) <- c("keggID", "name")
keggJoh[54,] <- c("C00329", "D-Glucosamine")

mslist <- merge(keggJoh, dfsilshuffl, by = "keggID")
#mslist[,3:8] <- sapply(mslist[, 3:8],as.numeric)
#write.xlsx(mslist, "C:/Users/faesslerd/Documents/mslistNormalizedAbundance.xlsx", row.names = FALSE)

######################  all metabolites species presence
dfpres <- df
dfpres[110:257] <- dfpres[110:257] %>% mutate_if(is.numeric, ~1 * (. > 0))

silmetlistpres <- vector(mode = "list", length = 55)
for(i in 2:55){
  silmetlistpres[[i]] <- data.frame(beta=double(),
                                    pval=double(), 
                                    rsq=double(),
                                    pcons=double(),
                                    pexp=double(),
                                    cohensk=double(),
                                    stringsAsFactors=FALSE)
  
  logvi <- log(dfpres[, i])
  logvi[which(is.infinite(logvi))]=NA
  
  metname <- colnames(dfpres)[i]
  abun <- micr
  origflx <- paste("EX_", colnames(dfpres)[i], sep="")
  # permutationtest
  for (u in 1:5000){
    lmV <- lapply(abun, function(x) coeftest(lm(logvi ~ dfpres[, x] + dfpres$Age + dfpres$BMI + dfpres$Stratification + dfpres$Gender), vcov = vcovHC(lm(logvi ~ dfpres[, x] + dfpres$Age + dfpres$BMI + dfpres$Stratification + dfpres$Gender), type="HC1"))[2,1]) 
    lmV <- as.numeric(lmV)/1000
    shflx <- dfpres[, origflx][shuffle(length(dfpres[, origflx]))]
    
    lmS <- lapply(abun, function(x) coeftest(lm(shflx  ~ dfpres[, x] + dfpres$Age + dfpres$BMI + dfpres$Stratification + dfpres$Gender), vcov = vcovHC(lm(shflx ~ dfpres[, x] + dfpres$Age + dfpres$BMI + dfpres$Stratification + dfpres$Gender), type="HC1"))[2,1]) 
    lmS <- as.numeric(lmS)/1000
    
    linmod <- lm(lmV ~ lmS)
    
    coef <-  summary(linmod)$coef[2,1]
    pval <- summary(linmod)$coef[2,4]
    rsq <- summary(linmod)$r.squared
    
    silmetlistpres[[i]][u,1:3] <- c(coef, pval, rsq)
    
    ratertab <- xtabs(~ sign(lmS) + sign(lmV))
    if(length(ratertab) <=2){
      mat <- matrix(rep(0,4), nrow=2)
      mat[1,] <- ratertab[1,]
      ratertab <- mat
    }
    
    pcons <- (ratertab[1,1]+ratertab[2,2])/sum(ratertab[,])
    pexp <- 1/sum(ratertab[,])^2*(sum(ratertab[1,])*sum(ratertab[,1])+sum(ratertab[2,])*sum(ratertab[,2]))
    silmetlistpres[[i]][u,4] <- pcons
    silmetlistpres[[i]][u,5] <- pexp
    silmetlistpres[[i]][u,6] <- CohenKappa(ratertab)
  }
  print(i)
}

dfsilshufflpres <- data.frame(keggID=character(),
                              NoSiViS=double(),
                              RSQ=double(),
                              RSQG=double(),
                              pcons=double(),
                              pexp=double(),
                              pG=double(),
                              stringsAsFactors=FALSE)

for(i in 2:55){
  logvi <- log(dfpres[, i])
  logvi[which(is.infinite(logvi))]=NA
  origflx <- paste("EX_", colnames(dfpres)[i], sep="")
  
  metname <- colnames(dfpres)[i]
  abun <- micr
  flx <- dfpres[, origflx]
  
  lmV <- lapply(abun, function(x) coeftest(lm(logvi ~ dfpres[, x] + dfpres$Age + dfpres$BMI + dfpres$Stratification + dfpres$Gender), vcov = vcovHC(lm(logvi ~ dfpres[, x] + dfpres$Age + dfpres$BMI + dfpres$Stratification + dfpres$Gender), type="HC1"))[2,1]) 
  lmV <- as.numeric(lmV)/1000
  
  lmS <- lapply(abun, function(x) coeftest(lm(flx  ~ dfpres[, x] + dfpres$Age + dfpres$BMI + dfpres$Stratification + dfpres$Gender), vcov = vcovHC(lm(flx ~ dfpres[, x] + dfpres$Age + dfpres$BMI + dfpres$Stratification + dfpres$Gender), type="HC1"))[2,1]) 
  lmS <- as.numeric(lmS)/1000
  
  linmod <- lm(lmV ~ lmS)
  
  coef <-  summary(linmod)$coef[2,1]
  pval <- summary(linmod)$coef[2,4]
  rsq <- summary(linmod)$r.squared

  ratertab <- xtabs(~ sign(lmS) + sign(lmV))
  if(length(ratertab) <=2){
    mat <- matrix(rep(0,4), nrow=2)
    mat[1,] <- ratertab[1,]
    ratertab <- mat
  }
  
  pcons <- (ratertab[1,1]+ratertab[2,2])/sum(ratertab[,])
  pexp <- 1/sum(ratertab[,])^2*(sum(ratertab[1,])*sum(ratertab[,1])+sum(ratertab[2,])*sum(ratertab[,2]))
  
  dfsilshufflpres[i,] = c(origflx, length(abun) ,rsq, sum(silmetlistpres[[i]][, 3]>rsq)/5000, pcons, pexp, sum(abs(silmetlistpres[[i]][, 4]-silmetlistpres[[i]][, 5])>abs(pcons-pexp))/5000)
  print(i)
}

dfsilshufflpres <- na.omit(dfsilshufflpres)
rownames(dfsilshufflpres) <- c()

sum(as.numeric(dfsilshufflpres$RSQG)<100)

dfsilshufflpres$keggID[(as.numeric(dfsilshufflpres$RSQG)<100)]
unique(sort(as.numeric(dfsilshufflpres$RSQG)))

dfsilshufflpres$keggID <- gsub("EX_", "", dfsilshufflpres$keggID)

mslistp <- merge(keggJoh, dfsilshufflpres, by = "keggID")
mslistp[,3:8] <- sapply(mslistp[, 3:8],as.numeric)
#write.xlsx(mslistp, "C:/Users/faesslerd/Documents/mslistPres.xlsx", row.names = FALSE)

### signsizetable abundance
metnames <- colnames(df)[2:55]
flxnames <- gsub("^", "EX_", metnames)

signsizetable <- data.frame(keggID=character(),
                            Name=double(),
                            NoSigVivSpec=double(),
                            Accuracy=double(),
                            ExpectedAcc=double(),
                            pvalsig=double(),
                            FDRsig=double(),
                            regrcoef=character(),
                            Rsq=double(),
                            pvalsize=double(),
                            FDRsize=double(),
                            stringsAsFactors=FALSE)


for(i in 1:length(metnames)){
  flxname <-gsub("^", "EX_", metnames[i])
  outp <- vivsil(metnames[i], flxname, df)
  signsizetable[i, 1] <- metnames[i]
  signsizetable[i,3:6] <- outp[1:4]
  signsizetable[i,8] <- paste(round(outp[5], digits=5), ' (', round(outp[6], digits=5), ',', round(outp[7], digits=5), ')', sep='')
  signsizetable[i,9] <- outp[8]
  signsizetable[i,10] <- outp[9]
  print(i)
}

signsizetable$FDRsig <- p.adjust(signsizetable$pvalsig, method = "fdr", n=length(signsizetable$pvalsig))
signsizetable$FDRsize <- p.adjust(signsizetable$pvalsize, method = "fdr", n=length(signsizetable$pvalsig))

### signsizetable presence
metnames <- colnames(dfpres)[2:55]
flxnames <- gsub("^", "EX_", metnames)

signsizetablep <- data.frame(keggID=character(),
                             Name=double(),
                             NoSigVivSpec=double(),
                             Accuracy=double(),
                             ExpectedAcc=double(),
                             pvalsig=double(),
                             FDRsig=double(),
                             regrcoef=character(),
                             Rsq=double(),
                             pvalsize=double(),
                             FDRsize=double(),
                             stringsAsFactors=FALSE)

for(i in 1:length(metnames)){
  flxname <-gsub("^", "EX_", metnames[i])
  outp <- vivsil(metnames[i], flxname, dfpres)
  signsizetablep[i, 1] <- metnames[i]
  signsizetablep[i,3:6] <- outp[1:4]
  signsizetablep[i,8] <- paste(round(outp[5], digits=5), ' (', round(outp[6], digits=5), ',', round(outp[7],digits=5), ')', sep='')
  signsizetablep[i,9] <- outp[8]
  signsizetablep[i,10] <- outp[9]
  print(i)
}

signsizetablep$FDRsig <- p.adjust(signsizetablep$pvalsig, method = "fdr", n=length(signsizetablep$pvalsig))
signsizetablep$FDRsize <- p.adjust(signsizetablep$pvalsize, method = "fdr", n=length(signsizetablep$pvalsize))


##############################################
# Put everything together abundance         ##
##############################################

g <- merge(x = signsizetable, y = mslist[ , c("keggID", "name")], by = "keggID", all.x=TRUE)
TableAbundance <- data.frame(keggID=g$keggID,
                            Name=g$name,
                            NoSigVivSpec=g$NoSigVivSpec,
                            Accuracy=g$Accuracy,
                            ExpectedAcc=g$ExpectedAcc,
                            pvalsig=g$pvalsig,
                            FDRsig=g$FDRsig,
                            pvalPermutation=as.numeric(mslist$pG),
                            regrcoef=g$regrcoef,
                            Rsq=g$Rsq,
                            pvalsize=g$pvalsize,
                            FDRsize=g$FDRsize,
                            pvalPermutation=mslist$RSQG,
                            stringsAsFactors=FALSE)

##############################################
# Put everything together presence          ##
##############################################

g2 <- merge(x = signsizetablep, y = mslistp[ , c("keggID", "name")], by = "keggID", all.x=TRUE)
TablePresence <- data.frame(keggID=g2$keggID,
                             Name=g2$name,
                             NoSigVivSpec=g2$NoSigVivSpec,
                             Accuracy=g2$Accuracy,
                             ExpectedAcc=g2$ExpectedAcc,
                             pvalsig=g2$pvalsig,
                             FDRsig=g2$FDRsig,
                             pvalPermutation=as.numeric(mslistp$pG),
                             regrcoef=g2$regrcoef,
                             Rsq=g2$Rsq,
                             pvalsize=g2$pvalsize,
                             FDRsize=g2$FDRsize,
                             pvalPermutation=mslistp$RSQG,
                             stringsAsFactors=FALSE)

##############################################
# Supplementary tables                      ##
##############################################

keggidname <- read.csv('keggidname.csv')

abun <- micr
InSilinVivStat <- data.frame(keggID=character(),
                             Name = character(),
                             Species=character(),
                             Regression.coefficientSI=character(),
                             Regression.coefficientVI=character(),
                             pvalInVivo=double(),
                             FDRinVivo=double(),
                             stringsAsFactors=FALSE)

for(i in 2:55){
  logvi <- log(df[, i])
  logvi[which(is.infinite(logvi))]=NA
  origflx <- paste("EX_", colnames(df)[i], sep="")
  
  metname <- colnames(df)[i]
  
  regrecoefs <- c()
  confind <- c()
  
  flx <- df[, origflx]
  
  lmV <- lapply(abun, function(x) coeftest(lm(logvi ~ scale(df[, x]) + df$Age + df$BMI + df$Stratification + df$Gender), vcov = vcovHC(lm(logvi ~ scale(df[, x]) + df$Age + df$BMI + df$Stratification + df$Gender), type="HC1")))
  lmS <- lapply(abun, function(x) coeftest(lm(flx  ~ scale(df[, x]) + df$Age + df$BMI + df$Stratification + df$Gender), vcov = vcovHC(lm(flx ~ scale(df[, x]) + df$Age + df$BMI + df$Stratification + df$Gender), type="HC1"))) 
  
  for(j in 1:length(abun))
  {
    regrcsil <- paste0(round(lmS[[j]][2,1], digits=6), '(', round(confint((lmS[[j]]))[2,1], digits=6), ',',  round(confint(lmS[[j]])[2,2], digits=6), ')')
    regrcviv <- paste0(round(lmV[[j]][2,1], digits=6), '(', round(confint((lmV[[j]]))[2,1], digits=6), ',',  round(confint(lmV[[j]])[2,2], digits=6), ')')
    keggidname$Name[which(metname == keggidname$keggID)]
    InSilinVivStat[(i-2)*148+j, ] <- c(metname, keggidname$name[which(metname == keggidname$keggID)], abun[j], regrcsil, regrcviv,  as.numeric(lmV[[j]][2,4]), NA)
  }
  print(i)
}

InSilinVivStat$FDRinVivo <- p.adjust(InSilinVivStat$pvalInVivo, method='fdr')
# SupplementaryTableS2
#write.csv(InSilinVivStat, "C:/Users/faesslerd/Documents/Table_SpeciesAbundance_new.csv", row.names=FALSE)

dfpres <- df
dfpres[110:257] <- dfpres[110:257] %>% mutate_if(is.numeric, ~1 * (. > 0))
abun <- micr

InSilinVivStat <- data.frame(keggID=character(),
                            Name = character(),
                            Species=character(),
                            Regression.coefficientSI=character(),
                            Regression.coefficientVI=character(),
                            pvalInVivo=double(),
                            FDRinVivo=double(),
                            stringsAsFactors=FALSE)


for(i in 2:55){
  logvi <- log(dfpres[, i])
  logvi[which(is.infinite(logvi))]=NA
  origflx <- paste("EX_", colnames(dfpres)[i], sep="")
  
  metname <- colnames(dfpres)[i]
  
  regrecoefs <- c()
  confind <- c()

  flx <- dfpres[, origflx]

  lmV <- lapply(abun, function(x) coeftest(lm(logvi ~ dfpres[, x] + dfpres$Age + dfpres$BMI + dfpres$Stratification + dfpres$Gender), vcov = vcovHC(lm(logvi ~ dfpres[, x] + dfpres$Age + dfpres$BMI + dfpres$Stratification + dfpres$Gender), type="HC1")))
  lmS <- lapply(abun, function(x) coeftest(lm(flx  ~ dfpres[, x] + dfpres$Age + dfpres$BMI + dfpres$Stratification + dfpres$Gender), vcov = vcovHC(lm(flx ~ dfpres[, x] + dfpres$Age + dfpres$BMI + dfpres$Stratification + dfpres$Gender), type="HC1")))

  for(j in 1:length(abun))
  {
    regrcsil <- paste0(round(lmS[[j]][2,1], digits=6), '(', round(confint((lmS[[j]]))[2,1], digits=6), ',',  round(confint(lmS[[j]])[2,2], digits=6), ')')
    regrcviv <- paste0(round(lmV[[j]][2,1], digits=6), '(', round(confint((lmV[[j]]))[2,1], digits=6), ',',  round(confint(lmV[[j]])[2,2], digits=6), ')')
    keggidname$Name[which(metname == keggidname$keggID)]
    InSilinVivStat[(i-2)*148+j, ] <- c(metname, keggidname$name[which(metname == keggidname$keggID)], abun[j], regrcsil, regrcviv,  as.numeric(lmV[[j]][2,4]), NA)
  }
  print(i)
}
  
InSilinVivStat$FDRinVivo <- p.adjust(InSilinVivStat$pvalInVivo, method='fdr')
# SupplementaryTableS1
#write.csv(InSilinVivStat, "C:/Users/faesslerd/Documents/Table_SpeciesPresence.csv", row.names=FALSE)

### plot
# nspecies = 148, p=2, k=1, R^2~F(1, 146)
rsq=as.numeric(unlist(silmetlist[[9]][3]))
dfplotbut <- data.frame(beta=as.numeric(unlist(silmetlist[[9]][1])), rsq=rsq, rsqt=rsq/(1-rsq)*146)

mean(dfplotbut$beta)

# nspecies = 148, p=2, k=1, R^2~F(1, 146)
rsq=as.numeric(unlist(silmetlist[[9]][3]))
dfplotbut <- data.frame(beta=as.numeric(unlist(silmetlist[[9]][1])), rsq=rsq, rsqt=rsq/(1-rsq)*146)

mean(dfplotbut$beta)
sum(dfplotbut$rsqt)*0.05
sort(dfplotbut$rsqt)[cumsum(sort(dfplotbut$rsqt)) > sum(dfplotbut$rsqt)*0.05]

#hist(dfplotbut$rsqt, freq=FALSE, xlim = c(0,150), breaks=50)
#curve(df(x, df1 = 1, df2 = 44), from = 0, to = 150, n = 1000, col= 'black', lwd=2, add = T)
cumsum(table(cut(dfplotbut$rsqt, breaks = seq.int(from = 0, to = 100, by = 1))))

quant <- qf(0.95, df1 = 1, df2 = 146, lower.tail = TRUE)
empquant <- quantile(dfplotbut$rsqt, probs=c(0.95))

#Regression coefficient in vivo butyrate
regrsilbut <- c()
regrcvivbut <- c()
logvi <- log(df[, 9])
logvi[which(is.infinite(logvi))]=NA
origflx <- paste("EX_", colnames(df)[9], sep="")

regrecoefs <- c()
confind <- c()

flx <- df[, origflx]

lmS <- lapply(abun, function(x) coeftest(lm(flx  ~ scale(df[, x]) + df$Age + df$BMI + df$Stratification + df$Gender), vcov = vcovHC(lm(flx ~ scale(df[, x]) + df$Age + df$BMI + df$Stratification + df$Gender), type="HC1"))) 
lmV <- lapply(abun, function(x) coeftest(lm(logvi ~ scale(df[, x]) + df$Age + df$BMI + df$Stratification + df$Gender), vcov = vcovHC(lm(logvi ~ scale(df[, x]) + df$Age + df$BMI + df$Stratification + df$Gender), type="HC1")))

for(j in 1:length(abun)){
  regrcvivbut <- c(regrcvivbut, lmV[[j]][2,1])
  regrsilbut <- c(regrsilbut, lmS[[j]][2,1])
}

lma <- lm(regrcvivbut~regrsilbut)
summary(lma)
k=length(lma$coefficients)-1
SSE=sum(lma$residuals**2)
n=length(lma$residuals)
RSE <- sqrt(SSE/(n-(1+k)))
sd_hat <- RSE/sqrt(sum((regrsilbut - mean(regrsilbut))^2))

#rm(df) for stat_function
plot1 <- ggplot(dfplotbut, aes(beta)) +
  geom_histogram(aes(y=after_stat(density)), binwidth=.0015, color="black", fill="lightblue") +
  theme_bw() +
  stat_function(fun = dnorm, args = list(mean = 0, sd = sd_hat), size=1) +
  theme(axis.text=element_text(colour="black"), text=element_text(size=15))

# Create the second plot
plot2 <- ggplot(dfplotbut, aes(x=rsqt)) +
  geom_histogram(aes(y=..density..), boundary=0, binwidth = 2.5, closed="left", color="black", fill="lightblue") +
  theme_bw() +
  theme(axis.text=element_text(colour="black"), text=element_text(size=15)) +
  xlab("transformed r-squared") +
  xlim(c(0,75)) +
  ylim(c(0,.21)) +
  stat_function(fun = df,
                geom = "line",
                size = 1,
                args = list(
                  df1 = 1,
                  df2 = 146
                )) +
  geom_vline(xintercept = quant, size=1, color='orange') +
  geom_vline(xintercept = empquant, size=1, color='red')


fig7 <- grid.arrange(plot1, plot2, ncol = 2)

ggsave(fig7, filename = '~/Figure_S1.tiff', width = 21, height = 6, scale=0.8, dpi = 300)

df$R <- rowSums(df[110:257]>0)

#Do:
#mean(df[df$Stratification == "CRC", "R"])
#sd(df[df$Stratification == "CRC", "R"])
#Do:
#t.test(df$Age~df$Stratification, var.equal = FALSE, two.sided=TRUE)
#fisher.test(table(df$Gender,df$Stratification))

# Reload the three commands of simres
simres[nrow(simres)+1, ] <- c("Metprod", colSums(simres[2:ncol(simres)]>0))
simres2 <- data.frame(sampName=colnames(simres)[2:ncol(simres)], Nmet=as.numeric(simres[392,2:ncol(simres)]))
simres2$sampName <-  gsub("^sample_", "", simres2$sampName)

dfn <- merge(df, simres2, by="sampName")
t.test(dfn$Nmet~dfn$Stratification, var.equal = FALSE, two.sided=TRUE)

rxnsprod <- data.table::fread("ModelStatistics.csv")

rxnsprod$sampName <- gsub("^sample_", "", rxnsprod$ModelIDs)
dfn <- merge(dfn, rxnsprod, by="sampName")

#t.test(dfn$Reactions~df$Stratification, var.equal = FALSE, two.sided=TRUE)
#sd(dfn[dfn$Stratification == "CRC", "Reactions"])