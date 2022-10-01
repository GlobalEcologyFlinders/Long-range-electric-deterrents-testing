## WHITE SHARK BOAT01 DETERRENT TRIAL ANALYSIS

## Remove everything
rm(list = ls())

## source functions & libraries
setwd("")
library(lme4)
library(boot)
library(Hmisc)
library(ggplot2)
library(plyr)
library(ggpubr)
library(dplyr)
library(tidyr)
library(tidyverse)
library(hrbrthemes)
library(viridis)


source("new_lmer_AIC_tables3.r")
source("r.squared.R")

## functions
delta.IC <- function(x) x - min(x) ## where x is a vector of an IC
weight.IC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dIC


# DISTANCE GLMM - stereo camera - no low passes ----------------------------------------------------
setwd("")

# import distance data
dist.dat <- read.csv("dist_data_boat01.csv", header=T, sep=",")
str(dist.dat)
dist.dat$trial <- as.integer(dist.dat$trial)
dist.dat$dist <- as.numeric(dist.dat$dist)
View(dist.dat)
dist.dat$trial <- ordered(dist.dat$trial)

# remove zero distances
dist.no0dist <- subset(dist.dat, dist > 0)
dist.no0dist$ldist <- log10(dist.no0dist$dist)
dist.no0dist$trialint <- as.integer(dist.no0dist$trial)
dist.no0dist$ID <- factor(dist.no0dist$ID)

## GLMM
# remove unknown sharks
dist.nounk <- subset(dist.no0dist, ID != "Unknown")


##Remove low passes
dist.EM.01 <- subset(dist.nounk, intent=="M" | intent=="H")

# model set
m1 <- "ldist ~ det + trip/trialint + (1|ID)"
m2 <- "ldist ~ trip/trialint + (1|ID)"
m3 <- "ldist ~ det + (1|ID)"
m4 <- "ldist ~ 1 + (1|ID)"

## Make model vector
mod.vec <- c(m1,m2,m3,m4)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]), data=dist.EM.01, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- round(as.numeric(logLik(fit)), 3)
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- round(r.squared(fit)$AIC, 3)
  Rm[i] <- round(100*r.squared(fit)$Marginal, 1) # marginal R-squared
  Rc[i] <- round(100*r.squared(fit)$Conditional, 1) # conditional R-squared
  print(i)
}

dAICc <- round(delta.IC(AICc.vec), 3)
wAICc <- round(weight.IC(dAICc), 3)

sumtable <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,10),round(Rm,4),round(Rc,4))
colnames(sumtable) <- c("model","k","LL","AICc","dAICc","wAICc","Rm","Rc")
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,5],decreasing=F),1:8]
summary.table

fit.sat <- lmer(as.formula(mod.vec[4]), data=dist.EM.01, na.action=na.omit)
summary(fit.sat)

UncertCoef(table(dist.EM.01$ID, dist.EM.01$trialint), direction = "row")

# DISTANCE GLMM - observer data - no passes 5m+ -------------------------------------------------

setwd("")

# import distance data
dist.dat.obs <- read.csv("dist_data_boat01_obs.csv", header=T, sep=",")
str(dist.dat.obs)
dist.dat.obs$trial <- as.integer(dist.dat.obs$trial)
str(dist.dat.obs)
dist.dat.obs$trial <- ordered(dist.dat.obs$trial)

# remove zero distances
dist.no0dist.obs <- subset(dist.dat.obs, dist > 0)
str(dist.no0dist.obs)
dist.no0dist.obs$ldist <- log10(dist.no0dist.obs$dist)
dist.no0dist.obs$trialint <- as.integer(dist.no0dist.obs$trial)
dist.no0dist.obs$ID <- factor(dist.no0dist.obs$ID)

## GLMM
# remove unknown sharks
dist.OBS.01 <- subset(dist.no0dist.obs, ID != "Unknown")
View(dist.OBS.01)

##Remove passes 5m+
dist.OBS.01 <- dist.OBS.01[dist.OBS.01$dist < 5.1, ]
View(dist.OBS.01)

# model set
m1 <- "ldist ~ det + trip/trialint + (1|ID)"
m2 <- "ldist ~ trip/trialint + (1|ID)"
m3 <- "ldist ~ det + (1|ID)"
m4 <- "ldist ~ 1 + (1|ID)"

## Make model vector
mod.vec <- c(m1,m2,m3,m4)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]), data=dist.OBS.01, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- round(as.numeric(logLik(fit)), 3)
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- round(r.squared(fit)$AIC, 3)
  Rm[i] <- round(100*r.squared(fit)$Marginal, 1) # marginal R-squared
  Rc[i] <- round(100*r.squared(fit)$Conditional, 1) # conditional R-squared
  print(i)
}

dAICc <- round(delta.IC(AICc.vec), 3)
wAICc <- round(weight.IC(dAICc), 3)

sumtable <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,10),round(Rm,4),round(Rc,4))
colnames(sumtable) <- c("model","k","LL","AICc","dAICc","wAICc","Rm","Rc")
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,5],decreasing=F),1:8]
summary.table

fit.sat <- lmer(as.formula(mod.vec[4]), data=dist.OBS.01, na.action=na.omit)
summary(fit.sat)

# DISTANCE GLM stereo camera data -----------------------------------------------------------------
# import data
setwd("")
dist.dat <- read.csv("dist_data_boat01.csv", header=T, sep=",")
str(dist.dat)
dist.dat$trial <- as.integer(dist.dat$trial)
dist.dat$dist <- as.numeric(dist.dat$dist)
View(dist.dat)
dist.dat$trial <- ordered(dist.dat$trial)

# remove zero distances
dist.no0dist <- subset(dist.dat, dist > 0)
dist.no0dist$ldist <- log10(dist.no0dist$dist)
dist.no0dist$trialint <- as.integer(dist.no0dist$trial)
dist.no0dist$ID <- factor(dist.no0dist$ID)

##Remove low passes
dist.EM.01 <- subset(dist.nounk, intent=="M" | intent=="H")
str(dist.EM.01)

# glm without ID
# model set
m1 <- "ldist ~ det + trip/trialint"
m2 <- "ldist ~ trip/trialint"
m3 <- "ldist ~ det"
m4 <- "ldist ~ 1"

## Make model vector
mod.vec <- c(m1,m2,m3,m4)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)
for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]), data=dist.EM.01, na.action=na.omit)
  assign(paste('fit',i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}
sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),]
summary.table

fit.sat <- glm(as.formula(mod.vec[2]), data=dist.EM.01, na.action=na.omit)
summary(fit.sat)


# DISTANCE GLM observer data -----------------------------------------------------------------
# import data
setwd("")

# import distance data
dist.dat.obs <- read.csv("dist_data_boat01_obs.csv", header=T, sep=",")
str(dist.dat.obs)
dist.dat.obs$trial <- as.integer(dist.dat.obs$trial)
str(dist.dat.obs)
dist.dat.obs$trial <- ordered(dist.dat.obs$trial)

# remove zero distances
dist.no0dist.obs <- subset(dist.dat.obs, dist > 0)
str(dist.no0dist.obs)
dist.no0dist.obs$ldist <- log10(dist.no0dist.obs$dist)
dist.no0dist.obs$trialint <- as.integer(dist.no0dist.obs$trial)
dist.no0dist.obs$ID <- factor(dist.no0dist.obs$ID)

##Remove passes 5m+
dist.OBS.01 <- dist.no0dist.obs[dist.no0dist.obs$dist < 5.1, ]
str(dist.OBS.01)
str(dist.EM.01)

# glm without ID
# model set
m1 <- "ldist ~ det + trip/trialint"
m2 <- "ldist ~ trip/trialint"
m3 <- "ldist ~ det"
m4 <- "ldist ~ 1"

## Make model vector
mod.vec <- c(m1,m2,m3,m4)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)
for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]), data=dist.OBS.01, na.action=na.omit)
  assign(paste('fit',i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}
sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),]
summary.table

fit.sat <- glm(as.formula(mod.vec[4]), data=dist.OBS.01, na.action=na.omit)
summary(fit.sat)

# Distance figures --------------------------------------------------------

#means
#stereo camera data means
View(dist.dat)
dist.EM.01 <- na.omit(dist.EM.01)
cdat <- ddply(dist.EM.01, "det", summarise, dist.mean=mean(dist))
cdat

#Observer data means
dist.OBS.01$dist <- dist.OBS.01$dist * 1000
View(dist.OBS.01)
dist.OBS.01 <- na.omit(dist.OBS.01)
cdat.obs <- ddply(dist.OBS.01, "det", summarise, dist.mean=mean(dist))
cdat.obs

cdat.obs$dist.mean / 1

#Stereo camera density distribution
p <- ggplot(dist.dat, aes(x=dist, fill=det, color=det)) +
  geom_density(alpha=0.05)+
  labs(x="Pass distance (mm)", y = "Density") +
  xlim(c(0,5000)) + 
  labs(color = "Deterrent") +
  geom_vline(data=cdat, aes(xintercept=dist.mean,  colour=det), linetype="dashed", size=1)
p + theme_classic()

#Observer data density distribution 
str(dist.OBS.01)
dist.obs <- ggplot(dist.OBS.01, aes(x=dist, fill=det, color=det)) +
  geom_density(alpha=0.05)+
  labs(x="Pass distance (mm)", y = "Density") +
  xlim(c(0,5000)) + 
  labs(color = "Deterrent") +
  geom_vline(data=cdat.obs, aes(xintercept=dist.mean,  colour=det), linetype="dashed", size=1) 
dist.obs + theme_classic()


### Event measure data density dist - no low passes
#means
#stereo camera data
View(dist.nounk.nolow)
dist.dat.EM.nolow <- na.omit(dist.nounk.nolow)
cdat.nolow <- ddply(dist.dat.EM.nolow, "det", summarise, dist.mean=mean(dist))
cdat.nolow

#Observer data
dist.dat.obs.naom$dist <- dist.dat.obs.naom$dist * 1000
View(dist.dat.obs.naom)


dist.dat.obs.naom.nolow <- na.omit(dist.nounk.obs.nolow)
cdat.obs.nolow <- ddply(dist.dat.obs.naom.nolow, "det", summarise, dist.mean=mean(dist))
cdat.obs.nolow
cdat.obs.nolow$dist.mean * 1000

#stereo camera data density distribution - no low passes
p <- ggplot(dist.nounk.nolow, aes(x=dist, fill=det, color=det)) +
  geom_density(alpha=0.05)+
  labs(x="Pass distance (mm)", y = "Density") +
  xlim(c(0,7500)) + 
  labs(color = "Deterrent") +
  geom_vline(data=cdat.nolow, aes(xintercept=dist.mean,  colour=det), linetype="dashed", size=1)
p + theme_classic()

# observer data density distribution 
View(dist.dat.obs)
dist.obs <- ggplot(dist.dat.obs.naom.nolow, aes(x=dist, fill=det, color=det)) +
  geom_density(alpha=0.05)+
  labs(x="Pass distance (mm)", y = "Density") +
  xlim(c(0,7.5)) + 
  labs(color = "Deterrent") +
  geom_vline(data=cdat.obs.nolow, aes(xintercept=dist.mean,  colour=det), linetype="dashed", size=1)
dist.obs + theme_classic()

# stereo camera vs observer data density dist 01 --------------------------------------------------
setwd("")
dist.data <- read.csv("01.dist.comb.csv")
dist.data$dist<-as.numeric(dist.data$dist)
str(dist.data)
em.obs.01 <- ggplot(dist.data, aes(x=dist, fill=data, color=data)) +
  geom_density(alpha=0.05) +
  labs(x="Pass distance (mm)", y = "Density") +
  xlim(c(0,5))  
em.obs.01 + theme_classic()


# PASSES GLMM - event measure data - no low passes --------------------------------------------------------
## Creating approach dataframe  
setwd("")
EM.dat.01 <- read.csv("Boat01.EM.dataset.csv", fill= TRUE, header=T, sep=",")
View(EM.dat.01)

EM.dat.01.subset <-select(EM.dat.01, trip, trial, Deterrent, Shark.ID, LOI)
View(EM.dat.01.subset)

EM.dat.01.subset.nolow <- subset(EM.dat.01.subset, LOI=="M" | LOI=="H")
View(EM.dat.01.subset.nolow)

EM.dat.01.subset.appr.nolow <- EM.dat.01.subset.nolow %>%
  group_by(trip, trial, Deterrent, Shark.ID)%>%
  summarise(count= n())
View(EM.dat.01.subset.appr.nolow)

write.csv(EM.dat.01.subset.appr.nolow,"", row.names = FALSE)

#model
setwd("")
A.dat.nolow <- read.csv("Appr.EM.01.nl.csv", header=T, sep=",")
A.dat.nolow$trial <- factor(A.dat.nolow$trial)
A.dat.nolow$trialint <- as.integer(A.dat.nolow$trial)

str(A.dat.nolow)

hist(A.dat.nolow$appr)
range(A.dat.nolow$appr, na.rm=T)
A.dat.nolow$lappr <- log10(scale((A.dat.nolow$appr), center=F, scale=T))

sum(A.dat.nolow$appr)

## GLMM
# remove unknown sharks
A.nounk.nolow <- subset(A.dat.nolow, ID != "Unknown")
A.nounk.nolow$ID <- factor(A.nounk.nolow$ID)
hist(A.nounk.nolow$appr)
A.nounk.nolow$lappr <- log10(scale((A.nounk.nolow$appr), center=F, scale=T))
hist((A.nounk.nolow$lappr))
str(A.nounk.nolow)

sum(A.nounk.nolow$appr)

# model set
m1 <- "lappr ~ det + trip/trialint + (1|ID)"
m2 <- "lappr ~ trip/trialint + (1|ID)"
m3 <- "lappr ~ det + (1|ID)"
m4 <- "lappr ~ 1 + (1|ID)"

## Make model vector
mod.vec <- c(m1,m2,m3,m4)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]), data=A.nounk.nolow, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- round(as.numeric(logLik(fit)), 3)
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- round(r.squared(fit)$AIC, 3)
  Rm[i] <- round(100*r.squared(fit)$Marginal, 1) # marginal R-squared
  Rc[i] <- round(100*r.squared(fit)$Conditional, 1) # conditional R-squared
  print(i)
}

dAICc <- round(delta.IC(AICc.vec), 3)
wAICc <- round(weight.IC(dAICc), 3)

sumtable <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,10),round(Rm,4),round(Rc,4))
colnames(sumtable) <- c("model","k","LL","AICc","dAICc","wAICc","Rm","Rc")
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,5],decreasing=F),1:8]
summary.table

fit.sat <- lmer(as.formula(mod.vec[3]), data=A.nounk.nolow, na.action=na.omit)
summary(fit.sat)

# PASSES GLMM - observer data - no passes 5m+ --------------------------------------------------
#OBS no passes 5m+ 
## Creating approach dataframe  
setwd("")
OBS.dat.01 <- read.csv("OBS.01.full.csv", fill= TRUE, header=T, sep=",")
str(OBS.dat.01)
View(OBS.dat.01)

##Remove passes 5m+
OBS.dat.01.nolow <- OBS.dat.01[OBS.dat.01$dist < 5.1, ]
View(OBS.dat.01.nolow)

OBS.dat.01.subset <-select(OBS.dat.01.nolow, trip, trial, det, ID)
OBS.dat.01.subset

OBS.dat.01.subset.appr <- OBS.dat.01.subset %>%
  group_by(trip, trial, det, ID)%>%
  summarise(count= n())
OBS.dat.01.subset.appr

write.csv(OBS.dat.01.subset.appr,"", row.names = FALSE)

#model
setwd("")
A.dat.obs <- read.csv("Appr.OBS.01.5m.csv", header=T, sep=",")
View(A.dat.obs)
A.dat.obs$trial <- factor(A.dat.obs$trial)
A.dat.obs$trialint <- as.integer(A.dat.obs$trial)

str(A.dat.obs)

hist(A.dat.obs$appr)
range(A.dat.obs$appr, na.rm=T)
A.dat.obs$lappr <- log10(scale((A.dat.obs$appr), center=F, scale=T))

## GLMM
# remove unknown sharks
A.nounk.obs <- subset(A.dat.obs, ID != "Unknown")
A.nounk.obs$ID <- factor(A.nounk.obs$ID)
hist(A.nounk.obs$appr)
A.nounk.obs$lappr <- log10(scale((A.nounk.obs$appr), center=F, scale=T))
hist((A.nounk.obs$lappr))

sum(A.nounk.obs$appr)

# model set
m1 <- "lappr ~ det + trip/trialint + (1|ID)"
m2 <- "lappr ~ trip/trialint + (1|ID)"
m3 <- "lappr ~ det + (1|ID)"
m4 <- "lappr ~ 1 + (1|ID)"

## Make model vector
mod.vec <- c(m1,m2,m3,m4)
## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]), data=A.nounk.obs, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- round(as.numeric(logLik(fit)), 3)
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- round(r.squared(fit)$AIC, 3)
  Rm[i] <- round(100*r.squared(fit)$Marginal, 1) # marginal R-squared
  Rc[i] <- round(100*r.squared(fit)$Conditional, 1) # conditional R-squared
  print(i)
}

dAICc <- round(delta.IC(AICc.vec), 3)
wAICc <- round(weight.IC(dAICc), 3)

sumtable <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,10),round(Rm,4),round(Rc,4))
colnames(sumtable) <- c("model","k","LL","AICc","dAICc","wAICc","Rm","Rc")
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,5],decreasing=F),1:8]
summary.table

fit.sat <- lmer(as.formula(mod.vec[4]), data=A.nounk.obs, na.action=na.omit)
summary(fit.sat)


# PASSES GLM stereo camera data -----------------------------------------------------------

#model
setwd("")
A.dat.nolow <- read.csv("Appr.EM.01.nl.csv", header=T, sep=",")
A.dat.nolow$trial <- factor(A.dat.nolow$trial)
A.dat.nolow$trialint <- as.integer(A.dat.nolow$trial)

str(A.dat.nolow)

hist(A.dat.nolow$appr)
range(A.dat.nolow$appr, na.rm=T)
A.dat.nolow$lappr <- log10(scale((A.dat.nolow$appr), center=F, scale=T))

sum(A.dat.nolow$appr)

# glm without ID
# model set
m1 <- "lappr ~ det + trip/trialint"
m2 <- "lappr ~ trip/trialint"
m3 <- "lappr ~ det"
m4 <- "lappr ~ 1"

## Make model vector
mod.vec <- c(m1,m2,m3,m4)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)
for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]), data=A.dat.nolow, na.action=na.omit)
  assign(paste('fit',i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}
sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),]
summary.table

fit.sat <- glm(as.formula(mod.vec[3]), data=A.dat.nolow, na.action=na.omit)
summary(fit.sat)

# PASSES GLM observer data -----------------------------------------------------------
#model
setwd("")
A.dat.obs <- read.csv("Appr.OBS.01.5m.csv", header=T, sep=",")
View(A.dat.obs)
A.dat.obs$trial <- factor(A.dat.obs$trial)
A.dat.obs$trialint <- as.integer(A.dat.obs$trial)

str(A.dat.obs)

hist(A.dat.obs$appr)
range(A.dat.obs$appr, na.rm=T)
A.dat.obs$lappr <- log10(scale((A.dat.obs$appr), center=F, scale=T))
hist((A.dat.obs$lappr))

sum(A.dat.obs$appr)

# glm without ID
# model set
m1 <- "lappr ~ det + trip/trialint"
m2 <- "lappr ~ trip/trialint"
m3 <- "lappr ~ det"
m4 <- "lappr ~ 1"

## Make model vector
mod.vec <- c(m1,m2,m3,m4)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)
for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]), data=A.dat.obs, na.action=na.omit)
  assign(paste('fit',i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}
sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),]
summary.table

fit.sat <- glm(as.formula(mod.vec[3]), data=A.dat.obs, na.action=na.omit)
summary(fit.sat)


# Passes figures ----------------------------------------------------------
### stereo camera data passes 
EM.passes.boxplot <- ggplot(data = A.nounk.nolow, aes(x=det, y=appr, fill=det)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Boat02" = "gray", "Control" = "white")) +
  #scale_fill_viridis(discrete = TRUE, alpha=.5) +
  geom_jitter(color="black", size=.8, alpha=.8) +
  xlab("Deterrent") +
  ylab("Number of passes") +
  labs(fill = "Deterrent") +
  theme_classic() 
EM.passes.boxplot

### observer data passes 
OBS.passes.boxplot <- ggplot(data = A.nounk.obs, aes(x=det, y=appr, fill=det)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Boat02" = "gray", "Control" = "white")) +
  #scale_fill_viridis(discrete = TRUE, alpha=.5) +
  geom_jitter(color="black", size=.8, alpha=.8) +
  xlab("Deterrent") +
  ylab("Number of passes") +
  labs(fill = "Deterrent") +
  theme_classic() 
OBS.passes.boxplot


# EAT GLMM event measure data -------------------------------------------------------------

#binomial model
# import data
setwd("")
E.dat <- read.csv("Binomial.EM.01.csv", header=T, sep=",")
#E.dat$trialset <- as.integer(E.dat$Trial_set)
E.dat$trial <- factor(E.dat$trial)
E.dat$trialint <- as.integer(E.dat$trial)

## GLMM
# remove unknown sharks
E.nounk <- subset(E.dat, ID != "Unknown")
E.dat$ID <- factor(E.dat$ID)
E.dat$eat <- factor(E.dat$eat)
str(E.dat)


# model set
m1 <- "eat ~ det + trip/trialint + (1|ID)"
m2 <- "eat ~ trip/trialint + (1|ID)"
m3 <- "eat ~ det + (1|ID)"
m4 <- "eat ~ 1 + (1|ID)"

## Make model vector
mod.vec <- c(m1,m2,m3,m4)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glmer(as.formula(mod.vec[i]), family=binomial(link="logit"),data=E.dat, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- round(as.numeric(logLik(fit)), 3)
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- round(r.squared(fit)$AIC, 3)
  Rm[i] <- round(100*r.squared(fit)$Marginal, 1) # marginal R-squared
  Rc[i] <- round(100*r.squared(fit)$Conditional, 1) # conditional R-squared
  print(i)
}

dAICc <- round(delta.IC(AICc.vec), 3)
wAICc <- round(weight.IC(dAICc), 3)

sumtable <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,10),round(Rm,4),round(Rc,4))
colnames(sumtable) <- c("model","k","LL","AICc","dAICc","wAICc","Rm","Rc")
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,5],decreasing=F),1:8]
summary.table

fit.sat <- glmer(as.formula(mod.vec[1]), family=binomial(link="logit"), data=E.dat, na.action=na.omit)
summary(fit.sat)


# EAT GLMM observer data --------------------------------------------------------
## Creating OBS eat data frame  
setwd("")
Eat.OBS.dat.01 <- read.csv("OBS.01.data.csv", fill= TRUE, header=T, sep=",")
str(Eat.OBS.dat.01)

Eat.OBS.dat.01.subset <-select(Eat.OBS.dat.01, trip, trial, det, ID, Bait.taken)
View(Eat.OBS.dat.01.subset)

OBS.dat.01.subset.eat <- Eat.OBS.dat.01.subset %>%
  group_by(trip, trial, det, Bait.taken, ID)%>%
  summarise(count= n())
OBS.dat.01.subset.eat
View(OBS.dat.01.subset.eat)

write.csv(OBS.dat.01.subset.eat,"", row.names = FALSE)

### observer data eat GLM
setwd("")
E.dat.obs <- read.csv("Binomial.data.OBS.01.csv", header=T, sep=",")
#E.dat$trialset <- as.integer(E.dat$Trial_set)
E.dat.obs$trial <- factor(E.dat.obs$trial)
E.dat.obs$trialint <- as.integer(E.dat.obs$trial)

## GLMM
# remove unknown sharks
E.dat.obs <- subset(E.dat.obs, ID != "Unknown")
E.dat.obs$ID <- factor(E.dat.obs$ID)
E.dat.obs$eat <- factor(E.dat.obs$eat)
str(E.dat.obs)

# model set
m1 <- "eat ~ det + trip/trialint + (1|ID)"
m2 <- "eat ~ trip/trialint + (1|ID)"
m3 <- "eat ~ det + (1|ID)"
m4 <- "eat ~ 1 + (1|ID)"

## Make model vector
mod.vec <- c(m1,m2,m3,m4)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glmer(as.formula(mod.vec[i]), family=binomial(link="logit"),data=E.dat.obs, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- round(as.numeric(logLik(fit)), 3)
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- round(r.squared(fit)$AIC, 3)
  Rm[i] <- round(100*r.squared(fit)$Marginal, 1) # marginal R-squared
  Rc[i] <- round(100*r.squared(fit)$Conditional, 1) # conditional R-squared
  print(i)
}

dAICc <- round(delta.IC(AICc.vec), 3)
wAICc <- round(weight.IC(dAICc), 3)

sumtable <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,10),round(Rm,4),round(Rc,4))
colnames(sumtable) <- c("model","k","LL","AICc","dAICc","wAICc","Rm","Rc")
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,5],decreasing=F),1:8]
summary.table

fit.sat <- glmer(as.formula(mod.vec[1]), family=binomial(link="logit"), data=E.dat.obs, na.action=na.omit)
summary(fit.sat)

# Eat GLM stereo camera data -----------------------------------------------------------------
# import data
setwd("")
E.dat <- read.csv("Binomial.EM.01.csv", header=T, sep=",")
#E.dat$trialset <- as.integer(E.dat$Trial_set)
E.dat$trial <- factor(E.dat$trial)
E.dat$trialint <- as.integer(E.dat$trial)

E.dat$eat <- factor(E.dat$eat)


# glm without ID
# model set
m1 <- "eat ~ det + trip/trialint"
m2 <- "eat ~ trip/trialint"
m3 <- "eat ~ det"
m4 <- "eat ~ 1"

## Make model vector
mod.vec <- c(m1,m2,m3,m4)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)
for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=binomial(link='logit'), data=E.dat, na.action=na.omit)
  assign(paste('fit',i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}
sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),]
summary.table

fit.sat <- glm(as.formula(mod.vec[1]), family=binomial(link="logit"), data=E.dat, na.action=na.omit)
summary(fit.sat)


# Eat GLM observer data -----------------------------------------------------------------
# import data
setwd("")
E.dat.obs <- read.csv("Binomial.data.OBS.01.csv", header=T, sep=",")
#E.dat$trialset <- as.integer(E.dat$Trial_set)
E.dat.obs$trial <- factor(E.dat.obs$trial)
E.dat.obs$trialint <- as.integer(E.dat.obs$trial)

str(E.dat.obs)
E.dat.obs$eat <- factor(E.dat.obs$eat)
View(E.dat.obs)

# glm without ID
# model set
m1 <- "eat ~ det + trip/trialint"
m2 <- "eat ~ trip/trialint"
m3 <- "eat ~ det"
m4 <- "eat ~ 1"

## Make model vector
mod.vec <- c(m1,m2,m3,m4)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)
for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=binomial(link='logit'), data=E.dat.obs, na.action=na.omit)
  assign(paste('fit',i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}
sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),]
summary.table

fit.sat <- glm(as.formula(mod.vec[1]), family=binomial(link="logit"), data=E.dat.obs, na.action=na.omit)
summary(fit.sat)

# Eat or not figures ------------------------------------------------------
#stereo camera figures
setwd("")
E.dat.bt <- read.csv("Binomial.EM.01.csv", header=T, sep=",")
str(E.dat.bt)

library(dplyr)
E.dat.bt$det <- factor(E.dat.bt$det)
E.dat.bt$eat <- factor(E.dat.bt$eat)
str(E.dat.bt)
Bait.taken. <- E.dat.bt %>%
  group_by(det, eat) %>%
  dplyr::summarise(count= n()) %>%
  mutate(percentage = count / sum(count))
Bait.taken.

is.num <- sapply(Bait.taken., is.numeric)
Bait.taken.[is.num] <- lapply(Bait.taken.[is.num], round, 2)

ggplot(Bait.taken., aes(x = det, y = percentage, fill = eat)) +
  geom_col(position = "fill") +
  xlab("Deterrent") +
  ylab("Proportion of trials") +
  geom_text(aes(label = percentage), position = position_stack(vjust = 0.5), colour = "Red") +
  labs(fill='Trial outcome') +
  theme_classic() + scale_fill_grey(start = 0.9, end = 0.1) 

#observer data figures
setwd("")
E.dat.bt.obs <- read.csv("Binomial.data.OBS.01.csv", header=T, sep=",")
str(E.dat.bt.obs)

library(dplyr)
E.dat.bt.obs$det <- factor(E.dat.bt.obs$det)
E.dat.bt.obs$eat <- factor(E.dat.bt.obs$eat)
str(E.dat.bt.obs)
Bait.taken.obs <- E.dat.bt.obs %>%
  group_by(det, eat) %>%
  dplyr::summarise(count= n()) %>%
  mutate(percentage = count / sum(count))
Bait.taken.obs

is.num <- sapply(Bait.taken.obs, is.numeric)
Bait.taken.obs[is.num] <- lapply(Bait.taken.obs[is.num], round, 2)

ggplot(Bait.taken.obs, aes(x = det, y = percentage, fill = eat)) +
  geom_col(position = "fill") +
  xlab("Deterrent") +
  ylab("Proportion of trials") +
  geom_text(aes(label = percentage), position = position_stack(vjust = 0.5), colour = "Red") +
  labs(fill='Trial outcome') +
  theme_classic() + scale_fill_grey(start = 0.9, end = 0.1) 

# Reaction data frame --------------------------------------------------
## Creating dataframe  
setwd("")
EM.dat.01 <- read.csv("EM.01.data.csv", fill= TRUE, header=T, sep=",")
str(EM.dat.01)

EM.dat.01.subset <-select(EM.dat.01, trip, trial, dist, AT, Deterrent, ID, react)
View(EM.dat.01.subset)

EM.react.01 <- EM.dat.01.subset %>%
  group_by(trip, trial, AT, Deterrent, ID, react)%>%
  summarise(count= n())
View(EM.react.01)

write.csv(EM.react.01,"", row.names = FALSE)

# REACTION stereo camera data GLMM --------------------------------------------------------
setwd("")
EM.react <- read.csv("EM.react.01.csv", header=T, sep=",")
EM.react$trial <- factor(EM.react$trial)
EM.react$trialint <- as.integer(EM.react$trial)

# remove unknown sharks
EM.react.nounk <- subset(EM.react, ID != "Unknown")
EM.react.nounk$ID <- factor(EM.react.nounk$ID)

# redo treating trialset as covariate (ignoring 'time' per se, but using trialset as temporal marker)
EM.react.nounk.tscov <- EM.react.nounk
EM.react.nounk.tscov$trial <- as.integer(EM.react.nounk.tscov$trial)

EM.react.nounk.tscov$Reaction2[EM.react.nounk.tscov$react == 'N' ] <- '0'
EM.react.nounk.tscov$Reaction2[EM.react.nounk.tscov$react == 'Y' ] <- '1'
EM.react.nounk.tscov$Reaction2<-as.numeric(EM.react.nounk.tscov$Reaction2)

str(EM.react.nounk.tscov)

sum(EM.react.nounk.tscov$Reaction2)

# model set
m1 <- "Reaction.2 ~ det + trip/trialint + (1|ID)"
m2 <- "Reaction.2 ~ det + trip/trialint + AT + (1|ID)"
m3 <- "Reaction.2 ~ trip/trialint + (1|ID)"
m4 <- "Reaction.2 ~ trip/trialint + AT + (1|ID)"
m5 <- "Reaction.2 ~ det + (1|ID)"
m6 <- "Reaction.2 ~ det + AT + (1|ID)"
m7 <- "Reaction.2 ~ det*AT + (1|ID)"
m8 <- "Reaction.2 ~ AT + (1|ID)"
m9 <- "Reaction.2 ~ 1 + (1|ID)"

## Make model vector
mod.vec <- c(m1,m2,m3,m4,m5,m6,m7,m8,m9)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glmer(as.formula(mod.vec[i]),family=binomial(link="logit"),data=EM.react.nounk.tscov, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- r.squared(fit)$AIC
  Rm[i] <- 100*r.squared(fit)$Marginal # marginal R-squared
  Rc[i] <- 100*r.squared(fit)$Conditional # conditional R-squared
  print(i)
}

dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)

sumtable <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,10),round(Rm,4),round(Rc,4))
colnames(sumtable) <- c("model","k","LL","AICc","dAICc","wAICc","Rc","Rm")
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,5],decreasing=F),1:8]
summary.table

fit.sat <- lmer(as.formula(mod.vec[9]), data=EM.react.nounk.tscov, na.action=na.omit)
summary(fit.sat)

plot(allEffects(fit.sat), ylab="Reaction")


# REACTION observer data GLMM -------------------------------------------------------
## Creating dataframe  
setwd("")
OBS.dat.01 <- read.csv("OBS.01.full.csv", fill= TRUE, header=T, sep=",")
str(OBS.dat.01)

OBS.dat.01.subset <-select(OBS.dat.01, trip, trial, dist, det, ID, reaction)
View(OBS.dat.01.subset)

OBS.react.01 <- OBS.dat.01.subset %>%
  group_by(trip, trial, det, ID, reaction)%>%
  summarise(count= n())
View(OBS.react.01)

write.csv(OBS.react.01,"", row.names = FALSE)

###
setwd("")
OBS.react <- read.csv("OBS.react.01.csv", header=T, sep=",")
OBS.react$trial <- factor(OBS.react$trial)
OBS.react$trialint <- as.integer(OBS.react$trial)

# remove unknown sharks
OBS.react.nounk <- subset(OBS.react, ID != "Unknown")
OBS.react.nounk$ID <- factor(OBS.react.nounk$ID)

# redo treating trialset as covariate (ignoring 'time' per se, but using trialset as temporal marker)
OBS.react.nounk.tscov <- OBS.react.nounk
OBS.react.nounk.tscov$trial <- as.integer(OBS.react.nounk.tscov$trial)

OBS.react.nounk.tscov$Reaction2[OBS.react.nounk.tscov$reaction == 'N' ] <- '0'
OBS.react.nounk.tscov$Reaction2[OBS.react.nounk.tscov$reaction == 'Y' ] <- '1'
OBS.react.nounk.tscov$Reaction2<-as.numeric(OBS.react.nounk.tscov$Reaction2)

str(OBS.react.nounk.tscov)

sum(OBS.react.nounk.tscov$Reaction2)

# model set
m1 <- "Reaction2 ~ det + trip/trialint + (1|ID)"
m2 <- "Reaction2 ~ trip/trialint + (1|ID)"
m3 <- "Reaction2 ~ det + (1|ID)"
m4 <- "Reaction2 ~ 1 + (1|ID)"

## Make model vector
mod.vec <- c(m1,m2,m3,m4)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glmer(as.formula(mod.vec[i]),family=binomial(link="logit"),data=OBS.react.nounk.tscov, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- r.squared(fit)$AIC
  Rm[i] <- 100*r.squared(fit)$Marginal # marginal R-squared
  Rc[i] <- 100*r.squared(fit)$Conditional # conditional R-squared
  print(i)
}

dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)

sumtable <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,10),round(Rm,4),round(Rc,4))
colnames(sumtable) <- c("model","k","LL","AICc","dAICc","wAICc","Rc","Rm")
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,5],decreasing=F),1:8]
summary.table

fit.sat <- lmer(as.formula(mod.vec[3]), data=OBS.react.nounk.tscov, na.action=na.omit)
summary(fit.sat)

# REACTION stereo camera data GLM ---------------------------------------------------------
setwd("")
EM.react <- read.csv("EM.react.01.csv", header=T, sep=",")
EM.react$trial <- factor(EM.react$trial)
EM.react$trialint <- as.integer(EM.react$trial)

# redo treating trialset as covariate (ignoring 'time' per se, but using trialset as temporal marker)
EM.react.tscov <- EM.react
EM.react.tscov$trial <- as.integer(EM.react.tscov$trial)

EM.react.tscov$Reaction2[EM.react.tscov$react == 'N' ] <- '0'
EM.react.tscov$Reaction2[EM.react.tscov$react == 'Y' ] <- '1'
EM.react.tscov$Reaction2<-as.numeric(EM.react.tscov$Reaction2)

str(EM.react.tscov)

sum(EM.react.tscov$Reaction2)

# model set
m1 <- "Reaction.2 ~ det + trip/trialint"
m2 <- "Reaction.2 ~ det + trip/trialint + AT"
m3 <- "Reaction.2 ~ trip/trialint"
m4 <- "Reaction.2 ~ trip/trialint + AT"
m5 <- "Reaction.2 ~ det"
m6 <- "Reaction.2 ~ det + AT"
m7 <- "Reaction.2 ~ det*AT"
m8 <- "Reaction.2 ~ AT"
m9 <- "Reaction.2 ~ 1"

## Make model vector
mod.vec <- c(m1,m2,m3,m4,m5,m6,m7,m8,m9)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)
for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=binomial(link='logit'), data=EM.react.tscov, na.action=na.omit)
  assign(paste('fit',i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}
sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),]
summary.table

fit.sat <- glm(as.formula(mod.vec[1]), family=binomial(link="logit"), data=EM.react.tscov, na.action=na.omit)
summary(fit.sat)


# REACTION observer data GLM --------------------------------------------------------
setwd("")
OBS.react <- read.csv("OBS.react.01.csv", header=T, sep=",")
OBS.react$trial <- factor(OBS.react$trial)
OBS.react$trialint <- as.integer(OBS.react$trial)

# redo treating trialset as covariate (ignoring 'time' per se, but using trialset as temporal marker)
OBS.react.tscov <- OBS.react
OBS.react.tscov$trial <- as.integer(OBS.react.tscov$trial)

OBS.react.tscov$Reaction2[OBS.react.tscov$react == 'N' ] <- '0'
OBS.react.tscov$Reaction2[OBS.react.tscov$react == 'Y' ] <- '1'
OBS.react.tscov$Reaction2<-as.numeric(OBS.react.tscov$Reaction2)

str(OBS.react.tscov)

sum(OBS.react.tscov$Reaction2)

# model set
m1 <- "Reaction2 ~ det + trip/trialint"
m2 <- "Reaction2 ~ trip/trialint"
m3 <- "Reaction2 ~ det"
m4 <- "Reaction2 ~ 1"

## Make model vector
mod.vec <- c(m1,m2,m3,m4)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)
for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=binomial(link='logit'), data=OBS.react.tscov, na.action=na.omit)
  assign(paste('fit',i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}
sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),]
summary.table

fit.sat <- glm(as.formula(mod.vec[3]), family=binomial(link="logit"), data=OBS.react.tscov, na.action=na.omit)
summary(fit.sat)


### Figure - likelihood of reaction
str(dist.dat)
ggplot(dist.dat, aes(x=dist, y=reaction)) +
  geom_point() +
  geom_smooth(method = "glm",
              method.args = list(family = "binomial"),col="black",
              se = FALSE) +
  facet_grid(det ~ .) +
  xlab("Distance to electrode (cm)") +
  ylab("Likelihood of reaction") +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### Figure - trial
dist.dat$trial <- as.numeric(dist.dat$trial)
str(dist.dat)
trial <- ggplot(dist.dat, aes(x=trial, y=reaction)) +
  geom_point() +
  geom_smooth(method = "glm",
              method.args = list(family = "binomial"),col="black",
              se = FALSE) +
  facet_grid(det ~ .) +
  xlab("Trial") +
  ylab("Likelihood of reaction") +
  theme_bw()
trial + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# V2 REACTION GLMM  --------------------------------------------------------------

## 
setwd("")

dist.dat <- read.csv("Dist.data.noV.csv", header=T, sep=",")
#dist.dat<-read.csv(file.choose())
na.omit(dist.dat)
#dist.dat$trialset <- ordered(dist.dat$Trial_set)
dist.dat$trial <- ordered(dist.dat$trial)

str(dist.dat)
## GLMM
# remove unknown sharks
dist.nounk <- subset(dist.dat, ID != "Unknown")
dist.nounk$ID <- factor(dist.nounk$ID)


# redo treating trialset as covariate (ignoring 'time' per se, but using trialset as temporal marker)
dist.nounk.tscov <- dist.nounk
#dist.nounk.tscov$trial <- as.integer(dist.nounk.tscov$trial)
dist.nounk.tscov$trialint <- as.integer(dist.nounk.tscov$trial)

dist.nounk.tscov$Reaction2[dist.nounk.tscov$Reaction == 'N' ] <- '0'
dist.nounk.tscov$Reaction2[dist.nounk.tscov$Reaction == 'Y' ] <- '1'
dist.nounk.tscov$Reaction2<-as.numeric(dist.nounk.tscov$Reaction2)
dist.nounk.tscov$dist<-as.numeric(dist.nounk.tscov$dist)

str(dist.nounk.tscov)
View(dist.nounk.tscov)

# model set
m1 <- "Reaction2 ~ dist + (1|ID)"
m2 <- "Reaction2 ~ det + (1|ID)"
m3 <- "Reaction2 ~ det + dist + (1|ID)"
m4 <- "Reaction2 ~ trip/trialint + dist + (1|ID)"
m5 <- "Reaction2 ~ det + det*App_type + dist + (1|ID) "
m6 <- "Reaction2 ~ trip/trialint + det*App_type + dist + (1|ID) "
m7 <- "Reaction2 ~ det + trip/trialint + dist + (1|ID) "
m8 <- "Reaction2 ~ det + trip/trialint + det*App_type + dist + (1|ID) "
m9 <- "Reaction2 ~ det + trip/trialint + det*trip/trialint + dist + (1|ID) "
m10 <- "Reaction2 ~ det + trip/trialint + det*trip/trialint + det*App_type + dist + (1|ID) "
m11 <- "Reaction2 ~ 1 + (1|ID)"


## Make model vector
mod.vec <- c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glmer(as.formula(mod.vec[i]), data=dist.nounk.tscov, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- round(as.numeric(logLik(fit)), 3)
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- round(r.squared(fit)$AIC, 3)
  Rm[i] <- round(100*r.squared(fit)$Marginal, 1) # marginal R-squared
  Rc[i] <- round(100*r.squared(fit)$Conditional, 1) # conditional R-squared
  print(i)
}

dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)

sumtable <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,10),round(Rm,4),round(Rc,4))
colnames(sumtable) <- c("model","k","LL","AICc","dAICc","wAICc","Rc","Rm")
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,5],decreasing=F),1:8]
summary.table

fit.sat <- glmer(as.formula(mod.vec[10]), data=dist.nounk.tscov, na.action=na.omit)
summary(fit.sat)

# V2 REACTION Stereo camera GLM ---------------------------------------------------------
setwd("")
dist.dat <- read.csv("Dist.data.noV.csv", header=T, sep=",")
na.omit(dist.dat)
#dist.dat$trialset <- ordered(dist.dat$Trial_set)
dist.dat$trial <- ordered(dist.dat$trial)
dist.dat$ID <- factor(dist.dat$ID)

# redo treating trialset as covariate (ignoring 'time' per se, but using trialset as temporal marker)
dist.dat.tscov <- dist.dat
dist.dat.tscov$trialint <- as.integer(dist.dat.tscov$trial)

dist.dat.tscov$Reaction2[dist.dat.tscov$Reaction == 'N' ] <- '0'
dist.dat.tscov$Reaction2[dist.dat.tscov$Reaction == 'Y' ] <- '1'
dist.dat.tscov$Reaction2<-as.numeric(dist.dat.tscov$Reaction2) 
dist.dat.tscov$dist<-as.numeric(dist.dat.tscov$dist)

str(dist.dat.tscov)

# model set
m1 <- "Reaction2 ~ dist"
m2 <- "Reaction2 ~ det"
m3 <- "Reaction2 ~ det + dist"
m4 <- "Reaction2 ~ trip/trialint + dist"
m5 <- "Reaction2 ~ det + det*App_type + dist"
m6 <- "Reaction2 ~ trip/trialint + det*App_type + dist"
m7 <- "Reaction2 ~ det + trip/trialint + dist "
m8 <- "Reaction2 ~ det + trip/trialint + det*App_type + dist"
m9 <- "Reaction2 ~ det + trip/trialint + det*trip/trialint + dist"
m10 <- "Reaction2 ~ det + trip/trialint + det*trip/trialint + det*App_type + dist"
m11 <- "Reaction2 ~ 1"

## Make model vector
mod.vec <- c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)
for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=binomial(link='logit'), data=dist.dat.tscov, na.action=na.omit)
  assign(paste('fit',i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}
sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),]
summary.table

fit.sat <- glm(as.formula(mod.vec[10]), family=binomial(link="logit"), data=dist.dat.tscov, na.action=na.omit)
summary(fit.sat)



# Reaction time series ----------------------------------------------------------

### react by dist

boat01.dist.dat <- dist.nounk.tscov 
str(boat01.dist.dat)
boat01.dist.dat$det <- factor(boat01.dist.dat$det)
str(boat01.dist.dat)

React01 <- ggplot(boat01.dist.dat, aes(x=dist, y=Reaction2)) +
  geom_point(aes(colour = factor(ID)), size = 2) +
  geom_smooth(method = "glm",
              method.args = list(family = "binomial"),col="black",
              se = FALSE) +
  facet_grid(det ~ .) +
  xlab("distance to bait (cm)") +
  ylab("likelihood of reaction") +
  xlim(0, 7500) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position="none") +
  theme(text = element_text(size=35))
React01 


### React by trial

boat01.dist.dat <- dist.nounk.tscov 
str(boat01.dist.dat)
boat01.dist.dat$det <- factor(boat01.dist.dat$det)
str(boat01.dist.dat)

React01.trial <- ggplot(boat01.dist.dat, aes(x=trialint, y=Reaction2)) +
  geom_point(aes(colour = factor(ID)), size = 3.5) +
  geom_smooth(method = "glm",
              method.args = list(family = "binomial"),col="black",
              se = FALSE) +
  facet_grid(det ~ .) +
  xlab("trial") +
  ylab("likelihood of reaction") +
  scale_x_continuous(limits = c(0, 65), breaks = seq(0, 65, by = 5)) + 
  geom_vline(xintercept=c(2.5, 4.5, 15.5, 65), linetype='dashed') +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position="none") +
  theme(text = element_text(size=20))
React01.trial 


# Bait taken time series ----------------------------------------------------------
setwd("")
E.dat01 <- read.csv("BT01.csv", header=T, sep=",")
E.dat01$trial <- factor(E.dat$trial)
E.dat01$trialint <- as.integer(E.dat$trial)

str(E.dat01)


# remove unknown sharks
E.nounk01 <- subset(E.dat, ID != "Unknown")
E.nounk01$ID <- factor(E.nounk01$ID)
str(E.nounk01)


BT01 <- ggplot(E.dat01, aes(x=trialint, y=eat)) +
  geom_point(aes(colour = factor(ID)), size = 3.5) +
  geom_vline(xintercept=c(2.5, 4.5, 15.5, 65), linetype='dashed') +
  geom_smooth(method = "glm",
              method.args = list(family = "binomial"),col="black",
              se = FALSE) +
  facet_grid(det ~ .) +
  xlab("trial") +
  ylab("likelihood of bait taken") +
  scale_x_continuous(limits = c(0, 65), breaks = seq(0, 65, by = 5)) + 
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position="none") +
  theme(text = element_text(size=35))
BT01


# Distance time series ----------------------------------------------------
setwd("")

# import distance data
dist.dat <- read.csv("dist_data_boat01.csv", header=T, sep=",")
str(dist.dat)
dist.dat$trial <- as.integer(dist.dat$trial)
dist.dat$dist <- as.numeric(dist.dat$dist)
str(dist.dat)
dist.dat$trial <- ordered(dist.dat$trial)

# remove zero distances
dist.no0dist <- subset(dist.dat, dist > 0)
dist.no0dist$ldist <- log10(dist.no0dist$dist)
dist.no0dist$trialint <- as.integer(dist.no0dist$trial)
dist.no0dist$ID <- factor(dist.no0dist$ID)

## GLMM
# remove unknown sharks
dist.nounk <- subset(dist.no0dist, ID != "Unknown")


##Remove low passes
dist.EM.01 <- subset(dist.nounk, intent=="M" | intent=="H")
str(dist.EM.01)

box <- ggplot(dist.EM.01, aes(x=trial, y=dist, fill=det)) + 
  labs(x="Trial", y = "Pass distance (mm)") +
  geom_boxplot(trim=FALSE) +
  guides(fill=guide_legend(title="Deterrent")) 
box + theme_classic() + stat_summary(fun=mean, geom="point", shape=23, size=1)


###########

str(dist.EM.01)

BT01 <- ggplot(dist.EM.01, aes(x=trialint, y=dist)) +
  geom_boxplot(aes(colour = factor(ID)), size = .5) +
  geom_vline(xintercept=c(2.5, 4.5, 15.5, 65), linetype='dashed') +
  geom_smooth(method = "glm",
              method.args = list(family = "binomial"),col="black",
              se = FALSE) +
  facet_grid(det ~ .) +
  xlab("trial") +
  ylab("distance (mm)") +
  scale_x_continuous(limits = c(0, 65), breaks = seq(0, 65, by = 5)) + 
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position="none") +
  theme(text = element_text(size=20))
BT01


######## Dist boxplot

str(dist.EM.01)
dist.EM.01$trial <- character(dist.EM.01$trial)
as.character(dist.EM.01(trial))
str(dist.EM.01)


dist.box01 <- ggplot(dist.EM.01, aes(x=factor(trial), y=dist)) +
  geom_boxplot(aes(colour = factor(ID)), size = .5) +
  geom_vline(xintercept=c(2.5, 4.5, 15.5, 65), linetype='dashed') +
  geom_smooth(aes(as.numeric(trial), dist), col="black",
              se = FALSE) +
  facet_grid(det ~ .) +
  scale_y_continuous(limits = c(0, 7500)) +
  xlab("") +
  ylab("distance (mm)") +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position="none") +
  theme(text = element_text(size=20)) +
  guides(x = guide_axis(angle = 90)) 
dist.box01 

str(dist.EM.01)

# dist scatter
dist.EM.01$trial <- as.numeric(dist.EM.01$trial)
str(dist.EM.01)
distscatter01 <- ggplot(dist.EM.01, aes(x=trial, y=dist)) +
  geom_point(aes(colour = factor(ID)), size = .5) +
  geom_vline(xintercept=c(2.5, 4.5, 15.5, 65), linetype='dashed') +
  geom_smooth (col="black",
               se = FALSE) +
  facet_grid(det ~ .) +
  xlab("") +
  ylab("Distance (mm)") +
  scale_x_continuous(limits = c(0, 65), breaks = seq(0, 65, by = 5)) + 
  scale_y_continuous(limits = c(0, 7500)) + 
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position="none") +
  theme(text = element_text(size=20)) 
distscatter01

