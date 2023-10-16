## InvaPact simulation
## Corey Bradshaw
## October 2023
## Vers-Pont-du-Gard

## load necessary libraries
library(jtools)
library(ggplot2)
library(ggpubr)

## functions
# beta distribution shape parameter estimator
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

## impact score vector (follows EICAT criteria)
impactScore <- seq(0.2,1,0.2) # (0.2=minimal, 0.4=minor, 0.6=moderate, 0.8=major, 1.0=massive)

## confidence score vector # (follows EICAT criteria)
confScore <- c(0.25,0.5,1) # (0.25=low, 0.5=medium, 1.0=high)


####################################
## define regional characteristics
####################################

## region 1
## characteristics: many invasive species, many assessments, high confidence, high impacts
sppR1 <- 200 # number of invasive species
assess.limR1 <- c(10,80) # range of number of assessments per species
impact.wR1 <- c(0.2,0.3,0.8,0.9,0.95) # impact probability weights
conf.wR1 <- c(0.5,0.8,0.9) # confidence probability weights
natsppaffR1 <- 500 # native species affected

mn.omega.sppR1 <- sd.omega.sppR1 <- rep(NA, sppR1)
for (s in 1:sppR1) {
  Nassess.spp <- round(runif(1,min=assess.limR1[1], max=assess.limR1[2]), 0)
  impactspp <- sample(impactScore, Nassess.spp, replace=T, prob=impact.wR1)
  confspp <- sample(confScore, Nassess.spp, replace=T, prob=conf.wR1)
  mn.omega.sppR1[s] <- weighted.mean(impactspp, w=confspp)     
  sd.omega.sppR1[s] <- ifelse(length(unique(impactspp))==1 | is.na(wtd.sd(impactspp, w=confspp))==T |
                                is.infinite(wtd.sd(impactspp, w=confspp))==T, 0, wtd.sd(impactspp, w=confspp))
}
par(mfrow=c(3,1))
hist(mn.omega.sppR1, main="Region 1", xlim=c(0,1))
abline(v=median(mn.omega.sppR1, na.rm=T), lty=2)


## region 2
## characteristics: many invasive species, few assessments, low confidence, high impacts
## invasive species
sppR2 <- 250 # number of species
assess.limR2 <- c(1,16) # range of number of assessments per species
impact.wR2 <- c(0.2,0.3,0.8,0.9,0.95) # impact probability weights
conf.wR2 <- c(0.95,0.5,0.2) # confidence probability weights
natsppaffR2 <- 600 # native species affected

mn.omega.sppR2 <- sd.omega.sppR2 <- rep(NA, sppR2)
for (s in 1:sppR2) {
  Nassess.spp <- round(runif(1,min=assess.limR2[1], max=assess.limR2[2]), 0)
  impactspp <- sample(impactScore, Nassess.spp, replace=T, prob=impact.wR2)
  confspp <- sample(confScore, Nassess.spp, replace=T, prob=conf.wR2)
  mn.omega.sppR2[s] <- weighted.mean(impactspp, w=confspp)     
  sd.omega.sppR2[s] <- ifelse(length(unique(impactspp))==1 | is.na(wtd.sd(impactspp, w=confspp))==T |
                                is.infinite(wtd.sd(impactspp, w=confspp))==T, 0, wtd.sd(impactspp, w=confspp))
}
hist(mn.omega.sppR2, main="Region 2", xlim=c(0,1))
abline(v=median(mn.omega.sppR2, na.rm=T), lty=2)


## region 3
## characteristics: few invasive species, few assessments, low confidence, low impacts
## invasive species
sppR3 <- 50 # number of species
assess.limR3 <- c(2,8) # range of number of assessments per species
impact.wR3 <- c(0.95,0.95,0.5,0.05,0.05) # impact probability weights
conf.wR3 <- c(0.95,0.5,0.05) # confidence probability weights
natsppaffR3 <- 200 # native species affected

mn.omega.sppR3 <- sd.omega.sppR3 <- rep(NA, sppR3)
for (s in 1:sppR3) {
  Nassess.spp <- round(runif(1,min=assess.limR3[1], max=assess.limR3[2]), 0)
  impactspp <- sample(impactScore, Nassess.spp, replace=T, prob=impact.wR3)
  confspp <- sample(confScore, Nassess.spp, replace=T, prob=conf.wR3)
  mn.omega.sppR3[s] <- weighted.mean(impactspp, w=confspp) 
  sd.omega.sppR3[s] <- ifelse(length(unique(impactspp))==1 | is.na(wtd.sd(impactspp, w=confspp))==T |
                                 is.infinite(wtd.sd(impactspp, w=confspp))==T, 0, wtd.sd(impactspp, w=confspp))
}
hist(mn.omega.sppR2, main="Region 3", xlim=c(0,1))
abline(v=median(mn.omega.sppR3, na.rm=T), lty=2)
par(mfrow=c(1,1))



#########################################################################################
## resample to create per-grouping averages (i.e., across invasive species per region)
#########################################################################################
reg.iter <- 10000
omegaR1md <- omegaR2md <- omegaR3md <- rep(NA,reg.iter)
for (r in 1:reg.iter) {
  # region 1
  bparams1 <- estBetaParams(mn.omega.sppR1, sd.omega.sppR1^2)$alpha
  bparams2 <- estBetaParams(mn.omega.sppR1, sd.omega.sppR1^2)$beta
  omegaR1 <- ifelse(is.infinite(bparams1)==T, mn.omega.sppR1, rbeta(sppR1, bparams1, bparams2))
  omegaR1md[r] <- median(omegaR1, na.rm=T)
  
  # region 2
  bparams1 <- estBetaParams(mn.omega.sppR2, sd.omega.sppR2^2)$alpha
  bparams2 <- estBetaParams(mn.omega.sppR2, sd.omega.sppR2^2)$beta
  omegaR2 <- ifelse(is.infinite(bparams1)==T, mn.omega.sppR2, rbeta(sppR2, bparams1, bparams2))
  omegaR2md[r] <- median(omegaR2, na.rm=T)

  # region 3
  bparams1 <- estBetaParams(mn.omega.sppR3, sd.omega.sppR3^2)$alpha
  bparams2 <- estBetaParams(mn.omega.sppR3, sd.omega.sppR3^2)$beta
  omegaR3 <- ifelse(is.infinite(bparams1)==T, mn.omega.sppR3, rbeta(sppR3, bparams1, bparams2))
  omegaR3md[r] <- median(omegaR3, na.rm=T)
  
}

# region 1
omegaR1globmd <- median(omegaR1md)
omegaR1globup <- quantile(omegaR1md, probs=0.975)
omegaR1globlo <- quantile(omegaR1md, probs=0.025)
print(round(c(omegaR1globlo,omegaR1globmd, omegaR1globup), 3))
      
# region 2
omegaR2globmd <- median(omegaR2md)
omegaR2globup <- quantile(omegaR2md, probs=0.975)
omegaR2globlo <- quantile(omegaR2md, probs=0.025)
print(round(c(omegaR2globlo,omegaR2globmd, omegaR2globup), 3))

# region 3
omegaR3globmd <- median(omegaR3md)
omegaR3globup <- quantile(omegaR3md, probs=0.975)
omegaR3globlo <- quantile(omegaR3md, probs=0.025)
print(round(c(omegaR3globlo,omegaR3globmd, omegaR3globup), 3))

## build output dataframe
regComp.dat <- data.frame(region=c("1","2","3"), omega=c(omegaR1globmd,omegaR2globmd,omegaR3globmd),
                          omegaup=c(omegaR1globup,omegaR2globup,omegaR3globup),
                          omegalo=c(omegaR1globlo,omegaR2globlo,omegaR3globlo))
regComp.dat

## plot
A <- ggplot() + 
  geom_errorbar(data=regComp.dat, mapping=aes(x=region, ymin=omegalo, ymax=omegaup), width=0.2, size=1, color="blue") + 
  geom_point(data=regComp.dat, mapping=aes(x=region, y=omega), size=4, shape=21, fill="white")
A


########################################################################################
## omega per invasive species assessed (average impact magnitude per invasive species)
########################################################################################
regCompMnOmega.dat <- regComp.dat
regCompMnOmega.dat[1,2:4] <- regComp.dat[1,2:4]/sppR1
regCompMnOmega.dat[2,2:4] <- regComp.dat[2,2:4]/sppR2
regCompMnOmega.dat[3,2:4] <- regComp.dat[3,2:4]/sppR3
colnames(regCompMnOmega.dat)[2] <- "MnOmega"
regCompMnOmega.dat

## plot
B <- ggplot() + 
  geom_errorbar(data=regCompMnOmega.dat, mapping=aes(x=region, ymin=omegaup, ymax=omegalo), width=0.2, size=1, color="blue") + 
  geom_point(data=regCompMnOmega.dat, mapping=aes(x=region, y=MnOmega), size=4, shape=21, fill="white") +
  labs(x = "region", y = "mean omega/invasive sp") 
B


#################################################################################################
## standardise by native species affected (native species-weighted per-invasive species impact)
#################################################################################################
regCompInvaPacts.dat <- regCompMnOmega.dat
regCompInvaPacts.dat[1,2:4] <- natsppaffR1*regCompMnOmega.dat[1,2:4]
regCompInvaPacts.dat[2,2:4] <- natsppaffR2*regCompMnOmega.dat[2,2:4]
regCompInvaPacts.dat[3,2:4] <- natsppaffR3*regCompMnOmega.dat[3,2:4]
colnames(regCompInvaPacts.dat)[2] <- "InvaPacts"
regCompInvaPacts.dat

## plot
C <- ggplot() + 
  geom_errorbar(data=regCompInvaPacts.dat, mapping=aes(x=region, ymin=omegaup, ymax=omegalo), width=0.2, size=1, color="blue") + 
  geom_point(data=regCompInvaPacts.dat, mapping=aes(x=region, y=InvaPacts), size=4, shape=21, fill="white") +
  labs(x = "region", y = "affected native species omega/invasive sp") 
C


####################################
## plot all three indices together
####################################
ggarrange(A, B, C,
          ncol=3, nrow=1)


###########
## export
###########

write.table(regComp.dat, "regComp.csv", row.names=F, sep=",")
write.table(regCompMnOmega.dat, "regCompMnOmega.csv", row.names=F, sep=",")
write.table(regCompInvaPacts.dat, "regCompInvaPacts.csv", row.names=F, sep=",")

