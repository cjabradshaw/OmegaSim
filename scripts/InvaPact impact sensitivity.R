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


########################
## sensitivity analyses
########################

reg.iter <- 10000

## slowly increase impact in region 3

## region 1
## characteristics: many invasive species, many assessments, high confidence, high impacts
sppR1 <- 200 # number of invasive species
assess.limR1 <- c(10,80) # range of number of assessments per species
impact.wR1 <- c(0.2,0.3,0.8,0.9,0.95) # impact probability weights
conf.wR1 <- c(0.5,0.8,0.9) # confidence probability weights
natsppaffR1 <- 500 # native species affected

## region 2
## characteristics: many invasive species, few assessments, low confidence, high impacts
sppR2 <- 250 # number of species
assess.limR2 <- c(1,16) # range of number of assessments per species
impact.wR2 <- c(0.2,0.3,0.8,0.9,0.95) # impact probability weights
conf.wR2 <- c(0.95,0.5,0.2) # confidence probability weights
natsppaffR2 <- 600 # native species affected

## region 3
## characteristics: few invasive species, few assessments, low confidence, medium impacts
sppR3 <- 50 # number of species
assess.limR3 <- c(2,8) # range of number of assessments per species
conf.wR3 <- c(0.95,0.5,0.05) # confidence probability weights
natsppaffR3 <- 200 # native species affected



R3lowimpact.vec <- seq(0.95,0.05,-0.05)

omegasSdiff.md <- omegasSdiff.up <- omegasSdiff.lo <- rep(NA, length(R3lowimpact.vec))

for (f in 1:length(R3lowimpact.vec)) {
  
  ####################################
  ## generate regional characteristics
  ####################################
  
  # region 1
  mn.omega.sppR1 <- sd.omega.sppR1 <- rep(NA, sppR1)
  for (s in 1:sppR1) {
    Nassess.spp <- round(runif(1,min=assess.limR1[1], max=assess.limR1[2]), 0)
    impactspp <- sample(impactScore, Nassess.spp, replace=T, prob=impact.wR1)
    confspp <- sample(confScore, Nassess.spp, replace=T, prob=conf.wR1)
    mn.omega.sppR1[s] <- weighted.mean(impactspp, w=confspp)     
    sd.omega.sppR1[s] <- ifelse(length(unique(impactspp))==1 | is.na(wtd.sd(impactspp, w=confspp))==T |
                                  is.infinite(wtd.sd(impactspp, w=confspp))==T, 0, wtd.sd(impactspp, w=confspp))
  }
  
  # region 2
  mn.omega.sppR2 <- sd.omega.sppR2 <- rep(NA, sppR2)
  for (s in 1:sppR2) {
    Nassess.spp <- round(runif(1,min=assess.limR2[1], max=assess.limR2[2]), 0)
    impactspp <- sample(impactScore, Nassess.spp, replace=T, prob=impact.wR2)
    confspp <- sample(confScore, Nassess.spp, replace=T, prob=conf.wR2)
    mn.omega.sppR2[s] <- weighted.mean(impactspp, w=confspp)     
    sd.omega.sppR2[s] <- ifelse(length(unique(impactspp))==1 | is.na(wtd.sd(impactspp, w=confspp))==T |
                                  is.infinite(wtd.sd(impactspp, w=confspp))==T, 0, wtd.sd(impactspp, w=confspp))
  }
  
  # region 3  
  impact.wR3 <- c(R3lowimpact.vec[f],R3lowimpact.vec[f],0.5,1-R3lowimpact.vec[f],1-R3lowimpact.vec[f]) # impact probability weights
  
  mn.omega.sppR3 <- sd.omega.sppR3 <- rep(NA, sppR3)
  for (s in 1:sppR3) {
    Nassess.spp <- round(runif(1,min=assess.limR3[1], max=assess.limR3[2]), 0)
    impactspp <- sample(impactScore, Nassess.spp, replace=T, prob=impact.wR3)
    confspp <- sample(confScore, Nassess.spp, replace=T, prob=conf.wR3)
    mn.omega.sppR3[s] <- weighted.mean(impactspp, w=confspp) 
    sd.omega.sppR3[s] <- ifelse(length(unique(impactspp))==1 | is.na(wtd.sd(impactspp, w=confspp))==T |
                                  is.infinite(wtd.sd(impactspp, w=confspp))==T, 0, wtd.sd(impactspp, w=confspp))
  }
  
  omegaR1md <- omegaR2md <- omegaR3md <- omegasSdiff <- rep(NA,reg.iter)
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
    
    omegas <- c(omegaR1md[r], omegaR2md[r], omegaR3md[r])
    omegasSdiff[r] <- sum(unique(as.vector(abs(outer(omegas, omegas, '-')))))
    
  } # end r
  
  # sum of regional differences
  omegasSdiff.md[f] <- median(omegasSdiff)
  omegasSdiff.up[f] <- quantile(omegasSdiff, probs=0.975)
  omegasSdiff.lo[f] <- quantile(omegasSdiff, probs=0.025)
  
  print(f)
  
} # end f

plot(R3lowimpact.vec, omegasSdiff.md, type="l", xlab="probability low impact", ylab="sum region pairwise diffs",
     ylim=c(min(omegasSdiff.lo), max(omegasSdiff.up)))
lines(R3lowimpact.vec, omegasSdiff.lo, lty=2, col="red")
lines(R3lowimpact.vec, omegasSdiff.up, lty=2, col="red")


R3lowimpact.dat <- data.frame(R3lowimpact.vec, omegasSdiff.md, omegasSdiff.up, omegasSdiff.lo)
write.table(R3lowimpact.dat, "R3lowimpact.csv", row.names=F, sep=",")

