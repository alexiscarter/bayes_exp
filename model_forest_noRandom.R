#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Zero altered model
# Find R2 for forest individually without the random effect
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Load libraries and data ####
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
library(plyr)
library(tidyverse)
library(rjags)
library(R2jags)
library(lattice) # using in the MyBUGSChains function
library(ggpubr)

load("trait_as.rdata")

## For modeling with survival data, create a column with O/1
trait_as$surv <- trait_as$total.w
trait_as$surv[which(trait_as$surv > 0)] <- 1

## Remove NIR
trait_as_ir <- trait_as %>%
  filter(!treatment == 'NIR') %>% 
  droplevels()

#-==-=-=-=-=-
# Model #####
#-==-=-=-=-=-

## Create a covariable matrix for the biomasse data
Xc<-model.matrix(~ forest , data= trait_as_ir)

## Create a covariable matrix for the survival
Xb<-model.matrix(~ forest, data= trait_as_ir)

## Define the number of parameter biomass part
Kc=ncol(Xc)

## Define the number of parameter survival part
Kb=ncol(Xb)

## List all the required variables for JAGS
JAGS.data_full<-list(
  Y     =  trait_as_ir$total.w,
  Xc    =  Xc,
  Xb    =  Xb,
  Kc    =  Kc,
  Kb    =  Kb,
  N     = nrow(trait_as_ir),
  Zeros = rep(0, nrow(trait_as_ir)))

## Write down the model
load.module("glm")
sink("ZA_LASSO.txt")
cat("
    model{
    #1A Priors beta and Gamma
    for(i in 1:Kc)  { beta[i]  ~ dnorm(0, 0.333333)}
    for(i in 1:Kb)  { gamma[i] ~ dnorm(0, 0.333333)}
    
    #1B Priors for r parameter of gamma distribution
    r~ dunif(0,5)
    
    #2 likelihood ( zero trick)
    C <- 1000
    for( i in 1:N){
    Zeros[i] ~ dpois(-ll[i] + C)
    z[i] <- step(Y[i] - 0.0001)
    l1[i] <- (1 - z[i]) * log(1-Pi[i])
    l2[i] <- z[i] * ( log(Pi[i]) - loggam(r) + r*log(r/mu[i])+
    (r-1)*log(Y[i]) - (Y[i] * r)/mu[i] )
    
    ll[i]<- l1[i] + l2[i]
    
    log(mu[i]) <- inprod( beta[], Xc[i,])
    logit(Pi[i]) <- inprod( gamma[], Xb[i,])
    }
    
    #3 Discrepancy measures
    for(i in 1:N){
    ExpY[i] <- Pi[i] * mu[i]
    VarY[i] <- (Pi[i] * r + Pi[i] - Pi[i] *
    Pi[i]*r) * (mu[i]+mu[i]/r)
    
    PRes[i] <- (Y[i] - ExpY[i])/ sqrt(VarY[i])
    }
    }", fill = TRUE)
sink()

## Specify the initiale values
inits<- function() {
  list(beta  = rnorm(Kc, 0, 0.1),
       gamma = rnorm(Kb, 0, 0.1),
       r = runif(1,0,5))
}

## Specify the parameters to store
params<-c("beta", "gamma", "PRes", "r", "ExpY")

## Run the model
ir.semi.info <- jags( data        = JAGS.data_full,
                      inits       = inits,
                      parameters  = params,
                      model       = "ZA_LASSO.txt",
                      n.thin      = 10,
                      n.chains    = 3,
                      n.burnin    = 4000,
                      n.iter      = 5000)

ir.semi.info2 <- update(ir.semi.info, n.iter = 10000, n.thin = 10) # number of iterations should be higher for final run, 100000 for example
out.ir <-ir.semi.info2$BUGSoutput

## check if the mixing is good ####
MyNames<- c(
  paste(c(colnames(Xc)) , "Gamma", sep = " "),
  paste(c(colnames(Xb)) , "Bern" , sep = " "),
  "r Gamma")

MyBUGSChains(out.ir, # /!\/!\/!\ before run MCMCSupportHighstatV4.R
             c(uNames("beta",Kc),
               uNames("gamma",Kb),
               "r"),
             PanelNames = MyNames)

## get the numerical output of JAGS  ####
OUT.ir <- MyBUGSOutput(out.ir,  c(uNames("beta",Kc), uNames("gamma",Kb), "r"), VarNames = MyNames)

## Get the posterior distribution ####
## Gamma part of the model
MyNamesG<- c(colnames(Xc),"r")

MyBUGSHist(out.ir, c(uNames("beta",Kc) ,"r"), PanelNames = MyNamesG)

## Bernoulli part of the model
MyNamesB<- c(colnames(Xb))

MyBUGSHist(out.ir, c(uNames("gamma",Kb)), PanelNames = MyNamesB)

## Get Pearson residuals and expected values of the model ####
E <-out.ir$mean$PRes
ExpY <-out.ir$mean$ExpY

## calculate de dispersion of the model ####
D<- sum(E^2)/(nrow(trait_as_ir) + (ncol(Xc)+ ncol(Xb)))
D #0.3518682

## Model validation of the Hurdle model ####
theme_set(theme_bw())

data1<-data.frame(E,ExpY)
Pres_ExpV<-ggplot(data= data1, aes( x= ExpY, y=E))+
  geom_point()+
  geom_hline(yintercept = 0, linetype="dotted")+
  ylab('Pearson residual')+
  xlab('Expected values')

data3<-data.frame(trait_as_ir$forest, E)
Pres_abio<-ggplot(data= data3, aes( x= trait_as_ir$forest  , y=E))+
  geom_boxplot()+
  geom_hline(yintercept = 0, linetype="dotted")+
  ylab('Pearson residual')+
  xlab('Forest')

data4<-data.frame(trait_as_ir$total.w, ExpY)
ObvsExp<-ggplot(data= data4, aes(y = ExpY, x=trait_as_ir$total.w ))+
  geom_point()+
  labs(x ='Observed', y='Expected')

model_val<-ggarrange(Pres_ExpV,ObvsExp, Pres_abio,
                     ncol=2, nrow=2, labels=c("a", "b", "c", "d"), 
                     common.legend = T, legend = c("right"))
model_val

## get the R2 of the observed value ~ the expected value
tt<-lm(ExpY~trait_as_ir$total.w)
summary(tt)$adj.r.squared #0.0008364409

## sketch the model fit
MyData <- expand.grid(Forest = levels(trait_as_ir$forest))

X<-model.matrix(~ Forest , data= MyData)
beta.mcmc  <- out.ir$sims.list$beta
gamma.mcmc <- out.ir$sims.list$gamma

mu.mcmc    <- exp(X %*% t(beta.mcmc))
Pi.mcmc    <- exp(X %*% t(gamma.mcmc)) / (1 + exp(X %*% t(gamma.mcmc)))
ExpY.mcmc <- Pi.mcmc * mu.mcmc

## function for confidence intervals
GetCIs <- function(x) {
  OUT <- matrix(nrow = nrow(x), ncol = 4) 
  for(i in 1:nrow(x)){
    xi <- x[i,]	
    OUT[i,3:4] <- quantile(xi, probs = c(0.15, 0.85)) # 70% confidence interval
    OUT[i,1] <- mean(xi)
    OUT[i,2] <- sd(xi)
  }
  colnames(OUT) <- c("mean", "se", "lo", "up")
  OUT
}

## Apply the function
L <- GetCIs(ExpY.mcmc)
L.Pi <- GetCIs(Pi.mcmc)
L.mu <- GetCIs(mu.mcmc)

## Combine L with MyData for ggplot2
expH <- cbind(MyData, L)
expB <- cbind(MyData, L.Pi)
expG <- cbind(MyData, L.mu)

## Plot the Hurdle results
Hurdle_O<- ggplot() +
  geom_point(data = expH, aes(x = Forest,  y = mean),position=position_dodge(width = 0.5))+ 
  geom_errorbar(data = expH,aes(x = Forest, ymax = up,ymin = lo),position=position_dodge(width = 0.5), width = 0.2)
Hurdle_O # not the same output than full model

## Calculate AIC
nrow(trait_as_ir)*(log(2*pi)+1+log((sum(E^2)/nrow(trait_as_ir))))+((length(L)+1)*2)
#205.4231

expH
#Forest      mean        se        lo        up
#1     FT 0.5016358 0.1594072 0.3488253 0.6542230
#2     FM 0.4421596 0.1388865 0.3085716 0.5758111
#3     FB 0.2859975 0.1050780 0.1889673 0.3865533
