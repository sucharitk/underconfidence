################################## Katyal & Fleming ################################
### Gender and anxiety reveal distinct computational sources of underconfidence ####
####################################################################################
###############################################################################
####### Fit model to Dataset 4 (Rouault et al, Exp 2) ##############
###############################################################################

rm(list=ls())

setwd("~/OneDrive - University of Copenhagen/Projects/Experiments/confidenceTime/analysis/model/meta-ddm/")

library(Rcpp) # to source, compile and run C++ functions
library(DEoptim) # optimization algorithm
library(dplyr)
library(tidyverse)

fit.vratio <- T
fit.tasks.separately <- F

sourceCpp("meta_DDM.cpp") # this will give R access to the DDM_with_confidence_slow function 

#Create a function that compares observed data to simulate data using quantile optimisation 
chi_square_optim <- function(params, observations, additional){
  
  if (length(additional)==5){
    # vratio is not specified as a fixed parameter, so estimate it
    #First, generate predictions:
    names(params) <- c('v','a','ter','z', 'conf_add','conf_mult', 'v_bias', 'vratio' )
    names(additional) <- c('returnFit', 'dt', 'sigma', 'ntrials','nconflevels')
    vratio <- params['vratio']
    
  } else {
    #First, generate predictions:
    names(params) <- c('v','a','ter','z', 'conf_add','conf_mult', 'v_bias' )
    names(additional) <- c('returnFit', 'dt', 'sigma', 'vratio', 'ntrials', 'nconflevels')
    vratio <- additional['vratio']
  }  
  
  
  returnFit <- additional['returnFit']
  ntrials <- additional['ntrials']
  sigma <- additional['sigma']
  dt <- additional['dt']
  nconflevels <- additional['nconflevels']
  
  ntrials.obs <- length(observations$rt2)
  
  if (ntrials.obs < ntrials){
    t2distribution <- rep(observations$rt2, ceiling(ntrials/ntrials.obs))
  } else {
    t2distribution <- observations$rt2
  }
  
  predictions <- data.frame(meta_DDM(
    v=params['v'], a=params['a'], ter=params['ter'], z=params['z'],
    ntrials=ntrials, s=sigma, dt=dt,
    t2distribution=t2distribution,  conf_add = params['conf_add'], 
    conf_mult = params['conf_mult'], 
    vratio = vratio, v_bias = params['v_bias']))
  
  names(predictions) <- c('rt', 'resp','accu', 'stim', 'evidence', 'conf')
  
  if (nconflevels > 1){
    # for discrete confidence values tile them
    predictions$conf <- as.numeric(cut(predictions$conf, 
                                       breaks=seq(0,1,1/nconflevels)))
  }
  
  predictions$rt2 <- t2distribution[1:ntrials]
  
  #if we're only simulating data, return the predictions
  if(returnFit==0){ 
    return(predictions[,c('rt','accu','conf','rt2','stim')])
    
    #If we are fitting the model, now compare these predictions to the observations 
  }else{ 
    
    dim.obs <- dim(observations)[1]
    dim.pred <- dim(predictions)[1]
    
    # again, separate the predections according to the response
    c_predicted <- predictions[predictions$accu == 1,]
    e_predicted <- predictions[predictions$accu == 0,]
    
    # to make the next step easier, lets sort the predictions for correct and errors
    c_predicted_rt <- sort(c_predicted$rt)
    e_predicted_rt <- sort(e_predicted$rt)
    
    # First, separate the data in correct and error trials
    c_observed <- observations[observations$accu == 1,]
    e_observed <- observations[observations$accu == 0,]
    
    # Now, get the quantile RTs on the "observed data" for correct and error distributions separately (for quantiles .1, .3, .5, .7, .9)
    c_quantiles <- quantile(c_observed$rt, probs = c(.1,.3,.5,.7,.9), names = FALSE)
    e_quantiles <- quantile(e_observed$rt, probs = c(.1,.3,.5,.7,.9), names = FALSE)
    
    # to combine correct and incorrect we scale the expected interquantile probability by the proportion of correct and incorect respectively
    prop_obs_c <- dim(c_observed)[1] / dim.obs
    prop_obs_e <- 1 - prop_obs_c #dim(e_observed)[1] / dim.obs
    
    c_obs_proportion = prop_obs_c * c(.1, .2, .2, .2, .2, .1)
    e_obs_proportion = prop_obs_e * c(.1, .2, .2, .2, .2, .1)
    obs_props <- c(c_obs_proportion,e_obs_proportion)
    
    # now, get the proportion of responses that fall between the observed quantiles when applied to the predicted data (scale by N?)
    c_pred_proportion <- c(
      sum(c_predicted_rt <= c_quantiles[1]),
      sum(c_predicted_rt <= c_quantiles[2]) - sum(c_predicted_rt <= c_quantiles[1]),
      sum(c_predicted_rt <= c_quantiles[3]) - sum(c_predicted_rt <= c_quantiles[2]),
      sum(c_predicted_rt <= c_quantiles[4]) - sum(c_predicted_rt <= c_quantiles[3]),
      sum(c_predicted_rt <= c_quantiles[5]) - sum(c_predicted_rt <= c_quantiles[4]),
      sum(c_predicted_rt > c_quantiles[5])
    ) / dim.pred
    
    e_pred_proportion <- c(
      sum(e_predicted_rt <= e_quantiles[1]),
      sum(e_predicted_rt <= e_quantiles[2]) - sum(e_predicted_rt <= e_quantiles[1]),
      sum(e_predicted_rt <= e_quantiles[3]) - sum(e_predicted_rt <= e_quantiles[2]),
      sum(e_predicted_rt <= e_quantiles[4]) - sum(e_predicted_rt <= e_quantiles[3]),
      sum(e_predicted_rt <= e_quantiles[5]) - sum(e_predicted_rt <= e_quantiles[4]),
      sum(e_predicted_rt > e_quantiles[5])
    ) / dim.pred
    pred_props <- c(c_pred_proportion,e_pred_proportion)
    
    ###################################
    
    # separate into first and second half rt2 then take the difference of confidence distributions between them
    # this should make different predictions for bias in drift rate vs. response criteria
    t2.mid <- median(observations$rt2)
    
    # First, separate the data in correct and error trials
    c_observed.t1 <- observations[observations$accu == 1 & observations$rt2 <= t2.mid,]
    e_observed.t1 <- observations[observations$accu == 0 & observations$rt2 <= t2.mid,]
    c_observed.t2 <- observations[observations$accu == 1 & observations$rt2 >  t2.mid,]
    e_observed.t2 <- observations[observations$accu == 0 & observations$rt2 >  t2.mid,]
    
    # again, separate the predections according to the response
    c_predicted.t1 <- predictions[predictions$accu == 1 & predictions$rt2 <= t2.mid,]
    e_predicted.t1 <- predictions[predictions$accu == 0 & predictions$rt2 <= t2.mid,]
    c_predicted.t2 <- predictions[predictions$accu == 1 & predictions$rt2 >  t2.mid,]
    e_predicted.t2 <- predictions[predictions$accu == 0 & predictions$rt2 >  t2.mid,]
    
    if (nconflevels > 1){
      ## for discrete confidence levels
      
      c_obs_props_conf.t1 <- c()
      e_obs_props_conf.t1 <- c()
      c_obs_props_conf.t2 <- c()
      e_obs_props_conf.t2 <- c()
      for (ncl in 1:nconflevels){
        c_obs_props_conf.t1[ncl] <- sum(c_observed.t1$conf==ncl)/dim.obs
        e_obs_props_conf.t1[ncl] <- sum(e_observed.t1$conf==ncl)/dim.obs
        c_obs_props_conf.t2[ncl] <- sum(c_observed.t2$conf==ncl)/dim.obs
        e_obs_props_conf.t2[ncl] <- sum(e_observed.t2$conf==ncl)/dim.obs
      }
      obs_props_conf <- c(c_obs_props_conf.t1,e_obs_props_conf.t1,c_obs_props_conf.t2,e_obs_props_conf.t2)
      
      c_pred_props_conf.t1 <- c()
      e_pred_props_conf.t1 <- c()
      c_pred_props_conf.t2 <- c()
      e_pred_props_conf.t2 <- c()
      for (ncl in 1:nconflevels){
        c_pred_props_conf.t1[ncl] <- sum(c_predicted.t1$conf==ncl)/dim.pred
        e_pred_props_conf.t1[ncl] <- sum(e_predicted.t1$conf==ncl)/dim.pred
        c_pred_props_conf.t2[ncl] <- sum(c_predicted.t2$conf==ncl)/dim.pred
        e_pred_props_conf.t2[ncl] <- sum(e_predicted.t2$conf==ncl)/dim.pred
      }
      pred_props_conf <- c(c_pred_props_conf.t1,e_pred_props_conf.t1,c_pred_props_conf.t2,e_pred_props_conf.t2)
      
    } else {
      ## for continuous confidence levels
      # to make the next step easier, lets sort the predictions for correct and errors
      c_predicted_conf.t1 <- sort(c_predicted.t1$conf)
      e_predicted_conf.t1 <- sort(e_predicted.t1$conf)
      c_predicted_conf.t2 <- sort(c_predicted.t2$conf)
      e_predicted_conf.t2 <- sort(e_predicted.t2$conf)
      
      # Now, get the quantile RTs on the "observed data" for correct and error distributions separately (for quantiles .1, .3, .5, .7, .9)
      c_quantiles.t1 <- quantile(c_observed.t1$conf, probs = c(.1,.3,.5,.7,.9), names = FALSE)
      e_quantiles.t1 <- quantile(e_observed.t1$conf, probs = c(.1,.3,.5,.7,.9), names = FALSE)
      c_quantiles.t2 <- quantile(c_observed.t2$conf, probs = c(.1,.3,.5,.7,.9), names = FALSE)
      e_quantiles.t2 <- quantile(e_observed.t2$conf, probs = c(.1,.3,.5,.7,.9), names = FALSE)
      
      # to combine correct and incorrect we scale the expected interquantile probability by the proportion of correct and incorect respectively
      prop_obs_c.t1 <- dim(c_observed.t1)[1] / dim.obs
      prop_obs_e.t1 <- dim(e_observed.t1)[1] / dim.obs
      prop_obs_c.t2 <- dim(c_observed.t2)[1] / dim.obs
      prop_obs_e.t2 <- dim(e_observed.t2)[1] / dim.obs
      
      c_obs_proportion.t1 = prop_obs_c.t1 * c(.1, .2, .2, .2, .2, .1)
      e_obs_proportion.t1 = prop_obs_e.t1 * c(.1, .2, .2, .2, .2, .1)
      c_obs_proportion.t2 = prop_obs_c.t2 * c(.1, .2, .2, .2, .2, .1)
      e_obs_proportion.t2 = prop_obs_e.t2 * c(.1, .2, .2, .2, .2, .1)
      obs_props_conf <- c(c_obs_proportion.t1,e_obs_proportion.t1,c_obs_proportion.t2,e_obs_proportion.t2)
      
      # now, get the proportion of responses that fall between the observed quantiles when applied to the predicted data (scale by N?)
      c_pred_proportion_conf.t1 <- c(
        sum(c_predicted_conf.t1 <= c_quantiles.t1[1]),
        sum(c_predicted_conf.t1 <= c_quantiles.t1[2]) - sum(c_predicted_conf.t1 <= c_quantiles.t1[1]),
        sum(c_predicted_conf.t1 <= c_quantiles.t1[3]) - sum(c_predicted_conf.t1 <= c_quantiles.t1[2]),
        sum(c_predicted_conf.t1 <= c_quantiles.t1[4]) - sum(c_predicted_conf.t1 <= c_quantiles.t1[3]),
        sum(c_predicted_conf.t1 <= c_quantiles.t1[5]) - sum(c_predicted_conf.t1 <= c_quantiles.t1[4]),
        sum(c_predicted_conf.t1 > c_quantiles.t1[5])
      ) / dim.pred
      
      e_pred_proportion_conf.t1 <- c(
        sum(e_predicted_conf.t1 <= e_quantiles.t1[1]),
        sum(e_predicted_conf.t1 <= e_quantiles.t1[2]) - sum(e_predicted_conf.t1 <= e_quantiles.t1[1]),
        sum(e_predicted_conf.t1 <= e_quantiles.t1[3]) - sum(e_predicted_conf.t1 <= e_quantiles.t1[2]),
        sum(e_predicted_conf.t1 <= e_quantiles.t1[4]) - sum(e_predicted_conf.t1 <= e_quantiles.t1[3]),
        sum(e_predicted_conf.t1 <= e_quantiles.t1[5]) - sum(e_predicted_conf.t1 <= e_quantiles.t1[4]),
        sum(e_predicted_conf.t1 > e_quantiles.t1[5])
      ) / dim.pred
      
      c_pred_proportion_conf.t2 <- c(
        sum(c_predicted_conf.t2 <= c_quantiles.t2[1]),
        sum(c_predicted_conf.t2 <= c_quantiles.t2[2]) - sum(c_predicted_conf.t2 <= c_quantiles.t2[1]),
        sum(c_predicted_conf.t2 <= c_quantiles.t2[3]) - sum(c_predicted_conf.t2 <= c_quantiles.t2[2]),
        sum(c_predicted_conf.t2 <= c_quantiles.t2[4]) - sum(c_predicted_conf.t2 <= c_quantiles.t2[3]),
        sum(c_predicted_conf.t2 <= c_quantiles.t2[5]) - sum(c_predicted_conf.t2 <= c_quantiles.t2[4]),
        sum(c_predicted_conf.t2 > c_quantiles.t2[5])
      ) / dim.pred
      
      e_pred_proportion_conf.t2 <- c(
        sum(e_predicted_conf.t2 <= e_quantiles.t2[1]),
        sum(e_predicted_conf.t2 <= e_quantiles.t2[2]) - sum(e_predicted_conf.t2 <= e_quantiles.t2[1]),
        sum(e_predicted_conf.t2 <= e_quantiles.t2[3]) - sum(e_predicted_conf.t2 <= e_quantiles.t2[2]),
        sum(e_predicted_conf.t2 <= e_quantiles.t2[4]) - sum(e_predicted_conf.t2 <= e_quantiles.t2[3]),
        sum(e_predicted_conf.t2 <= e_quantiles.t2[5]) - sum(e_predicted_conf.t2 <= e_quantiles.t2[4]),
        sum(e_predicted_conf.t2 > e_quantiles.t2[5])
      ) / dim.pred
      pred_props_conf <- c(c_pred_proportion_conf.t1,e_pred_proportion_conf.t1,
                           c_pred_proportion_conf.t2,e_pred_proportion_conf.t2)
    }
    
    # Combine the quantiles for rts and cj
    obs_props <- c(obs_props,obs_props_conf)
    pred_props <- c(pred_props,pred_props_conf)
    
    pred_props[pred_props==0] <- .001
    
    # calculate chi square
    chiSquare = sum( ( (obs_props - pred_props) ^ 2) / pred_props )
    
    #Return chiSquare
    return(chiSquare)
  }
}

min.accept.rt <- 100 # in ms
vratio.string <- c('novratio', 'vratio')

# ########################
# ########################
# #################### rouaoult  exp 2
## read data file
setwd("~/OneDrive - University of Copenhagen/Projects/Experiments/confidenceTime/data/Rouault_et_al_BiolPsychiat_2018//")

df <- read.csv('rouault_2018_exp2.csv', header = TRUE, sep = ',')

df <- df %>%
  dplyr::rename(CIT = Compul) %>%
  dplyr::rename(accu = correct) %>%
  mutate(gender = recode_factor(gender, '0' = 'Male', '1' = 'Female'))


max_RT_deviation = 3
rt1_max = median(df$rt1, na.rm=T)+ max_RT_deviation*mad(df$rt1, na.rm=T ) # max rt for perception task
df <- df %>%
  mutate_at(vars(contains("rt1")), ~replace(.,.>rt1_max[1],NA))
rt2_max = median(df$rt2, na.rm=T)+ max_RT_deviation*mad(df$rt2, na.rm=T ) # max rt for perception task
df <- df %>%
  mutate_at(vars(contains("rt2")), ~replace(.,.>rt2_max[1],NA))

df$rt1[df$rt1 < min.accept.rt] <- NA
df <- df %>% drop_na(rt1, rt2)

df.subj <- df %>% group_by(subj, gender) %>%
  summarise(AD = mean(AD),
            CIT = mean(CIT),
            SW = mean(SW),
            age = mean(age))

AD.score <- df.subj$AD
CIT.score <- df.subj$CIT
SW.score <- df.subj$SW
gender <- df.subj$gender
age <- df.subj$age
subjs <- df.subj$subj
nsubj <- length(subjs)

exp.name <- 'rouault_2'

nconflevels <- max(unique(df$conf))

if (fit.vratio){
  additional <- c(1, .002, 1, 5000, nconflevels)
  lower <- c(0, .5,  0.05, .2, -8, .01, -7, -3)
  upper <- c(3,  4,  .2,   2,  8,   5,  7,  3)
  
} else {
  additional <- c(1, .005, 1, 1, 3500, nconflevels)
  lower <- c(0, .5,  0.05, .2, -6, .02, -5)
  upper <- c(3,  4,  .2,   2,  6,   5,  5)
}
nparams <- length(lower)

subjs <- unique(df$subj)
nsubj <- length(subjs)


recover <- array(dim = c(nsubj, nparams))
bestfnval <- array(dim = c(1, nsubj))

save.file <- sprintf('%s-%s.Rds', exp.name, vratio.string[fit.vratio+1])
print(save.file)

next.subj <- 1


for (ns in next.subj:nsubj){
  
  print(paste('Running participant',ns))
  
  next.subj <- ns + 1
  
  tempDat <- filter(df, subj==subjs[ns]) %>% drop_na(c(rt2, rt1)) %>%
    select(rt1, rt2, accu, conf) %>%
    rename(rt = rt1) %>%
    mutate(rt = rt/1000) %>%
    mutate(rt2 = rt2/1000)
  
  mean.accu <- mean(tempDat$accu)
  
  if (nrow(tempDat) > 0 && mean.accu>.5 && mean.accu<.95){
    
    upper[3] <- min(tempDat$rt)
    
    optimal_params <- DEoptim(chi_square_optim, # function to optimize
                              # v,a,ter,z,conf_add,conf_mult,v_pd
                              lower = lower, 
                              upper = upper,
                              observations = tempDat %>% select(c(rt, rt2, accu, conf)), 
                              additional=additional, 
                              DEoptim.control(itermax = 1000, trace = 200))
    
    
    print(paste('Participant ',ns))
    results <- summary(optimal_params)
    
    recover[ns,] <- as.numeric(results$optim$bestmem)  
    bestfnval[ns] <- as.numeric(results$optim$bestval)
  }
  
}

### save the recovered parameters
setwd("~/OneDrive - University of Copenhagen/Projects/Experiments/confidenceTime/analysis/model/meta-ddm/")

k2 <- list(recover, bestfnval, AD.score, CIT.score, SW.score, gender, age, exp.name, subjs)
names(k2) <- c('recover', 'bestfnval', 'AD', 'CIT', 'SW', 'gender', 'age', 'expname',
               'subj')
saveRDS(k2, save.file)



