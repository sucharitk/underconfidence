
################################## Katyal & Fleming ################################
### Gender and anxiety reveal distinct computational sources of underconfidence ####
####################################################################################
################  This file reproduces Figure 3H and 3J ##################
###############################################################################
####### Plot and do stats on model estimated parameters  ##############
###############################################################################

rm(list=ls())

setwd("~/OneDrive - University College London/Projects/Experiments/confidenceTime/analysis/model/meta-ddm/")

library(Rcpp) # to source, compile and run C++ functions
library(DEoptim) # optimization algorithm
library(dplyr)
library(sjPlot)

fit.vratio <- T
# fit.tasks.separately <- F

sourceCpp("meta_DDM_sim.cpp") # this will give R access to the DDM_with_confidence_slow function 
sourceCpp("meta_DDM.cpp") # this will give R access to the DDM_with_confidence_slow function 

#Create a function that compares observed data to simulate data using quantile optimisation 
chi_square_optim <- function(params, observations, additional){
  
  if (length(additional)==4){
    # vratio is not specified as a fixed parameter, so estimate it
    #First, generate predictions:
    names(params) <- c('v','a','ter','z', 'conf_add','conf_mult', 'v_bias', 'vratio' )
    names(additional) <- c('returnFit', 'dt', 'sigma', 'ntrials')
    vratio <- params['vratio']
  } else {
    #First, generate predictions:
    names(params) <- c('v','a','ter','z', 'conf_add','conf_mult', 'v_bias' )
    names(additional) <- c('returnFit', 'dt', 'sigma', 'vratio', 'ntrials')
    vratio <- additional['vratio']
    
  }
  
  returnFit <- additional['returnFit']
  ntrials <- additional['ntrials']
  sigma <- additional['sigma']
  dt <- additional['dt']
  
  
  ntrials.obs <- length(observations$rt2)
  
  if (ntrials.obs < ntrials){
    t2distribution <- rep(observations$rt2, ceiling(ntrials/ntrials.obs))
  } else {
    t2distribution <- observations$rt2
  }
  
  if (returnFit==0){
    # simulate using the long form of post-decisional processing
    predictions <- data.frame(confTime_DDM_sim(
      v=params['v'], a=params['a'], ter=params['ter'], z=params['z'],
      ntrials=ntrials, s=sigma, dt=dt,
      t2distribution=t2distribution,  conf_add = params['conf_add'], 
      conf_mult = params['conf_mult'], 
      vratio = vratio, v_bias = params['v_bias']))
    
  } else{
    # recover using the short form
    predictions <- data.frame(confTime_DDM(
      v=params['v'], a=params['a'], ter=params['ter'], z=params['z'],
      ntrials=ntrials, s=sigma, dt=dt,
      t2distribution=t2distribution,  conf_add = params['conf_add'], 
      conf_mult = params['conf_mult'], 
      vratio = vratio, v_bias = params['v_bias']))
  }
  
  names(predictions) <- c('rt', 'resp','accu', 'stim', 'evidence', 'conf')
  
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
    
    # again, separate the predections according to the response
    c_predicted.t1 <- predictions[predictions$accu == 1 & predictions$rt2 <= t2.mid,]
    e_predicted.t1 <- predictions[predictions$accu == 0 & predictions$rt2 <= t2.mid,]
    c_predicted.t2 <- predictions[predictions$accu == 1 & predictions$rt2 >  t2.mid,]
    e_predicted.t2 <- predictions[predictions$accu == 0 & predictions$rt2 >  t2.mid,]
    
    # to make the next step easier, lets sort the predictions for correct and errors
    c_predicted_conf.t1 <- sort(c_predicted.t1$conf)
    e_predicted_conf.t1 <- sort(e_predicted.t1$conf)
    c_predicted_conf.t2 <- sort(c_predicted.t2$conf)
    e_predicted_conf.t2 <- sort(e_predicted.t2$conf)
    
    # First, separate the data in correct and error trials
    c_observed.t1 <- observations[observations$accu == 1 & observations$rt2 <= t2.mid,]
    e_observed.t1 <- observations[observations$accu == 0 & observations$rt2 <= t2.mid,]
    c_observed.t2 <- observations[observations$accu == 1 & observations$rt2 >  t2.mid,]
    e_observed.t2 <- observations[observations$accu == 0 & observations$rt2 >  t2.mid,]
    
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

#Create a function that compares observed data to simulate data using quantile optimisation 
chi_square_optim2 <- function(params, observations, additional){
  
  if (length(additional)==4){
    # vratio is not specified as a fixed parameter, so estimate it
    #First, generate predictions:
    names(params) <- c('v','a','ter','z', 'conf_add','conf_mult', 'v_bias', 'vratio' )
    names(additional) <- c('returnFit', 'dt', 'sigma', 'ntrials')
    vratio <- params['vratio']
  } else {
    #First, generate predictions:
    names(params) <- c('v','a','ter','z', 'conf_add','conf_mult', 'v_bias' )
    names(additional) <- c('returnFit', 'dt', 'sigma', 'vratio', 'ntrials')
    vratio <- additional['vratio']
    
  }
  
  returnFit <- additional['returnFit']
  ntrials <- additional['ntrials']
  sigma <- additional['sigma']
  dt <- additional['dt']
  
  
  ntrials.obs <- length(observations$rt2)
  
  if (ntrials.obs < ntrials){
    t2distribution <- rep(observations$rt2, ceiling(ntrials/ntrials.obs))
  } else {
    t2distribution <- observations$rt2
  }
  
  if (returnFit==0){
    # simulate using the long form of post-decisional processing
    predictions <- data.frame(confTime_DDM_sim(
      v=params['v'], a=params['a'], ter=params['ter'], z=params['z'],
      ntrials=ntrials, s=sigma, dt=dt,
      t2distribution=t2distribution,  conf_add = params['conf_add'], 
      conf_mult = params['conf_mult'], 
      vratio = vratio, v_bias = params['v_bias']))
    
  } else{
    # recover using the short form
    predictions <- data.frame(confTime_DDM_sim(
      v=params['v'], a=params['a'], ter=params['ter'], z=params['z'],
      ntrials=ntrials, s=sigma, dt=dt,
      t2distribution=t2distribution,  conf_add = params['conf_add'], 
      conf_mult = params['conf_mult'], 
      vratio = vratio, v_bias = params['v_bias']))
  }
  
  names(predictions) <- c('rt', 'resp','accu', 'stim', 'evidence', 'conf')
  
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
    
    # again, separate the predections according to the response
    c_predicted.t1 <- predictions[predictions$accu == 1 & predictions$rt2 <= t2.mid,]
    e_predicted.t1 <- predictions[predictions$accu == 0 & predictions$rt2 <= t2.mid,]
    c_predicted.t2 <- predictions[predictions$accu == 1 & predictions$rt2 >  t2.mid,]
    e_predicted.t2 <- predictions[predictions$accu == 0 & predictions$rt2 >  t2.mid,]
    
    # to make the next step easier, lets sort the predictions for correct and errors
    c_predicted_conf.t1 <- sort(c_predicted.t1$conf)
    e_predicted_conf.t1 <- sort(e_predicted.t1$conf)
    c_predicted_conf.t2 <- sort(c_predicted.t2$conf)
    e_predicted_conf.t2 <- sort(e_predicted.t2$conf)
    
    # First, separate the data in correct and error trials
    c_observed.t1 <- observations[observations$accu == 1 & observations$rt2 <= t2.mid,]
    e_observed.t1 <- observations[observations$accu == 0 & observations$rt2 <= t2.mid,]
    c_observed.t2 <- observations[observations$accu == 1 & observations$rt2 >  t2.mid,]
    e_observed.t2 <- observations[observations$accu == 0 & observations$rt2 >  t2.mid,]
    
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


###########
### Second try by varying v_bias, conf_add and conf_mult parameters

gam_t2  = c(3, 2) # simulate confidence response times
phq_dist = c(1.25, 4.37) # simulate phq scores (parameters obtained by fitting Katyal et al Exp 1 phq distribution)

vratio_dist <- c(1,.5)

ntrials <- 200
nsubj <- 100
AD <- c(scale(rgamma(nsubj, 1.25,.3)))
gender <- 2*rbinom(nsubj,1,.5)-1
beta.AD <- -1
beta.gender <- .5
model.num <- 1
v_bias <- AD
conf_add <- AD
conf_mult <- AD
v <- AD
a <- v
z <- v
ter <- v
returnFit <- 0
sigma <- 1
dt <- .005
vratio <- v

additional <- c(returnFit, dt, sigma, ntrials)

sim.data <- c()
for (ns in 1:nsubj){
  v[ns] <- max(rnorm(1, 1, .2), .1)
  a[ns] <- max(rnorm(1, 1.5, .5), .5)
  z[ns] <- a[ns]/2
  ter[ns] <- max(min(rnorm(1, .5, .1), .8), 0.05)

  vratio[ns] <- rnorm(1, vratio_dist[1], vratio_dist[2]) 
  
  v_bias[ns] <- beta.AD*AD[ns]
  conf_add[ns] <- rnorm(1,0,1.25) + beta.gender*gender[ns]
  conf_mult[ns] <- rgamma(1,6,6)
  
  observations <- data.frame(rgamma(ntrials, gam_t2[1], gam_t2[2]))
  colnames(observations) <- c('rt2')
  
  params <- c(v[ns], a[ns] , ter[ns], z[ns], conf_add[ns], conf_mult[ns], 
              v_bias[ns], vratio[ns])
  
  subj.data <- chi_square_optim(params, observations, additional)
  
  subj.data$subj <- ns
  # subj.data$conf_add <- conf_add[ns]
  # subj.data$conf_mult <- conf_mult[ns]
  # subj.data$v_bias <- v_bias[ns]
  subj.data$AD <- AD[ns]
  subj.data$gender <- gender[ns]
  sim.data <- rbind(sim.data, subj.data)
}

sim.data <- sim.data %>% mutate(accu = as.factor(accu)) %>%
  mutate(gender = as.factor(gender))

## check the simulation for the expected regression effects
m.reg <- lm(conf ~ accu, sim.data)
summary(m.reg)
plot_model(m.reg, type = ("pred"),line.size = 1.5, ci.lvl = .95,
           terms = c('accu'))

m.reg <- lm(conf ~ accu*rt2, sim.data)
summary(m.reg)
plot_model(m.reg, type = ("pred"),line.size = 1.5, ci.lvl = .95,
           terms = c('rt2', 'accu'))

m.reg <- lm(conf ~ AD + gender, sim.data)
summary(m.reg)
plot_model(m.reg, type = ("pred"),line.size = 1.5, ci.lvl = .95,
           terms = c('AD'))

m.reg <- lm(conf ~ AD*rt2*accu, sim.data)
summary(m.reg)
plot_model(m.reg, type = ("pred"),line.size = 1.5, ci.lvl = .95,
           terms = c("rt2", "AD", 'accu'),
           colors = c('green4', 'darkgoldenrod3', 'firebrick2'))


#######

# Model recovery

if (fit.vratio){
  additional <- c(1, .005, 1, 2500)
  lower <- c(0, .5,  0.05, .2, -6, .02, -5, -3)
  upper <- c(3,  4,  .2,   2,  6,   5,  5,  3)
  
} else {
  additional <- c(1, .005, 1, 1, 2500)
  lower <- c(0, .5,  0.05, .2, -6, .02, -5)
  upper <- c(3,  4,  .2,   2,  6,   5,  5)
}
nparams <- length(lower)

recover <- array(dim = c(nsubj, nparams))
bestfnval <- array(dim = c(1, nsubj))
for (ns in 1:nsubj){
  
  print(paste('Running participant',ns))
  print(paste('simulated               : v=', round(v[ns],2), 
              'a=', round(a[ns],2), 'ter=', round(ter[ns],2), 
              'z=', round(z[ns],2), 'conf_add=', round(conf_add[ns],2), 
              'conf_mult=', round(conf_mult[ns],2), 'v_bias=', round(v_bias[ns],2),
              'vratio=', round(vratio[ns],2)))
  
  subj.data <- subset(sim.data,subj==ns) %>% 
    select(c(rt, rt2, accu, conf))
  
  upper[3] <- min(subj.data$rt)
  
  optimal_params <- DEoptim(chi_square_optim, # function to optimize
                            # v,a,ter,z,conf_add,conf_mult,v_bias
                            lower = lower, 
                            upper = upper, #
                            observations = subj.data, 
                            additional=additional, 
                            DEoptim.control(itermax = 800, trace = 200))
  
  results <- summary(optimal_params)
  
  recover[ns,] <- as.numeric(results$optim$bestmem)
  bestfnval[ns] <- as.numeric(results$optim$bestval)
  
}


setwd("~/OneDrive - University College London/Projects/Experiments/confidenceTime/analysis/model/meta-ddm/")

### save data
sim.rec.ddm <- list(recover, additional, a, v, z, ter, conf_add, conf_mult, 
                    v_bias, vratio, beta.AD, beta.gender, AD, gender)
names(sim.rec.ddm) <- c('recover', 'additional', 'a', 'v', 'z', 'ter', 'conf_add', 
                        'conf_mult', 'v_bias', 'vratio', 'beta.AD', 'beta.gender',
                        'AD', 'gender')
saveRDS(sim.rec.ddm, 'simulate-recovery.Rds')


### load and plot data data

library(ggplot2)
library(ggpubr)

df <- data.frame(c(v_bias), c(conf_add), c(conf_mult), c(vratio),
                 c(recover[,5]), c(recover[,6]), c(recover[,7]), c(recover[,8]),
                 AD, gender)
colnames(df) <- c('vb.sim', 'add.sim', 'mult.sim', 'vr.sim',
                  'add.rec', 'mult.rec', 'vb.rec', 'vr.rec', 
                  'AD', 'gender')

ggscatter(df, x='vb.sim', y='vb.rec', add = 'reg.line', conf.int = T) +
  stat_cor(method = 'pearson')

ggscatter(df, x='add.sim', y='add.rec', add = 'reg.line', conf.int = T) +
  stat_cor(method = 'pearson')
  
ggscatter(df, x='mult.sim', y='mult.rec', add = 'reg.line', conf.int = T) +
  stat_cor(method = 'pearson')

ggscatter(df, x='vr.sim', y='vr.rec', add = 'reg.line', conf.int = T) +
  stat_cor(method = 'pearson')

mr <- lm(vb.rec ~ AD + gender, df)
summary(mr)

mr <- lm(add.sim ~ AD + gender, df)
summary(mr)

mr <- lm(vr.rec ~ AD + gender, df)
summary(mr)
