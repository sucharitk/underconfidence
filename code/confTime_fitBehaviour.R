
################################## Katyal & Fleming ################################
### Gender and anxiety reveal distinct computational sources of underconfidence ####
####################################################################################
################  This file reproduces Figure 3G & 3I ##################
###############################################################################
####### Plot model fits alongside data  ##############
###############################################################################

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(R.matlab)
library(Rcpp) # to source, compile and run C++ functions

setwd("~/OneDrive - University of Copenhagen/Projects/Experiments/confidenceTime/analysis/model/meta-ddm/")

sourceCpp("meta_DDM.cpp") # this will give R access to the DDM_with_confidence_slow function 

#Create a function that compares observed data to simulate data using quantile optimisation 
simulate_ddm <- function(params, observations, additional){
  
  # browser()
  vratio <- params['vratio']
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
  
  if (length(params)==8){
    # vratio is not specified as a fixed parameter, so estimate it
    #First, generate predictions:
    names(params) <- c('v','a','ter','z', 'conf_add','conf_mult', 'v_bias', 'vratio' )
    names(additional) <- c('returnFit', 'dt', 'sigma', 'ntrials','nconflevels')
    predictions <- data.frame(meta_DDM(
      v=params['v'], a=params['a'], ter=params['ter'], z=params['z'],
      ntrials=ntrials, s=sigma, dt=dt,
      t2distribution=t2distribution,  conf_add = params['conf_add'], 
      conf_mult = params['conf_mult'], 
      vratio = vratio, v_bias = params['v_bias']), c(1:ntrials))
    
    names(predictions) <- c('rt', 'resp','accu', 'stim', 'evidence', 'conf', 'trialnum')
    
  } else {
    #First, generate predictions:
    names(params) <- c('v','a','ter','z', 'conf_add','conf_mult', 'v_bias', 'vratio','ter2')
    names(additional) <- c('returnFit', 'dt', 'sigma', 'ntrials','nconflevels')
    predictions <- data.frame(meta_DDM_ter2(
      v=params['v'], a=params['a'], ter=params['ter'], z=params['z'],
      ntrials=ntrials, s=sigma, dt=dt,
      t2distribution=t2distribution,  conf_add = params['conf_add'], 
      conf_mult = params['conf_mult'], 
      vratio = vratio, v_bias = params['v_bias'], ter2 = params['ter2']))
    
    names(predictions) <- c('rt', 'resp','accu', 'stim', 'evidence', 'conf')
    
  }
  
  if (nconflevels > 1){
    # for discrete confidence values tile them
    predictions$conf <- as.numeric(cut(predictions$conf, 
                                       breaks=seq(0,1,1/nconflevels)))
  }
  
  predictions$rt2 <- t2distribution[1:ntrials]
  
  #if we're only simulating data, return the predictions
  if(returnFit==0){ 
    return(predictions[,c('rt','accu','conf','rt2','stim','trialnum')])
    
    #If we are fitting the model, now compare these predictions to the observations 
  }
}


###
font_size <- 16
tag_size <- 16
alpha <- .4
max_RT_deviation = 3
excl_insuf_var <- T

taskNames = c('Perception', 'Memory')
posnegNames = c('Positive', 'Negative')
posnegColours = c("chartreuse3", "coral3")
min.accept.rt <- 100

n_exp <- 8+1  

num.tiles.rt2 <- 6


####
#### Katyal et al - Exp 1
##

setwd("~/OneDrive - University of Copenhagen/Projects/Experiments/confidenceTime/data/Katyal_et_al_2023//")

## read data file
katyal_1.df = read.csv('mbsExp1_noperfexcl.csv', header = TRUE, sep = ',')
## recode columns as factors
katyal_1.df$subj <- as.factor(katyal_1.df$subj)
katyal_1.df$task <- as.factor(katyal_1.df$task)
katyal_1.df$gender[is.nan(katyal_1.df$gender)] <- NA
katyal_1.df$gender <- as.factor(katyal_1.df$gender)

katyal_1.df <- katyal_1.df %>%
  mutate_at(vars(contains("fb")), as.factor) %>%
  dplyr::rename(accu = corr) %>%
  dplyr::rename(stair = incdec) %>%
  dplyr::select(-starts_with('posword'))%>%
  dplyr::select(-starts_with('negword'))%>%
  dplyr::select(-starts_with('endors'))%>%
  dplyr::select(-starts_with('aware'))%>%
  dplyr::select(-starts_with('affect')) %>%
  mutate(AD = scale(gad)) %>%
  dplyr::rename(trialinrun = trialnum) %>%
  mutate(trialnum = trialinrun + (runnum-1)*40)

katyal_1.df <- katyal_1.df %>%
  mutate(fbblock = recode_factor(fbblock, "1" = posnegNames[1], 
                                 "2" = posnegNames[2], "0" = "neut"),
         task = recode_factor(task, "0" = taskNames[1], "1" = taskNames[2])  )

katyal_1.df$rt1[katyal_1.df$rt1 < min.accept.rt] <- NA

## code highly deviant RTs as NA
rt1_max = c(median(katyal_1.df$rt1[katyal_1.df$task==taskNames[1]], na.rm=T)+
              max_RT_deviation*mad(katyal_1.df$rt1[katyal_1.df$task==taskNames[1]], na.rm=T ), # max rt for perception task
            median(katyal_1.df$rt1[katyal_1.df$task==taskNames[2]], na.rm=T)+
              max_RT_deviation*mad(katyal_1.df$rt1[katyal_1.df$task==taskNames[2]], na.rm=T ))
katyal_1.df <- katyal_1.df %>%
  mutate_at(vars(contains("rt1")), ~replace(.,task==taskNames[1] & .>rt1_max[1],NA)) %>%
  mutate_at(vars(contains("rt1")), ~replace(.,task==taskNames[2] & .>rt1_max[2],NA))

# ## scale RTs
katyal_1.df$scale_rt1 <- katyal_1.df$rt1
katyal_1.df$scale_rt1[katyal_1.df$task==taskNames[1]] <- scale(katyal_1.df$rt1[katyal_1.df$task==taskNames[1]])
katyal_1.df$scale_rt1[katyal_1.df$task==taskNames[2]] <- scale(katyal_1.df$rt1[katyal_1.df$task==taskNames[2]])

rt2_max = c(median(katyal_1.df$rt2[katyal_1.df$task==taskNames[1]], na.rm=T)+ 
              max_RT_deviation*mad(katyal_1.df$rt2[katyal_1.df$task==taskNames[1]], na.rm=T ), # max rt for perception task
            median(katyal_1.df$rt2[katyal_1.df$task==taskNames[2]], na.rm=T)+ 
              max_RT_deviation*mad(katyal_1.df$rt2[katyal_1.df$task==taskNames[2]], na.rm=T ))
katyal_1.df <- katyal_1.df %>%
  mutate_at(vars(contains("rt2")), ~replace(.,task==taskNames[1] & .>rt2_max[1],NA)) %>%
  mutate_at(vars(contains("rt2")), ~replace(.,task==taskNames[2] & .>rt2_max[2],NA))

katyal_1.df <- katyal_1.df %>% drop_na(rt1, rt2)

katyal_1.df$scale_rt2 <- katyal_1.df$rt2
katyal_1.df$scale_rt2[katyal_1.df$task==taskNames[1]] <- scale(katyal_1.df$rt2[katyal_1.df$task==taskNames[1]])
katyal_1.df$scale_rt2[katyal_1.df$task==taskNames[2]] <- scale(katyal_1.df$rt2[katyal_1.df$task==taskNames[2]])

katyal_1.df$conf.sc[katyal_1.df$task==taskNames[1]] <- scale(katyal_1.df$conf[katyal_1.df$task==taskNames[1]])
katyal_1.df$conf.sc[katyal_1.df$task==taskNames[2]] <- scale(katyal_1.df$conf[katyal_1.df$task==taskNames[2]])

katyal_1.df <- katyal_1.df %>%mutate_at(c("group", 'gender'), as.factor) %>%
  mutate(gender = recode_factor(gender, '1' = 'Women', '2' = 'Men'))

expNum <- 1

katyal_1.df <- katyal_1.df %>% 
  dplyr::select(c(subj, task, rt2, scale_rt2, conf, conf.sc, age, gender, AD, accu))



setwd("~/OneDrive - University of Copenhagen/Projects/Experiments/confidenceTime/analysis/model/meta-ddm/")

k1.perc <- readRDS('katyal_1-Perception-vratio.Rds')

k1.perc <- data.frame(k1.perc$AD, k1.perc$gender, k1.perc$age, 
                      k1.perc$recover[,7], k1.perc$recover[,6], 
                      k1.perc$recover[,5], k1.perc$recover[,8],
                      c(k1.perc$bestfnval), k1.perc$subj, 
                      k1.perc$recover[,1], k1.perc$recover[,2], 
                      k1.perc$recover[,3], k1.perc$recover[,4])
colnames(k1.perc) <- c('AD', 'gender', 'age', 
                       'v.bias', 'mult.bias', 'add.bias', 'vratio', 
                       'fnval', 'subj', 'v', 'a', 'ter', 'z')
k1.perc$task <- 'Perception'


k1.mem <- readRDS('katyal_1-Memory-vratio.Rds')

k1.mem <- data.frame(k1.mem$AD, k1.mem$gender, k1.mem$age, 
                     k1.mem$recover[,7], k1.mem$recover[,6], 
                     k1.mem$recover[,5], k1.mem$recover[,8],
                     c(k1.mem$bestfnval), k1.mem$subj, 
                     k1.mem$recover[,1], k1.mem$recover[,2], 
                     k1.mem$recover[,3], k1.mem$recover[,4])
colnames(k1.mem) <- c('AD', 'gender', 'age', 
                      'v.bias', 'mult.bias', 'add.bias', 'vratio', 
                      'fnval', 'subj', 'v', 'a', 'ter', 'z')
k1.mem$task <- 'Memory'

AD <- k1.mem$AD
age <- k1.mem$age
gender <- k1.mem$gender

# subjs <- c(1,32,98,176,182)
additional <- c(0, .001, 1, 5000, 1)
names(additional) <- c('returnFit', 'dt', 'sigma', 'ntrials', 'nconflevels')


### for each study simulate all subjects
subjs <- unique(k1.perc$subj)
nsubj <- length(subjs)
ss <- c(1:nsubj)
allsubj <- c()
all.sim.data.perc <- c()
all.sim.data.mem <- c()
for (ns in ss){
  
  task.name <- 'Perception'
  
  ks <- filter(k1.perc, subj==subjs[ns])
  real_params <- c(ks$v,ks$a,ks$ter,ks$z, 
                   ks$add.bias, ks$mult.bias, ks$v.bias, ks$vratio)
  names(real_params) <- c("v","a","ter","z", "conf_add","conf_mult",'v_bias', 'vratio')
  
  subj.data <- katyal_1.df %>% filter(subj==ns, task==task.name)
  conf.ext <- subj.data$conf
  conf.ext[conf.ext<.5] <- 1-conf.ext[conf.ext<.5]
  
  if (nrow(subj.data)>80 & !is.na(real_params['v'])){
    if (mad(conf.ext)>0){
      # all RT2
      observations <-  (subj.data %>% dplyr::select(c('rt2')) %>% drop_na())/1000
      
      sim.data <- simulate_ddm(real_params, observations, additional)
      
      sim.data <- sim.data %>% dplyr::select(c(rt2, conf, accu, trialnum))
      sim.data$subj <- as.factor(subjs[ns])
      sim.data$AD <- AD[as.numeric(ns)]
      sim.data$age <- age[as.numeric(ns)]
      sim.data$gender <- gender[as.numeric(ns)]
      sim.data$task <- task.name
      
      all.sim.data.perc <- rbind(all.sim.data.perc, sim.data)
      
    }
    
  }
  
  task.name <- 'Memory'
  
  ks <- filter(k1.mem, subj==subjs[ns])
  real_params <- c(ks$v,ks$a,ks$ter,ks$z, 
                   ks$add.bias, ks$mult.bias, ks$v.bias, ks$vratio)
  names(real_params) <- c("v","a","ter","z", "conf_add","conf_mult",'v_bias', 'vratio')
  
  subj.data <- katyal_1.df %>% filter(subj==ns, task==task.name)
  conf.ext <- subj.data$conf
  conf.ext[conf.ext<.5] <- 1-conf.ext[conf.ext<.5]
  
  if (nrow(subj.data)>80 & !is.na(real_params['v'])){
    if (mad(conf.ext)>0){
      observations <-  (subj.data %>% dplyr::select(c('rt2')) %>% drop_na())/1000
      
      sim.data <- simulate_ddm(real_params, observations, additional)
      
      sim.data <- sim.data %>% dplyr::select(c(rt2, conf, accu, trialnum))
      sim.data$subj <- as.factor(subjs[ns])
      sim.data$AD <- AD[as.numeric(ns)]
      sim.data$age <- age[as.numeric(ns)]
      sim.data$gender <- gender[as.numeric(ns)]
      sim.data$task <- task.name
      
      all.sim.data.mem <- rbind(all.sim.data.mem, sim.data)
      
    }
    
  }
  
}

all.sim.data.perc$scale_rt2 <- c(scale(all.sim.data.perc$rt2))
all.sim.data.perc$conf.sc <- c(scale(all.sim.data.perc$conf))

all.sim.data.mem$scale_rt2 <- c(scale(all.sim.data.mem$rt2))
all.sim.data.mem$conf.sc <- c(scale(all.sim.data.mem$conf))

k1.sim <- rbind(all.sim.data.perc, all.sim.data.mem) %>%
  dplyr::select(c(conf.sc, AD, scale_rt2, accu, gender, age, subj)) %>%
  mutate(gender = recode_factor(gender, 'Male' = 'Men', 'Female' = 'Women')) %>%
  mutate(subj = as.factor(as.numeric(subj) + expNum*1000)) %>%
  mutate(data = 'Model')


katyal_1.sim <- k1.sim %>%     
  # group_by(subj) %>%
  mutate(srt2 = ntile(scale_rt2,num.tiles.rt2)) %>%
  group_by(srt2, subj, gender, data) %>%
  dplyr::summarise(conf = mean(conf.sc, na.rm = T))

katyal_1.tile <- k1.sim %>% dplyr::group_by(subj) %>%
  dplyr::summarise(AD = mean(AD),
            age = mean(age))  %>%
  mutate(ADt = ntile(AD,2)) %>%
  mutate(ADt = recode_factor(ADt, '1' = 'Low', '2' = 'High')) %>%
  mutate(ADt2 = ntile(AD,3)) %>%
  mutate(ADt2 = recode_factor(ADt2, '1' = 'Low', '2' = 'Mid', '3' = 'High')) %>%
  mutate(aget = ntile(age,2)) %>%
  mutate(aget = recode_factor(aget, '1' = 'Younger', '2' = 'Older'))

katyal_1.sim <- katyal_1.sim %>% dplyr::left_join(katyal_1.tile)


k1.obs <- katyal_1.df %>% 
  dplyr::select(c(conf.sc, AD, scale_rt2, accu, gender, age, subj)) %>%
  mutate(gender = recode_factor(gender, 'Male' = 'Men', 'Female' = 'Women')) %>%
  drop_na() %>%
  mutate(data = 'Observed') %>%
  mutate(subj = as.factor(as.numeric(subj) + expNum*1000))

katyal_1.obs <- k1.obs %>%
  # group_by(subj) %>%
  mutate(srt2 = ntile(scale_rt2,num.tiles.rt2)) %>%
  group_by(srt2, subj, gender, data) %>%
  dplyr::summarise(conf = mean(conf.sc))

katyal_1.tile <- k1.obs %>% group_by(subj) %>%
  dplyr::summarise(AD = mean(AD),
            age = mean(age))  %>%
  mutate(ADt = ntile(AD,2)) %>%
  mutate(ADt = recode_factor(ADt, '1' = 'Low', '2' = 'High')) %>%
  mutate(ADt2 = ntile(AD,3)) %>%
  mutate(ADt2 = recode_factor(ADt2, '1' = 'Low', '2' = 'Mid', '3' = 'High')) %>%
  mutate(aget = ntile(age,2)) %>%
  mutate(aget = recode_factor(aget, '1' = 'Younger', '2' = 'Older'))

katyal_1.obs <- katyal_1.obs %>% dplyr::left_join(katyal_1.tile)

obs.val.subj <- unique(katyal_1.obs$subj)
sim.val.subj <- unique(katyal_1.sim$subj)
katyal_1.both <- rbind(filter(katyal_1.obs), filter(katyal_1.sim))


########################################
########################################
######## Katyal 2

setwd("~/OneDrive - University of Copenhagen//Projects/Experiments/confidenceTime/data/Katyal_et_al_2023/")

## read data file
katyal_2.df = read.csv(paste('mbsExp2_noperfexcl.csv', sep=''), header = TRUE, sep = ',')

katyal_2.df$subj <- as.factor(katyal_2.df$subj)
katyal_2.df$task <- as.factor(katyal_2.df$task)
katyal_2.df$gender[is.nan(katyal_2.df$gender)] <- NA
katyal_2.df$gender <- as.factor(katyal_2.df$gender)

cumtrials <- c(0,40,80,120,140,180)
katyal_2.df <- katyal_2.df %>%
  mutate_at(vars(contains("fb")), as.factor) %>%
  dplyr::rename(accu = corr) %>%
  dplyr::rename(stair = incdec)%>%
  dplyr::select(-c(starts_with('word')))%>%
  dplyr::select(-c(starts_with('endors'))) %>%
  dplyr::rename(trialinrun = trialnum) %>%
  mutate(trialnum = trialinrun + cumtrials[runnum])

katyal_2.df <- katyal_2.df %>%
  mutate(fbblock = recode_factor(fbblock, "1" = posnegNames[1], 
                                 "2" = posnegNames[2], "0" = "None"),
         task = recode_factor(task, "0" = taskNames[1], "1" = taskNames[2])  )

## code highly deviant RTs as NA
max_RT_deviation = 3
rt1_max = c(median(katyal_2.df$rt1[katyal_2.df$task==taskNames[1]], na.rm=T)+ 
              max_RT_deviation*mad(katyal_2.df$rt1[katyal_2.df$task==taskNames[1]], na.rm=T ), # max rt for perception task
            median(katyal_2.df$rt1[katyal_2.df$task==taskNames[2]], na.rm=T)+ 
              max_RT_deviation*mad(katyal_2.df$rt1[katyal_2.df$task==taskNames[2]], na.rm=T ))
katyal_2.df <- katyal_2.df %>%
  mutate_at(vars(contains("rt1")), ~replace(.,task==taskNames[1] & .>rt1_max[1],NA)) %>%
  mutate_at(vars(contains("rt1")), ~replace(.,task==taskNames[2] & .>rt1_max[2],NA))

## scale RTs 
katyal_2.df$scale_rt1 <- katyal_2.df$rt1
katyal_2.df$scale_rt1[katyal_2.df$task==taskNames[1]] <- scale(katyal_2.df$rt1[katyal_2.df$task==taskNames[1]])
katyal_2.df$scale_rt1[katyal_2.df$task==taskNames[2]] <- scale(katyal_2.df$rt1[katyal_2.df$task==taskNames[2]])

rt2_max = c(median(katyal_2.df$rt2[katyal_2.df$task==taskNames[1]], na.rm=T)+ 
              max_RT_deviation*mad(katyal_2.df$rt2[katyal_2.df$task==taskNames[1]], na.rm=T ), # max rt for perception task
            median(katyal_2.df$rt2[katyal_2.df$task==taskNames[2]], na.rm=T)+ 
              max_RT_deviation*mad(katyal_2.df$rt2[katyal_2.df$task==taskNames[2]], na.rm=T ))

katyal_2.df <- katyal_2.df %>%
  mutate_at(vars(contains("rt2")), ~replace(.,task==taskNames[1] & .>rt2_max[1],NA)) %>%
  mutate_at(vars(contains("rt2")), ~replace(.,task==taskNames[2] & .>rt2_max[2],NA))

katyal_2.df$rt1[katyal_2.df$rt1 < min.accept.rt] <- NA
katyal_2.df <- katyal_2.df %>% drop_na(rt1, rt2)

katyal_2.df$scale_rt2 <- katyal_2.df$rt2
katyal_2.df$scale_rt2[katyal_2.df$task==taskNames[1]] <- scale(katyal_2.df$rt2[katyal_2.df$task==taskNames[1]])
katyal_2.df$scale_rt2[katyal_2.df$task==taskNames[2]] <- scale(katyal_2.df$rt2[katyal_2.df$task==taskNames[2]])

katyal_2.df$stimev[katyal_2.df$task==taskNames[1]] <- scale(katyal_2.df$stair[katyal_2.df$task==taskNames[1]])
katyal_2.df$stimev[katyal_2.df$task==taskNames[2]] <- scale(katyal_2.df$stair[katyal_2.df$task==taskNames[2]])

katyal_2.df$conf.sc[katyal_2.df$task==taskNames[1]] <- scale(katyal_2.df$conf[katyal_2.df$task==taskNames[1]])
katyal_2.df$conf.sc[katyal_2.df$task==taskNames[2]] <- scale(katyal_2.df$conf[katyal_2.df$task==taskNames[2]])

katyal_2.df <- katyal_2.df %>%
  mutate_at(c("group", 'gender'), as.factor) %>%
  mutate(group = recode_factor(group, "1" = 'Group 1', 
                               "2" = 'Group 2', "3" = 'Group 3', "4" = 'Group 4', 
                               "5" = 'Group 5', "6" = 'Group 6', "7" = 'Group 7', 
                               "8" = 'Group 8')) %>%
  mutate(gender = recode_factor(gender, '1' = 'Female', '2' = 'Male'))

# #####################
# #####################
# ## read in the psychiatric scores
psych = read.csv(paste('factor_scores_noperfexc.csv', sep=''),
                 header = TRUE, sep = ',')
psych <- psych %>% dplyr::rename(subj = subjIDs) %>% dplyr::select(-X) %>%
  mutate_at("subj", as.factor) %>%
  dplyr::rename(CIT = Compul)

katyal_2.df <- katyal_2.df %>% dplyr::left_join(psych)

katyal_2.df <- katyal_2.df %>% drop_na(scale_rt1, scale_rt2)


expNum <- 2

setwd("~/OneDrive - University of Copenhagen/Projects/Experiments/confidenceTime/analysis/model/meta-ddm/")

k2.perc <- readRDS('katyal_2-Perception-vratio.Rds')

k2.perc <- data.frame(k2.perc$AD, k2.perc$gender, k2.perc$age, k2.perc$CIT,
                      k2.perc$recover[,7], k2.perc$recover[,6], 
                      k2.perc$recover[,5], k2.perc$recover[,8],
                      c(k2.perc$bestfnval), k2.perc$subj, 
                      k2.perc$recover[,1], k2.perc$recover[,2], 
                      k2.perc$recover[,3], k2.perc$recover[,4])
colnames(k2.perc) <- c('AD', 'gender', 'age', 'CIT',
                       'v.bias', 'mult.bias', 'add.bias', 'vratio', 
                       'fnval', 'subj', 'v', 'a', 'ter', 'z')
k2.perc$task <- 'Perception'


k2.mem <- readRDS('katyal_2-Memory-vratio.Rds')

k2.mem <- data.frame(k2.mem$AD, k2.mem$gender, k2.mem$age, k2.mem$CIT,
                     k2.mem$recover[,7], k2.mem$recover[,6], 
                     k2.mem$recover[,5], k2.mem$recover[,8],
                     c(k2.mem$bestfnval), k2.mem$subj, 
                     k2.mem$recover[,1], k2.mem$recover[,2], 
                     k2.mem$recover[,3], k2.mem$recover[,4])
colnames(k2.mem) <- c('AD', 'gender', 'age', 'CIT',
                      'v.bias', 'mult.bias', 'add.bias', 'vratio', 
                      'fnval', 'subj', 'v', 'a', 'ter', 'z')
k2.mem$task <- 'Memory'

AD <- k2.mem$AD
CIT <- k2.mem$CIT
age <- k2.mem$age
gender <- k2.mem$gender

### for each study simulate all subjects
subjs <- unique(k2.perc$subj)
nsubj <- length(subjs)
allsubj <- c()
ss <- c(1:nsubj)
all.sim.data.perc <- c()
all.sim.data.mem <- c()
for (ns in ss){
  
  task.name <- 'Perception'
  
  ks <- filter(k2.perc, subj==subjs[ns])
  real_params <- c(ks$v,ks$a,ks$ter,ks$z, 
                   ks$add.bias, ks$mult.bias, ks$v.bias, ks$vratio)
  names(real_params) <- c("v","a","ter","z", "conf_add","conf_mult",'v_bias', 'vratio')
  
  subj.data <- katyal_2.df %>% filter(subj==ns, task==task.name)
  conf.ext <- subj.data$conf
  conf.ext[conf.ext<.5] <- 1-conf.ext[conf.ext<.5]
  
  if (nrow(subj.data)>80 & !is.na(real_params['v'])){
    if (mad(conf.ext)>0){
      # all RT2
      observations <-  (subj.data %>% dplyr::select(c('rt2')) %>% drop_na())/1000
      
      sim.data <- simulate_ddm(real_params, observations, additional)
      
      sim.data <- sim.data %>%dplyr::select(c(rt2, conf, accu, trialnum))
      sim.data$subj <- as.factor(subjs[ns])
      sim.data$AD <- AD[as.numeric(ns)]
      sim.data$CIT <- CIT[as.numeric(ns)]
      sim.data$age <- age[as.numeric(ns)]
      sim.data$gender <- gender[as.numeric(ns)]
      sim.data$task <- task.name
      
      all.sim.data.perc <- rbind(all.sim.data.perc, sim.data)
      
    }
    
  }
  
  task.name <- 'Memory'
  
  ks <- filter(k2.mem, subj==subjs[ns])
  real_params <- c(ks$v,ks$a,ks$ter,ks$z, 
                   ks$add.bias, ks$mult.bias, ks$v.bias, ks$vratio)
  names(real_params) <- c("v","a","ter","z", "conf_add","conf_mult",'v_bias', 'vratio')
  
  subj.data <- katyal_2.df %>% filter(subj==ns, task==task.name)
  conf.ext <- subj.data$conf
  conf.ext[conf.ext<.5] <- 1-conf.ext[conf.ext<.5]
  
  if (nrow(subj.data)>80 & !is.na(real_params['v'])){
    if (mad(conf.ext)>0){
      observations <-  (subj.data %>%dplyr::select(c('rt2')) %>% drop_na())/1000
      
      sim.data <- simulate_ddm(real_params, observations, additional)
      
      sim.data <- sim.data %>% dplyr::select(c(rt2, conf, accu, trialnum))
      sim.data$subj <- as.factor(subjs[ns])
      sim.data$AD <- AD[as.numeric(ns)]
      sim.data$CIT <- CIT[as.numeric(ns)]
      sim.data$age <- age[as.numeric(ns)]
      sim.data$gender <- gender[as.numeric(ns)]
      sim.data$task <- task.name
      
      all.sim.data.mem <- rbind(all.sim.data.mem, sim.data)
      
    }
    
  }
  
}

all.sim.data.perc$scale_rt2 <- as.numeric(scale(all.sim.data.perc$rt2))
all.sim.data.perc$conf.sc <- as.numeric(scale(all.sim.data.perc$conf))

all.sim.data.mem$scale_rt2 <- as.numeric(scale(all.sim.data.mem$rt2))
all.sim.data.mem$conf.sc <- as.numeric(scale(all.sim.data.mem$conf))

k2.sim <- rbind(all.sim.data.perc, all.sim.data.mem) %>%
  dplyr::select(c(conf.sc, AD, CIT, scale_rt2, accu, gender, age, subj, trialnum)) %>%
  mutate(gender = recode_factor(gender, 'Male' = 'Men', 'Female' = 'Women')) %>%
  mutate(data = 'Model') %>%
  mutate(subj = as.factor(as.numeric(subj) + expNum*1000))


katyal_2.sim <- k2.sim %>%
  # group_by(subj) %>%
  mutate(srt2 = ntile(scale_rt2,num.tiles.rt2)) %>%
  group_by(srt2, subj, gender, data) %>%
  dplyr::summarise(conf = mean(conf.sc, na.rm = T))

katyal_2.tile <- k2.sim %>% group_by(subj) %>%
  dplyr::summarise(AD = mean(AD),
            CIT = mean(CIT),
            age = mean(age)) %>%
  mutate(ADt = ntile(AD,2)) %>%
  mutate(ADt = recode_factor(ADt, '1' = 'Low', '2' = 'High')) %>%
  mutate(CITt = ntile(CIT,2)) %>%
  mutate(CITt = recode_factor(CITt, '1' = 'Low', '2' = 'High')) %>%
  mutate(ADt2 = ntile(AD,3)) %>%
  mutate(ADt2 = recode_factor(ADt2, '1' = 'Low', '2' = 'Mid', '3' = 'High')) %>%
  mutate(aget = ntile(age,2)) %>%
  mutate(aget = recode_factor(aget, '1' = 'Younger', '2' = 'Older'))

katyal_2.sim <- katyal_2.sim %>% dplyr::left_join(katyal_2.tile)

k2.obs <- katyal_2.df %>% 
  dplyr::select(c(conf.sc, AD, CIT, scale_rt2, gender, age, accu, subj)) %>%
  mutate(gender = recode_factor(gender, 'Male' = 'Men', 'Female' = 'Women')) %>%
  drop_na() %>%
  mutate(data = 'Observed') %>%
  mutate(subj = as.factor(as.numeric(subj) + expNum*1000))

katyal_2.obs <- k2.obs  %>%
  # group_by(subj) %>%
  mutate(srt2 = ntile(scale_rt2,num.tiles.rt2)) %>%
  group_by(srt2, subj, gender, data) %>%
  dplyr::summarise(conf = mean(conf.sc))

katyal_2.tile <- k2.obs %>% group_by(subj) %>%
  dplyr::summarise(AD = mean(AD),
            CIT = mean(CIT),
            age = mean(age))  %>%
  mutate(CITt = ntile(CIT,2)) %>%
  mutate(CITt = recode_factor(CITt, '1' = 'Low', '2' = 'High')) %>%
  mutate(ADt = ntile(AD,2)) %>%
  mutate(ADt = recode_factor(ADt, '1' = 'Low', '2' = 'High')) %>%
  mutate(ADt2 = ntile(AD,3)) %>%
  mutate(ADt2 = recode_factor(ADt2, '1' = 'Low', '2' = 'Mid', '3' = 'High')) %>%
  mutate(aget = ntile(age,2)) %>%
  mutate(aget = recode_factor(aget, '1' = 'Younger', '2' = 'Older'))

katyal_2.obs <- katyal_2.obs %>% dplyr::left_join(katyal_2.tile)

obs.val.subj <- unique(katyal_2.obs$subj)
sim.val.subj <- unique(katyal_2.sim$subj)
katyal_2.both <- rbind(filter(katyal_2.obs), filter(katyal_2.sim))


#### plot the figures

allsubj.both <- rbind(katyal_1.both, katyal_2.both)

allsubj.both <- allsubj.both %>% dplyr::rename(Data = data)

ggplot(allsubj.both, aes(x = srt2, y = conf, colour = ADt, shape = Data, alpha = Data)) +
  scale_shape_manual(breaks = c('Observed', 'Model'), values=c(19, 2)) +
  scale_alpha_manual(breaks = c('Observed', 'Model'), values=c(1, .5)) +
  scale_color_manual(breaks = c('Low', 'High'), values=c('green4', 'firebrick2')) +
  stat_summary(fun.data = mean_cl_boot, size = .6, position=position_dodge(width = .5)) +
  # geom_smooth(method='lm')+
  # coord_cartesian(ylim=c(-.1,.25)) +
  theme_classic2()+
  labs(colour = 'Anxious-\nDepression score')+
  xlab('Post-decision time (tiled)') + ylab('Confidence (z-scored)') +
  theme(text = element_text(size=font_size-4),
        legend.position = 'off',
        axis.title.y = element_text(hjust = .5)) + 
  guides(colour = guide_legend(title.position = "top", 
                               label.hjust = .5, direction = 'horizontal',
                               title.hjust = .5, ),
         shape = guide_legend(title.position = "top", 
                              label.hjust = .5, direction = 'vertical',
                              title.hjust = .5),
         alpha = 'none') 


ggplot(allsubj.both %>% drop_na(gender), 
       aes(x = srt2, y = conf, colour = gender, shape = Data, alpha = Data)) +
  scale_shape_manual(breaks = c('Observed', 'Model'), values=c(19, 2)) +
  scale_alpha_manual(breaks = c('Observed', 'Model'), values=c(1, .5)) +
  scale_color_manual(breaks = c('Men', 'Women'), values=c('green4', 'firebrick2')) +
  stat_summary(fun.data = mean_cl_boot, size = .6, position=position_dodge(width = .2)) +
  # geom_smooth(method='lm')+
  theme_classic2()+
  labs(colour = 'Gender')+
  xlab('Post-decision time (tiled)') + ylab('Confidence (z-scored)') +
  theme(text = element_text(size=font_size-4),
        legend.position = 'off') + 
  guides(colour = guide_legend(title.position = "top",
                               label.hjust = .5, direction = 'horizontal',
                               title.hjust = 0.5),
         shape = guide_legend(title.position = "top", 
                              label.hjust = .5,direction = 'vertical',
                              title.hjust = .5),
         alpha = 'none')# +
  # geom_segment(aes(x = 1.65, y = -.06, xend = 1.96, yend = -.06), color = 'grey40') +
  # geom_segment(aes(x = 1.65, y = .08, xend = 1.96, yend = .08), color = 'grey40') +
  # geom_segment(aes(x = 1.8, y = -.06, xend = 1.8, yend = .08), color = 'grey40') +
  # annotate('text', label = '****', x = 1.7, y = .01, size=8, angle=90)

ggplot(allsubj.both %>% drop_na(aget), 
       aes(x = srt2, y = conf, colour = aget, shape = Data, alpha = Data)) +
  scale_shape_manual(breaks = c('Observed', 'Model'), values=c(19, 2)) +
  scale_alpha_manual(breaks = c('Observed', 'Model'), values=c(1, .5)) +
  scale_color_manual(breaks = c('Younger', 'Older'), values=c('green4', 'firebrick2')) +
  stat_summary(fun.data = mean_cl_boot, size = .6, position=position_dodge(width = .2)) +
  theme_classic2()+
  labs(colour = 'Age')+
  xlab('Post-decision time (tiled)') + ylab('Confidence (z-scored)') +
  theme(text = element_text(size=font_size-4),
        legend.position = 'off') + 
  guides(colour = guide_legend(title.position = "top",
                               label.hjust = .5, direction = 'horizontal',
                               title.hjust = 0.5),
         shape = guide_legend(title.position = "top", 
                              label.hjust = .5,direction = 'vertical',
                              title.hjust = .5),
         alpha = F) #+
  # geom_segment(aes(x = 1.65, y = -.06, xend = 1.96, yend = -.06), color = 'grey40') +
  # geom_segment(aes(x = 1.65, y = .08, xend = 1.96, yend = .08), color = 'grey40') +
  # geom_segment(aes(x = 1.8, y = -.06, xend = 1.8, yend = .08), color = 'grey40') +
  # annotate('text', label = '****', x = 1.7, y = .01, size=8, angle=90)



########################################
########################################
######## Rouault 1

expNum <- 3

setwd("~/OneDrive - University of Copenhagen/Projects/Experiments/confidenceTime/analysis/model/meta-ddm/")

r1 <- readRDS('rouault_1-vratio.Rds')

r1 <- data.frame(r1$AD, r1$gender, r1$age, 
                 r1$recover[,7], r1$recover[,6], 
                 r1$recover[,5], r1$recover[,8],
                 c(r1$bestfnval), r1$subj, 
                 r1$recover[,1], r1$recover[,2], 
                 r1$recover[,3], r1$recover[,4])
colnames(r1) <- c('AD', 'gender', 'age', 
                  'v.bias', 'mult.bias', 'add.bias', 'vratio', 
                  'fnval', 'subj', 'v', 'a', 'ter', 'z')

excl_insuf_var <- T
## read data file
setwd("~/OneDrive - University of Copenhagen/Projects/Experiments/confidenceTime/data/Rouault_et_al_BiolPsychiat_2018/")

rouault_1.df <- read.csv('rouault_2018_exp1.csv', header = TRUE, sep = ',')

# for discrete ratings: use confidence resolution to restrict subjects with insufficient variation in confidence
conf.resolution <- 1/(length(unique(rouault_1.df$conf))-1)

rouault_1.df <- rouault_1.df %>%
  dplyr::rename(subj = id) %>%
  dplyr::rename(accu = correct) %>%
  mutate(gender = recode_factor(gender, '0' = 'Male', '1' = 'Female')) %>%
  mutate(conf = (conf-1)/10) %>%
  mutate(AD = c(scale(anxiety)))

mvsubj <- rouault_1.df %>% group_by(subj) %>%
  dplyr::summarise(sconf = sd(conf),
                   mconf = mean(conf))
insuf_var_subj <- mvsubj$subj[mvsubj$sconf<conf.resolution]

# min.accept.rt <- 100
rouault_1.df$rt1[rouault_1.df$rt1 < min.accept.rt] <- NA
# rouault_1.df$rt2[rouault_1.df$rt2 < min.accept.rt] <- NA
rouault_1.df <- rouault_1.df %>% drop_na(rt1, rt2)

rt1_max = median(rouault_1.df$rt1, na.rm=T)+ max_RT_deviation*mad(rouault_1.df$rt1, na.rm=T ) # max rt for perception task
rouault_1.df <- rouault_1.df %>%
  mutate_at(vars(contains("rt1")), ~replace(.,.>rt1_max[1],NA))
rt2_max = median(rouault_1.df$rt2, na.rm=T)+ max_RT_deviation*mad(rouault_1.df$rt2, na.rm=T ) # max rt for perception task
rouault_1.df <- rouault_1.df %>%
  mutate_at(vars(contains("rt2")), ~replace(.,.>rt2_max[1],NA))

rouault_1.df$scale_rt1 <- c(scale(rouault_1.df$rt1))
rouault_1.df$scale_rt2 <- c(scale(rouault_1.df$rt2))

rouault_1.df$stimev <- c(scale(rouault_1.df$stimev))
rouault_1.df$conf.sc <- c(scale(rouault_1.df$conf))

rouault_1.df <- rouault_1.df %>% drop_na(rt1, rt2)

if (excl_insuf_var){
  rouault_1.df <- filter(rouault_1.df, !(subj %in% insuf_var_subj))
}

AD <- r1$AD
age <- r1$age
gender <- r1$gender

# subjs <- c(1,11,12,32,93,98,103,135,153,176,182,204)

# subjs <- c(1,32,98,176,182)
ntrials <- 2000
nconflevels <- length(unique(rouault_1.df$conf))
additional <- c(0, .001, 1, ntrials, nconflevels)
names(additional) <- c('returnFit', 'dt', 'sigma', 'ntrials', 'nconflevels')

### for each study simulate all subjects
subjs <- unique(r1$subj)
nsubj <- length(subjs)
ss <- c(1:nsubj)
allsubj <- c()
all.sim.data <- c()
sim.subj <- c()
bad.subj <- c()
for (ns in ss){
  
  ks <- filter(r1, subj==subjs[ns])
  
  real_params <- c(ks$v,ks$a,ks$ter,ks$z, 
                   ks$add.bias, ks$mult.bias, ks$v.bias, ks$vratio)
  names(real_params) <- c("v","a","ter","z", "conf_add","conf_mult",'v_bias', 
                          'vratio')
  
  # real_params <- c(ks$v,ks$a,ks$ter,ks$z, 
  #                  ks$add.bias, ks$mult.bias, ks$v.bias, ks$vratio, ks$ter2)
  # names(real_params) <- c("v","a","ter","z", "conf_add","conf_mult",'v_bias', 
  #                         'vratio', 'ter2')
  
  subj.data <- rouault_1.df %>% filter(subj==subjs[ns])
  
  if (nrow(subj.data)>80 & 
      !is.na(real_params['v'] & !(subjs[ns] %in% insuf_var_subj))){
    # if (mad(conf.ext)>0){
    
    sim.subj <- rbind(sim.subj, subjs[ns])
    # all RT2
    observations <-  (subj.data %>% dplyr::select(c('rt2')) %>% drop_na())/1000
    
    sim.data <- simulate_ddm(real_params, observations, additional)
    
    sim.data$conf <- sim.data$conf/nconflevels
    
    sim.data <- sim.data %>% dplyr::select(c(rt2, conf, accu, trialnum))
    sim.data$subj <- as.factor(subjs[ns])
    sim.data$AD <- AD[as.numeric(ns)]
    sim.data$age <- age[as.numeric(ns)]
    sim.data$gender <- gender[as.numeric(ns)]
    
    all.sim.data <- rbind(all.sim.data, sim.data)
    
  }else{
    bad.subj <- rbind(bad.subj, subjs[ns])
  }
  
}

all.sim.data$scale_rt2 <- c(scale(all.sim.data$rt2))
all.sim.data$conf.sc <- c(scale(all.sim.data$conf))

r1.sim <- all.sim.data %>%
  dplyr::mutate(gender = recode_factor(gender, 'Male' = 'Men', 'Female' = 'Women')) %>%
  dplyr::select(c(conf.sc, AD, scale_rt2, accu, gender, age, subj, trialnum)) %>%
  mutate(data = 'Model') 

##
rouault_1.sim <- r1.sim %>%
  # group_by(subj) %>%
  mutate(srt2 = ntile(scale_rt2,num.tiles.rt2)) %>%
  group_by(srt2, subj, gender, data) %>%
  dplyr::summarise(conf = mean(conf.sc, na.rm = T))

rouault_1.tile <- r1.sim %>% group_by(subj) %>%
  dplyr::summarise(AD = mean(AD),
            age = mean(age))  %>%
  mutate(ADt = ntile(AD,2)) %>%
  mutate(ADt = recode_factor(ADt, '1' = 'Low', '2' = 'High')) %>%
  mutate(ADt2 = ntile(AD,3)) %>%
  mutate(ADt2 = recode_factor(ADt2, '1' = 'Low', '2' = 'Mid', '3' = 'High')) %>%
  mutate(aget = ntile(age,2)) %>%
  mutate(aget = recode_factor(aget, '1' = 'Younger', '2' = 'Older'))

rouault_1.sim <- rouault_1.sim %>% dplyr::left_join(rouault_1.tile)


r1.obs <-  filter(rouault_1.df) %>% 
  dplyr::select(c(conf.sc, AD, scale_rt2, gender, age, accu, subj)) %>%
  mutate(gender = recode_factor(gender, 'Male' = 'Men', 'Female' = 'Women')) %>%
  drop_na()  %>%
  mutate(data = 'Observed') %>%
  mutate_at(c('subj'), as.factor) 

rouault_1.obs <- r1.obs %>%
  # group_by(subj) %>%
  mutate(srt2 = ntile(scale_rt2, num.tiles.rt2)) %>%
  group_by(srt2, subj, gender, data) %>%
  dplyr::summarise(conf = mean(conf.sc))

rouault_1.tile <- r1.obs %>% group_by(subj) %>%
  dplyr::summarise(AD = mean(AD),
            age = mean(age))  %>%
  mutate(ADt = ntile(AD,2)) %>%
  mutate(ADt = recode_factor(ADt, '1' = 'Low', '2' = 'High')) %>%
  mutate(ADt2 = ntile(AD,3)) %>%
  mutate(ADt2 = recode_factor(ADt2, '1' = 'Low', '2' = 'Mid', '3' = 'High')) %>%
  mutate(aget = ntile(age,2)) %>%
  mutate(aget = recode_factor(aget, '1' = 'Younger', '2' = 'Older'))

rouault_1.obs <- rouault_1.obs %>% dplyr::left_join(rouault_1.tile)

obs.val.subj <- unique(rouault_1.obs$subj)
sim.val.subj <- unique(rouault_1.sim$subj)
# val.subj <- intersect(obs.val.subj, sim.val.subj)
# val.subj <- union(obs.val.subj, sim.val.subj)
rouault_1.both <- rbind(filter(rouault_1.obs), 
                        filter(rouault_1.sim))

########################################
########################################
######## Rouault 2

expNum <- 4

setwd("~/OneDrive - University of Copenhagen/Projects/Experiments/confidenceTime/analysis/model/meta-ddm/")

r2 <- readRDS('rouault_2-vratio.Rds') 

r2 <- data.frame(r2$AD, r2$gender, r2$age,
                 r2$recover[,7], r2$recover[,6],
                 r2$recover[,5], r2$recover[,8],
                 c(r2$bestfnval), r2$subj,
                 r2$recover[,1], r2$recover[,2],
                 r2$recover[,3], r2$recover[,4])
colnames(r2) <- c('AD', 'gender', 'age',
                  'v.bias', 'mult.bias', 'add.bias', 'vratio',
                  'fnval', 'subj', 'v', 'a', 'ter', 'z')

excl_insuf_var <- T
## read data file
setwd("~/OneDrive - University of Copenhagen/Projects/Experiments/confidenceTime/data/Rouault_et_al_BiolPsychiat_2018//")

rouault_2.df <- read.csv('rouault_2018_exp2.csv', header = TRUE, sep = ',') %>% drop_na()

nconflevels <- max(unique(rouault_2.df$conf), na.rm = T)
# for discrete ratings: use confidence resolution to restrict subjects with insufficient variation in confidence
conf.resolution <- 1/(length(unique(rouault_2.df$conf))-1)

rouault_2.df <- rouault_2.df %>%
  dplyr::rename(CIT = Compul) %>%
  dplyr::rename(accu = correct) %>%
  mutate(gender = recode_factor(gender, '0' = 'Male', '1' = 'Female')) %>% 
  drop_na() %>%
  mutate(conf = (conf-1)/5)

mvsubj <- rouault_2.df %>% group_by(subj) %>%
  dplyr::summarise(sconf = sd(conf),
                   mconf = mean(conf))
insuf_var_subj <- mvsubj$subj[mvsubj$sconf<conf.resolution]

rouault_2.df$rt1[rouault_2.df$rt1 < min.accept.rt] <- NA
# rouault_2.df$rt2[rouault_2.df$rt2 < min.accept.rt] <- NA
rouault_2.df <- rouault_2.df %>% drop_na(rt1, rt2)

max_RT_deviation = 3
rt1_max = median(rouault_2.df$rt1, na.rm=T) + max_RT_deviation*mad(rouault_2.df$rt1, na.rm=T ) # max rt for perception task
rouault_2.df <- rouault_2.df %>%
  mutate_at(vars(contains("rt1")), ~replace(.,.>rt1_max[1],NA))
rt2_max = median(rouault_2.df$rt2, na.rm=T)+ max_RT_deviation*mad(rouault_2.df$rt2, na.rm=T ) # max rt for perception task
rouault_2.df <- rouault_2.df %>%
  mutate_at(vars(contains("rt2")), ~replace(.,.>rt2_max[1],NA))

rouault_2.df <- rouault_2.df %>% drop_na(rt1, rt2, conf)

# rouault_2.df$scale_rt1 <- c(scale(rouault_2.df$rt1))
rouault_2.df$scale_rt2 <- c(scale(rouault_2.df$rt2))

rouault_2.df$conf.sc <- c(scale(rouault_2.df$conf))

if (excl_insuf_var){
  rouault_2.df <- filter(rouault_2.df, !(subj %in% insuf_var_subj))
  length(unique(rouault_2.df))
}


AD <- r2$AD
age <- r2$age
gender <- r2$gender

# subjs <- c(1,32,98,176,182)
additional <- c(0, .001, 1, ntrials, nconflevels)
names(additional) <- c('returnFit', 'dt', 'sigma', 'ntrials', 'nconflevels')

### for each study simulate all subjects
subjs <- unique(r2$subj)
nsubj <- length(subjs)
ss <- c(1:nsubj)
allsubj <- c()
sim.subj <- c()
all.sim.data <- c()
for (ns in ss){
  
  ks <- filter(r2, subj==subjs[ns])
  
  real_params <- c(ks$v,ks$a,ks$ter,ks$z,
                   ks$add.bias, ks$mult.bias, ks$v.bias, ks$vratio)
  names(real_params) <- c("v","a","ter","z", "conf_add","conf_mult",'v_bias', 'vratio')
  
  subj.data <- rouault_2.df %>% filter(subj==subjs[ns])
  
  if (nrow(subj.data)>80 & 
      !is.na(real_params['v'] & !(subjs[ns] %in% insuf_var_subj))){
    # if (mad(conf.ext)>0){
    # all RT2
    observations <-  (subj.data %>% dplyr::select(c('rt2')) %>% drop_na())/1000
    
    sim.subj <- rbind(sim.subj, subjs[ns])
    
    sim.data <- simulate_ddm(real_params, observations, additional)
    sim.data$conf <- sim.data$conf/nconflevels
    
    sim.data <- sim.data %>% dplyr::select(c(rt2, conf, accu, trialnum))
    sim.data$subj <- as.factor(subjs[ns])
    sim.data$AD <- AD[as.numeric(ns)]
    sim.data$age <- age[as.numeric(ns)]
    sim.data$gender <- gender[as.numeric(ns)]
    
    all.sim.data <- rbind(all.sim.data, sim.data)
    # }
  }
}

all.sim.data$scale_rt2 <- c(scale(all.sim.data$rt2))
all.sim.data$conf.sc <- c(scale(all.sim.data$conf))

r2.sim <- all.sim.data %>%
  mutate(gender = recode_factor(gender, 'Male' = 'Men', 'Female' = 'Women')) %>%
  dplyr::select(c(conf.sc, AD, scale_rt2, accu, gender, age, subj)) %>%
  mutate(data = 'Model')


rouault_2.sim <- r2.sim %>%     
  # group_by(subj) %>%
  mutate(srt2 = ntile(scale_rt2,num.tiles.rt2)) %>%
  group_by(srt2, subj, gender, data) %>%
  dplyr::summarise(conf = mean(conf.sc, na.rm = T))

rouault_2.tile <- r2.sim %>% group_by(subj) %>%
  dplyr::summarise(AD = mean(AD),
            age = mean(age))  %>%
  mutate(ADt = ntile(AD,2)) %>%
  mutate(ADt = recode_factor(ADt, '1' = 'Low', '2' = 'High'))%>%
  mutate(ADt2 = ntile(AD,3)) %>%
  mutate(ADt2 = recode_factor(ADt2, '1' = 'Low', '2' = 'Mid', '3' = 'High')) %>%
  mutate(aget = ntile(age,2)) %>%
  mutate(aget = recode_factor(aget, '1' = 'Younger', '2' = 'Older'))

rouault_2.sim <- rouault_2.sim %>% dplyr::left_join(rouault_2.tile)


r2.obs <- filter(rouault_2.df, subj %in% sim.subj) %>% 
  dplyr::select(c(conf.sc, AD, scale_rt2, gender, age, accu, subj)) %>%
  mutate(gender = recode_factor(gender, 'Male' = 'Men', 'Female' = 'Women')) %>%
  drop_na() %>%
  mutate(data = 'Observed') %>%
  mutate_at(c('subj'), as.factor)


rouault_2.obs <- r2.obs %>%
  # group_by(subj) %>%
  mutate(srt2 = ntile(scale_rt2,num.tiles.rt2)) %>%
  group_by(srt2, subj, gender, data) %>%
  dplyr::summarise(conf = mean(conf.sc))

rouault_2.tile <- r2.obs %>% group_by(subj) %>%
  dplyr::summarise(AD = mean(AD),
            age = mean(age))  %>%
  mutate(ADt = ntile(AD,2)) %>%
  mutate(ADt = recode_factor(ADt, '1' = 'Low', '2' = 'High')) %>%
  mutate(ADt2 = ntile(AD,3)) %>%
  mutate(ADt2 = recode_factor(ADt2, '1' = 'Low', '2' = 'Mid', '3' = 'High')) %>%
  mutate(aget = ntile(age,2)) %>%
  mutate(aget = recode_factor(aget, '1' = 'Younger', '2' = 'Older'))

rouault_2.obs <- rouault_2.obs %>% dplyr::left_join(rouault_2.tile)

obs.val.subj <- unique(rouault_2.obs$subj)
sim.val.subj <- unique(rouault_2.sim$subj)

# val.subj <- union(obs.val.subj, sim.val.subj)
rouault_2.both <- rbind(filter(rouault_2.obs),filter(rouault_2.sim))

##### Supplementary Figure 6 

allsubj.both <- rbind(rouault_1.both, rouault_2.both)

allsubj.both <- allsubj.both %>% dplyr::rename(Data = data)


ggplot(allsubj.both, 
       aes(x = srt2, y = conf, colour = ADt, shape = Data, alpha = Data)) +
  scale_shape_manual(breaks = c('Observed', 'Model'), values=c(19, 2)) +
  scale_alpha_manual(breaks = c('Observed', 'Model'), values=c(1, .4)) +
  scale_color_manual(breaks = c('Low', 'High'), values=c('green4', 'firebrick2')) +
  stat_summary(fun.data = mean_cl_boot, size = .6, position=position_dodge(width = .5)) +
  theme_classic2()+
  labs(colour = 'Anxious-Depression\nscore')+
  xlab('Post-decision time (tiled)') + ylab('Confidence (z-scored)') +
  theme(text = element_text(size=font_size-4),
        legend.position = 'off') + 
  guides(colour = guide_legend(title.position = "top",
                               label.hjust = .5, direction = 'horizontal',
                               title.hjust = 0.5),
         shape = guide_legend(title.position = "top", 
                              label.hjust = .5,direction = 'vertical',
                              title.hjust = .5),
         alpha = F) #+
  # geom_segment(aes(x = 4, y = -.065, xend = 4.8, yend = -.08), color = 'grey40') +
  # geom_segment(aes(x = 4, y = .033, xend = 4.8, yend = .033), color = 'grey40') +
  # geom_segment(aes(x = 4.5, y = -.075, xend = 4.5, yend = .033), color = 'grey40') +
  # annotate('text', label = '****', x = 4.8, y = -.02, size=8, angle=90)


ggplot(allsubj.both %>% drop_na(gender), 
       aes(x = srt2, y = conf, colour = gender, shape = Data, alpha = Data)) +
  scale_shape_manual(breaks = c('Observed', 'Model'), values=c(19, 2)) +
  scale_alpha_manual(breaks = c('Observed', 'Model'), values=c(1, .5)) +
  scale_color_manual(breaks = c('Men', 'Women'), values=c('green4', 'firebrick2')) +
  stat_summary(fun.data = mean_cl_boot, size = .6, position=position_dodge(width = .2)) +
  # geom_smooth(method='lm')+
  theme_classic2()+
  labs(colour = 'Gender')+
  xlab('Post-decision time (tiled)') + ylab('Confidence (z-scored)') +
  theme(text = element_text(size=font_size-4),
        legend.position = 'off') + 
  guides(colour = guide_legend(title.position = "top",
                               label.hjust = .5, direction = 'horizontal',
                               title.hjust = 0.5),
         shape = guide_legend(title.position = "top", 
                              label.hjust = .5,direction = 'vertical',
                              title.hjust = .5),
         alpha = F) #+
  # geom_segment(aes(x = 1.65, y = -.06, xend = 1.96, yend = -.06), color = 'grey40') +
  # geom_segment(aes(x = 1.65, y = .08, xend = 1.96, yend = .08), color = 'grey40') +
  # geom_segment(aes(x = 1.8, y = -.06, xend = 1.8, yend = .08), color = 'grey40') +
  # annotate('text', label = '****', x = 1.7, y = .01, size=8, angle=90)

ggplot(allsubj.both %>% drop_na(aget), 
       aes(x = srt2, y = conf, colour = aget, shape = Data, alpha = Data)) +
  scale_shape_manual(breaks = c('Observed', 'Model'), values=c(19, 2)) +
  scale_alpha_manual(breaks = c('Observed', 'Model'), values=c(1, .5)) +
  scale_color_manual(breaks = c('Younger', 'Older'), values=c('green4', 'firebrick2')) +
  stat_summary(fun.data = mean_cl_boot, size = .6, position=position_dodge(width = .2)) +
  # geom_smooth(method='lm')+
  theme_classic2()+
  labs(colour = 'Age')+
  xlab('Post-decision time (tiled)') + ylab('Confidence (z-scored)') +
  theme(text = element_text(size=font_size-4),
        legend.position = 'off') + 
  guides(colour = guide_legend(title.position = "top",
                               label.hjust = .5, direction = 'horizontal',
                               title.hjust = 0.5),
         shape = guide_legend(title.position = "top", 
                              label.hjust = .5,direction = 'vertical',
                              title.hjust = .5),
         alpha = F) #+
  # geom_segment(aes(x = 1.65, y = -.06, xend = 1.96, yend = -.06), color = 'grey40') +
  # geom_segment(aes(x = 1.65, y = .08, xend = 1.96, yend = .08), color = 'grey40') +
  # geom_segment(aes(x = 1.8, y = -.06, xend = 1.8, yend = .08), color = 'grey40') +
  # annotate('text', label = '****', x = 1.7, y = .01, size=8, angle=90)



allsubj.obs <- rbind(k1.obs, k2.obs %>% dplyr::select(-CIT), r1.obs, r2.obs) %>%
  mutate(srt2 = ntile(scale_rt2,num.tiles.rt2)) %>%
  mutate(ADt = ntile(AD,2)) %>%
  mutate(ADt = recode_factor(ADt, '1' = 'Low', '2' = 'High')) %>%
  mutate(aget = ntile(age,2)) %>%
  mutate(aget = recode_factor(aget, '1' = 'Low', '2' = 'High')) %>%
  group_by(srt2, subj, ADt, aget, gender) %>%
  dplyr::summarise(
    conf = mean(conf.sc))

ggplot(allsubj.both, aes(x = srt2, y = conf, colour = ADt)) +
  scale_color_manual(breaks = c('Low', 'High'), values=c('green4', 'firebrick2')) +
  stat_summary(fun.data = mean_se, size = .6, position=position_dodge(width = .2)) +
  # geom_smooth(method='lm')+
  # coord_cartesian(ylim=c(-.1,.25)) +
  theme_classic2()+
  labs(colour = 'AD score')+
  xlab('Time after decision (z-scored & tiled)') + ylab('Confidence (z-scored)') +
  theme(text = element_text(size=font_size-2),
        # legend.position = c(.65,.1),
        legend.position = 'off',
        legend.direction = 'horizontal')

