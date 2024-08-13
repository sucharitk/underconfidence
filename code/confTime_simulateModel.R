
################################## Katyal & Fleming ################################
### Gender and anxiety reveal distinct computational sources of underconfidence ####
####################################################################################
################  This file reproduces Figure 2 ##################
###############################################################################
####### Model simulations   ##############
###############################################################################


library(lmerTest) # anova and summary with p-values
library(emmeans) # Least squares means
library(sjPlot)
library(ggpubr)
library(tidyverse)
library(ggplot2)
library(pracma)
library(patchwork)

font_size <- 14
annotate_font_size <- 5

######################################
######################################
##### Define functions   ######
######################################

###################################
#### DDM exp simulation function ######
###################################

simDDMBiasExp <- function(exp.params){
  
  vars.in <- names(exp.params)
  
  if ('conf_add' %in% vars.in){ # when specified as exact value
    conf_add <- as.numeric(exp.params['conf_add'])
  } else{ # when specified as distribution
    conf_add <- as.numeric(exp.params['conf_add_mean'])
  }
  if ('conf_sd' %in% vars.in){ # when specified as exact value
    conf_sd <- as.numeric(exp.params['conf_sd'])
  } else{ # when specified as distribution
    conf_sd <- as.numeric(exp.params['conf_sd_mean'])
  }
  
  beta_v <-  as.numeric(exp.params['beta_v'])
  beta_add <- as.numeric(exp.params['beta_add'])
  beta_sd <- as.numeric(exp.params['beta_sd'])
  
  nTrials = exp.params['nTrials']
  nSubj   = exp.params['nSubj']
  
  dt      = as.numeric(exp.params['dt'])
  sqrt_dt = sqrt(dt)
  
  mu_a    = 1
  sigma_a = .05
  
  mu_z    = mu_a/2
  sigma_z = sigma_a/8
  
  mu_v    = 1
  sigma_v = .2
  
  sigma   = 1
  
  phq_dist = c(1.24972, 4.36897)
  
  exp.rt   = matrix(0, nSubj, nTrials)
  exp.accu = exp.rt
  exp.resp = exp.rt
  exp.stim = exp.rt
  exp.conf = exp.rt
  exp.rt2  = exp.rt
  
  phq <- matrix(1, nSubj)
  conf.sd <- phq
  conf.mean <- phq
  
  for (ns in 1:nSubj){
    # subject loop
    
    phq[ns] = round(rgamma(1, phq_dist[1], 1/phq_dist[2]))

    if ('conf_add_mean' %in% vars.in){
      conf.mean[ns] <- rnorm(1,conf_add,2)
    } else{
      conf.mean[ns] <- conf_add
    }
    if ('conf_sd_mean' %in% vars.in){
      conf.sd[ns] <- rgamma(1,conf_sd,5)
    } else{
      conf.sd[ns] <- conf_sd
    }
    
    v_subj = rnorm(1, mu_v, sigma_v)
    a_subj = rnorm(1, mu_a, sigma_a)
    z_subj = rnorm(1, mu_z, sigma_z)
    
    for (nt in 1:nTrials){
      
      trial.params <- c(v_subj, a_subj, z_subj, 
                        phq[ns], dt, sqrt_dt,
                        conf.mean[ns], conf.sd[ns], sigma,
                        beta_v, beta_add, beta_sd)
      
      names(trial.params) <- c('v', 'a', 'z', 'phq', 'dt', 'sqrt_dt',
                               'conf_add', 'conf_sd', 'sigma',
                               'beta_v', 'beta_add', 'beta_sd')
      simTrial <- simDDMBiasTrial(trial.params)
      
      exp.accu[ns,nt] <- simTrial['trial.accu']
      exp.resp[ns,nt] <- simTrial['trial.resp']
      exp.stim[ns,nt] <- simTrial['trial.stim']
      exp.rt[ns,nt]   <- simTrial['trial.rt']
      exp.rt2[ns,nt]  <- simTrial['trial.rt2']
      exp.conf[ns,nt] <- simTrial['trial.conf']
    }
  }
  
  subj <- kronecker(matrix(1,1,nTrials), c(1:nSubj))
  trialnum <- t(kronecker(matrix(1,1,nSubj), c(1:nTrials)))
  phq2 <- kronecker(matrix(1,1,nTrials), c(phq))
  
  df = data.frame(c(exp.rt), c(exp.accu), c(exp.resp), c(exp.stim), 
                  c(exp.rt2), c(exp.conf), c(subj),
                  c(trialnum), c(phq2))
  names(df) <- c('rt', 'accu', 'resp', 'stim', 'rt2', 'conf', 'subj', 
                 'trialnum', 'phq')
  
  return(df)
}

###################################
#### DDM trial simulation function ######
###################################

simDDMBiasTrial <- function(trial.params){
  
  z       <- trial.params['z']
  v       <- trial.params['v']
  a       <- trial.params['a']
  sigma   <- trial.params['sigma']
  dt      <- trial.params['dt']
  sqrt_dt <- trial.params['sqrt_dt']
  phq     <- trial.params['phq']
  
  conf_add <- trial.params['conf_add']
  conf_sd  <- trial.params['conf_sd']
  
  beta_v   <- trial.params['beta_v']
  beta_add <- trial.params['beta_add']
  beta_sd  <- trial.params['beta_sd']
  
  gam_t2  = c(5, .1) # for simulating confidence response times
  gam_t2  = c(7, .033) # for simulating confidence response times
  
  trial.stim = 2*round(runif(1))- 1 # stimulus can be [-1, 1]
  
  evidence = a*z
  choiceTime = 0
  trial.accu <- 0
  v_trial = trial.stim * v
  sigma_sqrt_dt <- sigma*sqrt_dt
  
  while(evidence>=0 && evidence<=a){
    dev = v_trial*dt + sigma_sqrt_dt*rnorm(1)
    
    choiceTime = choiceTime + dt
    evidence = evidence + dev
  }
  
  trial.rt = choiceTime;
  
  if (evidence >= a){
    trial.resp = 1;
    if (trial.stim > 0){
      trial.accu = 1;  
    }
    
  } else{
    trial.resp = -1;
    if (trial.stim < 0){
      trial.accu = 1;
    }
  }
  
  ## post-decisional processing for reporting confidence
  trial.rt2 = rgamma(1, gam_t2[1], 1/gam_t2[2])

  vpd = v_trial + trial.resp*beta_v*phq # assume vratio=1
  
  evidence <- evidence + rnorm(1, vpd*trial.rt2, sqrt(trial.rt2))
  
  if (trial.resp > 0){
    trial.conf = evidence - a;
  }    else{
    trial.conf = -evidence;
  }
  
  addfact = conf_add + beta_add*phq
  multfact = conf_sd + beta_sd*phq
  
  trial.conf = 1/(1+exp(-(trial.conf-addfact)*multfact)); # logistic transform to get ratings between 0 and 1
  
  df = c(trial.rt, trial.accu, trial.resp, trial.stim, trial.rt2, trial.conf)
  names(df) <- c('trial.rt', 'trial.accu', 'trial.resp', 'trial.stim', 'trial.rt2', 'trial.conf')
  
  return(df)
}


sigmoid = function(params, x) {
  params[1] / (1 + exp(-params[2] * (x - params[3])))
}


simDDMBiasTrialPlot <- function(trial.params){
  
  z       <- trial.params['z']
  v       <- trial.params['v']
  a       <- trial.params['a']
  sigma   <- trial.params['sigma']
  dt      <- trial.params['dt']
  cov1     <- trial.params['cov1']
  cov2     <- trial.params['cov2']
  
  conf_add <- trial.params['conf_add']
  conf_sd  <- trial.params['conf_sd']
  
  beta_v   <- trial.params['beta_v']
  beta_add <- trial.params['beta_add']
  beta_sd  <- trial.params['beta_sd']
  
  cov <- c(as.numeric(cov1),as.numeric(cov2))
  
  sqrt_dt <- sqrt(dt)
  gam_t2  = c(2, .1) # for simulating confidence response times
  
  trial.stim = 2*round(runif(1))- 1 # stimulus can be [-1, 1]
  
  sampNum <- 0
  evChoice <- c()
  
  evidence = z
  choiceTime = 0
  trial.accu <- 0
  v_trial = trial.stim * v
  sigma_sqrt_dt <- sigma*sqrt_dt
  
  while(evidence>=0 && evidence<=a){
    dev = v_trial*dt + sigma_sqrt_dt*rnorm(1)
    choiceTime = choiceTime + dt
    evidence = evidence + dev
    
    sampNum <- sampNum + 1
    
    evChoice[sampNum] <- evidence
  }
  
  trial.rt = choiceTime;
  
  if (evidence >= a){
    trial.resp = 1;
    evidence <- a
    if (trial.stim > 0){
      trial.accu = 1;  
    }
    
  } else{
    trial.resp = -1
    evidence <- 0
    if (trial.stim < 0){
      trial.accu = 1;
    }
  }
  evChoice[sampNum] <- evidence
  
  ## post-decisional processing for reporting confidence
  evConf <- array(0,c(2,1))
  evConf[,1] <- evidence
  # sampNum <- 2
  
  t2 = rgamma(1, gam_t2[1], 1/gam_t2[2])
  postDecisionTime = choiceTime
  
  vs <- trial.resp*beta_v*cov
  
  nloop <- round(t2/dt)
  for (tt in 1:nloop){
    dev = (v_trial + vs)*dt + sigma_sqrt_dt*rnorm(1)
    
    postDecisionTime = postDecisionTime + dt
    
    evConf <- cbind(evConf, t(last(t(evConf))) + dev)
  }
  
  trial.rt2 = postDecisionTime-choiceTime
  evidence <- last(t(evConf))
  if (trial.resp > 0){
    trial.conf = evidence - a;
  }    else{
    trial.conf = -evidence;
  }
  
  addfact = conf_add + beta_add*cov
  multfact = conf_sd + beta_sd*cov
  
  trial.conf = 1/(1+exp(-(trial.conf-addfact)*multfact)); # logistic transform to get ratings between 0 and 1
  
  df <- list('evChoice' = evChoice, 'evConf' = evConf, 
             'trial.conf' = trial.conf, 'choiceTime' = choiceTime,
             'confTime' = postDecisionTime)
  
  return(df)
}

# 
# ######################################
# ######################################
# ##### #####   Figure 1. (previous) #####  #####  
# ######################################
# 
# 
# nTrials  = 100
# nSubj    = 200
# dt       = .005
# ylabel <- 'How likely was I correct?\n(Confidence)\nUnlikely             Very likely'
# 
# ################
# ### A. Distortion increases over time
# exp.params <- c(-.2, 2, -.08, 0,0, nTrials, nSubj, dt)
# names(exp.params) <- c('conf_add', 'conf_sd', 
#                        'beta_v', 'beta_add', 'beta_sd',
#                        'nTrials', 'nSubj', 'dt')
# 
# simexp <- simDDMBiasExp(exp.params)
# 
# m.reg <- lm(conf ~ phq*rt2, simexp)
# emm <- emtrends(m.reg, var = "rt2" ,pairwise~phq, at=list(phq=c(11.6,.6)))
# cr1 <- summary(emm$contrasts)
# cr1$SE <- cr1$SE*2
# 
# plot_model(m.reg, type = ("pred"), terms = c("rt2", "phq"), title='',
#                   colors = c('green4', 'darkgoldenrod3', 'firebrick2'),
#                   line.size = 2, ci.lvl = 0, auto.label = T)+
#   theme_classic2()+
#   xlab('Time after decision (seconds)')+
#   ylab(ylabel)+
#   coord_cartesian(ylim = c(.42,.72)) +
#   labs(colour = 'AD')+
#   theme(text = element_text(size=font_size),
#         legend.position = 'none',
#         axis.ticks = element_blank(),
#         axis.text = element_blank())
# 
# 
# ################
# ### B. Distortion stays the same over time
# exp.params <- c(0, .3, 0, .1, 0, nTrials, nSubj, dt)
# 
# names(exp.params) <- c('conf_add', 'conf_sd', 
#                        'beta_v', 'beta_add', 'beta_sd',
#                        'nTrials', 'nSubj', 'dt')
# 
# simexp <- simDDMBiasExp(exp.params)
# 
# m.reg <- lm(conf ~ phq*rt2, simexp)
# emm <- emtrends(m.reg, var = "rt2" ,pairwise~phq, at=list(phq=c(11.6,.6)))
# cr2 <- summary(emm$contrasts)
# cr2$SE <- cr2$SE*12
# 
# plot_model(m.reg, type = ("pred"), terms = c("rt2", "phq"), title='',
#                   colors = c('green4', 'darkgoldenrod3', 'firebrick2'),
#                   line.size = 2, ci.lvl = 0, auto.label = T)+
#   theme_classic2()+
#   xlab('Time after decision (seconds)')+
#   ylab(ylabel)+
#   coord_cartesian(ylim = c(.41,.54)) +
#   labs(colour = 'AD')+
#   theme(text = element_text(size=font_size),
#         legend.position = 'none',
#         axis.ticks = element_blank(),
#         axis.text = element_blank())
# 
# ################
# ### C. Distortion decreases over time
# exp.params <- c(-1.2, 6, 0, .15,0, nTrials, nSubj, dt)
# 
# names(exp.params) <- c('conf_add', 'conf_sd', 
#                        'beta_v', 'beta_add', 'beta_sd',
#                        'nTrials', 'nSubj', 'dt')
# 
# simexp <- simDDMBiasExp(exp.params)
# 
# m.reg <- lm(conf ~ phq*rt2, simexp)
# emm <- emtrends(m.reg, var = "rt2" ,pairwise~phq, at=list(phq=c(11.6,.6)))
# cr3 <- summary(emm$contrasts)
# cr3$SE <- cr3$SE*3
# 
# plot_model(m.reg, type = ("pred"), terms = c("rt2", "phq"), title='',
#                   colors = c('green4', 'darkgoldenrod3', 'firebrick2'),
#                   line.size = 2, ci.lvl = 0, auto.label = T)+
#   theme_classic2()+
#   xlab('Time after decision (seconds)')+
#   ylab(ylabel)+
#   labs(colour = 'AD')+
#   coord_cartesian(ylim = c(.35,1.1)) +
#   theme(text = element_text(size=font_size),
#         legend.position = 'none',
#         axis.ticks = element_blank(),
#         axis.text = element_blank())
# 
# 
# 
# ggplot(cr1, aes(y = estimate, x = 0)) +
#   geom_bar(stat = 'identity', colour = 'black', fill = 'mistyrose')+
#   geom_errorbar(aes(ymin = estimate - SE, ymax = estimate + SE),
#                 linewidth = .75, width = .5, color = 'black')+
#   geom_hline(yintercept = 0, color = 'black') +
#   coord_cartesian(ylim = c(-.4,.4), xlim = c(-1,1))   +
#   scale_y_continuous(breaks = 0)+
#   theme_classic2()+
#   ylab('Slope difference\n(high - low)') +
#   theme(text = element_text(size=font_size),
#         legend.position = 'none',
#         axis.ticks = element_blank(),
#         axis.text.x = element_blank(),
#         axis.title.x = element_blank())
# 
# ggplot(cr2, aes(y = estimate, x = 0)) +
#   geom_bar(stat = 'identity', colour = 'black', fill = 'mistyrose')+
#   geom_errorbar(aes(ymin = estimate - SE, ymax = estimate + SE),
#                 linewidth = .75, width = .5, color = 'black')+
#   geom_hline(yintercept = 0, color = 'black') +
#   coord_cartesian(ylim = c(-.3,.3), xlim = c(-1,1))   +
#   scale_y_continuous(breaks = 0)+
#   theme_classic2()+
#   ylab('Slope difference\n(high - low)') +
#   theme(text = element_text(size=font_size),
#         legend.position = 'none',
#         axis.ticks = element_blank(),
#         axis.text.x = element_blank(),
#         axis.title.x = element_blank())
# 
# ggplot(cr3, aes(y = estimate, x = 0)) +
#   geom_bar(stat = 'identity', colour = 'black', fill = 'mistyrose')+
#   geom_errorbar(aes(ymin = estimate - SE, ymax = estimate + SE),
#                 linewidth = .75, width = .5, color = 'black')+
#   geom_hline(yintercept = 0, color = 'black') +
#   coord_cartesian(ylim = c(-.7,.7), xlim = c(-1,1)) +
#   scale_y_continuous(breaks = 0)+
#   theme_classic2()+
#   ylab('Slope difference\n(high - low)') +
#   theme(text = element_text(size=font_size),
#         legend.position = 'none',
#         axis.ticks = element_blank(),
#         axis.text.x = element_blank(),
#         axis.title.x = element_blank())
# 
# # create dummy data to make colour bars
# dum <- data.frame(c(1:1000), runif(1000)-.5)
# colnames(dum) <- c('x','y')
# legend.labs.hi <- c('High', 'High', 'Older')
# legend.labs.lo <- c('Low', 'Low', 'Younger')
# ind.fact <- 1
# legend.titles <- c('Individual difference measure related to underconfidence', 'Anxious-Depression', '-Compulsivity', 'Age')
# ggplot(dum, aes(y = y, x = x, colour = y)) +
#   geom_point(size = 2.2, shape = 21)+
#   theme_classic2()+
#   labs(colour = legend.titles[ind.fact]) +
#   scale_colour_gradient2(low='green4', high='firebrick2', mid = 'darkgoldenrod3',
#                        breaks = c(-.45, .45),
#                        labels = c(legend.labs.lo[ind.fact], 
#                                   legend.labs.hi[ind.fact])) +
#   theme(text = element_text(size=font_size+2),
#         legend.position = 'top') +
#   guides(colour = guide_colourbar(title.position="top", title.hjust = .5,
#                                   barwidth = 15, barheight = 1.5))
# 
# dum <- data.frame(c(1:1000), as.factor((round(runif(1000))-.5)))
# colnames(dum) <- c('x','y')
# ggplot(dum, aes(y = y, x = x, fill = y)) +
#   geom_tile() + theme_classic2()+
#   scale_fill_manual(values=c('green4', 'firebrick2'), 
#                        labels = c('Male', 'Female')) +
#   labs(fill = 'Gender') +
#   theme(text = element_text(size=font_size),
#         legend.position = 'top') +
#   guides(fill = guide_legend(title.position="top", title.hjust = 0.5,
#                              label.position = 'bottom', label.hjust = .5,
#                                   keywidth = 4, keyheight = 2))


######################################
######################################


n_models <- 3

nTrials  = 200
nSubj    = 200
dt       = .005

disc_conf <- 0

conf_add_0    <- seq(-2,2,1)
conf_div_0    <- c(.5, 1, 2, 4)
# conf_add_0    <- -1
# conf_div_0    <- 1 # gamma params 5,5 to mimic a mean of 1

# 
v_beta        <- -c(.05, .1, .15, .2)
addbias_beta  <- c(.05, .1, .15, .2)
multbias_beta <- -c(.05, .1, .15, .2)

n_add0        <- length(conf_add_0)
n_div0        <- length(conf_div_0)
n_beta        <- length(addbias_beta) # keeping number of conditions to test for the two models the same

underconf   = array(NA, c(n_models, n_add0, n_div0, n_beta,2)) # does this combination of parameters model underconfidence
underconf.se = underconf
pval.main   = underconf
conf.accu   = underconf
pval.accu   = underconf

slope.dif   = array(NA, c(n_models, n_add0, n_div0, n_beta, 2))
pval.inter  = slope.dif
slope.dif.accu   = array(NA, c(n_models, n_add0, n_div0, n_beta, 2))
slope.dif.se.accu <- slope.dif.accu
pval.inter.accu  = slope.dif.accu

conf.sim    = array(NA, c(n_models, n_add0, n_div0, n_beta, nSubj))
conf.sd.sim <- conf.sim
accu.sim <- conf.sim
rt1.sim  <- conf.sim
rt2.sim  <- conf.sim
phq.sim  <- conf.sim

model.sim <- array(c(1:n_models), c(n_models, n_add0, n_div0, n_beta,2))
ca0.sim   <- aperm(array(conf_add_0, c(n_add0, n_models, n_div0, n_beta, 2)), c(2,1,3:5))
cd0.sim   <- aperm(array(conf_div_0, c(n_div0, n_models, n_add0, n_beta, 2)), c(2,3,1,4,5))
beta.sim  <- aperm(array(c(1:n_beta), c(n_beta, n_models,  n_add0, n_div0,2)), c(2,3,4,1,5))
corr.sim  <- aperm(array(c(0:1), c(2, n_models,  n_add0, n_div0, n_beta)), c(2,3,4,5,1))

nm <- 1
ca0 <- 1
cd0 <- 1
nb <- 1

model.run <- 0

for (nm in 1:n_models){
  for (ca0 in 1:n_add0){
    for (cd0 in 1:n_div0){
      
      for (nb in 1:n_beta){
        
        model.run <- model.run + 1
        print(sprintf('%d/%d: nm=%d/%d, ca0=%d/%d, cd0=%d/%d, nb=%d/%d', 
                      model.run, n_models*n_beta*n_add0*n_div0, nm, n_models, 
                      ca0, n_add0, cd0, n_div0, nb, n_beta))
        
        modelNumber <- nm
        
        exp.params <- c(conf_add_0[ca0], conf_div_0[cd0],
                        switch(modelNumber, v_beta[nb], 0, 0),
                        switch(modelNumber, 0, addbias_beta[nb], 0),
                        switch(modelNumber, 0, 0, multbias_beta[nb]),
                        nTrials, nSubj, dt)
        names(exp.params) <- c('conf_add', 'conf_sd', 
                               'beta_v', 'beta_add', 'beta_sd',
                               'nTrials', 'nSubj', 'dt')
        
        simexp <- simDDMBiasExp(exp.params)
        
        if (disc_conf){
          simexp$conf <- (ceil(simexp$conf*disc_conf)-1)/(disc_conf-1)
        }
        
        m.reg <- lm(conf ~ phq, simexp)
        summary(m.reg)
        plot_model(m.reg, type = ("pred"),
                   terms = c("phq"),
                   title='')+
          theme_pubclean()+
          xlab('AD')+
          ylab('Confidence')+
          theme(text = element_text(size=font_size))

        s.reg <- summary(m.reg)$coefficients
        underconf[nm,ca0,cd0,nb,] <- s.reg[2,1] # save effect size
        underconf.se[nm,ca0,cd0,nb,] <- s.reg[2,2] # save effect size
        pval.main[nm,ca0,cd0,nb,] <- s.reg[2,4] # save p avlue
        
        m.reg <- lm(conf ~ rt2*accu, simexp)
        summary(m.reg)
        plot_model(m.reg, type = ("pred"),
                   terms = c("rt2", 'accu'),
                   title='')+
          theme_pubclean()
        
        emm <- emtrends(m.reg, var = "rt2" ,pairwise~accu)
        emm <- summary(emm$contrasts)
        conf.accu[nm,ca0,cd0,nb,] <- as.numeric(emm[2]) # save effect size
        pval.accu[nm,ca0,cd0,nb,] <- as.numeric(emm[6]) # save p avlue
        
        
        m.reg <- lm(conf ~ phq*rt2, simexp)
        emm <- emtrends(m.reg, var = "rt2" ,pairwise~phq, at=list(phq=c(11.6,.6)))
        emm <- test(emm)
        slope.dif[nm,ca0,cd0,nb,] <- emm$contrasts[,2]
        pval.inter[nm,ca0,cd0,nb,] <- emm$contrasts[,6]
        
        m.reg <- lm(conf ~ phq*rt2*accu, simexp)
        summary(m.reg)
        pm <- plot_model(m.reg, type = ("pred"), terms = c("rt2", "phq", 'accu'), title='',
                   colors = c('green4', 'darkgoldenrod3', 'firebrick2'),
                   line.size = 1.5, ci.lvl = .95, auto.label = T)+
          theme_pubclean()+
          xlab('Time')+
          ylab('Confidence')+
          labs(colour = 'AD')+
          theme(text = element_text(size=font_size))
        pm$data$facet <-  ifelse(pm$data$facet == "accu = 0", "Incorrect", "Correct")
        pm
        
        simexp.subj <- filter(simexp, accu==1) %>%
          mutate(srt2 = ntile(rt2,6)) %>%
          mutate(ADt = ntile(phq,2)) %>%
          mutate(ADt = recode_factor(ADt, '1' = 'Low', '2' = 'High')) %>%
          group_by(srt2, subj, ADt) %>%
          dplyr::summarise(
            conf = mean(conf))
        ggplot(simexp.subj, aes(x = srt2, y = conf, colour = ADt)) +
          scale_color_manual(breaks = c('Low', 'High'), values=c('green4', 'firebrick2')) +
          stat_summary(fun.data = mean_se, size = .6, position=position_dodge(width = .2)) + 
          # geom_smooth(method='lm')+
          # coord_cartesian(ylim=c(-.1,.25)) +
          theme_classic2()+
          labs(colour = 'AD score')+
          xlab('Time after decision (z-scored & tiled)') + ylab('Confidence (z-scored)') +
          theme(text = element_text(size=font_size-2),
                # legend.position = c(.65,.1),
                legend.position = 'top',
                legend.direction = 'horizontal')
        
        
        emm <- emtrends(m.reg, var = "rt2" ,pairwise~phq|accu, at=list(phq=c(11.6,.6)))
        emm <- test(emm)
        slope.dif.accu[nm,ca0,cd0,nb,]    <- emm$contrasts[,3]
        slope.dif.se.accu[nm,ca0,cd0,nb,] <- emm$contrasts[,4]
        pval.inter.accu[nm,ca0,cd0,nb,]   <- emm$contrasts[,7]

        cc <- simexp %>% group_by(subj) %>%
          dplyr::summarise(sconf = sd(conf),
            conf = mean(conf),
            accu = mean(accu),
            rt1 = mean(rt),
            rt2 = mean(rt2),
            phq = mean(phq))
        
        conf.sim[nm,ca0,cd0,nb,] <- cc$conf
        conf.sd.sim[nm,ca0,cd0,nb,] <- cc$sconf
        accu.sim[nm,ca0,cd0,nb,] <- cc$accu
        rt1.sim[nm,ca0,cd0,nb,] <- cc$rt1
        rt2.sim[nm,ca0,cd0,nb,] <- cc$rt2
        phq.sim[nm,ca0,cd0,nb,] <- cc$phq
      }
    }
  }
}

meanconf <- underconf
meanconf[,,,,1] <- apply(conf.sim, c(1:2), mean)
meanconf[,,,,2] <- meanconf[,,,,1]

df <- data.frame(c(underconf), c(slope.dif), c(pval.main), c(pval.inter), c(meanconf), 
                 c(conf.accu), c(pval.accu), c(slope.dif.accu), c(underconf.se),
                 c(slope.dif.se.accu), c(pval.inter.accu),
                 c(model.sim), c(ca0.sim), c(cd0.sim), c(beta.sim), c(corr.sim))
colnames(df) <- c('underconf', 'slope.dif', 'pval.main', 'pval.inter', 'conf', 
                  'conf.accu', 'pval.accu', 'slope.dif.accu', 'underconf.se',
                  'slope.dif.se.accu', 'pval.inter.accu',
                  'model.num', 'add0.sim', 'div0.sim', 'beta.num', 'accu')


setwd("~/OneDrive - University of Copenhagen/Projects/Experiments/confidenceTime/data/sim")

if (disc_conf){
  write.csv(df, 'simdata_popuinterc_disconf6.csv')
} else{
  write.csv(df, 'simdata_popuinterc.csv')
}


##############
### read in the model if it has already been run and saved
conf.range <- c(.3, .95)

disc_conf <- 0

conf_add_0    <- seq(-2,1,2)
conf_div_0    <- c(.5, 1, 2, 4)
v_beta        <- -c(.05, .1, .15, .2)
addbias_beta  <- c(.05, .1, .15, .2)
multbias_beta <- -c(.05, .1, .15, .2)

pcutoff <- .05
setwd("~/OneDrive - University of Copenhagen/Projects/Experiments/confidenceTime/data/sim/")

# accu.colors =  c('lightskyblue2','lightskyblue4')
accu.colors =  c('deepskyblue2','blue4')
ds.colors =  c('mistyrose','deeppink2')


if (disc_conf){
  df <- read.csv('simdata_popuinterc_disconf6.csv', header = TRUE, sep = ',')
} else{
  df <- read.csv('simdata_popuinterc.csv', header = TRUE, sep = ',')
}

df <- df %>% 
  mutate_at(c('model.num', 'beta.num', 'accu', 'add0.sim', 'div0.sim'),
            factor)

df <- df %>% 
  mutate(slopedif = slope.dif.accu) %>%
  mutate(slopedif.se = slope.dif.se.accu) %>%
  mutate(underc = underconf) %>%
  mutate_at('underc', ~replace(., pval.main > pcutoff, NA)) %>%
  mutate(accu = as.factor(accu)) %>%
  mutate(accu = recode_factor(accu, "1" = 'Correct', '0' = 'Incorrect'))

conf_add_0 <- unique(df$add0.sim)
conf_div_0 <- unique(df$div0.sim)

add0.labs <- paste('A0 = ', conf_add_0)
names(add0.labs) <- paste(conf_add_0)
div0.labs <- paste('M0 = ', conf_div_0)
names(div0.labs) <- paste(conf_div_0)


addinterc.toplot <- c(1:5)
multinterc.toplot <- c(1:4)
beta.toplot <- c(1:5)

ggplot(df %>% filter(beta.num %in% beta.toplot, add0.sim %in% conf_add_0[addinterc.toplot],
                     div0.sim %in% conf_div_0[multinterc.toplot], 
                     beta.num %in% beta.toplot,
                     model.num==1),
       aes(x = beta.num, y = underconf)) + 
  scale_x_discrete(labels = v_beta[beta.toplot]) +
  facet_grid(add0.sim ~ div0.sim,
             labeller = labeller(add0.sim = add0.labs,
                                 div0.sim = div0.labs)) +
  geom_bar(stat = 'identity', position = position_dodge(width=.75),
           width = .75, fill = 'khaki1', colour = 'black') +
  geom_errorbar(aes(ymin = underconf - underconf.se,
                    ymax = underconf + underconf.se),
                position = position_dodge(width=.75), width=.3)+
  # geom_point(shape=21, fill = 'white')+
  labs(colour = 'Accuracy') +
  # coord_cartesian(ylim = c(-.25,.5)) +
  geom_hline(yintercept = 0, color = 'grey40') +
  ylab('Regression slope of AD upon confidence') +
  xlab('Simulated regression of AD upon V-bias') +
  theme_classic()  +
  theme(legend.position = 'top',
        text = element_text(size=font_size))

ggplot(df %>% filter(beta.num %in% beta.toplot, add0.sim %in% conf_add_0[addinterc.toplot],
                     div0.sim %in% conf_div_0[multinterc.toplot], 
                     beta.num %in% beta.toplot,
                     model.num==2),
       aes(x = beta.num, y = underconf)) + 
  scale_x_discrete(labels = addbias_beta[beta.toplot]) +
  facet_grid(add0.sim ~ div0.sim,
             labeller = labeller(add0.sim = add0.labs,
                                 div0.sim = div0.labs)) +
  geom_bar(stat = 'identity', position = position_dodge(width=.75),
           width = .75, fill = 'khaki1', colour = 'black') +
  geom_errorbar(aes(ymin = underconf - underconf.se,
                    ymax = underconf + underconf.se),
                position = position_dodge(width=.75), width=.3)+
  # geom_point(shape=21, fill = 'white')+
  labs(colour = 'Accuracy') +
  # coord_cartesian(ylim = c(-.25,.5)) +
  geom_hline(yintercept = 0, color = 'grey40') +
  ylab('Regression slope of AD upon confidence') +
  xlab('Simulated regression of AD upon A-bias') +
  theme_classic()  +
  theme(legend.position = 'top',
        text = element_text(size=font_size))

ggplot(df %>% filter(beta.num %in% beta.toplot, add0.sim %in% conf_add_0[addinterc.toplot],
                     div0.sim %in% conf_div_0[multinterc.toplot], 
                     beta.num %in% beta.toplot,
                     model.num==3),
       aes(x = beta.num, y = underconf)) + 
  scale_x_discrete(labels = multbias_beta[beta.toplot]) +
  facet_grid(add0.sim ~ div0.sim,
             labeller = labeller(add0.sim = add0.labs,
                                 div0.sim = div0.labs)) +
  geom_bar(stat = 'identity', position = position_dodge(width=.75),
           width = .75, fill = 'khaki1', colour = 'black') +
  geom_errorbar(aes(ymin = underconf - underconf.se,
                    ymax = underconf + underconf.se),
                position = position_dodge(width=.75), width=.3)+
  labs(colour = 'Accuracy') +
  geom_hline(yintercept = 0, color = 'grey40') +
  ylab('Regression slope of AD upon confidence') +
  xlab('Simulated regression of AD upon M-bias') +
  theme_classic()  +
  theme(legend.position = 'top',
        text = element_text(size=font_size))


ggplot(df %>% filter(beta.num %in% beta.toplot, add0.sim %in% conf_add_0[addinterc.toplot],
                     div0.sim %in% conf_div_0[multinterc.toplot], 
                     beta.num %in% beta.toplot,
                     model.num==1),
       aes(x = beta.num, y = slopedif, colour = accu)) + 
  scale_x_discrete(labels = v_beta[beta.toplot]) +
  scale_colour_manual(values = accu.colors) +
  facet_grid(add0.sim ~ div0.sim,
             labeller = labeller(add0.sim = add0.labs,
                                 div0.sim = div0.labs)) +
  geom_bar(stat = 'identity', position = position_dodge(width=.75),
           width = .75, fill = ds.colors[1]) +
  geom_errorbar(aes(ymin = slopedif - slopedif.se,
                    ymax = slopedif + slopedif.se), 
                position = position_dodge(width=.75), width=.5)+
  geom_point(position = position_dodge(width = .77), shape=21, fill = 'white')+
  labs(colour = 'Accuracy') +
  # coord_cartesian(ylim = c(-.55,.2)) +
  geom_hline(yintercept = 0, color = 'grey40') +
  ylab('Slope difference (high - low)') +
  xlab('Simulated regression of AD upon V-bias') +
  theme_classic()  +
  theme(legend.position = 'none',
        text = element_text(size=font_size))


ggplot(df %>% filter(beta.num %in% beta.toplot, add0.sim %in% conf_add_0[addinterc.toplot],
                     div0.sim %in% conf_div_0[multinterc.toplot], 
                     beta.num %in% beta.toplot,
                     model.num==2),
       aes(x = beta.num, y = slopedif, colour = accu)) + 
  scale_x_discrete(labels = addbias_beta[beta.toplot]) +
  scale_colour_manual(values = accu.colors) +
  facet_grid(add0.sim ~ div0.sim,
             labeller = labeller(add0.sim = add0.labs,
                                 div0.sim = div0.labs)) +
  geom_bar(stat = 'identity', position = position_dodge(width=.75),
           width = .75, fill = ds.colors[1]) +
  geom_errorbar(aes(ymin = slopedif - slopedif.se,
                    ymax = slopedif + slopedif.se), 
                position = position_dodge(width=.75), width=.5)+
  geom_point(position = position_dodge(width = .77), shape=21, fill = 'white')+
  labs(colour = 'Accuracy') +
  # coord_cartesian(ylim = c(-.21,.42)) +
  geom_hline(yintercept = 0, color = 'grey40') +
  ylab('Slope difference (high - low)') +
  xlab('Simulated regression of AD upon A-bias') +
  theme_classic()  +
  theme(legend.position = 'none',
        text = element_text(size=font_size))

ggplot(df %>% filter(beta.num %in% beta.toplot, add0.sim %in% conf_add_0[addinterc.toplot],
                     div0.sim %in% conf_div_0[multinterc.toplot], 
                     beta.num %in% beta.toplot,
                     model.num==3),
       aes(x = beta.num, y = slopedif, colour = accu)) + 
  scale_x_discrete(labels = multbias_beta[beta.toplot]) +
  scale_colour_manual(values = accu.colors) +
  facet_grid(add0.sim ~ div0.sim,
             labeller = labeller(add0.sim = add0.labs,
                                 div0.sim = div0.labs)) +
  geom_bar(stat = 'identity', position = position_dodge(width=.75),
           width = .75, fill = ds.colors[1]) +
  geom_errorbar(aes(ymin = slopedif - slopedif.se,
                    ymax = slopedif + slopedif.se), 
                position = position_dodge(width=.75), width=.5)+
  geom_point(position = position_dodge(width = .77), shape=21, fill = 'white')+
  labs(colour = 'Accuracy') +
  # coord_cartesian(ylim = c(-.25,.55)) +
  geom_hline(yintercept = 0, color = 'grey40') +
  ylab('Slope difference (high - low)') +
  xlab('Simulated regression of AD upon M-bias') +
  theme_classic()  +
  theme(legend.position = 'none',
        text = element_text(size=font_size))


if (disc_conf){
  df <- read.csv('simdata_popuinterc_disconf6.csv', header = TRUE, sep = ',')
} else{
  df <- read.csv('simdata_popuinterc.csv', header = TRUE, sep = ',')
}

df <- df %>% 
  mutate_at(c('model.num', 'beta.num', 'accu', 'add0.sim', 'div0.sim'),
            factor)

df <- df %>% 
  mutate(slopedif = slope.dif.accu) %>%
  mutate(slopedif.se = slope.dif.se.accu) %>%
  mutate_at('slopedif', ~replace(., underconf>0 | pval.main > pcutoff
                                 , NA)) %>%
  mutate_at('slopedif.se', ~replace(., underconf>0 | pval.main > pcutoff
                                    , NA)) %>%
  mutate(underc = underconf) %>%
  mutate_at('underc', ~replace(., pval.main > pcutoff, NA)) %>%
  mutate(accu = as.factor(accu)) %>%
  mutate(accu = recode_factor(accu, "1" = 'Correct', '0' = 'Incorrect'))


addinterc.toplot <- c(2)
multinterc.toplot <- c(4)
beta.toplot <- c(1:5)
ggplot(df %>% filter(beta.num %in% beta.toplot, add0.sim %in% conf_add_0[addinterc.toplot],
                     div0.sim %in% conf_div_0[multinterc.toplot], 
                     beta.num %in% beta.toplot,
                     model.num==1),
       aes(x = beta.num, y = slopedif, colour = accu)) + 
  scale_x_discrete(labels = v_beta[beta.toplot]) +
  scale_colour_manual(values = accu.colors) +
  geom_bar(stat = 'identity', position = position_dodge(width=.75),
           width = .75, fill = ds.colors[1]) +
  geom_errorbar(aes(ymin = slopedif - slopedif.se,
                    ymax = slopedif + slopedif.se), 
                position = position_dodge(width=.75), width=.5)+
  geom_point(position = position_dodge(width = .77), shape=21, fill = 'white')+
  labs(colour = 'Accuracy') +
  coord_cartesian(ylim = c(-1.2,.2)) +
  geom_hline(yintercept = 0, color = 'grey40') +
  ylab('Slope difference\n(higher - lower)') +
  xlab('Simulated Beta upon V-bias') +
  theme_classic()  +
  theme(legend.position = 'off',
        # legend.box.just = 'left',
        text = element_text(size=font_size-2))
  


ggplot(df %>% filter(beta.num %in% beta.toplot, add0.sim %in% conf_add_0[addinterc.toplot],
                     div0.sim %in% conf_div_0[multinterc.toplot], 
                     beta.num %in% beta.toplot,
                     model.num==2),
       aes(x = beta.num, y = slopedif, colour = accu)) + 
  scale_x_discrete(labels = addbias_beta[beta.toplot]) +
  scale_colour_manual(values = accu.colors) +
  geom_bar(stat = 'identity', position = position_dodge(width=.75),
           width = .75, fill = ds.colors[1]) +
  geom_errorbar(aes(ymin = slopedif - slopedif.se,
                    ymax = slopedif + slopedif.se), 
                position = position_dodge(width=.75), width=.5)+
  geom_point(position = position_dodge(width = .77), shape=21, fill = 'white')+
  labs(colour = 'Accuracy') +
  # coord_cartesian(ylim = c(-.21,.42)) +
  geom_hline(yintercept = 0, color = 'grey40') +
  ylab('Slope difference\n(higher - lower)') +
  xlab('Simulated Beta upon A-bias') +
  theme_classic()  +
  theme(legend.position = 'off',
        text = element_text(size=font_size-2))

ggplot(df %>% filter(beta.num %in% beta.toplot, add0.sim %in% conf_add_0[addinterc.toplot],
                     div0.sim %in% conf_div_0[multinterc.toplot], 
                     beta.num %in% beta.toplot,
                     model.num==3),
       aes(x = beta.num, y = slopedif, colour = accu)) + 
  scale_x_discrete(labels = multbias_beta[beta.toplot]) +
  scale_colour_manual(values = accu.colors) +
  # facet_grid(add0.sim ~ div0.sim,
  #            labeller = labeller(add0.sim = add0.labs,
  #                                div0.sim = div0.labs)) +
  geom_bar(stat = 'identity', position = position_dodge(width=.75),
           width = .75, fill = ds.colors[1]) +
  geom_errorbar(aes(ymin = slopedif - slopedif.se,
                    ymax = slopedif + slopedif.se), 
                position = position_dodge(width=.75), width=.5)+
  geom_point(position = position_dodge(width = .77), shape=21, fill = 'white')+
  labs(colour = 'Accuracy') +
  coord_cartesian(ylim = c(-.05,.25)) +
  geom_hline(yintercept = 0, color = 'grey40') +
  ylab('Slope difference\n(higher - lower)') +
  xlab('Simulated Beta upon M-bias') +
  theme_classic()  +
  theme(legend.position = 'off',
        text = element_text(size=font_size-2))


#################
#### figure 3 for sample ddm trial
setwd("~/OneDrive - University of Copenhagen/Projects/Experiments/confidenceTime/data/sim")

a <- 1
z <- .5
v <- 1
sigma <- 1
dt <- .001
cov1 <- 0
cov2 <- 10
conf_add <- 0
conf_sd <- 1
beta_v <- -.2
beta_add <- 0
beta_sd <- 0

trial.params <- c(a, z, v, sigma, dt, cov1, cov2, 
                  conf_add, conf_sd, 
                  beta_v, beta_add, beta_sd)
names(trial.params) <- c('a', 'z', 'v', 'sigma', 'dt', 'cov1', 'cov2', 
                         'conf_add', 'conf_sd', 
                         'beta_v', 'beta_add', 'beta_sd')


# ######
# ## run this for simulating again (otherwise just load the saved trial)
# ddmTrial <- simDDMBiasTrialPlot(trial.params)
# save(ddmTrial, file = 'simDDMTrial.RData')
# #######
setwd("~/OneDrive - University of Copenhagen/Projects/Experiments/confidenceTime/data/sim/")

dd <- load('simDDMTrial.RData')

nSampChoice <- length(ddmTrial$evChoice)
nSampConf <- dim(ddmTrial$evConf)
choiceTime <- ddmTrial$choiceTime
confTime <- ddmTrial$confTime

x = seq(.5,1.5,.05)
xconf <- sigmoid(c(.1,12,1),x)

df <- data.frame(c(ddmTrial$evChoice, ddmTrial$evConf[1,], ddmTrial$evConf[2,],
                   x), 
                 c(c(1:nSampChoice, nSampChoice:(nSampChoice+nSampConf[2]-1), 
                     nSampChoice:(nSampChoice+nSampConf[2]-1))*dt,
                   c(xconf+confTime+.001)),
                 c(array('Predec', nSampChoice), array('Low AD', nSampConf[2]), 
                   array('High AD', nSampConf[2]), array('Sigm', length(x)))
                 # c(array('Low', nSampChoice+nSampConf[2]), array('High', nSampChoice+nSampConf[2]))
)
colnames(df) <- c('Evidence', 'Time', 'AD')
annotate_font_size <- 4

ggplot(df, aes(x=Time, y=Evidence, colour=AD)) + 
  scale_colour_manual(breaks = c('Predec', 'Low AD', 'High AD', 'Sigm'),
                      values = c('black', 'green4', 'firebrick2', 'black')) + 
  scale_linetype_manual(breaks = c('Predec', 'Low AD', 'High AD', 'Sigm'),
                        values=c("solid", "solid", "solid", "solid"))+
  geom_segment(x = choiceTime, y = 0, xend = choiceTime, yend = 1.12, color = 'gray') +
  geom_segment(x = 0, y = 0, xend = choiceTime, yend = 0, color = 'gray') +
  geom_segment(x = 0, y = a, xend = choiceTime, yend = a, color = 'gray') +
  geom_segment(x = confTime, y = 0.5, xend = confTime, yend = 1.5, color = 'gray') +
  geom_segment(x = confTime+.102, y = 0.5, xend = confTime+.102, yend = 1.5, color = 'gray') +  
  geom_line(size=.8, aes(linetype=AD)) +
  geom_segment(x = confTime, y = tail(ddmTrial$evConf[1,],1), 
               xend = confTime+.089, yend = tail(ddmTrial$evConf[1,],1), 
               color = 'green4', size=.8) +  
  geom_segment(x = confTime+.089, y = tail(ddmTrial$evConf[1,],1), 
               xend = confTime+.089, yend = 1.5, color = 'green4', size=1) +  
  geom_segment(x = confTime, y = tail(ddmTrial$evConf[2,],1), 
               xend = confTime+.0265, yend = tail(ddmTrial$evConf[2,],1), 
               color = 'firebrick2', size=1) +  
  geom_segment(x = confTime+.0265, y = tail(ddmTrial$evConf[2,],1), 
               xend = confTime+.0265, yend = 1.5, color = 'firebrick2', size=.8) +  
  annotate('text', x=.327, y=1.75, label="Confidence", size=annotate_font_size) +
  annotate('text', x=.276, y=1.59, label="0", size=annotate_font_size) +
  annotate('text', x=.377, y=1.59, label="1", size=annotate_font_size) +
  annotate('text', x=-.01, y=.5, label="z", size=annotate_font_size) +
  annotate('text', x=-.01, y=1, label="a", size=annotate_font_size) +
  annotate('text', x=-.01, y=0, label="0", size=annotate_font_size) +
  annotate('text', x=choiceTime-.01, y=1.25, label="Choice", size=annotate_font_size) +
  coord_cartesian(ylim = c(0,1.75))+
  theme_pubclean()+ 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.title=element_blank(),
        legend.text = element_text(size = font_size),
        legend.position = 'none',
        text = element_text(size=font_size))


x = seq(.5,1.5,.05)
xconfL <- sigmoid(c(.1,20,1.2),x)

df <- data.frame(c(ddmTrial$evChoice, ddmTrial$evConf[1,],
                   x, x), 
                 c(c(1:nSampChoice, nSampChoice:(nSampChoice+nSampConf[2]-1))*dt,
                   c(xconf+confTime+.001), c(xconfL+confTime+.001)),
                 c(array('Predec', nSampChoice), array('Postdec', nSampConf[2]), 
                   array('Low AD', length(x)), array('High AD', length(x)))
                 # c(array('Low', nSampChoice+nSampConf[2]), array('High', nSampChoice+nSampConf[2]))
)
colnames(df) <- c('Evidence', 'Time', 'AD')

ggplot(df, aes(x=Time, y=Evidence, colour=AD)) + 
  scale_colour_manual(breaks = c('Predec', 'Postdec', 'Low AD', 'High AD'),
                      values = c('black', 'black', 'green4', 'firebrick2')) + 
  geom_segment(x = choiceTime, y = 0, xend = choiceTime, yend = 1.12, color = 'gray') +
  geom_segment(x = 0, y = 0, xend = choiceTime, yend = 0, color = 'gray') +
  geom_segment(x = 0, y = a, xend = choiceTime, yend = a, color = 'gray') +
  geom_segment(x = confTime, y = 0.5, xend = confTime, yend = 1.5, color = 'gray') +
  geom_segment(x = confTime+.101, y = 0.5, xend = confTime+.101, yend = 1.5, color = 'gray') +
  geom_line(size=.8, aes(linetype=AD)) +
  scale_linetype_manual(breaks = c('Predec', 'Postdec', 'Low AD', 'High AD'),
                        values=c("solid", "solid", "solid", "solid"))+
  geom_segment(x = confTime, y = tail(ddmTrial$evConf[1,],1), 
               xend = confTime+.089, yend = tail(ddmTrial$evConf[1,],1), 
               color = 'black', size=.8) +  
  geom_segment(x = confTime+.089, y = tail(ddmTrial$evConf[1,],1), 
               xend = confTime+.089, yend = 1.5, color = 'green4', size=1) +  
  geom_segment(x = confTime, y = tail(ddmTrial$evConf[1,],1), 
               xend = confTime+.038, yend = tail(ddmTrial$evConf[1,],1), 
               color = 'black', size=.8) +  
  geom_segment(x = confTime+.038, y = tail(ddmTrial$evConf[1,],1), 
               xend = confTime+.038, yend = 1.5, color = 'firebrick2', size=1) +  
  annotate('text', x=.327, y=1.75, label="Confidence", size=annotate_font_size) +
  annotate('text', x=.276, y=1.59, label="0", size=annotate_font_size) +
  annotate('text', x=.377, y=1.59, label="1", size=annotate_font_size) +
  annotate('text', x=-.01, y=.5, label="z", size=annotate_font_size) +
  annotate('text', x=-.01, y=1, label="a", size=annotate_font_size) +
  annotate('text', x=-.01, y=0, label="0", size=annotate_font_size) +
  annotate('text', x=choiceTime-.01, y=1.25, label="Choice", size=annotate_font_size) +
  coord_cartesian(ylim = c(0,1.75))+
  theme_pubclean()+ 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.title=element_blank(),
        legend.text = element_text(size = font_size),
        legend.position = 'none',
        text = element_text(size=font_size))


x = seq(.5,1.5,.05)
xconf <- sigmoid(c(.1,14,1),x)
xconfL <- sigmoid(c(.1,7,1),x)

df <- data.frame(c(ddmTrial$evChoice, ddmTrial$evConf[1,],
                   x, x), 
                 c(c(1:nSampChoice, nSampChoice:(nSampChoice+nSampConf[2]-1))*dt,
                   c(xconf+confTime+.001), c(xconfL+confTime+.001)),
                 c(array('Predec', nSampChoice), array('Postdec', nSampConf[2]), 
                   array('Low AD', length(x)), array('High AD', length(x)))
                 # c(array('Low', nSampChoice+nSampConf[2]), array('High', nSampChoice+nSampConf[2]))
)
colnames(df) <- c('Evidence', 'Time', 'AD')

ggplot(df, aes(x=Time, y=Evidence, colour=AD)) + 
  scale_colour_manual(breaks = c('Predec', 'Postdec', 'Low AD', 'High AD'),
                      values = c('black', 'black', 'green4', 'firebrick2')) + 
  geom_segment(x = choiceTime, y = 0, xend = choiceTime, yend = 1.12, color = 'gray') +
  geom_segment(x = 0, y = 0, xend = choiceTime, yend = 0, color = 'gray') +
  geom_segment(x = 0, y = a, xend = choiceTime, yend = a, color = 'gray') +
  geom_segment(x = confTime, y = 0.5, xend = confTime, yend = 1.5, color = 'gray') +
  geom_segment(x = confTime+.101, y = 0.5, xend = confTime+.101, yend = 1.5, color = 'gray') +
  geom_line(size=.8, aes(linetype=AD)) +
  scale_linetype_manual(breaks = c('Predec', 'Postdec', 'Low AD', 'High AD'),
                        values=c("solid", "solid", "solid", "solid"))+
  geom_segment(x = confTime, y = tail(ddmTrial$evConf[1,],1), 
               xend = confTime+.0925, yend = tail(ddmTrial$evConf[1,],1), 
               color = 'black', size=.8) +  
  geom_segment(x = confTime+.0925, y = tail(ddmTrial$evConf[1,],1), 
               xend = confTime+.0925, yend = 1.5, color = 'green4', size=.8) +  
  geom_segment(x = confTime, y = tail(ddmTrial$evConf[1,],1), 
               xend = confTime+.078, yend = tail(ddmTrial$evConf[1,],1), 
               color = 'black', size=.8) +  
  geom_segment(x = confTime+.078, y = tail(ddmTrial$evConf[1,],1), 
               xend = confTime+.078, yend = 1.5, color = 'firebrick2', size=.8) +    
  annotate('text', x=.327, y=1.75, label="Confidence", size=annotate_font_size) +
  annotate('text', x=.276, y=1.59, label="0", size=annotate_font_size) +
  annotate('text', x=.377, y=1.59, label="1", size=annotate_font_size) +
  annotate('text', x=-.01, y=.5, label="z", size=annotate_font_size) +
  annotate('text', x=-.01, y=1, label="a", size=annotate_font_size) +
  annotate('text', x=-.01, y=0, label="0", size=annotate_font_size) +
  annotate('text', x=choiceTime-.01, y=1.2, label="Choice", size=annotate_font_size) +
  coord_cartesian(ylim = c(0,1.74))+
  theme_pubclean()+ 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.title=element_blank(),
        legend.text = element_text(size = font_size),
        legend.position = 'none',
        text = element_text(size=font_size))


