
################################## Katyal & Fleming ################################
### Gender and anxiety reveal distinct computational sources of underconfidence ####
####################################################################################
################  This file reproduces Figure 2 ##################
###############################################################################
####### Analyse and plot how confidence changes over time with   ##############
####### different individual differences                         ##############
###############################################################################


library(lmerTest) # anova and summary with p-values
library(emmeans) # Least squares means
library(ggpubr)
library(tidyverse)
library(ggplot2)
library(ggbeeswarm) # for geom_quasirandom()
library(ggpp)
library(sjPlot)
library(see)
library(R.matlab)


###
emm_options(lmer.df = "satterthwaite")
emm_options(lmerTest.limit = 35000)
font_size <- 16
tag_size <- 16
alpha <- .4
max_RT_deviation = 3
excl_insuf_var <- T
n.rt.tiles <- 6

taskNames = c('Perception', 'Memory')
posnegNames = c('Positive', 'Negative')
posnegColours = c("chartreuse3", "coral3")

min.accept.rt <- 100

n_exp <- 4+1  

AD.inter     <- array(NA, c(n_exp, 2, 5))
CIT.inter    <- array(NA, c(n_exp, 2, 5))
age.inter    <- array(NA, c(n_exp, 2, 5))
gender.inter <- array(NA, c(n_exp, 2, 5))

AD.main     <- array(NA, c(n_exp, 5))
CIT.main    <- array(NA, c(n_exp, 5))
age.main    <- array(NA, c(n_exp, 5))
gender.main <- array(NA, c(n_exp, 5))

#############################################################################
####### First anlayse for each of the four datasets separately ##############
#############################################################################

####
#### Katyal et al - Exp 1
##

setwd("~/OneDrive - University of Copenhagen/Projects/Experiments/confidenceTime/")

## read data file
katyal_1.df = read.csv('data/Katyal_et_al_2023/mbsExp1_noperfexcl.csv', header = TRUE, sep = ',')
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
                                 "2" = posnegNames[2], "0" = "None"),
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

katyal_1.df$scale_rt2 <- katyal_1.df$rt2
katyal_1.df$scale_rt2[katyal_1.df$task==taskNames[1]] <- scale(katyal_1.df$rt2[katyal_1.df$task==taskNames[1]])
katyal_1.df$scale_rt2[katyal_1.df$task==taskNames[2]] <- scale(katyal_1.df$rt2[katyal_1.df$task==taskNames[2]])

katyal_1.df$stimev[katyal_1.df$task==taskNames[1]] <- scale(katyal_1.df$stair[katyal_1.df$task==taskNames[1]])
katyal_1.df$stimev[katyal_1.df$task==taskNames[2]] <- scale(katyal_1.df$stair[katyal_1.df$task==taskNames[2]])

katyal_1.df$conf.sc[katyal_1.df$task==taskNames[1]] <- scale(katyal_1.df$conf[katyal_1.df$task==taskNames[1]])
katyal_1.df$conf.sc[katyal_1.df$task==taskNames[2]] <- scale(katyal_1.df$conf[katyal_1.df$task==taskNames[2]])


# give both interv and transf run the same value of fbblock (feedback block type)
katyal_1.df$fbblock[katyal_1.df$runnum==4] = katyal_1.df$fbblock[katyal_1.df$runnum==3]
katyal_1.df$fbblock[katyal_1.df$runnum==6] = katyal_1.df$fbblock[katyal_1.df$runnum==5]

katyal_1.df <- katyal_1.df %>%mutate_at(c("group", 'gender'), as.factor) %>%
  mutate(gender = recode_factor(gender, '1' = 'Female', '2' = 'Male'))

katyal_1.df.subj <- katyal_1.df %>% group_by(subj,task,gender) %>%
  dplyr::summarise(conf.sc = mean(conf.sc),
            accu = mean(accu),
            AD = mean(AD),
            age = mean(age),
            stimev = mean(stimev),
            scale_rt1 = mean(scale_rt1, na.rm = T),
            scale_rt2 = mean(scale_rt2, na.rm = T),)

expNum <- 1

m.reg <- lmer(conf.sc ~ AD + gender + age + stimev + accu + fbblock +
                (1 + accu+stimev|trialnum) + (1|task),
              control = lmerControl(optimizer = c('bobyqa'), calc.derivs = F),
              filter(katyal_1.df))
summary(m.reg)

sm <- summary(m.reg)
AD.main[expNum,c(1:3)] <- as.numeric(sm$coefficients[2,c(1,2,5)])
gender.main[expNum,c(1:3)] <- as.numeric(sm$coefficients[3,c(1,2,5)])
age.main[expNum,c(1:3)] <- as.numeric(sm$coefficients[4,c(1,2,5)])

## conf interv take really long with random slopes so remove them
m.reg <- lmer(conf.sc ~ AD + gender + age + stimev + accu + fbblock + (1|task),
              filter(katyal_1.df))
summary(m.reg)

ci <- confint(m.reg)
AD.main[expNum,c(4:5)] <- as.numeric(ci[4,])
gender.main[expNum,c(4:5)] <- as.numeric(ci[5,])
age.main[expNum,c(4:5)] <- as.numeric(ci[6,])


for (accuracy in c(0,1)){
  k1 <- filter(katyal_1.df, accu==accuracy) %>% 
    dplyr::select(c(conf.sc, AD, scale_rt1, scale_rt2, gender, age, stimev,
             runnum, trialinrun, task, accu, subj, trialnum)) %>% drop_na()
  
  m.reg <- lmer(conf.sc ~ AD*(scale_rt1 + scale_rt2) +
                  gender*(scale_rt1 + scale_rt2) +
                  age*(scale_rt1 + scale_rt2) +
                  stimev + (1|trialnum) + (1|task), k1)
  summary(m.reg)
  
  emm <- emtrends(m.reg, var = "scale_rt2" ,pairwise~AD, at=list(AD=c(1,-1)))
  test(emm)
  df <- summary(emm$contrasts)
  ci <- confint(emm)
  AD.inter[expNum,2-accuracy,] <- c(as.numeric(df[c(2,3,6)]), as.numeric(ci$contrasts[5:6]))
  
  emm <- emtrends(m.reg, var = "scale_rt2" ,pairwise~gender, at=list(gender=c('Female','Male')))
  test(emm)
  df <- summary(emm$contrasts)
  ci <- confint(emm)
  gender.inter[expNum,2-accuracy,] <- c(as.numeric(df[c(2,3,6)]), as.numeric(ci$contrasts[5:6]))
  
  
  emm <- emtrends(m.reg, var = "scale_rt2" ,pairwise~age, at=list(age=c(50,20)))
  test(emm)
  df <- summary(emm$contrasts)
  ci <- confint(emm)
  age.inter[expNum,2-accuracy,] <- c(as.numeric(df[c(2,3,6)]), as.numeric(ci$contrasts[5:6]))
  

}

###################
##########
#### Exp 2

setwd("~/OneDrive - University of Copenhagen/Projects/Experiments/confidenceTime/data/Katyal_et_al_2023//")

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

katyal_2.df$scale_rt2 <- katyal_2.df$rt2
katyal_2.df$scale_rt2[katyal_2.df$task==taskNames[1]] <- scale(katyal_2.df$rt2[katyal_2.df$task==taskNames[1]])
katyal_2.df$scale_rt2[katyal_2.df$task==taskNames[2]] <- scale(katyal_2.df$rt2[katyal_2.df$task==taskNames[2]])

katyal_2.df$stimev[katyal_2.df$task==taskNames[1]] <- scale(katyal_2.df$stair[katyal_2.df$task==taskNames[1]])
katyal_2.df$stimev[katyal_2.df$task==taskNames[2]] <- scale(katyal_2.df$stair[katyal_2.df$task==taskNames[2]])

katyal_2.df$conf.sc[katyal_2.df$task==taskNames[1]] <- scale(katyal_2.df$conf[katyal_2.df$task==taskNames[1]])
katyal_2.df$conf.sc[katyal_2.df$task==taskNames[2]] <- scale(katyal_2.df$conf[katyal_2.df$task==taskNames[2]])

katyal_2.df <- katyal_2.df %>%
  mutate(blockType = factor(runnum%%2)) %>%
  mutate(blockType = recode_factor(blockType, "0" = "transf", "1" = "interv"))
katyal_2.df <- katyal_2.df %>%
  mutate_at(vars(contains("runnum")), as.factor)
katyal_2.df$fbblock <- droplevels(katyal_2.df$fbblock)

# give both interv and transf run the same value of fbblock (feedback block type)
katyal_2.df$fbblock[katyal_2.df$runnum==4] = katyal_2.df$fbblock[katyal_2.df$runnum==3 &
                                                          katyal_2.df$trialinrun<=20]
katyal_2.df$fbblock[katyal_2.df$runnum==6] = katyal_2.df$fbblock[katyal_2.df$runnum==5
                                                        & katyal_2.df$trialinrun<=20]

katyal_2.df <- katyal_2.df %>%
  mutate_at(c("group", 'gender'), as.factor) %>%
  mutate(group = recode_factor(group, "1" = 'Group 1', 
                               "2" = 'Group 2', "3" = 'Group 3', "4" = 'Group 4', 
                               "5" = 'Group 5', "6" = 'Group 6', "7" = 'Group 7', 
                               "8" = 'Group 8')) %>%
  mutate(gender = recode_factor(gender, '1' = 'Female', '2' = 'Male'))

# #####################
# ## read in the psychiatric scores
psych = read.csv(paste('factor_scores_noperfexc.csv', sep=''),
                 header = TRUE, sep = ',')
psych <- psych %>% dplyr::rename(subj = subjIDs) %>% dplyr::select(-X) %>%
  mutate_at("subj", as.factor) %>%
  dplyr::rename(CIT = Compul)

katyal_2.df <- katyal_2.df %>% left_join(psych)

katyal_2.df <- katyal_2.df %>% drop_na(scale_rt1, scale_rt2)

katyal_2.df.subj <- katyal_2.df %>% group_by(subj,task) %>%
  dplyr::summarise(conf = mean(conf))

expNum <- 2

m.reg <- lmer(conf.sc ~ AD + CIT + gender + age + stimev + accu + fbblock + (1|task),
              control = lmerControl(optimizer = c('bobyqa'), calc.derivs = F),
              filter(katyal_2.df))
summary(m.reg)

sm <- summary(m.reg)

AD.main[expNum,c(1:3)] <- as.numeric(sm$coefficients[2,c(1,2,5)])
CIT.main[expNum,c(1:3)] <- as.numeric(sm$coefficients[3,c(1,2,5)])
gender.main[expNum,c(1:3)] <- as.numeric(sm$coefficients[4,c(1,2,5)])
age.main[expNum,c(1:3)] <- as.numeric(sm$coefficients[5,c(1,2,5)])

ci <- confint(m.reg)
AD.main[expNum,c(4:5)] <- as.numeric(ci[4,])
CIT.main[expNum,c(4:5)] <- as.numeric(ci[5,])
gender.main[expNum,c(4:5)] <- as.numeric(ci[6,])
age.main[expNum,c(4:5)] <- as.numeric(ci[7,])


for (accuracy in c(0:1)){
  
  k2 <- filter(katyal_2.df, accu==accuracy) %>% drop_na()
  
  m.reg <- lmer(conf.sc ~ AD*(scale_rt1 + scale_rt2) +
                  CIT*(scale_rt1 + scale_rt2) +
                  gender*(scale_rt1 + scale_rt2) +
                  age*(scale_rt1 + scale_rt2) +
                  stimev + (1|trialnum) + (1|task), k2)
  summary(m.reg)
  
  emm <- emtrends(m.reg, var = "scale_rt2" ,pairwise~AD, at=list(AD=c(1,-1)))
  test(emm)
  df <- summary(emm$contrasts)
  ci <- confint(emm)
  AD.inter[expNum,2-accuracy,] <- c(as.numeric(df[c(2,3,6)]), as.numeric(ci$contrasts[5:6]))
  
  emm <- emtrends(m.reg, var = "scale_rt2" ,pairwise~CIT, at=list(CIT=c(1,-1)))
  test(emm)
  df <- summary(emm$contrasts)
  ci <- confint(emm)
  CIT.inter[expNum,2-accuracy,] <- c(as.numeric(df[c(2,3,6)]), as.numeric(ci$contrasts[5:6]))
  
  emm <- emtrends(m.reg, var = "scale_rt2" ,pairwise~gender, at=list(gender=c('Female','Male')))
  test(emm)
  df <- summary(emm$contrasts)
  ci <- confint(emm)
  gender.inter[expNum,2-accuracy,] <- c(as.numeric(df[c(2,3,6)]), as.numeric(ci$contrasts[5:6]))
  
  emm <- emtrends(m.reg, var = "scale_rt2" ,pairwise~age, at=list(age=c(50,20)))
  test(emm)
  df <- summary(emm$contrasts)
  ci <- confint(emm)
  age.inter[expNum,2-accuracy,] <- c(as.numeric(df[c(2,3,6)]), as.numeric(ci$contrasts[5:6]))
  
}


########################
#################### rouaoult  exp 1

## read data file
setwd("~/OneDrive - University of Copenhagen/Projects/Experiments/confidenceTime/data/Rouault_et_al_BiolPsychiat_2018//")

rouault_1.df <- read.csv('rouault_2018_exp1.csv', header = TRUE, sep = ',')

# for discrete ratings: use confidence resolution to restrict subjects with insufficient variation in confidence
conf.resolution <- 1/(length(unique(rouault_1.df$conf))-1)

rouault_1.df <- rouault_1.df %>%
  dplyr::rename(subj = id) %>%
  dplyr::rename(accu = correct) %>%
  mutate(gender = recode_factor(gender, '0' = 'Male', '1' = 'Female')) %>%
  mutate(conf = (conf-1)/10) %>%
  mutate(AD = c(scale(anxiety+zung)))

mvsubj <- rouault_1.df %>% group_by(subj) %>%
  dplyr::summarise(sconf = sd(conf),
                   mconf = mean(conf))
insuf_var_subj <- mvsubj$subj[mvsubj$sconf<conf.resolution]

rouault_1.df$rt1[rouault_1.df$rt1 < min.accept.rt] <- NA
rouault_1.df <- rouault_1.df %>% drop_na(rt1, rt2)

max_RT_deviation = 3
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

if (excl_insuf_var){
  rouault_1.df <- filter(rouault_1.df, !(subj %in% insuf_var_subj))
}

length(unique(rouault_1.df$subj))

expNum <- 3

m.reg <- lmer(conf.sc ~ AD + gender + age + stimev + accu + 
                (accu + 0|trialnum), 
              # control = lmerControl(optimizer = c('Nelder_Mead'), calc.derivs = F),
              filter(rouault_1.df))
summary(m.reg)

sm <- summary(m.reg)

AD.main[expNum,c(1:3)] <- as.numeric(sm$coefficients[2,c(1,2,5)])
CIT.main[expNum,c(1:3)] <- NA
gender.main[expNum,c(1:3)] <- -as.numeric(sm$coefficients[3,c(1,2,5)])
age.main[expNum,c(1:3)] <- as.numeric(sm$coefficients[4,c(1,2,5)])

ci <- confint(m.reg)
AD.main[expNum,c(4:5)] <- as.numeric(ci[4,])
CIT.main[expNum,c(4:5)] <- NA
gender.main[expNum,c(4:5)] <- -as.numeric(ci[5,])
age.main[expNum,c(4:5)] <- as.numeric(ci[6,])

for (accuracy in c(0:1)){
  
  m.reg <- lmer(conf.sc ~ AD*(scale_rt1 + scale_rt2) +
                  gender*(scale_rt1 + scale_rt2) +
                  age*(scale_rt1 + scale_rt2) +
                  stimev + (1|trialnum),
                filter(rouault_1.df, accu==accuracy))
  summary(m.reg)
  
  emm <- emtrends(m.reg, var = "scale_rt2" ,pairwise~AD, at=list(AD=c(1,-1)))
  test(emm)
  df <- summary(emm$contrasts)
  ci <- confint(emm)
  AD.inter[expNum,2-accuracy,] <- c(as.numeric(df[c(2,3,6)]), as.numeric(ci$contrasts[5:6]))
  
  emm <- emtrends(m.reg, var = "scale_rt2" ,pairwise~gender, at=list(gender=c('Female','Male')))
  test(emm)
  df <- summary(emm$contrasts)
  ci <- confint(emm)
  gender.inter[expNum,2-accuracy,] <- c(as.numeric(df[c(2,3,6)]), as.numeric(ci$contrasts[5:6]))
  
  emm <- emtrends(m.reg, var = "scale_rt2" ,pairwise~age, at=list(age=c(50,20)))
  test(emm)
  df <- summary(emm$contrasts)
  ci <- confint(emm)
  age.inter[expNum,2-accuracy,] <- c(as.numeric(df[c(2,3,6)]), as.numeric(ci$contrasts[5:6]))
  
}


# 
# ########################
# ########################
# #################### rouaoult  exp 2

## read data file
setwd("~/OneDrive - University of Copenhagen/Projects/Experiments/confidenceTime/data/Rouault_et_al_BiolPsychiat_2018//")

rouault_2.df <- read.csv('rouault_2018_exp2.csv', header = TRUE, sep = ',') %>% drop_na()

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

max_RT_deviation = 3
rt1_max = median(rouault_2.df$rt1, na.rm=T)+ max_RT_deviation*mad(rouault_2.df$rt1, na.rm=T ) # max rt for perception task
rouault_2.df <- rouault_2.df %>%
  mutate_at(vars(contains("rt1")), ~replace(.,.>rt1_max[1],NA))
rt2_max = median(rouault_2.df$rt2, na.rm=T)+ max_RT_deviation*mad(rouault_2.df$rt2, na.rm=T ) # max rt for perception task
rouault_2.df <- rouault_2.df %>%
  mutate_at(vars(contains("rt2")), ~replace(.,.>rt2_max[1],NA))

rouault_2.df$scale_rt1 <- c(scale(rouault_2.df$rt1))
rouault_2.df$scale_rt2 <- c(scale(rouault_2.df$rt2))

rouault_2.df$stimev <- c(scale(rouault_2.df$stimev))
rouault_2.df$conf.sc <- c(scale(rouault_2.df$conf))

if (excl_insuf_var){
  rouault_2.df <- filter(rouault_2.df, !(subj %in% insuf_var_subj))
  length(unique(rouault_2.df))
}

expNum <- 4

m.reg <- lm(conf.sc ~ AD + CIT + gender + age + stimev + accu, 
              filter(rouault_2.df))
summary(m.reg)
sm <- summary(m.reg)

AD.main[expNum,c(1:3)] <- as.numeric(sm$coefficients[2,c(1,2,4)])
CIT.main[expNum,c(1:3)] <- as.numeric(sm$coefficients[3,c(1,2,4)])
gender.main[expNum,c(1:3)] <- as.numeric(summary(emmeans(m.reg, pairwise~gender)$contrasts)[c(2,3,6)])
age.main[expNum,c(1:3)] <- as.numeric(sm$coefficients[5,c(1,2,4)])

ci <- confint(m.reg)
AD.main[expNum,c(4:5)] <- as.numeric(ci[2,])
CIT.main[expNum,c(4:5)] <- as.numeric(ci[3,])
gender.main[expNum,c(4:5)] <- -as.numeric(ci[4,])
age.main[expNum,c(4:5)] <- as.numeric(ci[5,])


for (accuracy in c(0:1)){
  
  m.reg <- lmer(conf.sc ~ AD*(scale_rt1 + scale_rt2) +
                  CIT*(scale_rt1 + scale_rt2) +
                  gender*(scale_rt1 + scale_rt2) +
                  age*(scale_rt1 + scale_rt2) +
                  stimev + (1|trialnum),
                filter(rouault_2.df, accu==accuracy))
  summary(m.reg)
  
  emm <- emtrends(m.reg, var = "scale_rt2" ,pairwise~AD, at=list(AD=c(1,-1)))
  test(emm)
  df <- summary(emm$contrasts)
  ci <- confint(emm)
  AD.inter[expNum,2-accuracy,] <- c(as.numeric(df[c(2,3,6)]), as.numeric(ci$contrasts[5:6]))
  
  emm <- emtrends(m.reg, var = "scale_rt2" ,pairwise~CIT, at=list(CIT=c(1,-1)))
  test(emm)
  df <- summary(emm$contrasts)
  ci <- confint(emm)
  CIT.inter[expNum,2-accuracy,] <- c(as.numeric(df[c(2,3,6)]), as.numeric(ci$contrasts[5:6]))

  emm <- emtrends(m.reg, var = "scale_rt2" ,pairwise~gender, at=list(gender=c('Female','Male')))
  test(emm)
  df <- summary(emm$contrasts)
  ci <- confint(emm)
  gender.inter[expNum,2-accuracy,] <- c(as.numeric(df[c(2,3,6)]), as.numeric(ci$contrasts[5:6]))
  
  emm <- emtrends(m.reg, var = "scale_rt2" ,pairwise~age, at=list(age=c(50,20)))
  test(emm)
  df <- summary(emm$contrasts)
  ci <- confint(emm)
  age.inter[expNum,2-accuracy,] <- c(as.numeric(df[c(2,3,6)]), as.numeric(ci$contrasts[5:6]))
}


######################################################
####################################
################## Combined dataset
num.trials <- c()
katyal_1.df <- katyal_1.df %>%
  mutate(AD.sc = c(scale(gad))) %>%
  mutate(CIT.sc = NA) %>%
  mutate(SW.sc = NA) %>%
  mutate(expnum = as.factor(1)) %>%
  mutate(mid = as.factor(as.numeric(subj)+10000))

katyal_2.df <- katyal_2.df %>%
  mutate(AD.sc = c(scale(AD))) %>%
  mutate(CIT.sc = c(scale(CIT))) %>%
  mutate(SW.sc = c(scale(SW))) %>%
  mutate(expnum = as.factor(2)) %>%
  mutate(mid = as.factor(as.numeric(subj)+20000))

rouault_1.df <- rouault_1.df %>%
  mutate(AD.sc = c(scale(AD))) %>%
  mutate(CIT.sc = NA) %>%
  mutate(SW.sc = NA) %>%
  mutate(runnum = 1)%>%
  mutate(expnum = as.factor(3)) %>%
  mutate(task = 'Perception') %>%
  mutate(mid = as.factor(as.numeric(subj)+30000))

rouault_2.df <- rouault_2.df %>%
  mutate(AD.sc = c(scale(AD))) %>%
  mutate(CIT.sc = c(scale(CIT))) %>%
  mutate(SW.sc = c(scale(SW))) %>%
  mutate(runnum = 1)%>%
  mutate(expnum = as.factor(4)) %>%
  mutate(task = 'Perception') %>%
  mutate(mid = as.factor(as.numeric(subj)+40000))

expNum <- 5

cols2select <- c('conf.sc', 'accu', 'AD.sc', 'CIT.sc', 'SW.sc', 'age', 'gender', 
                 'trialnum', 'runnum', 'scale_rt1', 'scale_rt2', 'expnum', 
                 'stimev', 'task', 'mid')

meta.df <- rbind(katyal_1.df %>% dplyr::select(all_of(cols2select)),
                 katyal_2.df %>% dplyr::select(all_of(cols2select)),
                 rouault_1.df %>% dplyr::select(all_of(cols2select)),
                 rouault_2.df %>% dplyr::select(all_of(cols2select))) %>%
  dplyr::rename(conf = conf.sc)


### main effect of the basic confidence bias in combined dataset

meta.df <- meta.df %>%
  mutate(gender = factor(gender, levels = c('Male', 'Female'))) %>%
  mutate(accu = as.factor(accu))

m.reg <- lmer(conf ~ AD.sc + gender + age + stimev + accu + 
                (1|expnum), filter(meta.df %>% drop_na()))
sm <- summary(m.reg)
sm
AD.main[expNum,c(1:3)] <- as.numeric(sm$coefficients[2,c(1,2,5)])
gender.main[expNum,c(1:3)] <- as.numeric(summary(emmeans(
  m.reg, pairwise~gender)$contrasts, at=list(gender=c('Male', 'Female')))[c(2,3,6)])
age.main[expNum,c(1:3)] <- as.numeric(sm$coefficients[4,c(1,2,5)])

ci <- confint(m.reg)
AD.main[expNum,c(4:5)] <- as.numeric(ci[4,])
gender.main[expNum,c(4:5)] <- -as.numeric(ci[5,])
age.main[expNum,c(4:5)] <- as.numeric(ci[6,])


m2 <- update(m.reg, ~.-AD.sc)
anova(m.reg, m2)

m2 <- update(m.reg, ~.-gender)
anova(m.reg, m2)

m2 <- update(m.reg, ~.-age)
anova(m.reg, m2)

# calculate this separately for CIT because data is only available in Exp 2 and Exp 4
meta.df <- meta.df %>%
  mutate(gender = factor(gender, levels = c('Male', 'Female'))) %>%
  mutate(accu = as.factor(accu))

m.reg <- lmer(conf ~ AD.sc + CIT.sc + gender + age + stimev + accu + 
                (1|expnum), filter(meta.df %>% drop_na()))
sm <- summary(m.reg)
sm

CIT.main[expNum,c(1:3)] <- as.numeric(sm$coefficients[3,c(1,2,5)])

m2 <- update(m.reg, ~.-CIT.sc)
anova(m.reg, m2)

ci <- confint(m.reg)
CIT.main[expNum,c(4:5)] <- as.numeric(ci[5,])

n.rt.tiles <- 6
### confidence with postdecision time in combined dataset
for (accuracy in c(0,1)){
  
  meta <- filter(meta.df, accu==accuracy) %>% 
    dplyr::select(c(conf, AD.sc, scale_rt1, scale_rt2, gender, age, stimev,
             runnum, trialnum, expnum, accu, mid)) %>%
    mutate(gender = recode_factor(gender, 'Male' = 'Men', 'Female' = 'Women')) %>%
    drop_na()
  
  meta.subj <- meta %>% group_by(mid) %>%
    dplyr::summarise(AD.sc = mean(AD.sc),
                     age = mean(age)) %>%
    mutate(ADt = ntile(AD.sc,2)) %>%
    mutate(ADt = recode_factor(ADt, '1' = 'Lower', '2' = 'Higher')) %>%
    mutate(aget = ntile(age,2)) %>%
    mutate(aget = recode_factor(aget, '1' = 'Younger', '2' = 'Older')) %>%
    dplyr::select(mid, ADt, aget)
    
  meta.res <- meta %>% 
    mutate(srt2 = ntile(scale_rt2,n.rt.tiles)) %>% left_join(meta.subj) 
  
  srt.tilevals <- meta.res %>% group_by(srt2) %>%
    dplyr::summarise(scale_rt2 = mean(scale_rt2, na.rm=T))

  m.reg <- lmer(conf ~ AD.sc*scale_rt2 + gender*scale_rt2 + age*scale_rt2 +
                  stimev + (1|expnum), meta)
  
  summary(m.reg)
  
  meta.res2 <- meta.res %>% mutate(srt2 = as.numeric(srt2)) %>%
    mutate(ADt = recode_factor(ADt, "Higher" = '1', "Lower" = '-1')) %>%
    mutate(aget = recode_factor(aget, "Younger" = '24', "Older" = '44'))
  for (nrt in 1:n.rt.tiles){
    meta.res2$srt2[meta.res2$srt2==nrt] <- 
      srt.tilevals$scale_rt2[srt.tilevals$srt2==nrt]
  }

  if (accuracy==1){
    xlim = c(-1.5,1.7)
    ylim=c(0,.36)
  } else {
    xlim = c(-1.5,1.7)
    ylim=c(-.55,-.09)
  }
  plot_model(m.reg, type = ("pred"),line.size = 1, ci.lvl = .95,
             terms = c("scale_rt2", "AD.sc [-1,1]"),
             title='',
             colors = c('green4', 'firebrick2'))+
    theme_classic2()+
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(colour = 'AD')+
    ylab('Confidence (z-scored)') +
    xlab('Time (z-scored)')+
    theme(text = element_text(size=font_size),
          legend.position = 'off') +
    stat_summary(data = meta.res2,
                 fun.data = mean_cl_boot,
               aes(x=srt2,y=conf,colour=ADt),
               inherit.aes = FALSE)
  
  if (accuracy==1){
    xlim = c(-1.5,1.6)
    ylim=c(.03,.3)
  } else {
    xlim = c(-1.5,1.7)
    ylim=c(-.44,-.12)
  }
  plot_model(m.reg, type = ("pred"),line.size = 1, ci.lvl = .95,
             terms = c("scale_rt2", "gender"),
             title='', 
             colors = c('green4', 'firebrick2'))+
    theme_classic2()+
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(colour = 'AD')+
    ylab('Confidence (z-scored)') +
    xlab('Time (z-scored)')+
    theme(text = element_text(size=font_size),
          legend.position = 'off') +
    stat_summary(data = meta.res2,
                 fun.data = mean_cl_boot,
                 # fun = mean,
                 aes(x=srt2,y=conf,colour=gender),
                 inherit.aes = FALSE)
  
  if (accuracy==1){
    xlim = c(-1.5,1.6)
    ylim=c(0,.35)
  } else {
    xlim = c(-1.5,1.7)
    ylim=c(-.45,-.08)
  }
  plot_model(m.reg, type = ("pred"),line.size = 1, ci.lvl = .95,
             terms = c("scale_rt2", "age [24,44]"),
             title='', 
             colors = c('green4', 'firebrick2'))+
    theme_classic2()+
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(colour = 'AD')+
    ylab('Confidence (z-scored)') +
    xlab('Time (z-scored)') +
    theme(text = element_text(size=font_size),
          legend.position = 'off') +
    stat_summary(data = meta.res2,
                 # fun = mean,
                 fun.data = mean_cl_boot,
                 aes(x=srt2,y=conf,colour=aget),
                 inherit.aes = FALSE)
  
  m2 <- update(m.reg, ~.-AD.sc:scale_rt2)
  anova(m.reg, m2)

  m2 <- update(m.reg, ~.-gender:scale_rt2)
  anova(m.reg, m2)

  m2 <- update(m.reg, ~.-age:scale_rt2)
  anova(m.reg, m2)
  
  emm <- emtrends(m.reg, var = "scale_rt2" ,pairwise~AD.sc, at=list(AD.sc=c(1,-1)))
  test(emm)
  df <- summary(emm$contrasts)
  ci <- confint(emm)
  AD.inter[expNum,2-accuracy,] <- c(as.numeric(df[c(2,3,6)]), as.numeric(ci$contrasts[5:6]))

  emm <- emtrends(m.reg, var = "scale_rt2" ,pairwise~gender, at=list(gender=c('Women','Men')))
  test(emm)
  df <- summary(emm$contrasts)
  ci <- confint(emm)
  gender.inter[expNum,2-accuracy,] <- c(as.numeric(df[c(2,3,6)]), as.numeric(ci$contrasts[5:6]))
  
  emm <- emtrends(m.reg, var = "scale_rt2" ,pairwise~age, at=list(age=c(50,20)))
  test(emm)
  df <- summary(emm$contrasts)
  ci <- confint(emm)
  age.inter[expNum,2-accuracy,] <- c(as.numeric(df[c(2,3,6)]), as.numeric(ci$contrasts[5:6]))
  
}

corinc <- c('Incorrect', 'Correct')

for (accuracy in c(0,1)){
  
  meta <- filter(meta.df, accu==accuracy) %>% drop_na()
  
  m.reg <- lmer(conf ~ AD.sc*scale_rt2 + CIT.sc*scale_rt2 + gender*scale_rt2 +
                  age*scale_rt2 + stimev + (1|expnum), meta)
  summary(m.reg)
  
  m2 <- update(m.reg, ~.-CIT.sc:scale_rt2)
  anova(m.reg, m2)
  
  emm <- emtrends(m.reg, var = "scale_rt2" ,pairwise~CIT.sc, at=list(CIT.sc=c(1,-1)))
  test(emm)
  df <- summary(emm$contrasts)
  ci <- confint(emm)
  CIT.inter[expNum,2-accuracy,] <- c(as.numeric(df[c(2,3,6)]), as.numeric(ci$contrasts[5:6]))
  
  
  meta.subj <- meta %>% group_by(mid) %>%
    dplyr::summarise(CIT.sc = mean(CIT.sc),
                     age = mean(age)) %>%
    mutate(CITt = ntile(CIT.sc,2)) %>%
    mutate(CITt = recode_factor(CITt, '1' = 'Lower', '2' = 'Higher')) %>%
    dplyr::select(mid, CITt)
  
  meta.res <- meta %>% 
    mutate(srt2 = ntile(scale_rt2,n.rt.tiles)) %>% left_join(meta.subj) 
  
  srt.tilevals <- meta.res %>% group_by(srt2) %>%
    dplyr::summarise(scale_rt2 = mean(scale_rt2, na.rm=T))
  
  meta.res2 <- meta.res %>% mutate(srt2 = as.numeric(srt2)) %>%
    mutate(CITt = recode_factor(CITt, "Higher" = '1', "Lower" = '-1')) 
  for (nrt in 1:n.rt.tiles){
    meta.res2$srt2[meta.res2$srt2==nrt] <- 
      srt.tilevals$scale_rt2[srt.tilevals$srt2==nrt]
  }
  
  if (accuracy==1){
    xlim = c(-1.5,1.7)
    ylim=c(-.03,.42)
  } else {
    xlim = c(-1.5,1.7)
    ylim=c(-.58,-.05)
  }
  plot_model(m.reg, type = ("pred"),line.size = 1, ci.lvl = .95,
             terms = c("scale_rt2", "CIT.sc [-1,1]"),
             title='',
             colors = c('green4', 'firebrick2'))+
    theme_classic2()+
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(colour = 'CIT')+
    ylab('Confidence (z-scored)') +
    xlab('Time (z-scored)')+
    theme(text = element_text(size=font_size),
          legend.position = 'off') +
    stat_summary(data = meta.res2,
                 fun.data = mean_cl_boot,
                 aes(x=srt2,y=conf,colour=CITt),
                 inherit.aes = FALSE)
  
}

############
# library(ggpattern)

dim.all <- dim(AD.inter)

exp.labs <- c('Katyal\nE1', 'Katyal\nE2', 'Rouault\nE3', 'Rouault\nE4')
exp.labs[5] <- 'Combined'
# accu.colors =  c('lightskyblue2','lightskyblue4')
accu.colors =  c('deepskyblue2','blue4')
ds.colors.confmain =  c('white','khaki1')
ds.colors =  c('white','mistyrose')
dstype.colors <- c('')

exp.name <- array(exp.labs, c(dim.all))
exp.num <- array(c(1:5), c(dim.all))
ds.type <- array(c(rep(.5,4),1), c(dim.all))
accuracy <- aperm(array(c('Correct', 'Incorrect'), c(dim.all[c(2,1,3)])), c(2,1,3))
meas <- aperm(array(c('est', 'se', 'pval', 'cilo', 'cihi'), c(dim.all[c(3,1,2)])), c(2,3,1))

AD.main2 <- aperm(array(AD.main, c(dim(AD.main),2)), c(1,3,2))
CIT.main2 <- aperm(array(CIT.main, c(dim(CIT.main),2)), c(1,3,2))
gender.main2 <- aperm(array(gender.main, c(dim(CIT.main),2)), c(1,3,2))
age.main2 <- aperm(array(age.main, c(dim(CIT.main),2)), c(1,3,2))

df.all <- data.frame(c(AD.inter), c(CIT.inter), c(age.inter), c(gender.inter), 
                     c(exp.name), c(accuracy), c(meas), c(exp.num), c(ds.type),
                     c(AD.main2), c(CIT.main2), c(age.main2), c(gender.main2))
colnames(df.all) <- c('AD.inter', 'CIT.inter', 'age.inter', 'gender.inter', 
                     'exp.name', 'accuracy', 'meas', 'exp.num', 'ds.type',
                     'AD.main', 'CIT.main', 'age.main', 'gender.main')

df.all <- df.all %>%
  mutate(ds.type = as.factor(ds.type)) %>%
  pivot_wider(names_from = meas,
              values_from = c(AD.inter, CIT.inter, age.inter, gender.inter,
                              AD.main, CIT.main, age.main, gender.main)) %>%
  mutate(exp.name = factor(exp.name, levels = c(exp.labs[c(9, 1:8)]))) %>%
  mutate(exp.name = relevel(exp.name, 'Combined'))
  

studies <- c(1:4, 5)


ggplot(df.all %>% filter(exp.num %in% studies, accuracy=='Correct'), 
       aes(x = exp.name, y = AD.main_est, fill = ds.type)) +
  scale_fill_manual(values = c(ds.colors.confmain)) +
  geom_bar(stat = 'identity', position = position_dodge(width = .77), 
           width = .75, color='black', linewidth=.3) +
  geom_hline(yintercept = 0, color = 'black') +
  geom_errorbar(aes(ymin = AD.main_cilo, 
                    ymax = AD.main_cihi), width = .4,
                position = position_dodge(width = .77),
                linewidth=.75) +
  geom_point(shape=21, fill = 'white')+
  theme_classic2() +
  coord_cartesian(ylim = c(-.25,.04)) +
  labs(fill = 'Accuracy') + 
  xlab('Dataset') + 
  ylab('Beta upon confidence') +
  # ylab('Main effect upon confidence:\nAnxiety-Depression') +
  theme(legend.position = 'none',
        text = element_text(size=font_size-2),
        axis.title.y = element_text(hjust = 1.5)) + 
  geom_signif(y_position = -.18,
              xmin = c(1),
              xmax = c(1),
              annotation = '****',
              tip_length = 0, color = 'black', textsize = 5, size = 0)

ggplot(df.all %>% filter(exp.num %in% studies, accuracy=='Correct'), 
       aes(x = exp.name, y = CIT.main_est, fill = ds.type)) +
  scale_fill_manual(values = c(ds.colors.confmain)) +
  geom_bar(stat = 'identity', position = position_dodge(width = .77), 
           width = .75, color='black', linewidth=.3) +
  geom_hline(yintercept = 0, color = 'black') +
  geom_errorbar(aes(ymin = CIT.main_cilo, 
                    ymax = CIT.main_cihi), width = .4,
                position = position_dodge(width = .77),
                linewidth=.75) +
  geom_point(shape=21, fill = 'white')+
  theme_classic2() +
  coord_cartesian(ylim = c(-.05,.22)) +
  labs(fill = 'Accuracy') + 
  xlab('Dataset') + 
  ylab('Beta upon confidence') +
  # ylab('Main effect upon confidence:\nAnxiety-Depression') +
  theme(legend.position = 'none',
        text = element_text(size=font_size-2),
        axis.title.y = element_text(hjust = 1.5)) + 
  guides(colour="none") +
  geom_signif(y_position = .12,
              xmin = c(1),
              xmax = c(1),
              annotation = '****',
              tip_length = 0, color = 'black', textsize = 5, size = 0)

ggplot(df.all %>% filter(exp.num %in% studies, accuracy=='Correct'), 
       aes(x = exp.name, y = -gender.main_est, fill = ds.type)) +
  scale_fill_manual(values = c(ds.colors.confmain)) +
  geom_bar(stat = 'identity', position = position_dodge(width = .77), 
           width = .75, color='black', linewidth=.3) +
  geom_hline(yintercept = 0, color = 'black') +
  geom_errorbar(aes(ymin = -gender.main_cilo, 
                    ymax = -gender.main_cihi), width = .4,
                position = position_dodge(width = .77),
                linewidth=.75) +
  geom_point(shape=21, fill = 'white')+
  theme_classic2() +
  coord_cartesian(ylim = c(-.15,.05)) +
  labs(fill = 'Accuracy') + 
  xlab('Dataset') + 
  ylab('Beta upon confidence') +
  theme(legend.position = 'none',
        text = element_text(size=font_size-2),
        axis.title.y = element_text(hjust = 1.5)) + 
  guides(colour="none") +
  geom_signif(y_position = -.13,
              xmin = c(1),
              xmax = c(1),
              annotation = '****',
              tip_length = 0, color = 'black', textsize = 5, size = 0)

ggplot(df.all %>% filter(exp.num %in% studies, accuracy=='Correct'), 
       aes(x = exp.name, y = age.main_est, fill = ds.type)) +
  scale_fill_manual(values = c(ds.colors.confmain)) +
  geom_bar(stat = 'identity', position = position_dodge(width = .77), 
           width = .75, color='black', linewidth=.3) +
  geom_hline(yintercept = 0, color = 'black') +
  geom_errorbar(aes(ymin = age.main_cilo, 
                    ymax = age.main_cihi), width = .4,
                position = position_dodge(width = .77),
                linewidth=.75) +
  geom_point(shape=21, fill = 'white')+
  theme_classic2() +
  coord_cartesian(ylim = c(-.015,.004)) +
  labs(fill = 'Accuracy') + 
  xlab('Dataset') + 
  ylab('Beta upon confidence') +
  # ylab('Main effect upon confidence:\nAnxiety-Depression') +
  theme(legend.position = 'none',
        text = element_text(size=font_size-2),
        axis.title.y = element_text(hjust = 1.5)) + 
  guides(colour="none") +
  geom_signif(y_position = -.01,
              xmin = c(1),
              xmax = c(1),
              annotation = '****',
              tip_length = 0, color = 'black', textsize = 5, size = 0)


sig.y <- c(-.11, -.1)
ylim <- c(-.16,.05)

ggplot(df.all %>% filter(exp.num %in% studies), 
       aes(x = exp.name, y = AD.inter_est, colour = accuracy, fill = ds.type)) +
  scale_colour_manual(values = accu.colors) +
  scale_fill_manual(values = c(ds.colors)) +
  geom_bar(stat = 'identity', position = position_dodge(width = .77), width = .75) +
  geom_hline(yintercept = 0, color = 'darkgrey') +
  geom_errorbar(aes(ymin = AD.inter_cilo, 
                    ymax = AD.inter_cihi), width = .3,
                position = position_dodge(width = .77),
                linewidth=.75) +
  geom_point(position = position_dodge(width = .77), shape=21, fill = 'white')+
  # scale_color_manual(values = 'black')+
  theme_classic2() +
  coord_cartesian(ylim = ylim) +
  scale_y_continuous(breaks = seq(-.1, .1, .1))+
  labs(colour = 'Accuracy') + 
  xlab('Dataset') + ylab('Slope difference\n(higher - lower)') +
  theme(legend.position = 'top',
        text = element_text(size=font_size-2)) + 
  guides(fill="none",
         color=guide_legend(override.aes=list(fill=NA))) +
  geom_signif(y_position = sig.y,
              xmin = c(.79, 1.23),
              xmax = c(.79, 1.23),
              annotation = c('****', '****'),
              tip_length = 0, color = 'black', textsize = 5, size = 0)

sig.y <- c(.06, .07)
ylim <- c(-.05,.11)
anno <- c('****', '****')

ggplot(df.all %>% filter(exp.num %in% studies), 
       aes(x = exp.name, y = gender.inter_est, colour = accuracy, fill = ds.type)) +
  scale_colour_manual(values = accu.colors) +
  scale_fill_manual(values = c(ds.colors)) +
  geom_bar(stat = 'identity', position = position_dodge(width = .77), width = .75) +
  geom_hline(yintercept = 0, color = 'darkgrey') +
  geom_errorbar(aes(ymin = gender.inter_cilo, 
                    ymax = gender.inter_cihi), width = .4,
                position = position_dodge(width = .77),
                linewidth=.75) +
  geom_point(position = position_dodge(width = .77), shape=21, fill = 'white')+
  # scale_color_manual(values = 'black')+
  theme_classic2() +
  coord_cartesian(ylim = ylim) +
  labs(colour = 'Accuracy') + 
  xlab('Dataset') + ylab('Slope difference\n(women - men)') +
  theme(legend.position = 'top',
        text = element_text(size=font_size-2)) + 
  guides(fill="none",
         color=guide_legend(override.aes=list(fill=NA))) +
  geom_signif(y_position = sig.y,
              xmin = c(.79, 1.23),
              xmax = c(.79, 1.23),
              annotation = anno,
              tip_length = 0, color = 'black', textsize = 5, size = 0)

sig.y <- c(-.1, -.105)
ylim <- c(-.13,.06)

ggplot(df.all %>% filter(exp.num %in% studies), 
       aes(x = exp.name, y = -CIT.inter_est, colour = accuracy, fill = ds.type)) +
  scale_colour_manual(values = accu.colors) +
  scale_fill_manual(values = c(ds.colors)) +
  geom_bar(stat = 'identity', position = position_dodge(width = .77), width = .75) +
  geom_hline(yintercept = 0, color = 'darkgrey') +
  geom_errorbar(aes(ymin = -(CIT.inter_cilo), 
                    ymax = -(CIT.inter_cihi)), width = .3,
                position = position_dodge(width = .77),
                linewidth=.75) +
  geom_point(position = position_dodge(width = .77), shape=21, fill = 'white')+
  # scale_color_manual(values = 'black')+
  theme_classic2() +
  coord_cartesian(ylim = ylim) +
  scale_y_continuous(breaks = seq(-.1, .1, .1))+
  labs(colour = 'Accuracy') + 
  xlab('Dataset') + ylab('Slope difference\n(lower - higher)') +
  theme(legend.position = c(.5,.9),
        legend.direction = 'horizontal',
        text = element_text(size=font_size-2)) + 
  guides(fill="none",
         color=guide_legend(override.aes=list(fill=NA))) +
  geom_signif(y_position = sig.y,
              xmin = c(.79, 1.23),
              xmax = c(.79, 1.23),
              annotation = c('****', '****'),
              tip_length = 0, color = 'black', textsize = 5, size = 0)

sig.y <-  c(.03,.05)
ylim <- c(-.1,.19)
anno <- c('****', '****')

ggplot(df.all %>% filter(exp.num %in% studies), 
       aes(x = exp.name, y = age.inter_est, colour = accuracy, fill = ds.type)) +
  scale_colour_manual(values = accu.colors) +
  scale_fill_manual(values = c(ds.colors)) +
  geom_bar(stat = 'identity', position = position_dodge(width = .77), width = .75) +
  geom_hline(yintercept = 0, color = 'darkgrey') +
  geom_errorbar(aes(ymin = age.inter_est-age.inter_se, 
                    ymax = age.inter_est+age.inter_se), width = .4,
                position = position_dodge(width = .77),
                linewidth=.75) +
  geom_point(position = position_dodge(width = .77), shape=21, fill = 'white')+
  # scale_color_manual(values = 'black')+
  theme_classic2() +
  coord_cartesian(ylim = ylim) +
  labs(colour = 'Accuracy') + 
  xlab('Dataset') + ylab('Slope difference\n(older - younger)') +
  theme(legend.position = 'top',
        text = element_text(size=font_size-2)) + 
  guides(fill="none",
         color=guide_legend(override.aes=list(fill=NA))) +
  geom_signif(y_position = sig.y,
              xmin = c(.78, 1.2),
              xmax = c(.79, 1.2),
              annotation = anno,
              tip_length = 0, color = 'black', textsize = 5, size = 0)


