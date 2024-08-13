
################################## Katyal & Fleming ################################
### Gender and anxiety reveal distinct computational sources of underconfidence ####
####################################################################################
################  This file reproduces Figure 3H and 3J ##################
###############################################################################
####### Plot and do stats on model estimated parameters  ##############
###############################################################################

library(ggplot2)
library(ggpubr)
library(lmerTest)
library(corrplot)
library(psych)
library(tidyverse)
library(lmtest)

setwd("~/OneDrive - University of Copenhagen/Projects/Experiments/confidenceTime/analysis/model/meta-ddm/")

sim.rec <- readRDS('simulate-recovery.Rds')

df <- data.frame(sim.rec$v_bias, sim.rec$conf_add, sim.rec$conf_mult, 
                 sim.rec$vratio,
                 c(sim.rec$recover[,5]), c(sim.rec$recover[,6]), 
                 c(sim.rec$recover[,7]), c(sim.rec$recover[,8]))
colnames(df) <- c('vb.sim', 'add.sim', 'mult.sim', 'vr.sim',
                  'add.rec', 'mult.rec', 'vb.rec', 'vr.rec')

ggplot(df, aes(x=vb.sim,y=vb.rec)) + 
  geom_point() +
  geom_abline(slope=1, intercept=0) +
  theme_classic2() +
  coord_cartesian(ylim = c(-2.3,1), xlim = c(-2.3,1)) +
  xlab('Vbias simulated') +
  ylab('Vbias recovered') +
  stat_cor(label.y.npc=.5, label.x.npc = .6, 
           method = 'pearson', size = 5)+
  theme(text = element_text(size=16),
        legend.position = 'none')

ggplot(df, aes(x=add.sim,y=add.rec)) + geom_point() +
  geom_abline(slope=1, intercept=0) +
  theme_classic2() +
  coord_cartesian(ylim = c(-2.7,2.7), xlim = c(-2.7,2.7)) +
  xlab('Abias simulated') +
  ylab('Abias recovered') +
  stat_cor(label.y.npc=.15, label.x.npc = .3, 
           method = 'pearson', size = 5)+
  theme(text = element_text(size=16),
        legend.position = 'none')

ggplot(df, aes(x=mult.sim,y=mult.rec)) + geom_point() +
  geom_abline(slope=1, intercept=0) +
  theme_classic2() +
  coord_cartesian(ylim = c(0,2.5), xlim = c(0,2.5)) +
  xlab('Mbias simulated') +
  ylab('Mbias recovered') +
  stat_cor(label.y.npc=.1, label.x.npc = .3, 
           method = 'pearson', size = 5)+
  theme(text = element_text(size=16),
        legend.position = 'none')

ggplot(df, aes(x=vr.sim,y=vr.rec)) + geom_point() +
  geom_abline(slope=1, intercept=0) +
  theme_classic2() +
  coord_cartesian(ylim = c(-.5,3), xlim = c(-.5,3)) +
  xlab('Vratio simulated') +
  ylab('Vratio recovered') +
  stat_cor(label.y.npc=.1, label.x.npc = .4, 
           method = 'pearson', size = 5)+
  theme(text = element_text(size=16),
        legend.position = 'none')


cor.mat <- cor(drop_na(df))
cor.test <- corr.test(drop_na(df))
corrplot(cor.mat, method = 'ellipse', type='lower', diag = F,
         p.mat = cor.test$p)


########################
########### combined

setwd("~/OneDrive - University of Copenhagen/Projects/Experiments/confidenceTime/analysis/model/meta-ddm/")

ol.sd <- 5

k1.perc <- readRDS('katyal_1-Perception-vratio.Rds')

k1.perc <- data.frame(k1.perc$AD, k1.perc$gender, k1.perc$age, 
                      k1.perc$recover[,7], k1.perc$recover[,6], 
                      k1.perc$recover[,5], k1.perc$recover[,8],
                      c(k1.perc$bestfnval), k1.perc$subj, 
                      k1.perc$recover[,1])
colnames(k1.perc) <- c('AD', 'gender', 'age', 
                       'v.bias', 'mult.bias', 'add.bias', 'vratio', 
                       'fnval', 'subj', 'v')

k1.perc$task <- 'Perception'

mv <- mean(k1.perc$vratio,na.rm=T)
outlier <- ol.sd*sd(k1.perc$vratio,na.rm=T)
k1.perc$vratio[k1.perc$vratio<(mv-outlier) | k1.perc$vratio>(mv+outlier) ] <- NA
k1.perc$vratio <- c(scale(k1.perc$vratio))

mv <- mean(k1.perc$v.bias,na.rm=T)
outlier <- ol.sd*sd(k1.perc$v.bias,na.rm=T)
k1.perc$v.bias[k1.perc$v.bias<(mv-outlier) | k1.perc$v.bias>(mv+outlier) ] <- NA
k1.perc$v.bias <- c(scale(k1.perc$v.bias))

k1.perc$mult.bias <- c(scale(k1.perc$mult.bias))

mv <- mean(k1.perc$add.bias,na.rm=T)
outlier <- ol.sd*sd(k1.perc$add.bias,na.rm=T)
k1.perc$add.bias[k1.perc$add.bias<(mv-outlier) | k1.perc$add.bias>(mv+outlier) ] <- NA
k1.perc$add.bias <- c(scale(k1.perc$add.bias))
k1.perc$v <- c(scale(k1.perc$v))


k1.mem <- readRDS('katyal_1-Memory-vratio.Rds')

k1.mem <- data.frame(k1.mem$AD, k1.mem$gender, k1.mem$age, 
                     k1.mem$recover[,7], k1.mem$recover[,6], 
                     k1.mem$recover[,5], k1.mem$recover[,8],
                     c(k1.mem$bestfnval), k1.mem$subj, 
                     k1.mem$recover[,1])
colnames(k1.mem) <- c('AD', 'gender', 'age', 
                      'v.bias', 'mult.bias', 'add.bias', 'vratio', 
                      'fnval', 'subj', 'v')
k1.mem$task <- 'Memory'

mv <- mean(k1.mem$vratio,na.rm=T)
outlier <- ol.sd*sd(k1.mem$vratio,na.rm=T)
k1.mem$vratio[k1.mem$vratio<(mv-outlier) | k1.mem$vratio>(mv+outlier) ] <- NA
k1.mem$vratio <- c(scale(k1.mem$vratio))

mv <- mean(k1.mem$v.bias,na.rm=T)
outlier <- ol.sd*sd(k1.mem$v.bias,na.rm=T)
k1.mem$v.bias[k1.mem$v.bias<(mv-outlier) | k1.mem$v.bias>(mv+outlier) ] <- NA
k1.mem$v.bias <- c(scale(k1.mem$v.bias))

k1.mem$mult.bias <- c(scale(k1.mem$mult.bias))

mv <- mean(k1.mem$add.bias,na.rm=T)
outlier <- ol.sd*sd(k1.mem$add.bias,na.rm=T)
k1.mem$add.bias[k1.mem$add.bias<(mv-outlier) | k1.mem$add.bias>(mv+outlier) ] <- NA
k1.mem$add.bias <- c(scale(k1.mem$add.bias))
k1.mem$v <- c(scale(k1.mem$v))

k1 <- rbind(k1.perc, k1.mem) %>% group_by(subj, age, gender) %>%
  dplyr::summarise(AD = mean(AD, na.rm=T),
                   v.bias = mean(v.bias, na.rm=T),
                   mult.bias = mean(mult.bias, na.rm=T),
                   add.bias = mean(add.bias, na.rm=T),
                   vratio = mean(vratio, na.rm=T),
                   v = mean(v, na.rm=T)) 
k1$AD <- c(scale(k1$AD))
k1$CIT <- NA
k1$SW <- NA
k1$expnum <- 1
k1$subj <- as.factor(as.numeric(k1$subj)+10000)


k2.perc <- readRDS('katyal_2-Perception-vratio.Rds')

k2.perc <- data.frame(k2.perc$AD, k2.perc$CIT, k2.perc$SW, 
                      k2.perc$gender, k2.perc$age, 
                      k2.perc$recover[,7], k2.perc$recover[,6], 
                      k2.perc$recover[,5], k2.perc$recover[,8],
                      c(k2.perc$bestfnval), k2.perc$subj, 
                      k2.perc$recover[,1])
colnames(k2.perc) <- c('AD', 'CIT', 'SW', 'gender', 'age', 
                       'v.bias', 'mult.bias', 'add.bias', 'vratio', 
                       'fnval', 'subj', 'v')
k2.perc$task <- 'Perception'

mv <- mean(k2.perc$vratio,na.rm=T)
outlier <- ol.sd*sd(k2.perc$vratio,na.rm=T)
k2.perc$vratio[k2.perc$vratio<(mv-outlier) | k2.perc$vratio>(mv+outlier) ] <- NA
k2.perc$vratio <- c(scale(k2.perc$vratio))

mv <- mean(k2.perc$v.bias,na.rm=T)
outlier <- ol.sd*sd(k2.perc$v.bias,na.rm=T)
k2.perc$vratio[k2.perc$v.bias<(mv-outlier) | k2.perc$v.bias>(mv+outlier) ] <- NA
k2.perc$v.bias <- c(scale(k2.perc$v.bias))

k2.perc$mult.bias <- c(scale(k2.perc$mult.bias))

mv <- mean(k2.perc$add.bias,na.rm=T)
outlier <- ol.sd*sd(k2.perc$add.bias,na.rm=T)
k2.perc$add.bias[k2.perc$add.bias<(mv-outlier) | k2.perc$add.bias>(mv+outlier) ] <- NA
k2.perc$add.bias <- c(scale(k2.perc$add.bias))
k2.perc$v <- c(scale(k2.perc$v))


k2.mem <- readRDS('katyal_2-Memory-vratio.Rds')
k2.mem <- data.frame(k2.mem$AD, k2.mem$CIT, k2.perc$SW, 
                     k2.mem$gender, k2.mem$age, 
                     k2.mem$recover[,7], k2.mem$recover[,6], 
                     k2.mem$recover[,5], k2.mem$recover[,8],
                     c(k2.mem$bestfnval), k2.mem$subj, 
                     k2.mem$recover[,1])
colnames(k2.mem) <- c('AD', 'CIT', 'SW', 'gender', 'age', 
                      'v.bias', 'mult.bias', 'add.bias', 'vratio', 
                      'fnval', 'subj', 'v')
k2.mem$task <- 'Memory'
mv <- mean(k2.mem$vratio,na.rm=T)
outlier <- ol.sd*sd(k2.mem$vratio,na.rm=T)
k2.mem$vratio[k2.mem$vratio<(mv-outlier) | k2.mem$vratio>(mv+outlier) ] <- NA
k2.mem$vratio <- c(scale(k2.mem$vratio))

mv <- mean(k2.mem$v.bias,na.rm=T)
outlier <- ol.sd*sd(k2.mem$v.bias,na.rm=T)
k2.mem$v.bias[k2.mem$v.bias<(mv-outlier) | k2.mem$v.bias>(mv+outlier) ] <- NA
k2.mem$v.bias <- c(scale(k2.mem$v.bias))

k2.mem$mult.bias <- c(scale(k2.mem$mult.bias))

mv <- mean(k2.mem$add.bias,na.rm=T)
outlier <- ol.sd*sd(k2.mem$add.bias,na.rm=T)
k2.mem$add.bias[k2.mem$add.bias<(mv-outlier) | k2.mem$add.bias>(mv+outlier) ] <- NA
k2.mem$add.bias <- c(scale(k2.mem$add.bias))
k2.mem$v <- c(scale(k2.mem$v))

k2 <- rbind(k2.perc, k2.mem) %>% group_by(subj, age, gender) %>%
  dplyr::summarise(AD = mean(AD, na.rm=T),
                   CIT = mean(CIT, na.rm=T),
                   SW = mean(SW, na.rm=T),
                   v.bias = mean(v.bias, na.rm=T),
                   mult.bias = mean(mult.bias, na.rm=T),
                   add.bias = mean(add.bias, na.rm=T),
                   vratio = mean(vratio, na.rm=T),
                   v = mean(v, na.rm=T)) 
k2$AD <- c(scale(k2$AD))
k2$CIT <- c(scale(k2$CIT))
k2$SW <- c(scale(k2$SW))
k2$expnum <- 2


r1 <- readRDS('rouault_1-vratio.Rds')

r1 <- data.frame(r1$AD, r1$gender, r1$age, 
                 r1$recover[,7], r1$recover[,6], 
                 r1$recover[,5], r1$recover[,8],
                 r1$subj, 
                 r1$recover[,1])
colnames(r1) <- c('AD', 'gender', 'age', 
                  'v.bias', 'mult.bias', 'add.bias', 'vratio', 
                  'subj', 'v')
mv <- mean(r1$vratio,na.rm=T)
outlier <- ol.sd*sd(r1$vratio,na.rm=T)
r1$vratio[r1$vratio<(mv-outlier) | r1$vratio>(mv+outlier) ] <- NA
r1$vratio <- c(scale(r1$vratio))

mv <- mean(r1$v.bias,na.rm=T)
outlier <- ol.sd*sd(r1$v.bias,na.rm=T)
r1$v.bias[r1$v.bias<(mv-outlier) | r1$v.bias>(mv+outlier) ] <- NA
r1$v.bias <- c(scale(r1$v.bias))

r1$mult.bias <- c(scale(r1$mult.bias))

mv <- mean(r1$add.bias,na.rm=T)
outlier <- ol.sd*sd(r1$add.bias,na.rm=T)
r1$add.bias[r1$add.bias<(mv-outlier) | r1$add.bias>(mv+outlier) ] <- NA
r1$add.bias <- c(scale(r1$add.bias))
r1$v <- c(scale(r1$v))

r1$AD <- c(scale(r1$AD))
r1$CIT <- NA
r1$SW <- NA
r1$expnum <- 3
r1$subj <- as.factor(r1$subj)


r2 <- readRDS('rouault_2-vratio.Rds')

r2 <- data.frame(r2$AD, r2$CIT, r2$SW, r2$gender, r2$age, 
                 r2$recover[,7], r2$recover[,6], 
                 r2$recover[,5], r2$recover[,8],
                 r2$subj, 
                 r2$recover[,1])
colnames(r2) <- c('AD', 'CIT', 'SW', 'gender', 'age', 
                  'v.bias', 'mult.bias', 'add.bias', 'vratio', 
                  'subj', 'v')

mv <- mean(r2$vratio,na.rm=T)
outlier <- ol.sd*sd(r2$vratio,na.rm=T)
r2$vratio[r2$vratio<(mv-outlier) | r2$vratio>(mv+outlier) ] <- NA
r2$vratio <- c(scale(r2$vratio))

mv <- mean(r2$v.bias,na.rm=T)
outlier <- ol.sd*sd(r2$v.bias,na.rm=T)
r2$v.bias[r2$v.bias<(mv-outlier) | r2$v.bias>(mv+outlier) ] <- NA
r2$v.bias <- c(scale(r2$v.bias))

r2$mult.bias <- c(scale(r2$mult.bias))

mv <- mean(r2$add.bias,na.rm=T)
outlier <- ol.sd*sd(r2$add.bias,na.rm=T)
r2$add.bias[r2$add.bias<(mv-outlier) | r2$add.bias>(mv+outlier) ] <- NA
r2$add.bias <- c(scale(r2$add.bias))

r2$v <- c(scale(r2$v))
r2$AD <- c(scale(r2$AD))
r2$CIT <- c(scale(r2$CIT))
r2$SW <- c(scale(r2$SW))
r2$expnum <- 4
r2$subj <- as.factor(r2$subj)



df <- rbind(k1, k2, r1, r2) %>%
  mutate(gender = as.factor(gender)) %>%
  mutate(gender = factor(gender, levels = c('Male', 'Female')))


reg.coef     <- array(NA, c(4, 4, 4))

m.reg <- lm(AD ~ v.bias + add.bias + mult.bias + vratio + v, 
            df %>% drop_na(AD, v.bias, add.bias, mult.bias, vratio))
s.reg <- summary(m.reg)
s.reg
m2 <- update(m.reg, ~.-v.bias)
lrtest(m.reg,m2)
m2 <- update(m.reg, ~.-add.bias)
lrtest(m.reg,m2)
m2 <- update(m.reg, ~.-mult.bias)
lrtest(m.reg,m2)
m2 <- update(m.reg, ~.-vratio)
lrtest(m.reg,m2)

reg.coef[1,,c(1,4)] <- s.reg$coefficients[c(2:5),c(1,4)]
ci <- confint(m.reg)
reg.coef[1,,c(2:3)] <- ci[2:5,]

m.reg <- glmer(gender ~ v.bias + add.bias + mult.bias + vratio + v +
                 (1|expnum), family = binomial(link = 'logit'),
               df %>% drop_na(gender, v.bias, add.bias, mult.bias, vratio))
s.reg <- summary(m.reg)
s.reg
m2 <- update(m.reg, ~.-v.bias)
anova(m.reg,m2)
m2 <- update(m.reg, ~.-add.bias)
anova(m.reg,m2)
m2 <- update(m.reg, ~.-mult.bias)
anova(m.reg,m2)
m2 <- update(m.reg, ~.-vratio)
anova(m.reg,m2)

reg.coef[3,,c(1,4)] <- s.reg$coefficients[c(2:5),c(1,4)]
ci <- confint(m.reg)
reg.coef[3,,c(2:3)] <- ci[3:6,]


m.reg <- lm(CIT ~ v.bias + add.bias + mult.bias + vratio + AD + SW + v, 
            df %>% drop_na(CIT, v.bias, add.bias, mult.bias, vratio))
s.reg <- summary(m.reg)
s.reg
m2 <- update(m.reg, ~.-v.bias)
lrtest(m.reg,m2)
m2 <- update(m.reg, ~.-add.bias)
lrtest(m.reg,m2)
m2 <- update(m.reg, ~.-mult.bias)
lrtest(m.reg,m2)
m2 <- update(m.reg, ~.-vratio)
lrtest(m.reg,m2)

reg.coef[2,,c(1,4)] <- s.reg$coefficients[c(2:5),c(1,4)]
ci <- confint(m.reg)
reg.coef[2,,c(2:3)] <- ci[2:5,]

m.reg <- lmer(age ~ v.bias + add.bias + mult.bias + vratio + v +
                (1 |expnum), df %>% drop_na(age, v.bias, add.bias, mult.bias, vratio))
s.reg <- summary(m.reg)
s.reg
m2 <- update(m.reg, ~.-v.bias)
anova(m.reg,m2)
m2 <- update(m.reg, ~.-add.bias)
anova(m.reg,m2)
m2 <- update(m.reg, ~.-mult.bias)
anova(m.reg,m2)
m2 <- update(m.reg, ~.-vratio)
anova(m.reg,m2)
reg.coef[4,,c(1,4)] <- s.reg$coefficients[c(2:5),c(1,5)]
ci <- confint(m.reg)
reg.coef[4,,c(2:3)] <- ci[4:7,]



model.params <- c('V-bias', 'A-bias', 'M-bias', 'V-ratio');
indiv.fact <- array(c('Anxious-Depression', 'Compulsivity', 'Gender', 'Age'), 
                    c(4,4,4))
regressor <- aperm(array(model.params, 
                         c(4,4,4)), c(2,1,3))
stat <- aperm(array(c('est', 'cilo', 'cihi', 'pval'), c(4,4,4)), c(2,3,1))

pf <- data.frame(c(reg.coef), c(indiv.fact), c(regressor), c(stat))
colnames(pf) <- c('reg.coef', 'indiv.fact', 'regressor', 'stat')

pf <- pf %>% pivot_wider(names_from = stat,
                         values_from = c(reg.coef)) %>%
  mutate(regressor = factor(regressor, levels = c(model.params))) 

font_size <- 16

ggplot(pf %>% filter(indiv.fact=='Anxious-Depression'), 
       aes(x = regressor, y = est)) +
  geom_bar(stat = 'identity', position = position_dodge(width = .77), 
           width = .75, color='black', linewidth=.3, fill = 'deepskyblue1') +
  geom_hline(yintercept = 0, color = 'black') +
  geom_errorbar(aes(ymin = cilo, 
                    ymax = cihi), width = .3,
                position = position_dodge(width = .77),
                linewidth=.75) +
  geom_point(position = position_dodge(width = .77), shape=21, fill = 'white',
             size=2)+
  theme_classic2() +
  coord_cartesian(ylim = c(-.34,.35)) +
  labs(fill = 'Accuracy') + 
  xlab('Model parameters') + 
  ylab('Regression slope') +
  # ylab('Main effect upon confidence:\n') +
  theme(#legend.position = 'none',
    text = element_text(size=font_size-4)) + 
  # guides(colour="none") +
  geom_signif(y_position = c(-.31, .3, -.23),
              xmin = c(1.1, 2, 4),
              xmax = c(1.1, 2, 4),
              # annotation = c('****', '*'),
              annotation = c('p < .0001', 'p < .0001', 'p = .004'),
              tip_length = 0, color = 'black', textsize = 3, size = 0)

ggplot(pf %>% filter(indiv.fact=='Gender'), 
       aes(x = regressor, y = est)) +
  # scale_fill_manual(values = c(ds.colors)) +
  geom_bar(stat = 'identity', position = position_dodge(width = .77), 
           width = .75, color='black', linewidth=.3, fill = 'deepskyblue1') +
  geom_hline(yintercept = 0, color = 'black') +
  geom_errorbar(aes(ymin = cilo, 
                    ymax = cihi), width = .3,
                position = position_dodge(width = .77),
                linewidth=.75) +
  geom_point(position = position_dodge(width = .77), shape=21, fill = 'white',
             size=2)+
  theme_classic2() +
  coord_cartesian(ylim = c(-.35,.42)) +
  labs(fill = 'Accuracy') + 
  xlab('Model parameters') + 
  ylab('Regression slope') +
  # ylab('Main effect upon confidence:\nAnxiety-Depression') +
  theme(legend.position = 'none',
        text = element_text(size=font_size-4)) + 
  guides(colour="none") +
  geom_signif(y_position = c(.36),
              xmin = c(2),
              xmax = c(2),
              # annotation = c('****', '****'),
              annotation = c('p = .012'),
              tip_length = 0, color = 'black', textsize = 3, size = 0)

ggplot(pf %>% filter(indiv.fact=='Age'), 
       aes(x = regressor, y = est)) +
  geom_bar(stat = 'identity', position = position_dodge(width = .77), 
           width = .75, color='black', linewidth=.3, fill = 'deepskyblue1') +
  geom_hline(yintercept = 0, color = 'black') +
  geom_errorbar(aes(ymin = cilo, 
                    ymax = cihi), width = .3,
                position = position_dodge(width = .77),
                linewidth=.75) +
  geom_point(position = position_dodge(width = .77), shape=21, fill = 'white',
             size=2)+
  theme_classic2() +
  coord_cartesian(ylim = c(-2.5,1.9)) +
  labs(fill = 'Accuracy') + 
  xlab('Model parameters') + 
  ylab('Regression slope') +
  theme(legend.position = 'none',
        text = element_text(size=font_size-4)) + 
  guides(colour="none") +
  geom_signif(y_position = c(1.6, -2.5),
              xmin = c(2,3),
              xmax = c(2,3),
              # annotation = c('****', '****'),
              annotation = c('p = .048', 'p < .0001'),
              tip_length = 0, color = 'black', textsize = 3, size = 0)

ggplot(pf %>% filter(indiv.fact=='Compulsivity'), 
       aes(x = regressor, y = est)) +
  geom_bar(stat = 'identity', position = position_dodge(width = .77), 
           width = .75, color='black', linewidth=.3, fill = 'deepskyblue1') +
  geom_hline(yintercept = 0, color = 'black') +
  geom_errorbar(aes(ymin = cilo, 
                    ymax = cihi), width = .3,
                position = position_dodge(width = .77),
                linewidth=.75) +
  geom_point(position = position_dodge(width = .77), shape=21, fill = 'white',
             size=2)+
  theme_classic2() +
  coord_cartesian(ylim = c(-.36,.29)) +
  labs(fill = 'Accuracy') + 
  xlab('Model parameters') + 
  ylab('Regression slope') +
  theme(legend.position = 'none',
        text = element_text(size=font_size-4)) + 
  # guides(colour="none") +
  geom_signif(y_position = c(.2, -.31),
              xmin = c(1.1, 2),
              xmax = c(1.1, 2),
              annotation = c('p = .032', 'p = .0009'),
              tip_length = 0, color = 'black', textsize = 3, size = 0)

#########
#######
### Descriptve stats on model params

setwd("~/OneDrive - University of Copenhagen/Projects/Experiments/confidenceTime/analysis/model/meta-ddm/")

k1.perc <- readRDS('katyal_1-Perception-vratio.Rds')

k1.perc <- data.frame(k1.perc$AD, k1.perc$gender, k1.perc$age, 
                      k1.perc$recover[,7], k1.perc$recover[,6], 
                      k1.perc$recover[,5], k1.perc$recover[,8],
                      c(k1.perc$bestfnval), k1.perc$subj, 
                      k1.perc$recover[,1])
colnames(k1.perc) <- c('AD', 'gender', 'age', 
                       'v.bias', 'mult.bias', 'add.bias', 'vratio', 
                       'fnval', 'subj', 'v')

k1.perc$task <- 'Perception'
mv <- mean(k1.perc$vratio,na.rm=T)
outlier <- ol.sd*sd(k1.perc$vratio,na.rm=T)
k1.perc$vratio[k1.perc$vratio<(mv-outlier) | k1.perc$vratio>(mv+outlier) ] <- NA

mv <- mean(k1.perc$v.bias,na.rm=T)
outlier <- ol.sd*sd(k1.perc$v.bias,na.rm=T)
k1.perc$v.bias[k1.perc$v.bias<(mv-outlier) | k1.perc$v.bias>(mv+outlier) ] <- NA

mv <- mean(k1.perc$add.bias,na.rm=T)
outlier <- ol.sd*sd(k1.perc$add.bias,na.rm=T)
k1.perc$add.bias[k1.perc$add.bias<(mv-outlier) | k1.perc$add.bias>(mv+outlier) ] <- NA


k1.mem <- readRDS('katyal_1-Memory-vratio.Rds')

k1.mem <- data.frame(k1.mem$AD, k1.mem$gender, k1.mem$age, 
                     k1.mem$recover[,7], k1.mem$recover[,6], 
                     k1.mem$recover[,5], k1.mem$recover[,8],
                     c(k1.mem$bestfnval), k1.mem$subj, 
                     k1.mem$recover[,1])
colnames(k1.mem) <- c('AD', 'gender', 'age', 
                      'v.bias', 'mult.bias', 'add.bias', 'vratio', 
                      'fnval', 'subj', 'v')
k1.mem$task <- 'Memory'

mv <- mean(k1.mem$vratio,na.rm=T)
outlier <- ol.sd*sd(k1.mem$vratio,na.rm=T)
k1.mem$vratio[k1.mem$vratio<(mv-outlier) | k1.mem$vratio>(mv+outlier) ] <- NA

mv <- mean(k1.mem$v.bias,na.rm=T)
outlier <- ol.sd*sd(k1.mem$v.bias,na.rm=T)
k1.mem$v.bias[k1.mem$v.bias<(mv-outlier) | k1.mem$v.bias>(mv+outlier) ] <- NA

mv <- mean(k1.mem$add.bias,na.rm=T)
outlier <- ol.sd*sd(k1.mem$add.bias,na.rm=T)
k1.mem$add.bias[k1.mem$add.bias<(mv-outlier) | k1.mem$add.bias>(mv+outlier) ] <- NA


k2.perc <- readRDS('katyal_2-Perception-vratio.Rds')

k2.perc <- data.frame(k2.perc$AD, k2.perc$CIT, k2.perc$SW, 
                      k2.perc$gender, k2.perc$age, 
                      k2.perc$recover[,7], k2.perc$recover[,6], 
                      k2.perc$recover[,5], k2.perc$recover[,8],
                      c(k2.perc$bestfnval), k2.perc$subj, 
                      k2.perc$recover[,1])
colnames(k2.perc) <- c('AD', 'CIT', 'SW', 'gender', 'age', 
                       'v.bias', 'mult.bias', 'add.bias', 'vratio', 
                       'fnval', 'subj', 'v')
k2.perc$task <- 'Perception'

mv <- mean(k2.perc$vratio,na.rm=T)
outlier <- ol.sd*sd(k2.perc$vratio,na.rm=T)
k2.perc$vratio[k2.perc$vratio<(mv-outlier) | k2.perc$vratio>(mv+outlier) ] <- NA

mv <- mean(k2.perc$v.bias,na.rm=T)
outlier <- ol.sd*sd(k2.perc$v.bias,na.rm=T)
k2.perc$vratio[k2.perc$v.bias<(mv-outlier) | k2.perc$v.bias>(mv+outlier) ] <- NA

mv <- mean(k2.perc$add.bias,na.rm=T)
outlier <- ol.sd*sd(k2.perc$add.bias,na.rm=T)
k2.perc$add.bias[k2.perc$add.bias<(mv-outlier) | k2.perc$add.bias>(mv+outlier) ] <- NA


k2.mem <- readRDS('katyal_2-Memory-vratio.Rds')
k2.mem <- data.frame(k2.mem$AD, k2.mem$CIT, k2.perc$SW, 
                     k2.mem$gender, k2.mem$age, 
                     k2.mem$recover[,7], k2.mem$recover[,6], 
                     k2.mem$recover[,5], k2.mem$recover[,8],
                     c(k2.mem$bestfnval), k2.mem$subj, 
                     k2.mem$recover[,1])
colnames(k2.mem) <- c('AD', 'CIT', 'SW', 'gender', 'age', 
                      'v.bias', 'mult.bias', 'add.bias', 'vratio', 
                      'fnval', 'subj', 'v')
k2.mem$task <- 'Memory'
mv <- mean(k2.mem$vratio,na.rm=T)
outlier <- ol.sd*sd(k2.mem$vratio,na.rm=T)
k2.mem$vratio[k2.mem$vratio<(mv-outlier) | k2.mem$vratio>(mv+outlier) ] <- NA

mv <- mean(k2.mem$v.bias,na.rm=T)
outlier <- ol.sd*sd(k2.mem$v.bias,na.rm=T)
k2.mem$v.bias[k2.mem$v.bias<(mv-outlier) | k2.mem$v.bias>(mv+outlier) ] <- NA

mv <- mean(k2.mem$add.bias,na.rm=T)
outlier <- ol.sd*sd(k2.mem$add.bias,na.rm=T)
k2.mem$add.bias[k2.mem$add.bias<(mv-outlier) | k2.mem$add.bias>(mv+outlier) ] <- NA


r1 <- readRDS('rouault_1-vratio.Rds')

r1 <- data.frame(r1$AD, r1$gender, r1$age, 
                 r1$recover[,7], r1$recover[,6], 
                 r1$recover[,5], r1$recover[,8],
                 r1$subj, 
                 r1$recover[,1])
colnames(r1) <- c('AD', 'gender', 'age', 
                  'v.bias', 'mult.bias', 'add.bias', 'vratio', 
                  'subj', 'v')
mv <- mean(r1$vratio,na.rm=T)
outlier <- ol.sd*sd(r1$vratio,na.rm=T)
r1$vratio[r1$vratio<(mv-outlier) | r1$vratio>(mv+outlier) ] <- NA

mv <- mean(r1$v.bias,na.rm=T)
outlier <- ol.sd*sd(r1$v.bias,na.rm=T)
r1$v.bias[r1$v.bias<(mv-outlier) | r1$v.bias>(mv+outlier) ] <- NA

mv <- mean(r1$add.bias,na.rm=T)
outlier <- ol.sd*sd(r1$add.bias,na.rm=T)
r1$add.bias[r1$add.bias<(mv-outlier) | r1$add.bias>(mv+outlier) ] <- NA


r2 <- readRDS('rouault_2-vratio.Rds')

r2 <- data.frame(r2$AD, r2$CIT, r2$SW, r2$gender, r2$age, 
                 r2$recover[,7], r2$recover[,6], 
                 r2$recover[,5], r2$recover[,8],
                 r2$subj, 
                 r2$recover[,1])
colnames(r2) <- c('AD', 'CIT', 'SW', 'gender', 'age', 
                  'v.bias', 'mult.bias', 'add.bias', 'vratio', 
                  'subj', 'v')

mv <- mean(r2$vratio,na.rm=T)
outlier <- ol.sd*sd(r2$vratio,na.rm=T)
r2$vratio[r2$vratio<(mv-outlier) | r2$vratio>(mv+outlier) ] <- NA

mv <- mean(r2$v.bias,na.rm=T)
outlier <- ol.sd*sd(r2$v.bias,na.rm=T)
r2$v.bias[r2$v.bias<(mv-outlier) | r2$v.bias>(mv+outlier) ] <- NA

mv <- mean(r2$add.bias,na.rm=T)
outlier <- ol.sd*sd(r2$add.bias,na.rm=T)
r2$add.bias[r2$add.bias<(mv-outlier) | r2$add.bias>(mv+outlier) ] <- NA


# test if vratio > 0
t.test(k1.perc$vratio)
t.test(k1.mem$vratio)
t.test(k2.perc$vratio)
t.test(k2.mem$vratio)
t.test(r1$vratio)
t.test(r2$vratio)

# test if vratio(mem) > vratio(perc)
t.test(k1.perc$vratio,k1.mem$vratio)
t.test(k2.perc$vratio,k2.mem$vratio, paired = T)

# test if vratio(mem) is correlated with vratio(perc)
e1.bothtasks.subj <- intersect((k1.perc %>% drop_na())$subj, 
                               (k1.mem %>% drop_na())$subj)

cor.test(x = filter(k1.perc, subj %in% e1.bothtasks.subj)$vratio,
         y = filter(k1.mem, subj %in% e1.bothtasks.subj)$vratio, 
         method = c( 'pearson'))

cor.test(x = k2.perc$vratio, y = k2.mem$vratio)

df.pm <- data.frame(
  perc = c(filter(k1.perc, subj %in% e1.bothtasks.subj)$vratio, k2.perc$vratio),
  mem = c(filter(k1.mem, subj %in% e1.bothtasks.subj)$vratio, k2.perc$vratio))
df.pm <- df.pm %>% drop_na()
library(BayesFactor)
correlationBF(y =df.pm$mem, x = df.pm$perc)

# test if A-bias < 0
t.test(k1.perc$add.bias)
t.test(k1.mem$add.bias)
t.test(k2.perc$add.bias)
t.test(k2.mem$add.bias)
t.test(r1$add.bias)
t.test(r2$add.bias)

t.test(k1.perc$add.bias,k1.mem$add.bias)
t.test(k2.perc$add.bias,k2.mem$add.bias, paired = T)

cor.test(x = filter(k1.perc, subj %in% e1.bothtasks.subj)$add.bias,
         y = filter(k1.mem, subj %in% e1.bothtasks.subj)$add.bias, 
         method = c( 'pearson'))

cor.test(x = k2.perc$add.bias, y = k2.mem$add.bias, method = c( 'pearson'))


t.test(k1.perc$v.bias)
t.test(k1.mem$v.bias)
t.test(k2.perc$v.bias)
t.test(k2.mem$v.bias)
t.test(r1$v.bias)
t.test(r2$v.bias)

pvbias <- c(k1.perc$v.bias, k2.perc$v.bias, r1$v.bias, r2$v.bias)
pvbias <- pvbias[!is.na(pvbias)]
ttestBF(pvbias)

mvbias <- c(k1.mem$v.bias, k2.mem$v.bias)
mvbias <- mvbias[!is.na(mvbias)]
ttestBF(mvbias)

