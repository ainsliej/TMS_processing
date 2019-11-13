#Script for analysis of CE data
#First do an analysis of the T1 data

setwd('~/../../Volumes/Ainslie_USB/VibData/PreProcessedData')
library(ggplot2)

#Sanity checks! 

#Load ALL data
ALL_CE <- read.delim('ALL_tableCE.txt', header=TRUE, sep=',' )
ALL_CE$muscle <- factor(ALL_CE$muscle)
ALL_CE$tDCS <- factor(ALL_CE$tDCS)
ALL_CE$ptp <- factor(ALL_CE$ptp)
ALL_CE$timept <- factor(ALL_CE$timept)

ALL_CE_FDI <- subset(ALL_CE,muscle=="FDI")
T1_CE <- subset(ALL_CE,timept=="T1")

ggplot(T1_CE, aes(x=muscle, y=DATA, group=interaction(tDCS,muscle)))+
  geom_violin(position="dodge",aes(fill=factor(tDCS))) +theme_minimal()+
  geom_jitter(inherit.aes=TRUE)+
  labs(title='Baseline CE, split by muscle',
       x='time point', y='MEP size in mV')
 
modelAOV_T1CE <- aov(DATA~muscle*tDCS+Error(ptp), data=T1_CE)
print("Check whether there is any CE differences at T1")
print(summary(modelAOV_T1CE))



#Now check the effect of tDCS on CE 

means_ALL_CE<- summarySEwithin(ALL_CE, measurevar="DATA", 
                               withinvars=c( "tDCS", "timept", "muscle"), idvar="ptp")

ggplot(means_ALL_CE, aes(x=timept, y=DATA, group=tDCS, color=tDCS, shape=muscle)) + 
  geom_errorbar(aes(ymin=DATA-se, ymax=DATA+se), width=.1, 
                position=position_dodge(0.05)) +
  geom_line() + geom_point()+
  scale_color_brewer(palette="Paired")+theme_minimal()+
  facet_wrap(facets= vars(muscle))+ ylim(0,1.5)+
  labs(title='The effect of tDCS on coritcospinal excitability, split by muscle',
       x='time point', y='MEP size in mV')


modelAOV_ALLCE<-aov(DATA~timept*tDCS*muscle+Error(ptp), data=ALL_CE)
print("Check whether there is any effect of tDCS on CE")
print(summary(modelAOV_ALLCE))

modelAOV_ALLCE_FDI<-aov(DATA~timept*tDCS+Error(ptp), data=ALL_CE_FDI)
print("Check whether there is any effect of tDCS on CE in FDI only")
print(summary(modelAOV_ALLCE_FDI))



#Analysis of the BL data
BL_CE <- read.delim('BL_tableCE.txt', header=TRUE, sep=',' )

BL_CE$muscle<-factor(BL_CE$muscle)
BL_CE$tDCS<- factor(BL_CE$tDCS)
BL_CE$ptp<- factor(BL_CE$ptp)
BL_CE$timept<- factor(BL_CE$timept)

BL_CE_FDI=subset(BL_CE,muscle=="FDI")

means_BL_CE<- summarySEwithin(BL_CE, measurevar="DATA", 
                               withinvars=c( "tDCS", "timept", "muscle"), idvar="ptp")

ggplot(means_BL_CE, aes(x=timept, y=DATA, group=tDCS, color=tDCS, shape=muscle)) + 
  geom_errorbar(aes(ymin=DATA-se, ymax=DATA+se), width=.1, 
                position=position_dodge(0.05)) +
  geom_line() + geom_point()+
  scale_color_brewer(palette="Paired")+theme_minimal()+
  facet_wrap(facets= vars(muscle))+ ylim(0.5,2.7)+
  labs(title='The effect of tDCS on coritcospinal excitability change, split by muscle',
       x='time point', y='MEP size in mV')

modelAOV_BLCE<-aov(DATA~timept*tDCS+Error(ptp), data=BL_CE_FDI)
print("Check whether there is any effect of baselined CE")
print(summary(modelAOV_BLCE))



