#Script for analysis of SAI data
#First do an analysis of the T1 data

setwd('~/../../Volumes/Ainslie_USB/VibData/PreProcessedData')
library(ggplot2)

#Sanity checks! 

#Load all data and check for baseline differences

ALL_SAI <- read.delim('ALL_tableSAI.txt', header=TRUE, sep=',' )
ALL_SAI$muscle<-factor(ALL_SAI$muscle)
ALL_SAI$tDCS<- factor(ALL_SAI$tDCS)
ALL_SAI$ptp<- factor(ALL_SAI$ptp)
ALL_SAI$timept<- factor(ALL_SAI$timept)

T1_SAI<- subset(ALL_SAI,timept=="T1")
ALL_SAI_FDI<-subset(ALL_CE,muscle=="FDI")

ggplot(T1_CE, aes(x=muscle, y=DATA, group=interaction(tDCS,muscle)))+
  geom_violin(position="dodge",aes(fill=factor(tDCS))) +theme_minimal()+
  geom_jitter(inherit.aes=TRUE)

modelAOV_T1SAI<-aov(DATA~muscle*tDCS+Error(ptp), data=T1_SAI)
print("Check whether there is any SAI differences at T1")
print(summary(modelAOV_T1SAI))


#All data just the single pulses 

ALL_SAIwoMNS <- read.delim('ALL_tableSAIwoMNS.txt', header=TRUE, sep=',' )
ALL_SAIwoMNS$muscle<-factor(ALL_SAIwoMNS$muscle)
ALL_SAIwoMNS$tDCS<- factor(ALL_SAIwoMNS$tDCS)
ALL_SAIwoMNS$ptp<- factor(ALL_SAIwoMNS$ptp)
ALL_SAIwoMNS$timept<- factor(ALL_SAIwoMNS$timept)

means_ALL_SAIwoMNS<- summarySEwithin(ALL_SAIwoMNS, measurevar="DATA", 
                                     withinvars=c( "tDCS", "timept", "muscle"), idvar="ptp")

ggplot(means_ALL_SAIwoMNS, aes(x=timept, y=DATA, group=tDCS, color=tDCS, shape=muscle)) + 
  geom_errorbar(aes(ymin=DATA-se, ymax=DATA+se), width=.1, 
                position=position_dodge(0.05)) +
  geom_line() + geom_point()+
  scale_color_brewer(palette="Paired")+theme_minimal()+
  facet_wrap(facets= vars(muscle))+ ylim(0,1.5)+
  labs(title='The effect of tDCS on single pulses in the SAI protocol, split by muscle',
       x='time point', y='MEP in mV')

ALL_SAIwoMNS_FDI=subset(ALL_SAIwoMNS,muscle=="FDI")
modelAOV_ALLSAIwoMNS<-aov(DATA~timept*tDCS+Error(ptp), data=ALL_SAIwoMNS_FDI)
print("Check whether there is the SAI single pulses change over time")
print(summary(modelAOV_ALLSAIwoMNS))






# Analysis of all the data 

means_ALL_SAI<- summarySEwithin(ALL_SAI, measurevar="DATA", 
                                withinvars=c( "tDCS", "timept", "muscle"), idvar="ptp")

ggplot(means_ALL_SAI, aes(x=timept, y=DATA, group=tDCS, color=tDCS, shape=muscle)) + 
  geom_errorbar(aes(ymin=DATA-se, ymax=DATA+se), width=.1, 
                position=position_dodge(0.05)) +
  geom_line() + geom_point()+
  scale_color_brewer(palette="Paired")+theme_minimal()+
  facet_wrap(facets= vars(muscle))+ ylim(0,1.5)+
  labs(title='The effect of tDCS on the SAI effect, split by muscle',
       x='time point', y='MEP wMNS/ MEP woMNS')

ALL_SAI_FDI=subset(ALL_SAI,muscle=="FDI")
modelAOV_ALLSAI<-aov(DATA~timept*tDCS*muscle+Error(ptp), data=ALL_SAI)
print("Check whether there is any chance in SAI with tDCS")
print(summary(modelAOV_ALLSAI))

modelAOV_ALLSAI<-aov(DATA~timept*tDCS+Error(ptp), data=ALL_SAI_FDI)
print("Check whether there is any chance in SAI with tDCS in just FDI")
print(summary(modelAOV_ALLSAI))



# Analysis of the baselined data
BL_SAI <- read.delim('BL_tableSAI.txt', header=TRUE, sep=',' )

BL_SAI$muscle<-factor(BL_SAI$muscle)
BL_SAI$tDCS<- factor(BL_SAI$tDCS)
BL_SAI$ptp<- factor(BL_SAI$ptp)
BL_SAI$timept<- factor(BL_SAI$timept)
BL_SAI_FDI=subset(BL_SAI,muscle=="FDI")

means_BL_SAI<- summarySEwithin(BL_SAI, measurevar="DATA", 
                                withinvars=c( "tDCS", "timept", "muscle"), idvar="ptp")

ggplot(means_BL_SAI, aes(x=timept, y=DATA, group=tDCS, color=tDCS, shape=muscle)) + 
  geom_errorbar(aes(ymin=DATA-se, ymax=DATA+se), width=.1, 
                position=position_dodge(0.05)) +
  geom_line() + geom_point()+
  scale_color_brewer(palette="Paired")+theme_minimal()+
  facet_wrap(facets= vars(muscle))+ ylim(0,2)+
  labs(title='The effect of tDCS on the change in SAI effect, split by muscle',
       x='time point', y='Change in MEP wMNS/ MEP woMNS from BL')

modelAOV_BLSAI<-aov(DATA~timept*tDCS+Error(ptp), data=BL_SAI_FDI)
print("Check whether there is any chance in baselined SAI with tDCS in just FDI")
print(summary(modelAOV_BLSAI))




