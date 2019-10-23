#Script for analysis of VIB data

#Sanity checks! Single Pulse

setwd('~/../../Volumes/Ainslie_USB/VibData/PreProcessedData')
library(ggplot2)

#Analysis of full SP data set with all timepoints

ALL_SP_VIB <- read.delim('ALL_SP_tableVIB.txt', header=TRUE, sep=',' )

ALL_SP_VIB$muscle<-factor(ALL_SP_VIB$muscle)
ALL_SP_VIB$tDCS<- factor(ALL_SP_VIB$tDCS)
ALL_SP_VIB$ptp<- factor(ALL_SP_VIB$ptp)
ALL_SP_VIB$vibCond<- factor(ALL_SP_VIB$vibCond)
ALL_SP_VIB$timept<- factor(ALL_SP_VIB$timept)

#Check whether sp MEPs change across the session
ALL_SP_noVIB_FDI=subset(ALL_SP_VIB,vibCond=="vibNO" & muscle=="FDI")
modelAOV_ALL_SP_noVIB_FDI<-aov(DATA~timept*tDCS+Error(ptp), data=ALL_SP_noVIB_FDI)
print("A check that sp MEPs do not change across the session")
print(summary(modelAOV_ALL_SP_noVIB_FDI))

means_ALL_SP_noVIB_FDI<- summarySEwithin(ALL_SP_noVIB_FDI, measurevar="DATA", 
                               withinvars=c( "tDCS", "timept"), idvar="ptp")

ggplot(means_ALL_SP_noVIB_FDI, aes(x=timept, y=DATA, group=tDCS, color=tDCS)) + 
  geom_errorbar(aes(ymin=DATA-se, ymax=DATA+se), width=.1, 
                position=position_dodge(0.05)) +
  geom_line() + geom_point()+
  scale_color_brewer(palette="Paired")+theme_minimal()+ 
  labs(title='The effect of tDCS on coritcospinal excitability, split by muscle',
       x='time point', y='MEP size in mV')

#Check whether we get vibration effect at the T1 timepoint

T1_SP_VIB<-subset(ALL_SP_VIB,timept=="T1")
modelAOV_T1_SP_VIB<-aov(DATA~muscle*vibCond*tDCS+Error(ptp), data=T1_SP_VIB)
print("A check that we are getting the VIB effect at baseline")
print(summary(modelAOV_T1_SP_VIB))

means_ALL_SP_VIB<- summarySEwithin(ALL_SP_VIB, measurevar="DATA", 
                                   withinvars=c( "tDCS", "timept", "muscle", "vibCond"), idvar="ptp")
ggplot(subset(means_ALL_SP_VIB,timept==("T1")), aes(x=muscle, y=DATA, fill=muscle)) + 
  geom_errorbar(aes(ymin=DATA-se, ymax=DATA+se)) +
  geom_bar(stat='identity')+ ylim(0,2)+
  scale_color_brewer(palette="Paired")+theme_minimal()+
  facet_grid(rows=vars(vibCond),cols=vars(tDCS)) + 
  labs(title='MEPs at baseline, split by stim condition and vib condition',
       x='time point', y='MEP size in mV')

#Analysis of VIBeffect SP data at the T1 timepoint

ALL_SP_VIBeffect <- read.delim('ALL_SP_tableVIBeffect.txt', header=TRUE, sep=',' )

ALL_SP_VIBeffect$muscle<-factor(ALL_SP_VIBeffect$muscle)
ALL_SP_VIBeffect$tDCS<- factor(ALL_SP_VIBeffect$tDCS)
ALL_SP_VIBeffect$ptp<- factor(ALL_SP_VIBeffect$ptp)
ALL_SP_VIBeffect$vibCond<- factor(ALL_SP_VIBeffect$vibCond)
ALL_SP_VIBeffect$timept<- factor(ALL_SP_VIBeffect$timept)

means_ALL_SP_VIBeffect<- summarySEwithin(ALL_SP_VIBeffect, measurevar="DATA", 
                                         withinvars=c( "tDCS", "timept", "muscle", "vibCond"), idvar="ptp")

ggplot(subset(means_ALL_SP_VIBeffect,timept==("T1")), aes(x=muscle, y=DATA, fill=muscle)) + 
  geom_errorbar(aes(ymin=DATA-se, ymax=DATA+se)) +
  geom_bar(stat='identity')+ ylim(0,2)+
  scale_color_brewer(palette="Paired")+theme_minimal()+
  facet_grid(rows=vars(vibCond),cols=vars(tDCS)) + 
  labs(title='Vib effect on MEPs at baseline, split by stim condition',
       x='Muscle', y='MEP size in mV')




#Sanity Checks! SICI measures 

#The effect of SICI measures at T1

ALL_SICI_VIB <- read.delim('ALL_SICI_tableVIB.txt', header=TRUE, sep=',' )

ALL_SICI_VIB$muscle<-factor(ALL_SICI_VIB$muscle)
ALL_SICI_VIB$tDCS<- factor(ALL_SICI_VIB$tDCS)
ALL_SICI_VIB$ptp<- factor(ALL_SICI_VIB$ptp)
ALL_SICI_VIB$vibCond<- factor(ALL_SICI_VIB$vibCond)
ALL_SICI_VIB$timept<- factor(ALL_SICI_VIB$timept)

T1_SP_VIB<-subset(ALL_SP_VIB,timept=="T1")
modelAOV_T1_SP_VIB<-aov(DATA~muscle*vibCond*tDCS+Error(ptp), data=T1_SP_VIB)
print("A check that we are getting the VIB on SICI at baseline")
print(summary(modelAOV_T1_SP_VIB))

means_ALL_SICI_VIB<- summarySEwithin(ALL_SICI_VIB, measurevar="DATA", 
                                           withinvars=c( "tDCS", "timept", "vibCond","muscle"), idvar="ptp")

ggplot(subset(means_ALL_SICI_VIB,timept==("T1")), aes(x=muscle, y=DATA, fill=muscle)) + 
  geom_errorbar(aes(ymin=DATA-se, ymax=DATA+se)) +
  geom_bar(stat='identity')+ 
  scale_color_brewer(palette="Paired")+theme_minimal()+
  facet_grid(rows=vars(vibCond),cols=vars(tDCS)) + 
  labs(title='SICI baseline, split by stim condition and vib condition',
       x='time point', y='ppMEP/spMEP')

#Visualisng the vibration effect of SICI measures at T1

ALL_SP_VIBeffect <- read.delim('ALL_SP_tableVIBeffect.txt', header=TRUE, sep=',' )

ALL_SP_VIBeffect$muscle<-factor(ALL_SP_VIBeffect$muscle)
ALL_SP_VIBeffect$tDCS<- factor(ALL_SP_VIBeffect$tDCS)
ALL_SP_VIBeffect$ptp<- factor(ALL_SP_VIBeffect$ptp)
ALL_SP_VIBeffect$vibCond<- factor(ALL_SP_VIBeffect$vibCond)
ALL_SP_VIBeffect$timept<- factor(ALL_SP_VIBeffect$timept)

means_ALL_SP_VIBeffect<- summarySEwithin(ALL_SP_VIBeffect, measurevar="DATA", 
                                         withinvars=c( "tDCS", "timept", "muscle", "vibCond"), idvar="ptp")

ggplot(subset(means_ALL_SICI_VIB,timept==("T1")), aes(x=muscle, y=DATA, fill=muscle)) + 
  geom_errorbar(aes(ymin=DATA-se, ymax=DATA+se)) +
  geom_bar(stat='identity')+ 
  scale_color_brewer(palette="Paired")+theme_minimal()+
  facet_grid(rows=vars(vibCond),cols=vars(tDCS)) + 
  labs(title='SICI vib effect at baseline, split by stim condition',
       x='time point', y='vib effect on ppMEP/spMEP')






#Testing whether there is any effect of tDCS on SICI

ALL_SICI_noVIB=subset(ALL_SICI_VIB,vibCond=="vibNO" )

means_ALL_SICI_noVIB<- summarySEwithin(ALL_SICI_noVIB, measurevar="DATA", 
                                       withinvars=c( "tDCS", "timept","muscle"), idvar="ptp")

ggplot(means_ALL_SICI_noVIB, aes(x=timept, y=DATA, group=tDCS, color=tDCS, shape=muscle)) + 
  geom_errorbar(aes(ymin=DATA-se, ymax=DATA+se), width=.1, 
                position=position_dodge(0.05)) +
  geom_line() + geom_point()+
  scale_color_brewer(palette="Paired")+theme_minimal()+
  facet_wrap(facets= vars(muscle))+ ylim(0,1.5)+
  labs(title='The effect of tDCS on SICI, split by muscle',
       x='time point', y='ppMEP/spMEP')

modelAOV_ALL_SICI_noVIB<-aov(DATA~timept*tDCS*muscle+Error(ptp), data=ALL_SICI_noVIB)
print("Testing whether tDCS influences SICI")
print(summary(modelAOV_ALL_SICI_noVIB))







#Analysis of SP MEPs across the session (full model)
means_ALL_SP_VIB<- summarySEwithin(ALL_SP_VIB, measurevar="DATA", 
                                   withinvars=c( "tDCS", "timept", "muscle", "vibCond"), idvar="ptp")
ggplot(means_ALL_SP_VIB, aes(x=muscle, y=DATA, fill=muscle)) + 
  geom_errorbar(aes(ymin=DATA-se, ymax=DATA+se)) +
  geom_bar(stat='identity')+ ylim(0,2)+
  scale_color_brewer(palette="Paired")+theme_minimal()+
  facet_grid(rows=vars(vibCond),cols=vars(tDCS,timept)) + 
  labs(title='The effect of tDCS on MEPs',
       x='time point', y='MEP size in mV')

modelAOV_ALL_VIB_SP<-aov(DATA~timept*tDCS*vibCond*muscle+Error(ptp), data=ALL_SP_VIB)
print("Does tDCS influence the vib effect on MEPs")
print(summary(modelAOV_ALL_VIB_SP))



#Analysis of VIBeffect SP data set with all timepoints

means_ALL_SP_VIBeffect<- summarySEwithin(ALL_SP_VIBeffect, measurevar="DATA", 
                                         withinvars=c( "tDCS", "timept", "muscle", "vibCond"), idvar="ptp")

ggplot(means_ALL_SP_VIBeffect, aes(x=muscle, y=DATA, fill=muscle)) + 
  geom_errorbar(aes(ymin=DATA-se, ymax=DATA+se)) +
  geom_bar(stat='identity')+ ylim(0,2)+
  scale_color_brewer(palette="Paired")+theme_minimal() +
  facet_grid(rows=vars(vibCond),cols=vars(tDCS,timept)) + 
  labs(title='The effect of tDCS on vibration effect on SP MEP',
       x='time point', y='VIB MEP size/ noVIB MEP size')

modelAOV_ALL_VIBeffect_SP<-aov(DATA~timept*tDCS*vibCond*muscle+Error(ptp), data=ALL_SP_VIBeffect)
print(summary(modelAOV_ALL_VIBeffect_SP))


#Analysis of SPvibEffect data set baselined

BL_SP_VIBeffect <- read.delim('BL_SP_tableVIBeffect.txt', header=TRUE, sep=',' )

BL_SP_VIBeffect$muscle<-factor(BL_SP_VIBeffect$muscle)
BL_SP_VIBeffect$tDCS<- factor(BL_SP_VIBeffect$tDCS)
BL_SP_VIBeffect$ptp<- factor(BL_SP_VIBeffect$ptp)
BL_SP_VIBeffect$VIBeffectCond<- factor(BL_SP_VIBeffect$VIBeffectCond)
BL_SP_VIBeffect$timept<- factor(BL_SP_VIBeffect$timept)

means_BL_SP_VIBeffect<- summarySEwithin(BL_SP_VIBeffect, measurevar="DATA", 
                                        withinvars=c( "tDCS", "timept", "muscle", "vibCond"), idvar="ptp")

ggplot(means_BL_SP_VIBeffect, aes(x=muscle, y=DATA, fill=muscle)) + 
  geom_errorbar(aes(ymin=DATA-se, ymax=DATA+se)) +
  geom_bar(stat='identity')+ 
  scale_color_brewer(palette="Paired")+theme_minimal()+
  facet_grid(rows=vars(vibCond),cols=vars(tDCS,timept)) + 
  labs(title='The effect of tDCS on SP effect',
       x='time point', y='MEP size in mV')

modelAOV_BL_VIBeffect_SP<-aov(DATA~timept*tDCS*vibCond*muscle+Error(ptp), data=BL_SP_VIBeffect)
print("Does tDCS influence the vib effect once we baseline?")
print(summary(modelAOV_BL_VIBeffect_SP))





#Analysis of full SICI data set with all timepoints

modelAOV_ALL_VIB_SICI<-aov(DATA~timept*tDCS*vibCond*muscle+Error(ptp), data=ALL_SICI_VIB)
print("Testing for an effect of tDCS on the VIB SICI effect ")
print(summary(modelAOV_ALL_VIB_SICI))

means_ALL_SICI_VIB<- summarySEwithin(ALL_SICI_VIB, measurevar="DATA", 
                                     withinvars=c( "tDCS", "timept", "muscle", "vibCond"), idvar="ptp")

ggplot(means_ALL_SICI_VIB, aes(x=muscle, y=DATA, fill=muscle)) + 
  geom_errorbar(aes(ymin=DATA-se, ymax=DATA+se)) +
  geom_bar(stat='identity')+ ylim(0,1)+
  scale_color_brewer(palette="Paired")+theme_minimal()+
  facet_grid(rows=vars(vibCond),cols=vars(tDCS,timept)) + 
  labs(title='The effect of tDCS on SICI vibration',
       x='time point', y='MEP size in mV')

modelAOV_ALL_VIB_SICI<-aov(DATA~timept*tDCS*vibCond*muscle+Error(ptp), data=ALL_SICI_VIB)
print(summary(modelAOV_ALL_VIB_SICI))

#Analysis of SICIvibEffect data set with all timepoints

means_ALL_SICI_VIBeffect<- summarySEwithin(ALL_SICI_VIBeffect, measurevar="DATA", 
                                           withinvars=c( "tDCS", "timept", "muscle", "vibCond"), idvar="ptp")

ggplot(means_ALL_SICI_VIBeffect, aes(x=muscle, y=DATA, fill=muscle)) + 
  geom_errorbar(aes(ymin=DATA-se, ymax=DATA+se)) +
  geom_bar(stat='identity')+ 
  scale_color_brewer(palette="Paired")+theme_minimal()+
  facet_grid(rows=vars(vibCond),cols=vars(tDCS,timept)) + 
  labs(title='The effect of tDCS on vibration effect on SICI MEP',
       x='time point', y='VIB SICI size/ noVIB SICI size')

modelAOV_ALL_VIBeffect_SICI<-aov(DATA~timept*tDCS*vibCond*muscle+Error(ptp), data=ALL_SICI_VIBeffect)
print(summary(modelAOV_ALL_VIBeffect_SICI))


#Analysis of SICIvibEffect data set baselined

BL_SICI_VIBeffect <- read.delim('BL_SICI_tableVIBeffect.txt', header=TRUE, sep=',' )

BL_SICI_VIBeffect$muscle<-factor(BL_SICI_VIBeffect$muscle)
BL_SICI_VIBeffect$tDCS<- factor(BL_SICI_VIBeffect$tDCS)
BL_SICI_VIBeffect$ptp<- factor(BL_SICI_VIBeffect$ptp)
BL_SICI_VIBeffect$vibCond<- factor(BL_SICI_VIBeffect$vibCond)
BL_SICI_VIBeffect$timept<- factor(BL_SICI_VIBeffect$timept)

means_BL_SICI_VIBeffect<- summarySEwithin(BL_SICI_VIBeffect, measurevar="DATA", 
                                          withinvars=c( "tDCS", "timept", "muscle", "vibCond"), idvar="ptp")

ggplot(means_BL_SICI_VIBeffect, aes(x=muscle, y=DATA, fill=muscle)) + 
  geom_errorbar(aes(ymin=DATA-se, ymax=DATA+se)) +
  geom_bar(stat='identity')+ 
  scale_color_brewer(palette="Paired")+theme_minimal()+
  facet_grid(rows=vars(vibCond),cols=vars(tDCS,timept)) + 
  labs(title='The effect of tDCS on vibration effect on SICI MEP',
       x='time point', y='VIB SICI size/ noVIB SICI size')

modelAOV_BL_VIBeffect_SICI<-aov(DATA~timept*tDCS*vibCond*muscle+Error(ptp), data=BL_SICI_VIBeffect)
print("Does tDCS influence the SICI effect once we have baselined")
print(summary(modelAOV_BL_VIBeffect_SICI))

