# TMS_processing
TMS processing scripts

Pipeline
1. cfs2mat.m for transforming signal script to a 3D matrix (time, frame, muscle)
2. pk2pk_[datatype].m for extracting average pk2pk mV values for each participant, for each timepoint, TMS condition, TDCS session, and muscle
3.[datatype]ANOVA.R for running the statistical analysis and creating graphs

N.B. to create bar/line graphs using ggplot2 you require the summarySE_within.R script.   
  You must run this to load in the functions which create the means, SD and SE for each measure type

Ignore SignalScripts folder
