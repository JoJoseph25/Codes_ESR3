# Clear plots
if(!is.null(dev.list())) dev.off()
# Clean workspace
rm(list=ls())
# Clear console
cat("\014") 

library(car)
library(ggplot2)
library(reshape2)
library(readxl)
library(gplots)
library(DescTools)
library(emmeans)
library(multcomp)
library(ggpubr)
library(ggsignif)

library(Rmisc) 
library(EnvStats)
library(ggplot2)
library(reshape2)
library(readxl)
library(psych)
library(dplyr)
library(permuco)

library(coin)
library(FSA)
library(afex)
library(rcompanion)
library(multcompView)


##### FETCH DATA ################
getwd()
setwd("G:/Socrates/Exp_Data/Henrique_data/Analysis")

data = read_excel("Exp_data.xlsx")
data = data.frame(data)


data$subject_id = factor(data$subject_id)

data$emo_num = ordered(data$emo_num, labels = c("Neutral","Relaxed","Sad","Happy","Fear"))
data$emo_num = ordered(data$emo_num, levels = c("Neutral","Relaxed","Sad","Happy","Fear"))
data$Gender = ordered(data$Gender, labels=c("Male","Female"))

###### OG ANOVA #################
col_name = "COP_clust_rel_std_COP_X"


model1=aov_ez(id = "subject_id", dv = col_name, data,  within = "emo_num")

One_Way_ANOVA = summary(model1)$univariate.tests
One_Way_Emo = One_Way_ANOVA[2,5]


model2=aov_ez(id = "subject_id", dv = col_name, data, between = "Gender", within = "emo_num")

Two_Way_ANOVA = summary(model2)$univariate.tests
Two_Way_Gen = Two_Way_ANOVA[2,5]
Two_Way_Emo = Two_Way_ANOVA[3,5]
Two_Way_Int = Two_Way_ANOVA[4,5]


#### PERMUTATION TEST ###########
set.seed(1979)

nreps = 10000

OW_E = numeric(nreps)
OW_E[1] = One_Way_Emo

TW_G = numeric(nreps)  
TW_E = numeric(nreps)
TW_GE = numeric(nreps)
TW_G[1] = Two_Way_Gen  
TW_E[1] = Two_Way_Emo
TW_GE[1] = Two_Way_Int


