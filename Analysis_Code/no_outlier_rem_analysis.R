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

create_dir= function(current_Dir,Dir_name){
  if (dir.exists(paste(current_Dir, Dir_name, sep = "/", collapse = "/"))) {
    cat("Dir_name exists in current_Dir and is a directory\n")
  } else {
    cat("Dir_name does not exist in current_Dir - creating\n")
    dir.create(file.path(current_Dir, Dir_name))
  }
  setwd(paste(current_Dir, Dir_name, sep = "/", collapse = "/"))
}

create_txt_file= function(current_Dir,file_name){
  if (file.exists(paste(current_Dir, file_name, sep = "/", collapse = "/"))) {
    cat("file_name exists in current_Dir and is a cleaned\n")
    file.remove(file.path(current_Dir, file_name))
    file.create(file.path(current_Dir, file_name))
  } else {
    cat("file_name does not exist in current_Dir - creating\n")
    file.create(file.path(current_Dir, file_name))
  }
}

getwd()
setwd("G:/Socrates/Exp_Data/Henrique_data/Analysis")

data= read_excel("Exp_data.xlsx")
data=data.frame(data)

data$subject_id= factor(data$subject_id)
data$emo_name= factor(data$emo_name)
data$emo_num= factor(data$emo_num, labels = c("Neutral","Relaxed","Sad","Happy","Fear"))
data$emo_num= ordered(data$emo_num, levels = c("Neutral","Relaxed","Sad","Happy","Fear"))
data$Gender= factor(data$Gender, labels=c("Male","Female"))

describeBy(data[1:15],data$emo_num)


# Create Main Folder
mainDir = getwd()
subDir = "No_Outlier_removed"
create_dir(mainDir,subDir)
main_wd = getwd()


all_columns=colnames(data) # all columns in excel file
# All column names seperated by different sections of analysis
# and stored in a nested list
all_col_names = list(list(all_columns[16:56]),list(all_columns[57:80]),list(all_columns[81:98]),list(all_columns[99:192]),list(all_columns[192:230]))
# All different sections of analysis stored in list
dir_names= list("Stand_COP","Stand_Angles","First_Step","Walking_Steps","Walking_Angles") 
dir_names= unlist(dir_names)
i=1

for (cur_dir in dir_names){  
  setwd(main_wd)
  
  cat(cur_dir,'\n') # Current directory name (section of analysis)
  
  cur_col_names = unlist(all_col_names[[i]]) # all columns for current analysis section
  
  # Create Section of analysis 
  mainDir = getwd()
  subDir = cur_dir
  create_dir(mainDir,subDir)
  section_wd = getwd()
  
  # Create Outlier Plot Folder 
  mainDir = getwd()
  subDir = "Outlier_Plots"
  create_dir(mainDir,subDir)
  outlier_wd = getwd()
  
  setwd("../") # Back to section folder
  
  # Create Main Effects Plot Folder 
  mainDir = getwd()
  subDir = "Main_Effect_Plots"
  create_dir(mainDir,subDir)
  main_eff_plots_wd = getwd()
  
  setwd("../") # Back to section folder
  
  # Create Gender Effects Plot Folder 
  mainDir = getwd()
  subDir = "Gender_Effect_Plots"
  create_dir(mainDir,subDir)
  gender_eff_plots_wd = getwd()
  
  setwd("../") # Back to section folder
  
  # Create Main Effects txt file
  create_txt_file(getwd(),"Main_effect.txt")
  # Create Pairwise Effects txt file
  create_txt_file(getwd(),"Pairwise.txt")
  # Create Gender Main Effects txt file
  create_txt_file(getwd(),"Gender_effect.txt")
  # Create Gender Pairwise Effects txt file
  create_txt_file(getwd(),"Gender_pairwise.txt")
  
  
  j=101
  
  for (col_name in cur_col_names) {
  
    # OUTLIER
    setwd(outlier_wd)
    jpeg(paste(c(j,"_",col_name,"_boxplot.jpg"),collapse = ""), width = 1366, height = 700)
    boxplot(data[[col_name]] ~ data$emo_num, main = 'Boxplot', xlab = 'Emotions', ylab = col_name)
    dev.off()
    setwd("../") # Back to section folder
    
    
    # One-Way ANOVA - Emotion
    model1= aov(data[[col_name]] ~ data[["emo_num"]] + Error(data[["subject_id"]]/data[["emo_num"]]))
    
    out1=capture.output(col_name)
    cat("\nVariable Name", out1, file="Main_effect.txt", sep="\n", append=TRUE)
    out2=capture.output(summary(model1)) # main effect
    cat(out2, file="Main_effect.txt", sep="\n", append=TRUE)
    
    # Pair-wise Main Effect
    emms = emmeans(model1, specs = ~ emo_num)
    pair_comp = summary(as.glht(pairs(emms)), test=adjusted("bonferroni")) # pair-wise comparison
    
    out3=capture.output(col_name)
    cat("\nVariable Name", out3, file="Pairwise.txt", sep="\n", append=TRUE)
    out4=capture.output(pair_comp)
    cat(out4, file="Pairwise.txt", sep="\n", append=TRUE)
    
    # Main Effect Plot
    data_summ = summarySEwithin(data, measurevar=col_name, withinvars="emo_num",idvar="subject_id", na.rm=TRUE, conf.interval=.95)
    names(data_summ)[3] = "var_value" 
    
    names(data)[names(data) == col_name] = "cur_var"


    setwd(main_eff_plots_wd)
    jpeg(paste(c(j,"_",col_name,".jpg"),collapse = ""), width = 1366, height = 700)
    
    
    anno=pair_comp$test$pvalues
    y.max = max(data_summ$var_value+data_summ$ci)
    p.value.y.coord = rep(y.max, length(anno))
    
    step.increase = (1:length(anno))*(y.max*.01)
    p.value.y.coord = p.value.y.coord + step.increase
    
    
    group_names = matrix(unlist(strsplit(names(anno)," - ")), ncol=2, byrow=TRUE)
    pvalues = data.frame(cbind(group_names[1:10,1],group_names[1:10,2],anno, p.value.y.coord))
    pvalues$anno=format.pval(anno, digits = 1)
    colnames(pvalues) = c("group1","group2","p.adj","y.coord")
    pvalues$group1= ordered(pvalues$group1, levels = c("Neutral","Relaxed","Sad","Happy"))
    pvalues$group1= factor(pvalues$group1, labels = c("Neutral","Relaxed","Sad","Happy"))
    pvalues$group2= ordered(pvalues$group2, levels = c("Relaxed","Sad","Happy","Fear"))
    pvalues$group2= factor(pvalues$group2, labels = c("Relaxed","Sad","Happy","Fear"))
    pvalues$y.coord = p.value.y.coord
    rownames(pvalues)= c(1:10)
    pvalues=pvalues[pvalues$p.adj <= 0.05,]
    
    
    
    plot(ggplot(data_summ, aes(x=emo_num, y=var_value, ymin=var_value-ci, ymax=var_value+ci, group=1)) + geom_line() + geom_errorbar(width=.1) + geom_point(shape=21, size=3, fill='white') + labs(y=col_name, x="Emotions"))
    
    p = ggplot(data_summ, aes(x=emo_num, y=var_value, ymin=var_value-ci, ymax=var_value+ci, group=1)) + geom_line() + geom_errorbar(width=.1) + geom_point(shape=21, size=3, fill='white') + labs(y=col_name, x="Emotions")
    p = p +geom_signif(inherit.aes = FALSE, data = pvalues, aes(xmin = group1, xmax = group2, annotations = p.adj, y_position = y.coord), manual= TRUE,  map_signif_level = TRUE, textsize=5, tip_length = 0.001)
    p
    
    names(data)[names(data) == "cur_var"] = col_name
    
    compare_means(cur_val ~ emo_num,  data = data_summ)
    p=plot(ggplot(data_summ, aes(x=emo_num, y=var_value, ymin=var_value-ci, ymax=var_value+ci, group=1)) + geom_line() + geom_errorbar(width=.1) + geom_point(shape=21, size=3, fill='white') + labs(y=col_name, x="Emotions"))
    
    
    p=p+stat_compare_means(aes(label=..p.adj..)) 
    p
    
    dev.off()
    setwd("../") # Back to section folder
    
    
    # Two-Way ANOVA - Gender * Emotion
    model2= aov(data[[col_name]] ~ data[["emo_num"]] * data[["Gender"]] + Error(data[["subject_id"]]/(data[["emo_num"]] * data[["Gender"]])))
    
    out5=capture.output(col_name)
    cat("\nVariable Name", out5, file="Gender_effect.txt", sep="\n", append=TRUE)
    out6=capture.output(summary(model2)) # main effect
    cat(out6, file="Gender_effect.txt", sep="\n", append=TRUE)
    
    # Gender Pairwise Effect
    emms2 = emmeans(model2, ~ emo_num | Gender)
    pair_comp2 = summary(as.glht(pairs(emms2)), test=adjusted("bonferroni")) # pair-wise comparison
    
    out7=capture.output(col_name)
    cat("\nVariable Name", out7, file="Gender_pairwise.txt", sep="\n", append=TRUE)
    out8=capture.output(pair_comp2)
    cat(out8, file="Gender_pairwise.txt", sep="\n", append=TRUE)
    
    # Gender Effect Plot
    data_summ2 = summarySEwithin(data, measurevar=col_name, betweenvars ="Gender", withinvars="emo_num",idvar="subject_id", na.rm=TRUE, conf.interval=.95)
    names(data_summ2)[4] = "var_value" 
    
    setwd(gender_eff_plots_wd)
    jpeg(paste(c(j,"_",col_name,".jpg"),collapse = ""), width = 1366, height = 700)
    plot(ggplot(data_summ2, aes(x=emo_num, y=var_value, ymin=var_value-ci, ymax=var_value+ci, group=Gender, color=Gender)) + geom_line() + geom_errorbar(width=.1) + geom_point(shape=21, size=3, fill='white') + labs(y=col_name, x="Emotions"))
    dev.off()
    setwd("../") # Back to section folder
    
    j= j+1
  }
  i= i+1
}
