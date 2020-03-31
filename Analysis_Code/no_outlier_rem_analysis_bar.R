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

# Create Main Folder
mainDir = getwd()
subDir = "Bar_Plots_ppr"
create_dir(mainDir,subDir)
main_wd = getwd()


all_columns = colnames(data) # all columns in excel file

# Manipulation Check Variables
manipul_check = all_columns[8:13]

# Create Manipulation txt file

for (cur_manipul_col in manipul_check){
  mani_model= aov(data[[col_name]] ~ data[["emo_num"]] * data[["Gender"]] + Error(data[["subject_id"]]/(data[["emo_num"]] * data[["Gender"]])))
  mani_emms = emmeans(mani_model, ~ Gender | emo_num)
  mani_pair_comp = summary(as.glht(pairs(mani_emms)), test=adjusted("bonferroni")) # pair-wise comparison
  mani_out1=capture.output(cur_manipul_col)
}


# Behavioural Measures
psy_col_names = all_columns[14:17]

for (cur_psy_col in psy_col_names){
  Male <- subset(data,  Gender == "Male", cur_psy_col,drop = TRUE)
  Female <- subset(data,  Gender == "Female", cur_psy_col,drop = TRUE)
  
  homogenous_test=var.test(Male, Female)
  if (homogenous_test$p.value>0.05){
    t_out=t.test(Male,Female,var.equal = TRUE)
  } else {
    t_out=t.test(Male,Female,var.equal = FALSE)
  }  
}  


# All column names for biomechanic parameters seperated by different sections of analysis
# and stored in a nested list
all_col_names = list(list(all_columns[18:58]),list(all_columns[59:82]),list(all_columns[83:100]),list(all_columns[101:194]),list(all_columns[195:232]))
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
  
    
    # One-Way ANOVA - Emotion
    model1= aov(data[[col_name]] ~ data[["emo_num"]] + Error(data[["subject_id"]]/data[["emo_num"]]))
    

    # Pair-wise Main Effect
    emms = emmeans(model1, ~ emo_num)
    pair_comp = summary(as.glht(pairs(emms)), test=adjusted("bonferroni")) # pair-wise comparison
    pair_comp

    # Main Effect Plot
    data_summ = summarySEwithin(data, measurevar=col_name, withinvars="emo_num",idvar="subject_id", na.rm=TRUE, conf.interval=.95)
    names(data_summ)[3] = "var_value" 
  
    pair_comp_df=pair_comp$test$pvalues
    y.max = max(data_summ$var_value+data_summ$ci)
    p.value.y.coord = rep(y.max, length(pair_comp_df))
    step.increase = (1:length(pair_comp_df))*(y.max*.01)
    p.value.y.coord = p.value.y.coord + step.increase
    
    group_names = matrix(unlist(strsplit(names(pair_comp_df)," - ")), ncol=2, byrow=TRUE)
    pvalues = data.frame(cbind(group_names[1:10,1],group_names[1:10,2],pair_comp_df, p.value.y.coord))
    pvalues$pair_comp_df=format.pval(pair_comp_df, digits = 2)
    colnames(pvalues) = c("group1","group2","p.adj","y.coord")
    pvalues$group1= ordered(pvalues$group1, levels = c("Neutral","Relaxed","Sad","Happy"))
    pvalues$group1= factor(pvalues$group1, labels = c("Neutral","Relaxed","Sad","Happy"))
    pvalues$group2= ordered(pvalues$group2, levels = c("Relaxed","Sad","Happy","Fear"))
    pvalues$group2= factor(pvalues$group2, labels = c("Relaxed","Sad","Happy","Fear"))
    pvalues$y.coord = p.value.y.coord
    rownames(pvalues)= c(1:10)
    pvalues=pvalues[pvalues$p.adj <= 0.05,]
    
    p = ggplot(data_summ, aes(x=emo_num, y=var_value)) + geom_bar(position=position_dodge(), stat="identity",colour="black",size=.3) + geom_errorbar(aes(ymin=var_value-se, ymax=var_value+se), size=.3, width=.2, position=position_dodge(.9)) + xlab("Emotions") + ylab(col_name) + theme_bw()
    if (dim(pvalues)[1] != 0){
      p = p +geom_signif(inherit.aes = FALSE, data = pvalues, aes(xmin = group1, xmax = group2, annotations = "*", y_position = y.coord), manual= TRUE,  map_signif_level = TRUE, textsize=5, tip_length = 0.001)
    }
    p
   
    
    # Two-Way ANOVA - Gender * Emotion
    model2= aov(data[[col_name]] ~ data[["emo_num"]] * data[["Gender"]] + Error(data[["subject_id"]]/(data[["emo_num"]] * data[["Gender"]])))
    
    # Gender Pairwise Effect
    emms2 = emmeans(model2, ~ emo_num | Gender)
    pair_comp2 = summary(as.glht(pairs(emms2)), test=adjusted("bonferroni")) # pair-wise comparison
    
    # Gender Effect Plot
    data_summ2 = summarySEwithin(data, measurevar=col_name, betweenvars ="Gender", withinvars="emo_num",idvar="subject_id", na.rm=TRUE, conf.interval=.95)
    names(data_summ2)[4] = "var_value" 
    
    pair_comp_df=cbind(pair_comp2[[1]]$test$pvalues,pair_comp2[[2]]$test$pvalues)
    y.max = max(data_summ2$var_value+data_summ2$ci)
    p.value.y.coord = rep(y.max, length(pair_comp_df))

    step.increase = (1:length(pair_comp_df))*(y.max*.01)
    p.value.y.coord = p.value.y.coord + step.increase

    group_names = matrix(unlist(strsplit(names(pair_comp2[[1]]$test$pvalues)," - ")), ncol=2, byrow=TRUE)
    pvalues=data.frame(pair_comp_df)
    pvalues$group1 = group_names[1:10,1]
    pvalues$group2 = group_names[1:10,2]
    pvalues$group1= ordered(pvalues$group1, levels = c("Neutral","Relaxed","Sad","Happy"))
    pvalues$group2= ordered(pvalues$group2, levels = c("Relaxed","Sad","Happy","Fear"))
    pvalues$male_y.coord = p.value.y.coord[11:20]
    pvalues$female_y.coord = p.value.y.coord[1:10]
    rownames(pvalues)= c(1:10)
    colnames(pvalues) = c("male_p.adj","female_p.adj","group1","group2","male_y.coord","female_y.coord")
    pvalues$male_p.adj=format.pval(pvalues$male_p.adj, digits = 2)
    pvalues$female_p.adj=format.pval(pvalues$female_p.adj, digits = 2)
    pvalues=pvalues[pvalues$male_p.adj <= 0.05 | pvalues$female_p.adj <= 0.05 ,]


    p = ggplot(data_summ2, aes(x=emo_num, y=var_value, fill=Gender)) + geom_bar(position=position_dodge(), stat="identity",colour="black",size=.3) + geom_errorbar(aes(ymin=var_value-se, ymax=var_value+se), size=.3, width=.2, position=position_dodge(.9)) + xlab("Emotions") + ylab(col_name) + scale_fill_hue(name="Gender",labels=c("Male", "Female")) + theme_bw() + scale_fill_grey()
    if (dim(pvalues)[1] != 0){
      if (pvalues$female_p.adj <= 0.05){
        p = p +geom_signif(inherit.aes = FALSE, data = pvalues, aes(xmin = group1, xmax = group2, annotations = "*", y_position = female_y.coord), manual= TRUE,  map_signif_level = TRUE, textsize=5, tip_length = 0.001)
      }
      if (pvalues$male_p.adj <= 0.05){
      p = p +geom_signif(inherit.aes = FALSE, data = pvalues, aes(xmin = group1, xmax = group2, annotations = "*", y_position = male_y.coord), manual= TRUE,  map_signif_level = TRUE, textsize=5, tip_length = 0.001)
      }
    }
    p
    
    # Gender Pairwise Effect
    emms3 = emmeans(model2, ~ Gender | emo_num)
    pair_comp3 = summary(as.glht(pairs(emms3)), test=adjusted("bonferroni")) # pair-wise comparison
    
    # Gender Effect Plot
    data_summ3 = summarySEwithin(data, measurevar=col_name, betweenvars ="Gender", withinvars="emo_num",idvar="subject_id", na.rm=TRUE, conf.interval=.95)
    names(data_summ3)[4] = "var_value" 
    
    pair_comp_df=rbind(pair_comp3[[1]]$test$pvalues,pair_comp3[[2]]$test$pvalues,pair_comp3[[3]]$test$pvalues,pair_comp3[[4]]$test$pvalues,pair_comp3[[5]]$test$pvalues)
    y.max = max(data_summ3$var_value+data_summ3$ci)
    p.value.y.coord = rep(y.max, length(pair_comp_df))
    
    step.increase = (1:length(pair_comp_df))*(y.max*.01)
    p.value.y.coord = p.value.y.coord + step.increase
    
    group_names = c("Neutral","Relaxed","Sad","Happy","Fear")
    pvalues=data.frame(pair_comp_df)
    pvalues$group = group_names
    pvalues$y.coord = p.value.y.coord
    pvalues$male_cord = c(0.75,1.75,2.75,3.75,4.75)
    pvalues$female_cord = c(1.25,2.25,3.25,4.25,5.25)
    colnames(pvalues) = c("p.adj","group","y.coord","male_cord","female_cord")
    pvalues=pvalues[pvalues$p.adj <= 0.05 ,]
    
    p = ggplot(data_summ3, aes(x=emo_num, y=var_value, fill=Gender)) + geom_bar(position=position_dodge(), stat="identity",colour="black",size=.3) + geom_errorbar(aes(ymin=var_value-se, ymax=var_value+se), size=.3, width=.2, position=position_dodge(.9)) + xlab("Emotions") + ylab(col_name) + scale_fill_hue(name="Gender",labels=c("Male", "Female")) + theme_bw() + scale_fill_grey()
    if (dim(pvalues)[1] != 0){
      p = p + geom_signif(inherit.aes = FALSE, data = pvalues, aes(xmin = male_cord, xmax = female_cord, annotations = "*", y_position = y.coord), manual= TRUE,  map_signif_level = TRUE, textsize=5, tip_length = 0.0005)
    }
    p
    
   # pair_comp
   # pair_comp2
   # pair_comp3
    summary(model1)
   
    
    j= j+1
  }
  i= i+1
}
