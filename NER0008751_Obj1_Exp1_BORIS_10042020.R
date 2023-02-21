#DavidHCollins
#Created: 6 April 2020

#Aim

#The aim of this script is to analyse data from our main Bombus experiment. This is a draft script which I am using in order to work out how to analyse the movement data generated from monitoring the BORIS analysis.

########################################################################################################################

#First you must always reset R
rm(list=ls())

#make these packages (and associated functions) available to use in the script

library(readr) #read in files
library(dplyr) #use dplyr functions for data manipulation
library(ggplot2) #graphics package

#Might be usefult to draw a survival curve for the longevity data later on, install the ggplot package fortify which has suvival curve functions
#(if on a new computer use:
#install.packages('ggfortify')

#these packages all help for modelling
library(ggfortify)
library(survival)
library("survminer")
library(bbmle) #produces AIC values which allows you to compare multiple models
library(lme4) ## for glmer
library(grid)
library(gridExtra)
library(glmmTMB) ##for zero inflated glmms

#here are some useful shortcuts to get to know
#ctrl+a highlights the whole script
#ctrl+enter runs the highlighted script or just the selected line of code if nothing is highlighted
#ctrl+1 moves cursor to the console
#ctrl+2 moves it back to the script editor


#set file directory

setwd("Documents/UEA Work documents/Post Doc_LHT/Objective 1 Experiment 1/Analysis")

#read in data for BORIS analysis
boris_data <- read_csv("BORIS_video.csv")
View(boris_data)
df1<-boris_data

#Turn treat into a factor and reorder so R comes first
df1$treat <- as.factor(df1$treat)
is.factor(df1$treat)
df1$treat <- factor(df1$treat, levels = c("R", "C"))

#Round the time columns to integers
df1$queen_move <- round(df1$queen_move)
df1$queen_stationary <- round(df1$queen_stationary)

#Make a new column with the total amount of time being filmed (i.e. stationary + moving)
df1 <- df1 %>% 
  mutate(total = queen_stationary + queen_move)

#make time a factor
df1$time<-as.factor(df1$time)

#Remove colony 65 as it doesn't contain data in the second monitoring period
df2<- filter(df1, colony!= "65")


####QUEEN ACTIVITY####

#make a figure of activity

plot_boris<-ggplot(df2, aes(x = time, y = queen_move)) +
  geom_boxplot(aes(fill=treat))+
theme_classic()+
  scale_fill_manual(values = c("orange", "gray68"))+
  labs(title="Amount of time each queen spent active over two time points",
       x ="Time point", y = "Time spent active (s)", fill="Treatment type")+ 
  scale_y_continuous(breaks=seq(0,2500,500))


plot_boris

#make a figure of inactivity (just to look, won't be used)

#plot_activity2<-ggplot(df2, aes(x = time, y = queen_stationary, colour = treat)) +
#  geom_boxplot(aes())
#plot_activity2

#curiously there is a higher time spent stationary in the treatment than control, this could be due to uneven video footage. So keep figures proportional as below. This will make a supplemental figure of activity

supp.labs <- c("Film period 1 (10% queen mortality)", "Film period 2 (60% queen mortality)")
names(supp.labs) <- c("1", "2")

plot_activity3<-ggplot(df2, aes(x = treat, y = prop_moving, fill=treat)) +
  facet_wrap(~time, labeller = labeller(time=supp.labs))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(), width = 0.2, colour = "black")+
  theme(plot.title.position = "plot", 
        text = element_text(size=8),
        axis.text.x = element_text(),
        axis.line= element_line(),
        legend.position = "none",
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.spacing.x = unit(0.3,"line"))+
  scale_fill_manual(values = c("orange", "gray68"))+
  scale_x_discrete(breaks=c("R","C"),
                   labels=c("R", "C"))+
  labs(x ="Treatment", y = "Proportion of time spent active", fill="Treatment")
  
plot_activity3

ggsave("Treatment_boris_activity.png", plot_activity3, width = 12, height = 12, units = "cm")



#Plot for inactivity (for curiosity)
#plot_activity4<-ggplot(df2, aes(x = time, y = prop_stationary)) +
#  geom_boxplot(aes(fill=treat))+
#  theme_classic()+
#  scale_fill_manual(values = c("gray68", "orange"))+
#  labs(title="Amount of time each queen spent stationary for 1 hour over two time points",
#       x ="Time point", y = "Proportion of time spent stationary", fill="Treatment type")
#plot_activity4

# Test effect of treatment on queen activity level
## I want to estimate the effect of treatment on queen activity for each of the two time points. Normally I would use an Anova, but because the data are proportional, and because of low sample sizes, I will use a binomial glmm model.

model<-glmer(prop_moving ~ treat * time + (1|colony), data = df2, family = "binomial")
summary(model)#This model is incomplete because it only models the raw proportions, rather than the number of success relative to failure. I.e. 5/10 is not the same as 500/1000 even though they produce the same proportion. Instead it is better to model the number of target effects relative to non-target effects

model2 <- glmer(cbind(queen_move,total-queen_move) ~ treat * time + (1|colony), data = df2, family = "binomial") #Note the response is the same as cbind(queen_move,queen_stationary)

#Drop the interaction term
model3 <- glmer(cbind(queen_move,total-queen_move) ~ treat + time + (1|colony), data = df2, family = "binomial") #Note the response is the same as cbind(queen_move,queen_stationary)

#Drop the time component (which is having the strongest effect on the overall result)
model4 <- glmer(cbind(queen_move,total-queen_move) ~ treat + (1|colony), data = df2, family = "binomial") #Note the response is the same as cbind(queen_move,queen_stationary)


#Make a null model
null.model2 <- glmer(cbind(queen_move,total-queen_move) ~ 1 + (1|colony), data = df2, family = "binomial")

anova(model2,model3, model4, null.model2)#The model with time and treatment performs much better than the nulll model. The model with just treatment performs much worse than the null model. Best model is model 2 (which includes the interaction term)

anova(model2, null.model2) #Use this to get the model selection stats for the stats table

summary(model2) #Model shows no effect of treatment on queen movement (but a very strong effect of time (i.e., film period))

summary (model4) #Model4 (the worst model) shows no effect of treatment on queen movement


# Also make figures of aggression and egg laying (though the power will be very low for these)

### Extra note: I worked out that the problem is that geom_point was summarising all of the datapoints within each day as a single point, while jitter was plotting points from individual data points. Therefore the figure I need to summarise these data is simply geom_jitter but not geom_point. As below:

######worker egg laying #####

(plot_boris.worker_egg<-ggplot(df2, aes(x = treat, y = worker_egg, fill = treat)) +
  facet_wrap(~time, labeller = labeller(time=supp.labs))+
  scale_fill_manual(values = c("orange", "gray68"))+
  geom_jitter(shape=21,width = 0.2, height = 0.2)+
  theme(plot.title.position = "plot", 
        text = element_text(size=10),
        axis.text.x = element_text(),
        axis.line= element_line(),
        legend.position = "none",
        panel.border = element_blank(), 
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.spacing.x = unit(0.3,"line"))+
  scale_x_discrete(breaks=c("R","C"),
                   labels=c("R", "C"))+
  scale_y_continuous(breaks=seq(0,2,1))+
  labs(x ="Treatment", y = "Number of egg-laying events"))

ggsave("Treatment_boris_egg.png", plot_boris.worker_egg, width = 12, height = 12, units = "cm")

df2 <-  df2 %>% 
  mutate(worker_egg_tf = if_else(worker_egg>0,TRUE,FALSE))

model.b.we.null1 <- glmer(worker_egg_tf ~ 1 + (1|colony), data = df2, family = "binomial")
model.b.we1 <- glmer(worker_egg_tf ~ treat*time + (1|colony), data = df2, family = "binomial")
anova(model.b.we.null1,model.b.we1) #model that includes treatment and time are not significantly different from the null model and fit the data less well
summary(model.b.we1) #very strange result shown by the model, not quite sure what's going on here. Make time a ranndom effect instead and just look for an effect of treatment across the whole model
model.b.we2 <- glmer(worker_egg_tf ~ treat + (1|time) + (1|colony), data = df2, family = "binomial")
summary(model.b.we2) 



#####worker aggression#####

#plot

(plot_boris.worker_agg2<-ggplot(df2, aes(x = treat, y = worker_agg, fill = treat)) +
  facet_wrap(~time, labeller = labeller(time=supp.labs))+
   scale_fill_manual(values = c("orange", "gray68"))+
   geom_jitter(shape=21,width = 0.2, height = 0.2)+
  theme(plot.title.position = "plot", 
        text = element_text(size=10),
        axis.text.x = element_text(),
        axis.line= element_line(),
        legend.position = "none",
        panel.border = element_blank(), 
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.spacing.x = unit(0.3,"line"))+
  scale_x_discrete(breaks=c("R","C"),
                   labels=c("R", "C"))+
  scale_y_continuous(breaks=seq(0,5,1))+
  labs(x ="Treatment", y = "Number of aggressive incidents"))
  

  ggsave("Treatment_boris_aggression.png", plot_boris.worker_agg2, width = 12, height = 12, units = "cm")
  
  
#models
  
#Test for worker aggression, I initially tried friedman test, but data is not set up for it. 
  friedman.test(worker_agg ~ treat | time,
                data = df2)
  
#glmer poisson model
  
  model.p.wa1<-glmer(worker_agg ~ treat*time + (1|colony), data = df2, family = "poisson")
  summary(model)
  
  
#glmer negative binomial model
  
  model.nb.wa1 <- glmer.nb(worker_agg ~ treat + (1|colony), data = df2)
  summary(model.wa1)

#I suspect that aggression coded as a count is not really appropriate here given the low numbers of each count (and high number of zeros). Might be more informative if we turn aggression into a binary model or if we we use a zero inflated model
  
#Make binary index
 df2 <-  df2 %>% 
    mutate(worker_agg_tf = if_else(worker_agg>0,TRUE,FALSE))

#Model the binary index with a null model
  model.b.wa.null1 <- glmer(worker_agg_tf ~ 1 + (1|colony), data = df2, family = "binomial")
  
  # the binary index with a treatment and time (and their interaction) model
  model.b.wa1 <- glmer(worker_agg_tf ~ treat*time + (1|colony), data = df2, family = "binomial")
  
  # the binary index with a treatment and time (without their interaction) model
  model.b.wa2 <- glmer(worker_agg_tf ~ treat + time + (1|colony), data = df2, family = "binomial")

  # the binary index with just treatment  model
  model.b.wa3 <- glmer(worker_agg_tf ~ treat + (1|colony), data = df2, family = "binomial")
  
  anova(model.b.wa1,model.b.wa2,model.b.wa3) # Best model appears to be model 3 but none of the models are that different from each other. Therefore stick to model 2 as it is more informative (retains the effect of time period)
  
  AIC(model.b.wa1,model.b.wa2,model.b.wa3,model.b.wa.null1) #Best model appears to be the null model.
  
  summary (model.b.wa2) # No effect of treatment on filmed worker aggression
   
  anova (model.b.wa.null1,model.b.wa2) #model that includes treatment and time are not significantly different from the null model and fit the data less well
  
  
# Zero inflated models
  model.nb.wa3 <- glmmTMB(worker_agg ~ treat*time + (1|colony), data = df2, ziformula=~1, family = poisson)
  summary (model.nb.wa3)
  
  
  #Overall there was no effect of treatment on worker aggression using any of the models. The models are not significantly different from each other (or from null models) using a liklihood ratio test. The most informative model is model.b.wa1 (which includes time), so this model is included in the manuscript.
  
######queen egg laying, probably no point in putting this in the paper as very few queen egg laying events were observed
  
  plot_boris.queen_egg<-ggplot(df1, aes(x = time, y = queen_egg, colour = treat)) +
    scale_colour_manual(values = c("orange", "gray68"),
                        breaks = c("R", "C"),
                        labels=c("R", "C"))+
    theme_classic()+
    labs(title="Number of queen egg laying events per colony",
         x ="Day", y = "Number of egg laying events", colour="Treatment")+
    geom_jitter(width = 0.1, height = 0.2)+
    scale_y_continuous(breaks=seq(0,5,1))
  plot_boris.queen_egg
  
  #Can't model queen egg laying, because the sample size is far too low.
