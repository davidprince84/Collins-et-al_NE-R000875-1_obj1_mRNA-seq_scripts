#DavidHCollins
#Created: 3 November 2020

#Aim

#The aim of this script is to analyse data from our main Bombus experiment. This script is to establish whether the numbers of workers that were added or removed or that died were different for treatment and control colonies and to establish a means of controlling for these factors.

###################################################################################################################################################

#First you must always reset R
rm(list=ls())

#make these packages (and associated functions) available to use in the script

library(readr) #read in files
library(tidyverse)
library(bbmle) #produces AIC values which allows you to compare multiple models
library(lme4) ## for glmer
library(glmmTMB) ##for zero inflated glmms
library(effects) ## All effects function which gives confidence intervals of actual values when using a glm
library(jtools) ## summ function which is a much nicer easier to read version of summary with some better info
library(insight)
library(performance)
library(AER) # for dispersion test
library(MASS) # for negative binomial glms

#here are some useful shortcuts to get to know
#ctrl+a highlights the whole script
#ctrl+enter runs the highlighted script or just the selected line of code if nothing is highlighted
#ctrl+1 moves cursor to the console
#ctrl+2 moves it back to the script editor


#set file directory

setwd("Documents/UEA Work documents/Post Doc_LHT/Objective 1 Experiment 1/Analysis")

#This is for the hp: setwd("C:/Users/dhcollins500/OneDrive - University of East Anglia/Documents/Post Doc_LHT/Objective 1 Experiment 1/Analysis/")

#Read in file which contains the numbers of workers added to each colony, the numbers removed due to aggression, the numbers removed due to egg laying, the numbers removed for other reasons (mostly to compensate for newly eclosed workers taking overall number over 20), total removed, workers that died, and number of workers in the colony at the end of the monitoring period (should be 20 but was less than this at the start of the experiment as the numbers were still increasing in each colony)

df1 <- read_csv("numbers.csv")
View(df1)

#Read in survival data to check effect of worker addition/removal/death on queen survival
df2 <- read_csv("survival.csv") %>% 
  dplyr::select(individual,treat,age,censor) %>% 
  rename(col_no = individual)

#Work out total numbers of workers added/removed/died to each colony
df3 <- df1 %>%
  group_by(col_no) %>% 
  summarise(w_add = sum(w_add,na.rm = TRUE),
            removed_aggression = sum(removed_aggression, na.rm = TRUE),
            removed_egglaying = sum(removed_egglaying, na.rm = TRUE),
            removed_other = sum (removed_other, na.rm = TRUE),
            removed_total = sum(removed_total, na.rm = TRUE),
            w_dead = sum (w_dead, na.rm = TRUE))

#Join the summarised numbers and the survival data together
df3 <- inner_join(df2,df3,"col_no")

#Remove the colony (col_no 62) that wasn't used in the experiment (keep the censored individuals, as they still provide useful data)
df3 <- df3 %>% 
  filter (col_no != "62")


#Make treat a factor and reorder it so treatment is the first level
df3$treat <- as.factor(df3$treat)
is.factor(df3$treat)
df3$treat <- factor(df3$treat, levels = c("R", "C"))

#Make a boxplot of each total to estimate whether these metrics differ by treatment

####workers added
df3 %>% 
  ggplot(aes(treat,w_add)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter()

#Get basic measures of how much each group differs
df3 %>% 
  group_by(treat) %>% 
summarize(mean=mean(w_add), n=n(), sd=sd(w_add), se=sd/sqrt(n))


worker_added_model1 <- glm(w_add ~ treat,family="poisson",data = df3)
dispersiontest(worker_added_model1,trafo=1)#Data are overdisperse so use a negative binomial model
worker_added_model2 <- glm.nb(w_add ~ treat,data = df3)
worker_added_model3 <- glm.nb(w_add ~ treat*age,data = df3)
worker_added_model4 <- glm.nb(w_add ~ treat + age,data = df3)
worker_added_null <- glm.nb(w_add ~ 1,data = df3)

AIC(worker_added_model2,worker_added_model3,worker_added_model4,worker_added_null)#Model4 performs significantly better
anova(worker_added_null,worker_added_model4)#Model4 has significantly more predictive value than the null model
summ(worker_added_model4, digits =4) #There was no effect of treatment on the number of workers added to each colony


###workers removed due to aggression
df3 %>% 
  ggplot(aes(treat,removed_aggression)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter()

#Get basic measures of how much each group differs
df3 %>% 
  group_by(treat) %>% 
  summarize(median = median (removed_aggression),mean=mean(removed_aggression), n=n(), sd=sd(removed_aggression), se=sd/sqrt(n))


worker_rem_agg_model1 <- glm(removed_aggression ~ treat,family="poisson",data = df3)
dispersiontest(worker_rem_agg_model1,trafo=1)#Data are overdisperse so use a negative binomial model
worker_rem_agg_model2 <- glm.nb(removed_aggression ~ treat,data = df3)
worker_rem_agg_model3 <- glm.nb(removed_aggression ~ treat*age,data = df3)
worker_rem_agg_model4 <- glm.nb(removed_aggression ~ treat + age,data = df3)
worker_rem_agg_nullmodel <- glm.nb(removed_aggression ~ 1,data = df3)

AIC(worker_rem_agg_model2,worker_rem_agg_model3,worker_rem_agg_model4, worker_rem_agg_nullmodel)#Model4 performs significantly better
anova(worker_rem_agg_nullmodel,worker_rem_agg_model4)#Model4 is better than the null model

summ(worker_rem_agg_model4, digits = 4) #No effect of treatment on the number of workers removed due to aggression

###workers removed due to egg_laying
df3 %>% 
  ggplot(aes(treat,removed_egglaying)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter()


#Get basic measures of how much each group differs
df3 %>% 
  group_by(treat) %>% 
  summarize(median= median(removed_egglaying), mean=mean(removed_egglaying), n=n(), sd=sd(removed_egglaying), se=sd/sqrt(n)) %>% 
  as.data.frame()

worker_rem_egg_model1 <- glm(removed_egglaying ~ treat,family="poisson",data = df3)
dispersiontest(worker_rem_egg_model1,trafo=1)#Data are overdisperse so use a negative binomial model
worker_rem_egg_model2 <- glm.nb(removed_egglaying ~ treat,data = df3)
worker_rem_egg_model3 <- glm.nb(removed_egglaying ~ treat*age,data = df3)
worker_rem_egg_model4 <- glm.nb(removed_egglaying ~ treat + age,data = df3)
worker_rem_egg_nullmodel <- glm.nb(removed_egglaying ~ 1,data = df3)


AIC(worker_rem_egg_model2,worker_rem_egg_model3, worker_rem_egg_model4,worker_rem_egg_nullmodel)#Model4 performs better
anova(worker_rem_egg_nullmodel,worker_rem_egg_model4) #Model 4 predicts significantly more than the null model

summ(worker_rem_egg_model4, digits = 4) #No effect of treatment on the number of workers removed due to egg laying
plot(worker_rem_egg_model4)
plot(allEffects(worker_rem_egg_model4))



#To summarise, significantly more workers were removed due to egg laying from R than from C colonies 


#workers removed due to other reasons
df3 %>% 
  ggplot(aes(treat,removed_other)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter()

#Get basic measures of how much each group differs
df3 %>% 
  group_by(treat) %>% 
  summarize(mean=mean(removed_other), n=n(), sd=sd(removed_other), se=sd/sqrt(n))

worker_rem_other_model1 <- glm(removed_other ~ treat,family="poisson",data = df3)
dispersiontest(worker_rem_other_model1,trafo=1)#Data are overdisperse so use a negative binomial model
worker_rem_other_model2 <- glm.nb(removed_other ~ treat,data = df3)
summary(worker_rem_other_model2)
worker_rem_other_model3 <- glm.nb(removed_other ~ treat*age,data = df3)
summary(worker_rem_other_model3)
anova(worker_rem_other_model2,worker_rem_other_model3)
AIC(worker_rem_other_model2,worker_rem_other_model3)#Model3 performs significantly better

#Model age as a random effect
worker_rem_other_model4<- glm.nb(removed_other ~ treat + (1|age),data = df3)
#Model age as an offset
worker_rem_other_model5<- glm.nb(removed_other ~ treat + offset(log(age)),data = df3)
#Model age as an additional factor
worker_rem_other_model6<- glm.nb(removed_other ~ treat + age,data = df3)
#Null model
worker_rem_other_nullmodel <- glm.nb(removed_other ~ 1,data = df3)

AIC(worker_rem_other_model2,worker_rem_other_model3,worker_rem_other_model4,worker_rem_other_model5,worker_rem_other_model6,worker_rem_other_nullmodel) #Model 6 (with age as a random factor) is the best model
anova(worker_rem_other_model6,worker_rem_other_nullmodel)

summary (worker_rem_other_model6) #Shows that both age and treatment have a significant effect on the number of workers removed for reasons other than egg laying or aggression

summ(worker_rem_other_model6, digits = 4)

##workers removed in total
#make boxplot
df3 %>% 
  ggplot(aes(treat,removed_total)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter()

#Boxplot shows that more were removed from the C than the R colonies

df3 %>% 
  ggplot(aes(age,removed_total,colour=treat)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)

#Scatterplot shows that the colonies that lived longer had more removed from them than the colonies that did not live as long. This indicates that age should be accoutned for in the model 




#make models
worker_removed_tot_model1 <- glm(removed_total ~ treat,family="poisson",data = df3)
dispersiontest(worker_removed_tot_model1,trafo=1)#Data are overdisperse so use a negative binomial model

worker_removed_tot_model2 <- glm.nb(removed_total ~ treat,data = df3)
summary(worker_removed_tot_model2) #significant effect of treatment

worker_removed_tot_model3 <- glm.nb(removed_total ~ treat*age,data = df3)
summary(worker_removed_tot_model3) #significant effect of teatment disappears if age is included
anova(worker_removed_tot_model2, worker_removed_tot_model3) # the two models produce significantly different effects

AIC(worker_removed_tot_model2, worker_removed_tot_model3) #model 3 performs significantly better
#Perhaps model age as a random effect

worker_removed_tot_model4 <- glm.nb(removed_total ~ treat + (1|age),data = df3)
summary(worker_removed_tot_model4)
anova(worker_removed_tot_model2,worker_removed_tot_model4)
anova(worker_removed_tot_model3,worker_removed_tot_model4)
AIC(worker_removed_tot_model2,worker_removed_tot_model3,worker_removed_tot_model4)#model 3 still better

plot(allEffects(worker_removed_tot_model3))
plot(allEffects(worker_removed_tot_model4)) # can see that treatment might still have an effect if age is accounted for



#When accounting for age - cannot use the index (number removed / age) with poisson or negative binomial because these models assume integer values. Instead could use age as an offset, which should account for the increased longevity of C queens
worker_removed_tot_model5 <- glm.nb(removed_total ~ treat + offset(log(age)),data = df3)
summary(worker_removed_tot_model5)#There was no significant effect of treatment on worker removal (for all reasons)
AIC(worker_removed_tot_model2, worker_removed_tot_model3, worker_removed_tot_model5) #Model 5 better than model 2 but not 3
anova(worker_removed_tot_model2,worker_removed_tot_model5)

#Try plotting with age as an additional fixed factor
worker_removed_tot_model6 <- glm.nb(removed_total ~ treat + age,data = df3)
summary(worker_removed_tot_model6)
plot(worker_removed_tot_model6) # diagnostics look better
AIC(worker_removed_tot_model2,worker_removed_tot_model3, worker_removed_tot_model5,worker_removed_tot_model6) #Model 6 provides the best fit to the data

#As model 6 is our best model, it seems apparent that treatment really does have an effect on the number of workers removed. Ideally the number of workers removed needs to be accounted for in my other models (e.g. on fecundity and longevity)



df3 %>% 
  ggplot(aes(treat,mean_removed_total)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter()


#As age has a significant effect might be better to plot these data as coloured scatterplot
df3 %>% 
  ggplot(aes(x=age,y=removed_total, colour=treat)) + 
  geom_point()+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)

#Due to the above we know that workers removed for reasons of aggression and egg laying are not significantly affecting the results, so will focus on workers removed in the paper



##workers that died in total
#Make boxplot
df3 %>% 
  ggplot(aes(treat,w_dead)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter()

#Get basic measures of how much each group differs
df3 %>% 
  group_by(treat) %>% 
  summarize(mean=mean(w_dead), n=n(), sd=sd(w_dead), se=sd/sqrt(n))

#Make models
dead.model1 <- glm(w_dead ~ treat,family="poisson",data = df3)
dispersiontest(dead.model1,trafo=1) #Data are overdispersed so use a negative binomial model

dead.model2 <- glm.nb(w_dead ~ treat,data = df3)
summary(dead.model2)

dead.model3 <- glm.nb(w_dead ~ treat*age,data = df3)
summary(dead.model3) 
anova(dead.model2, dead.model3) # the two models produce significantly different effects
AIC(dead.model2, dead.model3) #model 3 performs significantly better, now test if treatment still has an effect after accounting for age using glmm

dead.model4 <- glm.nb(w_dead ~ treat + (1|age),data = df3)

dead.model5 <- glm.nb(w_dead ~ treat + age,data = df3)

dead.model.null <- glm.nb(w_dead ~ 1,data = df3)

AIC(dead.model2, dead.model3, dead.model4, dead.model5,dead.model.null) #model 5 performs very slightly better than model 4 and both models perform better than the models without age included
anova(dead.model4, dead.model5) # LRS shows there is a difference in the fit of the two models
anova(dead.model5,dead.model.null) # LRS shows there is a difference in the fit of the two models


summ(dead.model5, digits = 4) # no effect of treatment on the number of worker deaths when treating age as a random factor
plot(allEffects(dead.model5))


#As age has a significant effect might be better to plot these data as coloured scatterplot
#Make boxplot
df3 %>% 
  ggplot(aes(x=age,y=w_dead, colour=treat)) + 
  geom_point()+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)

worker_removed_tot_model5 <- glm.nb(mean_removed_total ~ treat,data = df3)
summary(worker_removed_tot_model5)#There was no significant effect of treatment on worker death rate (number of dead workers as a function of age at which the queen died)
