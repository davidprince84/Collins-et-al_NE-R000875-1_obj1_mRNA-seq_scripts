#DavidHCollins
#Created: 6 April 2020

#Aim

#The aim of this script is to analyse data from our main Bombus experiment. The aim of this script was to analyse the data generated from monitoring the worker aggression to the queen, and worker egg-laying.

########################################################################################################################

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
library(survival)
library("survminer")


#here are some useful shortcuts to get to know
#ctrl+a highlights the whole script
#ctrl+enter runs the highlighted script or just the selected line of code if nothing is highlighted
#ctrl+1 moves cursor to the console
#ctrl+2 moves it back to the script editor


#set file directory

setwd("Documents/UEA Work documents/Post Doc_LHT/Objective 1 Experiment 1/Analysis")

#If using HP computer
#setwd("C:/Users/dhcollins500/OneDrive - University of East Anglia/Documents/Post Doc_LHT/Objective 1 Experiment 1/Analysis/")


#Read in individual level data for aggression

aggression_data <- read_csv("aggression.csv")
View(aggression_data)
df1<-aggression_data

df1 %>% count(aggression)


#Make treat a factor and reorder it so treatment is the first level
df1$treat <- as.factor(df1$treat)
is.factor(df1$treat)
df1$treat <- factor(df1$treat, levels = c("R", "C"))


#May Want to check whether models are overdispersed

overdisp_fun <- function(model) {
  ## number of variance parameters in 
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}


####Worker Egg Laying######

#The aim of the following code is to check whether treatment and time (day) have an effect on the worker egg laying recorded in the experiment

### Egg laying figure

#Make a figure showing the change in aggression with time 

(
  plot_sp.time.workereggs<-ggplot(df1, aes(x = day, y = eggs, fill = treat)) +
   geom_jitter(aes(), width = 0.5, shape = 21, colour = "black")+
  scale_fill_manual(values = c("orange", "gray68"))+
   theme(plot.title.position = "plot", 
         text = element_text(size=11),
         axis.text.x = element_text(),
         axis.line= element_line(),
         legend.position = "none",
         panel.border = element_blank(), 
         panel.background = element_blank(), 
         panel.grid = element_blank(), 
         panel.spacing.x = unit(0.3,"line"))+
  labs(x ="Day", y = "Number of worker egg laying events", 
       colour="Treatment")+
  scale_y_continuous(breaks=seq(0,5,1))+
  scale_x_continuous(breaks=seq(0,154,8))+
    geom_vline(xintercept=25, size = 0.5, linetype='dashed', col = 'black'))


ggsave("Treatment_time_workereggs.png", plot_sp.time.workereggs, width = 16, height = 10, units = "cm")
ggsave("Treatment_time_workereggs.svg", plot_sp.time.workereggs, width = 16, height = 6, units = "cm")




#Looks like there is more egg laying in the treatment early on and more in the control later on



### Egg laying models

##Poisson models


#First make the null models

null.model.e1 <- glmer( eggs ~ 1 + (1|ID) + (1|day), family=poisson, data=df1)
null.model.e2 <- glmer( eggs ~ 1 + (1|day), family=poisson, data=df1)
null.model.e3 <- glmer( eggs ~ 1 + (1|ID), family=poisson, data=df1)

anova(null.model.e1,null.model.e2,null.model.e3) ##nullmodel2 fails to converge, however as nullmodel1 is alright(and has the best AIC) I will keep going with this model

overdisp_fun(null.model.e1)
#model is not overdispersed

#full models with fixed effects
model.e1 <- glmer( eggs ~ day + (1|ID) + (1|day), family=poisson, data=df1)
model.e2 <- glmer( eggs ~ treat + (1|ID) + (1|day), family=poisson, data=df1)
model.e3 <- glmer( eggs ~ treat*day + (1|ID) + (1|day), family=poisson, data=df1)
model.e4 <- glmer( eggs ~ treat*poly(day,2) + (1|ID) + (1|day), family=poisson, data=df1) #just to check if a quadratic term makes the  model fit the data better

anova(model.e1, model.e2, model.e3)
#results suggest model 3 is the best, and therefore treatment and day should be incorportated as fixed effects

anova(model.e3,model.e4)
#results suggest model 4 is better than model 3, and therefore time (day) should be fitted as a quadratic effect

overdisp_fun(model.e4)
#model4 is not overdispersed

summary(model.e4)
#all of the terms, except treatment:poly(day,2)2 are significant

plot(model.e4)

plot (allEffects(model.e4))

###
#Note, these models might not be appropriate because of the huge number of zeros, instead we shall use Ben Bolkers glmmTMB package. We can fit a zero inflated poisson, then zi neg binom, then hurdle models

#zero inflated poisson
fit_zipoisson <- glmmTMB(eggs~ treat*day + (1|ID) + (1|day),data=df1,ziformula=~1,family = poisson)
summary(fit_zipoisson)
plot(allEffects(fit_zipoisson))

# Works as a poly model
fit_zipoisson2 <- glmmTMB(eggs~ treat*poly(day,2) + (1|ID) + (1|day),data=df1,ziformula=~1,family = poisson)
summary(fit_zipoisson2)
plot(allEffects(fit_zipoisson2))

anova(fit_zipoisson,fit_zipoisson2) #zi poisson model performs better as a poly model, however could not work out how to check for overdispersion



#zero inflated negative binomial
fit_zinbinom2 <- glmmTMB(eggs~ treat*poly(day,2) + (1|ID) + (1|day),data=df1,ziformula=~1,family = nbinom1)

#This model shows clear effect of treatment, time, and treatment:time interaction
summary(fit_zinbinom2)
plot(allEffects(fit_zinbinom2))



#The negative binomial model is slightly better
anova(fit_zinbinom2,fit_zipoisson2) 

#Try a hurdle model
fit_hurdlepoisson <- glmmTMB(eggs~ treat*day + (1|ID) + (1|day),data=df1,ziformula=~1,family=list(family="truncated_poisson",link="log"))
summary(fit_hurdlepoisson)

AIC(fit_zinbinom2,fit_hurdlepoisson) #zi nb model performs better, no need to check for overdispersion with this model

#Make null model
null_zinbinom2 <- glmmTMB(eggs~ 1 + (1|ID) + (1|day),data=df1,ziformula=~1,family = nbinom1)

anova(null_zinbinom2,fit_zinbinom2)
#the full model provides a significantly better fit to the data than the null model

confint(fit_zinbinom2) 

              ## Compared the time when egg laying was first observed in each colony between treatments ##

#If the first time we observe egg laying in a given colony is reflective of the competition point in that colony then this could be a useful analysis.


#If the first time we observe aggression in a given colony is reflective of the competition point in that colony then this could be a useful analysis.

#Create a new df which summarises the first day when worker egg laying was observed (or the date of death if egg laying was not observed) for each colony, a censor column which indicates whether egg laying was observed or not, and a treatment column that indicates the treatment for each colony
df2 <- df1 %>%
  filter(eggs != 0) %>% #removes all instances where zero aggression was observed
  group_by(col_no) %>%
  filter (row_number(day)==1) %>% #removes all rows following the first value where eggs were seen (i.e., if egglaying was seen on multiple different days)
  arrange(col_no) %>% #arrange in colony order
  mutate (censor  = 1) %>% #assign a value of 1 to a new censor column for all colonies when egg laying was observed
  select (col_no, day, eggs,censor) %>% #keep the values for colony number, the day that eggs were first seen, the treatment, and the number of eggs that were seen on the first observation.
  right_join(read_csv("survival.csv") %>% #combine the new df with the 'survival df' so values for individuals longevity are included
               select (individual, age, treat)
             , by = c("col_no" = "individual")) %>% 
  mutate(censor = if_else(is.na(censor), 0,censor)) %>% #assign a value of 0 in the censor column to colonies where aggression was not observed
  mutate (day = if_else(is.na(day),age,day)) %>%  #assign the longevity of individuals where aggression was not observed to the day column
  mutate (eggs = if_else(is.na(eggs),0,eggs)) %>% #assign a value of 0 aggression to colonies where no aggression was recorded (so it's easy to track that these are the colonies that were censored)
  select (col_no,treat,day,eggs,censor) %>%
  filter(col_no != "62") #exclude colony 62 because this colony was treated differently from the others

#The resulting df has a col_no column (for all colonies where egg laying was monitored), a day column (the first day egg laying was monitored OR the day the colony died if egg laying was not observed), a censor column (indicating if egg laying was observed (=1) or not observed (=0), an age column (the age the colony reached, should be the same as the day column for colonies where no egg laying was observed), and a treat column (indicating the treatment each colony was assigned to), note. colony 62 was also censored because it was treated in a different way to the other colonies, even though egg laying was seen in this colony)

#Use this new dataframe to determine whether egg laying was seen earlier in treatment compared to control colonies (remember to exclude the censored colonies)

df2 %>% 
  filter (censor == 1) %>% 
  group_by(treat) %>% 
  summarise(mean = mean(day), median = median(day), sd = sd(day), n = n()) #Eggs were seen earlier in R compared to C colonies

# View this using a survival curve

fit<- survfit(Surv(day,censor) ~ treat, data = df2)

(plot_worker_eggs_survival.1 <- ggsurvplot(fit,
                               data=df2, 
                               censor=TRUE, 
                               legend="right",
                               conf.int = FALSE,
                               palette = c("orange", "gray68"),
                               legend.title = "Treatment",
                               legend.labs= c("R","C"))+
    labs(x = "Time (days)",
         y = "Proportion"))

#Save the survival plot with the following function
grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}

#Then save the plot
ggsave("Treatment_first_worker_eggs.png", plot_worker_eggs_survival.1, width = 16, height = 10, units = "cm")
ggsave("Treatment_first_worker_eggs.svg", print(plot_worker_eggs_survival.1), width = 16, height = 10, units = "cm")


#Make a coxph model to test the effect of treatment on age of 1st reproduction
egg_coxmodel <- coxph(Surv(day,censor) ~ treat, data = df2)

#Test assumption that regression coefficients are constant over time
Test1<-cox.zph(egg_coxmodel)

Test1 ##we cannot accept that the regression coefficients are constant over time (as p<0.05) and that the slope is significantly different from zero for some treatment levels (this is because the censors are included)

ggcoxzph(Test1) ##Graphical illustration of the fact the decline hazard with time, i.e. hazards are non-proportional

#Check using dfbeta method
ggcoxdiagnostics(egg_coxmodel, type = "dfbeta",linear.predictions = FALSE, ggtheme = theme_bw()) ##Graphical illustration that there are three data points affecting the model

#Check using residual deviance method
ggcoxdiagnostics(egg_coxmodel, type = "deviance",linear.predictions = FALSE, ggtheme = theme_bw()) ##Graphical illustration that there are three data points affecting the model


####

##Create new model with outlying datapoints filtered

#Filter the data points
df2b <- df2 %>% 
  filter(col_no != 5) %>% 
  filter(col_no != 56) %>% 
  filter(col_no != 70)
#Note, these individuals are all censored in the analysis anyway, so removing them should not be a problem.


egg_coxmodel2 <- coxph(Surv(day,censor) ~ treat, data = df2b)

#Test assumption that regression coefficients are constant over time
Test2<-cox.zph(egg_coxmodel2)

Test2 #Regression coefficients are now constant

ggcoxzph(Test2) #

ggcoxdiagnostics(egg_coxmodel2, type = "dfbeta",linear.predictions = FALSE, ggtheme = theme_bw()) ##Graphical illustration that there is a datapoint affecting the model (removing point doesn't change the result)

####

##Make new model without censors
df2c <- df2 %>% 
  filter(censor == 1) 

egg_coxmodel3 <- coxph(Surv(day,censor) ~ treat, data = df2c)

Test3<-cox.zph(egg_coxmodel3)

Test3 ##we can now accept that the regression coefficients are constant over time (as p>0.05) and that the slope is significantly different from zero for some treatment levels (removing the censors improved the model)

ggcoxzph(Test3)

ggcoxdiagnostics(egg_coxmodel3, type = "dfbeta",linear.predictions = FALSE, ggtheme = theme_bw()) ##Graphical illustration that there is a datapoint affecting the model (removing point doesn't change the result)

#Summary of model
summary(egg_coxmodel) 
summary(egg_coxmodel2) 
summary(egg_coxmodel3) # Workers in R colonies were observed laying eggs earlier than workers in C colonies. However, please be aware that the hazards are not proportional. The proportionality assumption *is* met when either the outlier individuals or all the censored individuals are removed from the analysis altogether. Removing either group has no effect on the results.



                                                        ####Aggression######

#The aim of the following code is to check whether treatment and time (day) have an effect on the aggression recorded in the experiment

### Aggression figure:

#Make a figure showing the change in aggression with time 

(plot_sp.time.aggression<-ggplot(df1, aes(x = day, y = aggression, fill = treat)) +
   geom_jitter(aes(), width = 0.5, shape = 21, colour = "black")+
   scale_fill_manual(values = c("orange","gray68"))+
   theme(plot.title.position = "plot", 
         text = element_text(size=11),
         axis.text.x = element_text(),
         axis.line= element_line(),
         legend.position = "none",
         panel.border = element_blank(), 
         panel.background = element_blank(), 
         panel.grid = element_blank(), 
         panel.spacing.x = unit(0.3,"line"))+
   labs(x ="Day", 
        y = "Number of aggressive incidents", 
        colour="Treatment")+
   scale_y_continuous(breaks=seq(0,5,1))+
   scale_x_continuous(breaks=seq(0,154,8)))

(plot_sp.time.aggression2<-ggplot(df1, aes(x = day, y = aggression, fill = treat)) +
    geom_jitter(aes(), width = 0.5, shape = 21, colour = "black")+
    scale_fill_manual(values = c("#FC8E0B","gray68"))+
    theme(plot.title.position = "plot", 
          text = element_text(size=10),
          axis.text.x = element_text(),
          axis.line= element_line(),
          legend.position = "none",
          panel.border = element_blank(), 
          panel.background = element_blank(), 
          panel.grid = element_blank(), 
          panel.spacing.x = unit(0.3,"line"))+
    labs(x ="Day", 
         y = "Number of aggressive incidents", 
         colour="Treatment")+
    scale_y_continuous(breaks=seq(0,5,1))+
    scale_x_continuous(breaks=seq(0,154,8)))

#FC8E0B

ggsave("Treatment_time_aggression.svg", plot_sp.time.aggression, width = 16, height = 6, units = "cm")
ggsave("Treatment_time_aggression2.svg", plot_sp.time.aggression2, width = 16, height = 10, units = "cm")

### Aggression models:

##Poisson model # please note that the negative binomial models are *not* appropriate for these data due to the high number of zeros

#First make the null models

null.model.a1 <- glmer( aggression ~ 1 + (1|ID) + (1|day), family=poisson, data=df1)
null.model.a2 <- glmer( aggression ~ 1 + (1|day), family=poisson, data=df1)
null.model.a3 <- glmer( aggression ~ 1 + (1|ID), family=poisson, data=df1)

anova(null.model.a1,null.model.a2,null.model.a3) #Note model 3 failed to converge, but at any rate the null model with the larger number of random effects was fine, and performed the best


overdisp_fun(null.model.a1)
#model is not overdispersed

#full models with fixed effects
model.a1 <- glmer( aggression ~ day + (1|ID) + (1|day), family=poisson, data=df1)
model.a2 <- glmer( aggression ~ treat + (1|ID) + (1|day), family=poisson, data=df1)
model.a3 <- glmer( aggression ~ treat*day + (1|ID) + (1|day), family=poisson, data=df1)
model.a4 <- glmer( aggression ~ treat*poly(day,2) + (1|ID) + (1|day), family=poisson, data=df1) #just to check if a quadratic term makes the  model fit the data better



anova(null.model.a1,model.a1) #day increases the power of the model
anova(null.model.a1,model.a2) #treatment does not appear to increase the power of the model
anova(null.model.a1,model.a3) #day increases the power of the model

anova(model.a1,model.a3,model.a4) #the three models aren't significantly different from each other, and model 1 provides a slightly better fit to the data. From this I dont't think the quadratic terms in model 4 are necessary

overdisp_fun(model.a3)
#model is not overdispersed

summary(model.a3)

#drop treatment term

model.a5 <- glmer( aggression ~ treat*day + (1|ID), family=poisson, data=df1)
model.a6 <- glmer( aggression ~ treat*day + (1|day), family=poisson, data=df1)
model.a7 <- glm( aggression ~ treat*day, family=poisson, data=df1)


anova(model.a3,model.a5,model.a6,model.a7) #model.a3 is still the best fit, it would be an even better fit to the data without the treatment fixed effect. Therefore it seems that treatment had no effect on aggressive incidents

#zero inflated poisson: use the program glmmTMB and refit prefered poisson model (model3). Then fit it with a zero inflation formula
model.a3 <- glmmTMB(aggression~ treat*day + (1|ID) + (1|day),data=df1,family = poisson)
model.a8 <- glmmTMB(aggression~ treat*day + (1|ID) + (1|day),data=df1,ziformula=~1,family = poisson)

anova(model.a3, model.a8) #The zero inflated model provides a significantly better fit to the data

model.a9 <- glmmTMB(aggression~ treat*day + (1|ID) + (1|day),data=df1,ziformula=~1,family = nbinom1)

anova(model.a8, model.a9) #The zero inflated model poisson provides a slightly better fit to the data than the zero inflated negative binom, but not much in it. Stick to negative binom because of possible overdispersion

summary(model.a9) #Still shows significant effect of time but not of treatment

null.model.a9 <- glmmTMB(aggression~ 1 + (1|ID) + (1|day),data=df1,ziformula=~1,family = nbinom1)

anova(null.model.a9, model.a9) #The zero inflated model negative binomial with treatment and day as fixed effects provides a significantly better fit to the data than the null model


plot(allEffects(model.a8))

          ## Compared the time when aggression was first observed in each colony between treatments ##

#If the first time we observe aggression in a given colony is reflective of the competition point in that colony then this could be a useful analysis.

#Create a new df which summarises the first day when aggression was observed (or the date of death if aggression was not observed) for each colony, a censor column which indicates whether aggression was observed or not, and a treatment column that indicates the treatment for each colony
df3 <- df1 %>%
  filter(aggression != 0) %>% #removes all instances where zero aggression was observed
  group_by(col_no) %>%
  filter (row_number(day)==1) %>% #removes all rows following the first value where eggs were seen (i.e., if egglaying was seen on multiple different days)
  arrange(col_no) %>% #arrange in colony order
  mutate (censor  = 1) %>% #assign a value of 1 to a new censor column for all colonies when aggression was observed
  select (col_no, day, aggression,censor) %>% #keep the values for colony number, the day that eggs were first seen, the treatment, and the number of eggs that were seen on the first observation.
right_join(read_csv("survival.csv") %>% #combine the new df with the 'survival df' so values for individuals longevity are included
                   select (individual, age, treat)
                 , by = c("col_no" = "individual")) %>% 
  mutate(censor = if_else(is.na(censor), 0,censor)) %>% #assign a value of 0 in the censor column to colonies where aggression was not observed
  mutate (day = if_else(is.na(day),age,day)) %>%  #assign the longevity of individuals where aggression was not observed to the day column
  mutate (aggression = if_else(is.na(aggression),0,aggression)) %>% #assign a value of 0 aggression to colonies where no aggression was recorded (so it's easy to track that these are the colonies that were censored)
  select (col_no,treat,day,aggression,censor) %>%
  filter(col_no != "62") #exclude colony 62 because this colony was treated differently from the others

#The resulting df has a col_no column (for all colonies where aggression was monitored), a day column (the first day aggression was monitored OR the day the colony died if aggression was not observed), a censor column (indicating if aggression was observed (=1) or not observed (=0), an age column (the age the colony reached, should be the same as the day column for colonies where no aggression was observed), and a treat column (indicating the treatment each colony was assigned to), note. colony 62 was also censored because it was treated in a different way to the other colonies, even though aggression was seen in this colony)

#Use this new dataframe to determine whether aggression was seen earlier in treatment compared to control colonies (remember to exclude the censored colonies)

df3 %>% 
  filter(censor==1)%>% 
  group_by(treat) %>% 
  summarise(mean = mean(day), sd = sd(day), n = n()) #Aggression was seen slightly earlier in R compared to C colonies in the colonies where aggression was observed

# View this using a survival curve

fit2<- survfit(Surv(day,censor) ~ treat, data = df3)

(plot_aggression_survival <- ggsurvplot(fit2,
                               data=df3, 
                               censor=TRUE, 
                               legend="right",
                               conf.int = FALSE,
                               palette = c("orange", "gray68"),
                               legend.title = "Treatment",
                               legend.labs= c("R","C"))+
    labs(x = "Time (days)",
         y = "Proportion"))

ggsave("Treatment_first_aggression.png", plot_aggression_survival, width = 16, height = 10, units = "cm")

#Make a coxph model to test the effect of treatment on age of 1st aggression
aggression_coxmodel <- coxph(Surv(day,censor) ~ treat, data = df3)

#Test assumption that regression coefficients are constant over time
Test2<-cox.zph(aggression_coxmodel)

Test2 ##we can accept that the regression coefficients are constant over time (p>0.05) and that the slope is not significantly different from zero, for all treatment levels

ggcoxzph(Test2) ##Graphical illustration of the fact there is no pattern of hazzard with time, i.e. hazards are proportional, not non-proportional

ggcoxdiagnostics(aggression_coxmodel, type = "dfbeta",linear.predictions = FALSE, ggtheme = theme_bw()) ##Graphical illustration that there are no overly influential datapoints affecting the model

#Summary of model
summary(aggression_coxmodel) # There was no effect of treatment on the first day aggression was observed.



