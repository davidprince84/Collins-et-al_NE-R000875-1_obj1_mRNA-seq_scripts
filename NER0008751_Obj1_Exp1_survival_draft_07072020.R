#DavidHCollins
#Created: 7 July 2020

#Aim

#The aim of this draft script is to analyse the survival data from the main Bombus experiment. This script covers figure 2, and is the main result from the experiment.


#First you must always reset R
rm(list=ls())

#make these packages (and associated functions) available to use in the script

library(readr) #read in files
library(tidyverse)#ggplot and dplyr
#these packages all help build the survival curves and the coxp models
library(ggfortify)
library(survival)
library("survminer")
library(bbmle) #produces AIC values which allows you to compare multiple models
library(coxme) #for mixed effects coxph models

#here are some useful shortcuts to get to know
#ctrl+a highlights the whole script
#ctrl+enter runs the highlighted script or just the selected line of code if nothing is highlighted
#ctrl+1 moves cursor to the console
#ctrl+2 moves it back to the script editor


#individual level data for fecundity and cell counts

setwd("~/Documents/UEA Work documents/Post Doc_LHT/Objective 1 Experiment 1/Analysis/")
cells <- read_csv("cells.csv")
View(cells)
df1<-cells


#do you have the right data?
glimpse(df1)

#turn the count and day columns into a factors

df1$count<-as.factor(df1$count)
df1$day<-as.factor(df1$day)

is.factor(df1$count)
is.factor(df1$day)


######################### DATAFRAMES CONTAINING VARIOUS SUMMARIES OF THE EGG COUNT DATA #########################

# Make a new dataframe which excludes all of the egg counts that occurred after each queens death (i.e. Only use data that says 'T' in the queenright column)
df2 <- filter(df1, queenright == "TRUE")

#And a dataframe which only includes only the egg counts from colonies where the queen has died i.e. data that says 'F' in the queenright column 
df3 <- filter(df1, queenright == "FALSE")

# Make a new dataframe which excludes all of the egg counts that occurred before treatment (i.e. Only use data that says 'T' in the after_treatment column); and only includes the egg counts from before each queen death
df4 <- filter(df1, queenright == "TRUE", after_treatment == "TRUE")

#And a dataframe which only includes the egg counts from colonies before treatment had started i.e. data that says 'F' in the after_treatment column; and only includes the egg counts from before each queen death. Note that this removes colony 62 which did not have any egg counts for the first few days. This is fine as this individual should be censored anyway
df5 <- filter(df1, queenright == "TRUE", after_treatment == "FALSE")

# Make a new dataframe which excludes all of the egg counts that occurred from count 14 onwards when the first worker egg laying was observed (i.e. Only use data that says 'T' in the before_worker_eggs column), and which excludes egg counts from before treatment started, and which excludes egg counts from after the queen died (i.e. data that says 'T' in the queenright column and 'T' in the after_treatment column)
df6 <- filter(df1, queenright == "TRUE", after_treatment == "TRUE", before_worker_eggs == "TRUE")

########################################TABLE OF SUMMARY STATISTICS #############################################

# Make a new summary dataframe from df2 which summarises the mean and total egg counts of each queen across the experiment (except for the last two counts after the queen had died)

tot.summary.stats <- df2 %>%
  group_by(individual) %>%
  summarise (
    mean_eggs_across_experiment = mean(eggs),
    total_eggs_across_experiment = sum(eggs))

# Make a new summary dataframe from df4 which summarises the mean egg counts of each queen excluding the first two counts from before treatments had started and the last two counts after the queen had died)

after_treat.summary.stats <- df4 %>%
  group_by(individual) %>%
  summarise (
    mean_eggs_after_treatment = mean(eggs),
    total_eggs_after_treatment = sum(eggs))

# Make a new summary dataframe from df5 which summarises the mean and total egg counts of each queen from before treatment had started

before_treat.summary.stats <- df5 %>%
  group_by(individual) %>%
  summarise (
    mean_eggs_before_treatment = mean(eggs),
    total_eggs_before_treatment = sum(eggs))

# Make a new summary dataframe from df6 which summarises the mean and total egg counts of each queen in the experiment before workers were observed laying eggs and excluding the first two counts from before treatment had started (and except for the last two counts after the queen had died)


before_worker_eggs.summary.stats <- df6 %>%
  group_by(individual) %>%
  summarise (
    mean_eggs_before_worker_egg_laying = mean(eggs),
    total_eggs_before_worker_egg_laying = sum(eggs))


#Import survival data and add these summary tables to it. Survival is used as this censors all of the queens removed for RNA extraction, and the following queens: queens 57 (censored due to workers killing her after she fell in fairy liquid), 62 (replacement of 57, but had no measurements taken of her fecundity at the beginning of the experiment, and she spent several weeks with many more workers than 20).


survival <- read_csv("survival.csv")
View(survival)

#What was the mean age of queens
mean(survival$age)
#89.7d

#What was the mean age of queens for each treatment
survival %>% 
  group_by(treat) %>% 
  summarise(mean = mean(age))
#C = 106d ; R = 73.6d - these values might be different to the manuscript as the include censored individuals

# Add the metrics from the three summary stats dataframes to new df7

df7 <- inner_join(survival,tot.summary.stats,by="individual")
df7 <- inner_join(df7,after_treat.summary.stats,by="individual")
df7 <- inner_join(df7,before_treat.summary.stats,by="individual")
df7 <- inner_join(df7,before_worker_eggs.summary.stats,by="individual")

# Make a new dataframe from survival4 which excludes all of the censored individuals which were removed either for RNA extraction or accidently killed (i.e. exclude all entries that have 0 in the censor column)

#Make treat a factor
df7$treat <- as.factor(df7$treat)
is.factor(df7$treat)
df7$treat <- factor(df7$treat, levels = c("R", "C"))

df8 <- filter(df7, censor == "1")

# 1.6 Suvival curves for each treatment

# make the survival curve
fit<- survfit(Surv(age, censor) ~ treat, data = df7) #it is important to use df7 not df8 as the censored individuals are useful for survival curves
print(fit) 

#There are 36 R colonies, 35 C colonies, and 24 R deaths to 22 C deaths. The median survival is 85 days for R colonies and 119 days for C colonies.

#Plot it in ggplot3
(plot_survival.1 <- ggsurvplot(fit,
                              data=df7, 
                              censor=TRUE, 
                              legend="right",
                              conf.int = TRUE,
                              palette = c("gold", "gray68"),
                              legend.title = "Treatment",
                              legend.labs= c("R","C"))+
  labs(title= "Survival of Bombus terrestris females reared under two egg removal regimes", 
       x = "Time (days)",
       y = "Proportion alive"))

#Plot without confidence intervals which looks nicer
(plot_survival.2 <- ggsurvplot(fit,
                              data=df7, 
                              censor=FALSE, 
                              legend="none",
                              conf.int = FALSE,
                              palette = c("orange", "black"),
                              legend.title = "Treatment",
                              linetype = c( "solid", "dashed"),
                              legend.labs= c("R", "C"),
                              ggtheme = theme_classic2(base_size=11),
                              font.family = "Arial")+
  labs(x = "Time (days)",
       y = "Proportion alive"))
#Plot with larger font and 'a' removed for use in talks/other documents
(plot_survival.3 <- ggsurvplot(fit,
                               data=df7, 
                               censor=FALSE, 
                               legend="none",
                               conf.int = FALSE,
                               palette = c("#FC8E0B", "black"),
                               legend.title = "Treatment",
                               linetype = c( "solid", "dashed"),
                               legend.labs= c("R", "C"),
                               ggtheme = theme_classic2(base_size=14),
                               font.family = "Arial")+
    labs(x = "Time (days)",
         y = "Proportion alive"))


ggsave("plot_survival.svg", print(plot_survival.2), width = 8, height = 6, units = "cm")
ggsave("plot_survival2.svg", print(plot_survival.3), width = 16, height = 10, units = "cm")

svg("plot_survival.svg")
print(plot_survival.2)
dev.off()

#Also useful to visualise as a boxplot, use df8 as this does not include the censors
(plot_bp.survival<-ggplot(df8, aes(x = treat, y = age, fill = treat)) +
  geom_boxplot()+
  scale_fill_manual(values = c("orange", "gray68"))+ 
  labs(x="Treatment", y = "Longevity")+
  theme(plot.title.position = "plot", 
          text = element_text(size=11),
          axis.text.x = element_text(),
          axis.line= element_line(),
          legend.position = "none",
          panel.border = element_blank(), 
          panel.background = element_blank(), 
          panel.grid = element_blank(), 
          panel.spacing.x = unit(0.3,"line"))+ 
  geom_jitter(width = 0.2)+
  annotate("segment", x = 1, xend = 2, y = 170, yend = 170, colour = "black", size=0.5)+
  annotate("segment", x = 2, xend = 2, y = 165, yend = 170, colour = "black", size=0.5)+
  annotate("segment", x = 1, xend = 1, y = 165, yend = 170, colour = "black", size=0.5)+
    geom_text(aes( x=1.5, y=175, label="**"), color="black",size=5))

#This suitable for posters
(plot_bp.survival2<-ggplot(df8, aes(x = treat, y = age, fill = treat)) +
    geom_boxplot()+
    scale_fill_manual(values = c("#FC8E0B", "gray68"))+ 
    labs(x="Treatment", y = "Longevity")+
    theme(plot.title.position = "plot", 
          text = element_text(size=8),
          axis.text.x = element_text(),
          axis.line= element_line(),
          legend.position = "none",
          panel.border = element_blank(), 
          panel.background = element_blank(), 
          panel.grid = element_blank(), 
          panel.spacing.x = unit(0.3,"line"))+ 
    geom_jitter(width = 0.2)+
    annotate("segment", x = 1, xend = 2, y = 170, yend = 170, colour = "black", size=0.5)+
    annotate("segment", x = 2, xend = 2, y = 165, yend = 170, colour = "black", size=0.5)+
    annotate("segment", x = 1, xend = 1, y = 165, yend = 170, colour = "black", size=0.5)+
    geom_text(aes( x=1.5, y=175, label="**"), color="black",size=5))

ggsave("boxplot_longevity.svg", plot_bp.survival, width = 8, height = 6, units = "cm")
ggsave("boxplot_longevity2.svg", plot_bp.survival2, width = 16, height = 10, units = "cm")


##Test for effect of treatment on overall longevity

##Build model
model1 <- lm(age ~ treat, data = df7)

##test assumptions
autoplot(model1 , smooth.colour = NA) #these data aren't perfect but will do for now

## use anova() to test if there is an effect of Treatment
anova(model1)

##use summary() to test what the nature of the effect of treatment is, and see if it makes sense with the boxplot itself
summary(model1)


##Cox Proportional Hazzards test for effect of treatment on longevity (CoxPH is a better model as it takes into account the effect of deathrate with time)

#Make the model
model2 <- coxph(Surv(age, censor) ~ treat, data = df7)

#Test assumption that regression coefficients are constant over time
Test1<-cox.zph(model2)

Test1 ##we can accept that the regression coefficients are constant over time (p>0.05) and that the slope is not significantly different from zero, for all treatment levels

ggcoxzph(Test1) ##Graphical illustration of the fact there is no pattern of hazzard with time, i.e. hazzards are proportional, not non-proportional

ggcoxdiagnostics(model2, type = "dfbeta",linear.predictions = FALSE, ggtheme = theme_bw()) ##Graphical illustration that there are no overly influential datapoints affecting the model



#Summary of model
summary(model2) # There is a significant difference in survival between the two treatments, from the figure we can see that the treatment queens lived significantly less long than the control queens (coxph: z = -2.977, p = 0.00291). From summary we can see that the coefficient of C has a negative sign (indicating reduced hazard). We see that the exponentiated coefficient (-0.9192) is 0.3988 (i.e., the hazzards ratio). Therefore being in the control treatment reduces the hazzard by a factor of 0.3988, or 60%
null.model2 <- coxph(Surv(age, censor) ~ 1, data = df7)
summary(null.model2)
anova(model2,null.model2)


model3 <- coxph(Surv(age, censor) ~ treat + (1|individual), data = df7)

AIC(model2, model3)# No point in including a random effect for individual, because it's the same model

#Model 2 is the main model for these tests. From this point we will determine if any additional fixed factors should be included. As described in the supplementary methods, the factors to be tested are: marginal cell width, queen fertility, observed worker egg-laying, observed worker aggression, excess-workers removed, and egg-laying workers removed


#Queen fertility
model4<- coxph(Surv(age, censor) ~ treat*mean_eggs_after_treatment, data = df7)
summary(model4)

AIC(model2, model4)# model4 performs worse


#Test if size has an effect on survival


#Import size measurements of queens. Note we have the thorax width and marginal cell lengths (referred to as wing for convenience) of every queen except for queen 57 (which was destroyed by the workers after she fell into fairyliquid), and 3 of the 4 "extra" colonies (the one extra colony we measured, 62 was to be a replacement colony for colony 57 so was treated as a control for the second half of the experiment. However, we opted to censor her in the end because she had a slightly different treatment regime at the start off the experiment). Although thorax width was measured, it was hard to measure accurately (due to the undefined edge regions), whereas marginal cell length was easy to measure accurately.
df9 <- read_csv("queen_measurements.csv")

df9 <- rename(df9, individual = "Queen No", thorax = "Av thorax", wing = "Av wing") #For renaming tibble column

df10 <- inner_join(df7, df9, by = "individual", copy = FALSE, suffix = c(".x", ".y"))

is.numeric(df10$thorax)

#Quick check that thorax and wing are correlated
(plot_sp_size <- ggplot(df10, aes(x=thorax, y=wing)) +
    geom_point(size=2)+
    geom_smooth(method=lm))

lin.mod.size <-  lm(thorax ~ wing, data = df10)
summary(lin.mod.size)
#The two size measurements are significantly correlated


#Next check that the size measurements affect longevity
#due to missing value, df10 is a smaller dataset so first need to check that the effect of treatment is still significant
(plot_bp.survival2<-ggplot(df10, aes(x = treat, y = age)) +
  geom_boxplot()+
  labs(title="The effect of egg removal on survival in Bombus terrestris females",x="Treatment", y = "Age at death")+
  theme_classic()+
  geom_jitter(width = 0.2))
lin.mod0<-  lm(age ~ treat, data = df10)
summary(lin.mod0)
#results and boxplot indicate that treatment still has a significant effect on age at death

#Check if thorax or wing has any relationship with age
(plot_sp_wing <- ggplot(df10, aes(x=thorax, y=age, shape = treat, colour = treat)) +
    geom_point(size=2)+
    geom_smooth(method=lm))

(plot_sp_wing <- ggplot(df10, aes(x=wing, y=age, shape = treat, colour = treat)) +
    geom_point(size=2)+
    geom_smooth(method=lm))

#check whether there is an effect of thorax size on age
lin.mod1 <-  lm(age ~ thorax*treat, data = df10)
lin.mod2 <-  lm(age ~ thorax, data = df10)
summary(lin.mod1)
summary(lin.mod2)

#there is no significant effect of thorax on age
lin.mod3 <-  lm(age ~ wing*treat, data = df10)
lin.mod4 <-  lm(age ~ wing, data = df10)
summary(lin.mod3)
summary(lin.mod4)
#there is no effect of wing size on age

#as Pierre previously found that marginal cell length (as a proxy for body size) was a significant predictor for survival in workers, try a coxph model with wing as a predictor

#1st look at the effect without wing
model5<- coxph(Surv(age, censor) ~ treat, data = df10)
summary(model5)
#treatment still a significant predictor of survival

#Add wing as a predictor
model6<- coxph(Surv(age, censor) ~ treat + wing, data = df10)
print(model6)
summary(model6)

model7<- coxph(Surv(age, censor) ~ treat * wing, data = df10)
summary(model7)

AIC(model5,model6, model7)
anova(model5,model6)
anova(model5,model7)
#Adding wing as a predictor has no significant impact on the quality of the model and does not change its predictions

#Test effect of wing by itself
model8 <- coxph(Surv(age, censor) ~ wing, data = df10)
summary(model8)
#Wing/Marginal cell width (and hence body size) does not influence longevity


####Make a model where egg production is added as a significant predictor of longevity
model9 <- coxph(Surv(age, censor) ~ treat+mean_eggs_after_treatment, data = df10)
summary(model9)
AIC(model5,model9)
anova(model5,model9)
#Again no significant effect of adding eggs to the model

####Make a model where aggression is accounted for

#Import data
aggression <- read_csv("aggression.csv")%>% 
  group_by(col_no) %>% 
  summarise(total_agg = sum(aggression),
            behavioural_observations = n()) %>% 
  rename(individual = col_no)



df11 <- inner_join(df10,aggression,by="individual")

#Turn aggression into an aggression rate (i.e., number of incidences of aggression relative to the number of times an aggression observation was carried out. This is to account for the fact that colonies that live for longer have had more aggression data collected for them) 
df11 <-   df11 %>% 
  mutate(agg_index = total_agg/behavioural_observations)

#aggression as a fixed effect
model10 <- coxph(Surv(age, censor) ~ treat+agg_index, data = df11)
summary(model10)

#Total aggression, again, total aggression may correlate with age, but the most reasonable explannation for this is that longer lived colonies had higher numbers of observations. Ideally it's best to account for this by turning aggression into a rate as above for model10.
model10b <- coxph(Surv(age, censor) ~ treat+total_agg, data = df11)
summary(model10b)

AIC(model5,model10,model10b)
anova(model5,model10)
anova(model5,model10b)#No significant effect of aggression

####Make a model where the number of workers removed is accounted for

numbers <- read_csv("numbers.csv") %>% 
  group_by(col_no) %>% 
  summarise(w_add = sum(w_add,na.rm = TRUE),
            removed_aggression = sum(removed_aggression, na.rm = TRUE),
            removed_egglaying = sum(removed_egglaying, na.rm = TRUE),
            removed_other = sum (removed_other, na.rm = TRUE),
            removed_total = sum(removed_total, na.rm = TRUE),
            w_dead = sum (w_dead, na.rm = TRUE))%>% 
  ungroup() %>% 
  rename(individual = col_no)

df12 <- inner_join(df11,numbers,"individual")




###egg laying workers removed (i.e. workers that were removed because they were seen laying eggs during the observation periods)
model11 <- coxph(Surv(age, censor) ~ treat + removed_egglaying, data = df12)
summary(model11)
AIC(model5,model11)
anova(model5,model11) #These results show that adding "workers removed due to egg laying" as an additional fixed effect improves the fit of the model, however this is a total rather than a rate so this might be due to longer lived colonies having more observations

#Turn the number of workers removed due to egg laying into a rate rather than a totall (i.e., number of workers removed for observed egg laying behaviour relative to the number of times observation was carried out - observations were carried out every four days on the same days as aggression and worker egg laying were measured per se, so we can use the number of observations as calculated from the aggression datasheet above. We do this to account for the fact that colonies that live for longer have had more aggression data collected for them) 
df12 <-   df12 %>% 
  mutate(removed_egglaying_index = removed_egglaying/behavioural_observations)

model11b <- coxph(Surv(age, censor) ~ treat + removed_egglaying_index, data = df12)
summary(model11b)
anova(model5,model11b)
AIC(model5,model11,model11b) #Note that using totals improves the AIC, while using rate does not. However, my intuition is that rate is a better measure for this purpose because it more accurately captures what we want to control for.

#Check the two measures using a scattterplot
(plot_sp_egg_removed <- ggplot(df12, aes(x=removed_egglaying, y=age)) +
    geom_point(size=2)+
    geom_smooth(method = lm)+
    labs(x="Number of workers removed", y = "Longevity")+
    theme_classic())

(plot_sp_egglaying_removed_index <- ggplot(df12, aes(x=removed_egglaying_index, y=age)) +
    geom_point(size=2)+
    geom_smooth(method=lm)+
    labs(x="Rate of worker removal", y = "Longevity")+
    theme_classic()) 
#These plots again highlight the problem with using totals rather than rates. The total number of workers removed increases with age (because long lived colonies have more removed from them); whereas the rate of worker removal is not significantly associated with age.

#Test that new model still fulfills cox assumptions
Test2<-cox.zph(model11b)
Test2 #There is no evidence that the model fails the assumption of hazards being constant

ggcoxdiagnostics(model11b, type = "dfbeta",linear.predictions = FALSE, ggtheme = theme_bw()) #the residuals are huge for the egg laying worker removal rate data, probably an additional indication that these data have no explannatory value.


#Test whether including an interaction term between egg laying workers removed and treatment improves the fit
model11c <- coxph(Surv(age, censor) ~ treat * removed_egglaying_index, data = df12)
summary(model11c)
anova(model11,model11c)
AIC(model5,model11,model11b,model11c) #The best model excludes the interaction term

#As there was no effect of worker removal rate on longevity, it can be removed as a factor in the model.

###aggressive workers removed
#Worth testing if workers removed for aggression reasons had an effect on the survivorship curve

model12 <- coxph(Surv(age, censor) ~ treat + removed_aggression, data = df12)
summary(model12) #No effect of the total workers for aggression reasons (as above totals are spurious, and it is better to convert workers removed into a rate)

#Calculate rate of workers removed due to aggression
df12 <-   df12 %>% 
  mutate(removed_aggression_index = removed_aggression/behavioural_observations)

#Make model with removed_aggression_index as a predictor
model13 <- coxph(Surv(age, censor) ~ treat + removed_aggression_index, data = df12)
summary(model13)
AIC(model5,model12,model13) #Model12 is slightly better (but only slightly), as discussed this model makes less sense so better to stick to model13 to control for workers removed due to aggression. As this performs worse than the simple model, stick to the simple model.
anova(model11,model12)

###total workers removed as an additional fixed effect
#Lots of workers were removed from the colonies for reasons other than egglaying and aggression - so instead of looking at these factors seperately, look at the total number of workers removed for all reasons.
model13 <- coxph(Surv(age, censor) ~ treat + removed_total, data = df12)
summary(model13)

#Turn total workers removed into an index
df12 <-   df12 %>% 
  mutate(removed_total_index = removed_total/behavioural_observations)

#Run model
model14 <- coxph(Surv(age, censor) ~ treat + removed_total_index, data = df12)
summary(model14)


AIC(model5,model13,model14)
anova(model5,model14)
#The rate of total workers removed provides a better fit to the data than the models that includes total_workers_removed as a total, and the model that only included treatment.



###excess workers removed (i.e. workers removed because there were more than 20 workers in the colony) fixed effect
#Excess workers removed was the main reason workers were removed from the colony. As the rate of total workers removed (removed_total_index) was a significant fixed effect, and as egg-laying, aggression, removed_egglaying_index, removed_aggression_index were not significant effects, then it seems reasonable to test whether excess workers removed (as both a rate and total) are worth including in the model. As these are highly co-linear with total workers removed it does not make sense to have both fixed effects in the models (i.e., to avoid pseudoreplication).
model14 <- coxph(Surv(age, censor) ~ treat + removed_other, data = df12)
summary(model14) #Removed_other as a total is not significant

#Turn removed_other into an index.
df12 <-   df12 %>% 
  mutate(removed_other_index = removed_other/behavioural_observations)

#Create model
model15 <- coxph(Surv(age, censor) ~ treat + removed_other_index, data = df12)
summary(model15) #Both treatment and removed_other_index are significant
AIC(model5,model14,model15) #Model15 performs much better than either the simple model (model 5)

#Check it the two fixed effects have an interaction
model16 <- coxph(Surv(age, censor) ~ treat * removed_other_index, data = df12)
summary(model16) #There is no significant interaction between the two variables

AIC(model5,model14,model15,model16) #Model15 performs much better than the simple model (model 5) and the totals model (model 14), it performs slightly worse than the interaction effects model
anova(model16,model15) #There is no significant difference between the predictions of the two models. Therefore opt for the simple model.

#Check the two measures using a scattterplot
(plot_sp_removed_other <- ggplot(df12, aes(x=removed_other, y=age)) +
    geom_point(size=2)+
    geom_smooth(method = lm)+
    labs(x="Number of workers removed", y = "Longevity")+
    theme_classic())

#Check the two measures using a scattterplot
(plot_sp_removed_other_index <- ggplot(df12, aes(x=removed_other_index, y=age)) +
    geom_point(size=2)+
    geom_smooth(method = lm)+
    labs(x="Rate of workers removed (per event)", y = "Longevity")+
    theme_classic())

#The scatterplot shows the result quite clearly, longer lived colonies had a lower rate of worker removal than short-lived colonies. This could imply that increased disturbance (i.e., more workers removed) could have led to a decrease in colony longevity. It could also imply that longer lived colonies had a gradual reduction in worker output, meaning that fewer of them were removed with time as they lived longer. Either way, removed_other_index might be something we should include in the final model.

#Test assumption that regression coefficients are constant over time
Test_model15<-cox.zph(model15)

Test_model15 ##we can accept that the regression coefficients are constant over time (p>0.05) and that the slope is not significantly different from zero, for all treatment levels

ggcoxzph(Test_model15) ##Graphical illustration of the fact there is no pattern of hazzard with time, i.e. hazards are proportional, not non-proportional

ggcoxdiagnostics(model15, type = "dfbeta",linear.predictions = FALSE, ggtheme = theme_bw(), ox.scale = "observation.id") ##Graphical illustration of residuals. One of the treatment residuals has a very high value compared to the others, so this is a potential outlier. 

#Check the raw values of the residuals
residuals(model15, type = "dfbeta") #ResidualID 16 has a value of 0.11 which is quite high. This is a potential outlier so test if removing it affects the model

#Create new dataframe with outlier removed
df12b <- df12 %>% 
  filter(individual != "16")

#Recreate model
model15b <- coxph(Surv(age, censor) ~ treat + removed_other_index, data = df12b)
summary(model15b) #Summary shows that the effects are stronger (if anything) without the outlier

#Retest with ggcoxdiagnostics
ggcoxdiagnostics(model15b, type = "dfbeta",linear.predictions = FALSE, ggtheme = theme_bw(), ox.scale = "observation.id") ##Graphical illustration that there are no overly influential datapoints affecting the model - we can be confident with or without the outlier excluded.

### null model
model17 <- coxph(Surv(age, censor) ~ 1, data = df12)
AIC(model15,model17)
anova(model15,model17)

#Use LRT to test whether there is a significant effect (summary function for Cox's models already makes this comparison with the null model and gives the same result as for anova function)
