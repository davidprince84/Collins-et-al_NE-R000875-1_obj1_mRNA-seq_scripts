#DavidHCollins
#Created: 15 April 2020

#Aim

#The aim of this draft script is to analyse the feccundity data from the main Bombus experiment. This first script covers figure 1a and 1b, and the corresponding model fitting for each figure. For the models the analysis will be used to analyse the total egg count before and after treatment began for both the R and C colonies. Note we even include the censored individuals in these analyses because the censoring took place late on in the exeriment, and should not affect the egg laying rates in any case (censoring *does* matter for the longevity data however).


#First you must always reset R
rm(list=ls())

#make these packages (and associated functions) available to use in the script

library(readr) #read in files
library(tidyverse)#ggplot and dplyr
library(car) # for levenes test
library(MASS) # for negative binomial glms
library(AER) # for other glms on count data and dispersison test
library(MuMIn) # for analysing how well the model fits the data
library(bbmle) #produces AIC values which allows you to compare multiple models
library(lme4) ## for glmer
library(effects) ## All effects function which gives confidence intervals of actual values when using a glm
library(jtools) ## summ function which is a much nicer easier to read version of summary with some better info
library(ggsignif) #For adding significance bars to plot
library(svglite)

#here are some useful shortcuts to get to know
#ctrl+a highlights the whole script
#ctrl+enter runs the highlighted script or just the selected line of code if nothing is highlighted
#ctrl+1 moves cursor to the console
#ctrl+2 moves it back to the script editor


#individual level data for fecundity and cell counts

setwd("Documents/UEA Work documents/Post Doc_LHT/Objective 1 Experiment 1/Analysis")

#If on HP computer
#setwd("C:/Users/dhcollins500/OneDrive - University of East Anglia/Documents/Post Doc_LHT/Objective 1 Experiment 1/Analysis/")
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

#Make treat a factor and reorder it so R is the first level
df1$treat <- as.factor(df1$treat)
is.factor(df1$treat)
df1$treat <- factor(df1$treat, levels = c("R", "C"))


######################### DATAFRAMES CONTAINING VARIOUS SUMMARIES OF THE EGG COUNT DATA #########################

# Make a new dataframe which excludes all of the egg counts that occurred after each queens death (i.e. Only use data that says 'T' in the queenright column)
df2 <- filter(df1, queenright == "TRUE")

#And a dataframe which only includes only the egg counts from colonies where the queen has died i.e. data that says 'F' in the queenright column 
df3 <- filter(df1, queenright == "FALSE")

# Make a new dataframe which excludes all of the egg counts that occurred before treatment (i.e. Only use data that says 'T' in the after_treatment column); and only includes the egg counts from before each queen death
df4 <- filter(df1, queenright == "TRUE", after_treatment == "TRUE")

#And a dataframe which only includes the egg counts from colonies before treatment had started i.e. data that says 'F' in the after_treatment column; and only includes the egg counts from before each queen death
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


#Import survival data (use survival, as this is the data which censors the correct queens) and add these summary tables to it. Survival censors all of the queens removed for RNA extraction, and the following queens: queens 57 (censored due to workers killing her after she fell in fairy liquid), 62 (replacement of 57, but had no measurements taken of her fecundity at the beginning of the experiment, and she spent several weeks with many more workers than 20).


survival <- read_csv("survival.csv")
View(survival)

#Turn treat into a factor and reorder so R comes first
survival$treat <- as.factor(survival$treat)
is.factor(survival$treat)
survival$treat <- factor(survival$treat, levels = c("R", "C"))

# Add the metrics from the three summary stats dataframes to df3

df7 <- merge(survival,tot.summary.stats,by="individual")
df7 <- merge(df7,after_treat.summary.stats,by="individual")
df7 <- merge(df7,before_treat.summary.stats,by="individual")
df7 <- merge(df7,before_worker_eggs.summary.stats,by="individual")

# Make a new dataframe from survival4 which excludes all of the censored individuals which were removed either for RNA extraction or accidently killed (i.e. exclude all entries that have 0 in the censor column)

df8 <- filter(df7, censor == "1")

##################################################################################################################


###########1.0 Effect of Treatment on Longevity and Fecundity############

#These figures will be figures 1a and 1b

## boxplot: effect of treatment on total fecundity before treatment started
plot_bp.egg.before.treat <- ggplot(df7,aes(x=treat, y=total_eggs_before_treatment)) +
  geom_boxplot()+
  labs(title="The effect of egg removal on egg production before treatment began in Bombus terrestris females",x="Treatment", y = "Total number of eggs")+
  theme_classic()+
  geom_jitter(width = 0.2)

plot_bp.egg.before.treat

# 1.1 boxplot: effect of treatment on total egg production after treatment began when R is compared to C across all counts
plot_bp.egg.before.worker.eggs <- ggplot(df7,aes(x=treat, y=total_eggs_before_worker_egg_laying)) +
  geom_boxplot()+
  labs(title="The effect of egg removal on mean egg production after treatment began
       in Bombus terrestris females (before workers were recorded laying)",
       x="Treatment", 
       y = "Total number of eggs")+
  theme_classic()+
  geom_jitter(width = 0.2)

plot_bp.egg.before.worker.eggs


#ANOVA tests
#We want to do three difference tests, one that shows the effect of treatment on TOTAL fecundity before treatment, one that shows the effect on TOTAL fecundity after treatment, and a two way test that measures the effect of treatment, and treatment timing on mean number of eggs. In all cases we are only looking at the data from before the workers started to lay eggs. First I thought to use ANOVA for these tests. However because we are analysing count data, then glms with poission distribution (if data not overdispersed) or glm with. Note we use df3 not df9 because we want to include the censors (which were only censored *after* worker egg laying was recorded, so censoring was not relevent to these analyses)

##ANOVA 1: 

#Build model
model1 <- lm(total_eggs_before_treatment ~ treat, data = df7)
#test assumptions
autoplot(model1 , smooth.colour = NA) #these data aren't perfect, have done further tests below
# use anova() to test if there is an effect of Treatment
anova(model1 )
#use summary() to test what the nature of the effect of treatment is, and see if it makes sense with the boxplot itself
summary(model1)

#Are data normal?
resid.model1<-resid(model1)
shapiro.test(resid.model1) # Data are normal
#Is variance homogenous
bartlett.test(total_eggs_before_treatment ~ treat, df7)
leveneTest(total_eggs_before_treatment ~ treat, data = df7) # Variance is homogenous

## ANOVA 2

#Build model
model2 <- lm(total_eggs_before_worker_egg_laying ~ treat, data = df7)
#test assumptions
autoplot(model2 , smooth.colour = NA) #these data aren't perfect, have done further tests below
# use anova() to test if there is an effect of Treatment
anova(model2 )
#use summary() to test what the nature of the effect of treatment is, and see if it makes sense with the boxplot itself
summary(model2)

#Are data normal?
resid.model2<-resid(model2)
shapiro.test(resid.model2) # Data are normal
#Is variance homogenous
bartlett.test(total_eggs_before_worker_egg_laying ~ treat, df7)
leveneTest(total_eggs_before_worker_egg_laying ~ treat, data = df7) # Variance is not homogenous




#This might be better modelled using poisson or negative binomial distribution (ANOVA not really appropriate for count data)

model3 <- glm(total_eggs_before_treatment ~ treat,family="poisson",data = df7)
dispersiontest(model3,trafo=1) #data are overdispersed (c = 6.94, when it should be zero)

#As data are overdispersed a negative binomial glm is probably more appropriate
model4 <- glm.nb(total_eggs_before_treatment ~ treat,data = df7)
summary(model4)

#From these models the critical 5% qchisq is (this function calculates the critical value based on the degrees of freedom of a given model)
qchisq(0.95, df.residual(model3))
qchisq(0.95, df.residual(model4))

# And the actual deviance values are
deviance(model3)
deviance(model4)
#We can see that model 3 is above the critical value of 89.4 (it is 573.8), and model 4 is below the critical value (it is 73.8). Model 3 is therefore not a good fit to the data, whereas model 4 is.

#Make a similar model for total eggs before worker egg laying started



model5 <- glm.nb(total_eggs_before_worker_egg_laying ~ treat,data = df7)
summary(model5) #again deviance is 72.723 and is below the critical chisq value of 89.4
plot(model5) #plots look okay
summ(model5, digits=3)
summary(allEffects(model5))#gives confidence intervals of the actual values

#Use round to find percentage change between treatment and reference
round((model5$family$linkinv(0.71487)-1)*100)
#Answer is a 104% increase in egg number.

#We can report for the paper: 

#A negative binomial glm was used to examine the effects of treatment on total egg production after treatment began and before workers were recorded laying eggs. Treatment was found to have a significant effect on egg production (b = 0.715, SEb = 0.09, z = 7.901, p <0.001). Overall the R colonies (mean [95% CI] = 574 [506.5, 651.5] eggs) produced 104% more eggs than the C colonies (mean [95% CI] = 281 eggs [246.9, 319.6]).




#### Model to determine the effect of treatment on egg production before and after treatment

#These look at average rather than total egg laying
df9 <- df7 %>% 
  pivot_longer(c('mean_eggs_before_treatment','mean_eggs_before_worker_egg_laying'),names_to = "x", values_to = "eggs")%>% 
  mutate(treatment_started=ifelse(x=="mean_eggs_before_treatment","Before","After")) %>% 
  dplyr::select(1:7,treatment_started,eggs)

df9$treatment_started <- factor(df9$treatment_started, levels = c("Before", "After"))

# Calculate mean and SD before and during treatment for each treatment

df9 %>% 
  group_by(treat,treatment_started) %>% 
  summarise(mean = mean(eggs),
            sd = sd (eggs))


## boxplot: effect of treatment on total fecundity before treatment started, this will now be figure 1 for the paper, rather than the above boxplots

# New facet label names for supp variable
supp.labs <- c("Before treatment", "During treatment")
names(supp.labs) <- c("Before", "After")

(Av_bp <- ggplot(df9,aes(x=treat, y=eggs, fill = treat)) +
  facet_wrap(~treatment_started, labeller = labeller(treatment_started=supp.labs))+
  geom_boxplot(outlier.shape=NA)+
  labs(x="Treatment", y = "Mean queen fertility (eggs produced)")+
  scale_fill_manual(values = c("orange", "gray68"))+ 
  ylim(0,110)+
  geom_jitter(width = 0.2)+
  theme(plot.title.position = "plot", 
        text = element_text(size=11),
        axis.text.x = element_text(),
        axis.line= element_line(),
        legend.position = "none",
        panel.border = element_blank(), 
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.spacing.x = unit(0.3,"line"))+ 
  geom_signif(comparisons = list(c("R", "C")),
              map_signif_level = TRUE, textsize=4))

#Make one for poster 

(Av_bp2 <- ggplot(df9,aes(x=treat, y=eggs, fill = treat)) +
    facet_wrap(~treatment_started, labeller = labeller(treatment_started=supp.labs))+
    geom_boxplot(outlier.shape=NA)+
    labs(x="Treatment", y = "Mean queen fertility (eggs produced)")+
    scale_fill_manual(values = c("#FC8E0B", "gray68"))+ 
    ylim(0,110)+
    geom_jitter(width = 0.2)+
    theme(plot.title.position = "plot", 
          text = element_text(size=8),
          axis.text.x = element_text(),
          axis.line= element_line(),
          legend.position = "none",
          panel.border = element_blank(), 
          panel.background = element_blank(), 
          panel.grid = element_blank(), 
          panel.spacing.x = unit(0.3,"line"))+ 
    geom_signif(comparisons = list(c("R", "C")),
                map_signif_level = TRUE, textsize=4))

#FC8E0B is my preferred poster colour

#Note the added significance bars calculate significance based on a wilcox test but this is fine for display purposes. When editing this figure be sure to remove the fullstop from the NS

ggsave("Treatment_egg_boxplot.png", Av_bp, width = 8, height = 10, units = "cm")
ggsave("Treatment_egg_boxplot.svg", Av_bp, width = 16, height = 6, units = "cm") #For making up multipanel figure in the paper

ggsave("Treatment_egg_boxplot2.png", Av_bp2, width = 16, height = 10, units = "cm")
ggsave("Treatment_egg_boxplot2.svg", Av_bp2, width = 16, height = 10, units = "cm")

##########################################################################################


model6 <- lmer(eggs ~ treat*treatment_started + (1|individual),data = df9)
summary(model6)
plot(allEffects(model6))
summ(model6)
summary(allEffects(model6))

model7 <- glm.nb(eggs ~ treat*treatment_started + (1|individual),data = df9)
summary(model7)
plot(allEffects(model7))
summ(model7, digits=3)
summary(allEffects(model7))#gives confidence intervals of the actual values

AIC(model6,model7) #model 7 is better


#Compare model 7 with the intercept only model
model8 <- glm.nb(eggs ~ 1 + (1|individual),data = df9)

AIC(model7, model8)
anova(model7,model8) #Seems the null model comparison is not possible using likelihood ratio test

#A negative binomial glmm was used to examine the effects of treatment on average daily egg production before and after treatment began and before workers were recorded laying eggs.There was no difference between the two treatment groups before egg production had begun (b = 0.057, SEb = 0.099, z = 0.578, p = 0.563). Following the start of treatment both treatment groups showed a significant increase in egg laying rate (b = 0.243, SEb = 0.099, z= 2.457, p = 0.014) and there was a significant interaction between the egg production before and after treatment had begun and the treatment given (b = 0.650, SEb = 0.136, z = 4.766, p < 0.0001). Overall before treatment the R colonies (mean [95% CI] = 21.51 [18.76, 24.68] eggs) had the same average daily egg production rate as the C colonies (mean [95% CI] = 20.34 eggs [17.66, 23.37]), and after treatment the R colonies (mean [95% CI] = 52.54 [46.33, 59.58]) produced eggs at an overall rate of approximately twice that of the C colonies (mean [95% CI] = 25.89 [22.60, 29.66].

