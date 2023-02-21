##DavidHCollins
#Created: 31 August 2020

#Aim

#The aim of this draft script is to fecundity data from the main Bombus experiment.

#First you must always reset R
rm(list=ls())

#make these packages (and associated functions) available to use in the script

library(readr) #read in files
library(tidyverse)#ggplot and dplyr
library(rstatix)#pipe_friendly_stats


#here are some useful shortcuts to get to know
#ctrl+a highlights the whole script
#ctrl+enter runs the highlighted script or just the selected line of code if nothing is highlighted
#ctrl+1 moves cursor to the console
#ctrl+2 moves it back to the script editor


#individual level data for fecundity and cell counts

setwd("Documents/UEA Work documents/Post Doc_LHT/Objective 1 Experiment 1/Analysis")

cells <- read_csv("cells.csv")
View(cells)
df1<-cells


#do you have the right data?
glimpse(df1)


######################### DATAFRAMES CONTAINING VARIOUS SUMMARIES OF THE EGG COUNT DATA #########################

# Make a new dataframe which excludes all of the egg counts that occurred after each queens death (i.e. Only use data that says 'T' in the queenright column)
df2 <- filter(df1, queenright == "TRUE")

#And a dataframe which only includes only the egg counts from colonies where the queen has died i.e. data that says 'F' in the queenright column 
df3 <- filter(df1, queenright == "FALSE")

# Make a new dataframe which excludes all of the egg counts that occurred before treatment (i.e. Only use data that says 'T' in the after_treatment column); and only includes the egg counts from before each queen death
df4 <- filter(df1, queenright == "TRUE", after_treatment == "TRUE")

#And a dataframe which only includes the egg counts from colonies before treatment had started i.e. data that says 'F' in the after_treatment column; and only includes the egg counts from before each queen death
df5 <- filter(df1, queenright == "TRUE", after_treatment == "FALSE")

# Make a new dataframe which excludes all of the egg counts that occurred from count 14 onwards when the first worker egg laying was observed (i.e. Only use data that says 'T' in the before_worker_eggs column), and which excludes egg counts from before treatment started (i.e. data that says 'T' in the queenright column), and which excludes egg counts from after the queen died (i.e. data that says'T' in the after_treatment column)
df6 <- filter(df1, queenright == "TRUE", after_treatment == "TRUE", before_worker_eggs == "TRUE")



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

#Import survival data

survival <- read_csv("survival.csv")

#Turn treat into a factor and reorder so R comes first
survival$treat <- as.factor(survival$treat)
is.factor(survival$treat)
survival$treat <- factor(survival$treat, levels = c("R", "C"))


# Add the metrics from the three summary stats dataframes to new df7 and survival data

df7 <- inner_join(survival,tot.summary.stats,by="individual")
df7 <- inner_join(df7,after_treat.summary.stats,by="individual")
df7 <- inner_join(df7,before_treat.summary.stats,by="individual")
df7 <- inner_join(df7,before_worker_eggs.summary.stats,by="individual")

# Make a new dataframe from df7 which excludes all of the censored individuals which were removed either for RNA extraction or accidently killed (i.e. exclude all entries that have 0 in the censor column)

df8 <- filter(df7, censor == "1")

#Relationship between queen fecundity and longevity before the first worker egg laying was observed. Note. See Kramer et al. 2015. PLOSONE for further validation of this approach to the relationships (i.e. that egg production is the response variable, treatment the predictor, and longevity a covariate)




## Scatterplot of the relationship (before worker egg laying was detected)

(scatterplot1<- ggplot(df8,aes(x=age, y= mean_eggs_before_worker_egg_laying, colour = treat)) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  labs(title="The relationship between longevity and average egg production in Bombus terrestris queens
        reared under two egg removal treatments before workers were observed laying eggs",x="Longevity (days)", y = "Average egg production per count", colour = 'Treatment') +
  theme_classic())



(scatterplot2<- ggplot(df8,aes(x=age, y= total_eggs_before_worker_egg_laying, colour = treat)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(title="The relationship between longevity and average egg production in Bombus terrestris queens
        reared under two egg removal treatments before workers were observed laying eggs",x="Longevity (days)", y = "Total egg production per count", colour = 'Treatment') +
  theme_classic())


## Linear model of the relationship

model1 <- lm(mean_eggs_before_worker_egg_laying ~ treat*age, data = df8)


#graphic test of assumptions
autoplot(model1)

#Shapiro-wilks test of normality of residuals
model.metrics <- augment(model1) %>%
  select(-.hat, -.sigma, -.fitted, -.se.fit) # Remove details
head(model.metrics, 3)
shapiro_test(model.metrics$.resid)
#Shows the residuals were normal (Shapiro-Wilk, w = 0.975, p = 0.410)

model.metrics %>% levene_test(.resid ~ treat)


anova(model1) 

df8 %>%
  group_by(treat) %>%
  anova_test(mean_eggs_before_worker_egg_laying ~ age) #Anova test of longevity-fecundity after accounting for treatment

summary(model1) #Shows that there was a highly significiant effect of treatment on the average number of eggs produced, however there was no evidence of a relationship between age and average egg production in either treatment, and there was no evidence of an interaction between treatment and age. While not signficant, it is of interest that the C group had a negative relationship between longevity and average egg production, and the R group had a positive relationship. We would have expected this under the H1 scenario (i.e. queen costs of reproduction are present but latent)




###Add the model to the figure

# Make 10 new 'age' values at which we need some predictions. These will go from 17 to 157 (the range of values we have for age in the dataset). To these add all of the values in the treat 'column' so there is a full factorial representation of the randomly generated numbers and the values for treat. To do this treat needs to be a factor 
df8$treat<-as.factor(df8$treat)
new3.x<-expand.grid(age= seq (from = 17, to = 157, length.out = 10),
                    treat = levels(df8$treat))
#Did it work?
head (new3.x)

#use predict() to generate new y values
new3.y <- predict(model1,newdata=new3.x,interval = 'confidence')

#Did it work?
head (new3.y)

#housekeeping: collect the new x and new y for plotting
addThese3<-data.frame(new3.x,new3.y)
#change name from fit to match the original data (in this case EGGS)
addThese3<-rename(addThese3, mean_eggs_before_worker_egg_laying = fit)
#check it worked
head(addThese3)

#remake figure and add lines and CI
(scatterplot3<- 
  ggplot (df8,aes(x=age, y= mean_eggs_before_worker_egg_laying, fill = treat)) +
  geom_point(aes(), size = 2, shape = 21, colour = "black")+
  geom_smooth(data = addThese3,
              aes(ymin = lwr, ymax = upr, fill = treat, colour = treat),
              stat = 'identity',show.legend = FALSE) +
    theme(plot.title.position = "plot", 
          text = element_text(size=10),
          axis.text.x = element_text(),
          axis.line= element_line(),
          legend.position = "none",
          panel.border = element_blank(), 
          panel.background = element_blank(), 
          panel.grid = element_blank(), 
          panel.spacing.x = unit(0.3,"line"))+
      labs(title = NULL, x="Queen longevity (days)", y = "During treatment mean queen fertility (eggs produced)")+
  scale_colour_manual (values = c(R = "orangered", C = "gray45")) +
  scale_fill_manual (values = c(R = "orange", C = "gray68")))

ggsave("Longevity_fecundity_relationship.png", scatterplot3, width = 16, height = 10, units = "cm")
ggsave("Longevity_fecundity_relationship.svg", scatterplot3, width = 16, height = 10, units = "cm" )

### Relationship after worker egg laying was detected but before queens died

model2 <- lm(mean_eggs_after_treatment ~ treat*age, data = df8)

autoplot(model2)

anova(model2) 

df8 %>%
  group_by(treat) %>%
  anova_test(mean_eggs_after_treatment ~ age) #Anova test of longevity-fecundity after accounting for treatment

summary(model2) 

