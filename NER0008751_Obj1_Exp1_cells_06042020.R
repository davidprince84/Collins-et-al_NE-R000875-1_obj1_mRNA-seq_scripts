#DavidHCollins
#Created: 28 August 2020



#This script covers the effect of treatment on average cell count over time throughout the course of the experiment. Censoring queens should not matter for these analyses, because these analyses are about using the data for every individual up until the point it died/was censored). It complements the scripts dealing with egg count, and uses the same models as they do.


#First you must always reset R
rm(list=ls())

#make these packages (and associated functions) available to use in the script
library(rlang)
library(readr) #read in files
library(tidyverse)#ggplot and dplyr
library(car) # for levenes test
library(MASS) # for negative binomial glms
library(AER) # for other glms on count data
library(MuMIn) # for analysing how well the model fits the data
library(stats) #for scaling data
library(bbmle) #produces AIC values which allows you to compare multiple models
library(lme4) ## for glmer
library(effects) ## All effects function which gives confidence intervals of actual values when using a glm
library(jtools) ## summ function which is a much nicer easier to read version of summary with some better info

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

#Make treat a factor and reorder it so R is the first level
df1$treat <- as.factor(df1$treat)
is.factor(df1$treat)
df1$treat <- factor(df1$treat, levels = c("R", "C"))

#Make day a factor
df1$day <- as.factor(df1$day)
is.factor(df1$day)

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

#Before modeling, rescale the variables in case the models are unidentifiable


df4 <- df4 %>% 
  mutate(treat.logic = ifelse(treat=="R",T,F),
         treat.scale = scale(treat.logic),
         day.scale= scale(day),
         eggs.scale=scale(eggs))

#Use these to count highest values of cells and eggs per cell - relevant to y axis
#x <- df2 %>% arrange(desc(cells)) 
#x <- df2 %>% arrange(desc(eggs_per_cell)) 

#Number of egg cells (note, only up to count 68 is included as there was only 1 treatment colony after that count)
(plot_bp.cells<-df2 %>% 
    filter(queenright == T) %>% 
    group_by(day, treat) %>% 
    filter(count < 68)%>% 
    filter(cells<20)%>% 
    ungroup()%>% 
    ggplot(aes(x = day, y = cells)) +
    geom_boxplot(aes(fill=treat), outlier.shape=NA)+
    theme(plot.title.position = "plot", 
          text = element_text(size=8),
          axis.text.x = element_text(),
          axis.line= element_line(),
          legend.position = "none",
          panel.border = element_blank(), 
          panel.background = element_blank(), 
          panel.grid = element_blank(), 
          panel.spacing.x = unit(0.3,"line"))+ 
    scale_fill_manual(values = c("orange", "gray68"))+
    labs(title=NULL,
         x ="Day", y = "Number of egg cells", fill="Treatment")+ 
    scale_x_discrete(breaks=seq(5,155,5))+
    scale_y_continuous(breaks=seq(0,20,5))+
    geom_vline(xintercept=2.5, size = 0.5, linetype='dashed', col = 'red')+
    geom_vline(xintercept=13.5, size = 0.5, linetype='dashed', col = 'black'))

#Number of eggs per cell
(plot_bp.eggspercell<-df2 %>% 
    filter(queenright == T) %>% 
    group_by(day, treat) %>% 
    filter(count < 68)%>% 
    filter(eggs_per_cell < 21)%>% 
    ungroup()%>% 
    ggplot(aes(x = day, y = eggs_per_cell)) +
    geom_boxplot(aes(fill=treat), outlier.shape=NA)+
    theme(plot.title.position = "plot", 
          text = element_text(size=8),
          axis.text.x = element_text(),
          axis.line= element_line(),
          legend.position = "none",
          panel.border = element_blank(), 
          panel.background = element_blank(), 
          panel.grid = element_blank(), 
          panel.spacing.x = unit(0.3,"line"))+ 
    scale_fill_manual(values = c("orange", "gray68"))+
    labs(title=NULL,
         x ="Day", y = "Number of eggs per cell", fill="Treatment")+ 
    scale_x_discrete(breaks=seq(5,155,5))+
    scale_y_continuous(breaks=seq(0,20,5))+
    geom_vline(xintercept=2.5, size = 0.5, linetype='dashed', col = 'red')+
    geom_vline(xintercept=13.5, size = 0.5, linetype='dashed', col = 'black'))

ggsave("Figure_Scells.png", plot_bp.cells, width = 16, height = 10, units = "cm")
ggsave("Figure_Seggspercell.png", plot_bp.eggspercell, width = 16, height = 10, units = "cm")


####For cells: 

#Make null models

null.p.model1 <- glmer( cells ~ 1 + (1|ID), family=poisson, data=df4)
summary(null.p.model1)

##test for overdispersion using the following useful function
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

overdisp_fun(null.p.model1)
#Data are very overdispersed, so will use the negative binomial regression instead

#Start with the best null model and full model for eggs, but with cells instead

null.nb.model1 <- glmer.nb(cells ~ 1 +(1|ID)+(1|day), family=poisson, data=df4)
nb.model1 <-  glmer.nb(cells ~ treat*poly(day,2) + (1|ID) + (1|day), data=df4)
nb.model2 <-  glmer.nb(cells ~ treat*poly(day,2) + (1|ID), data=df4)
nb.model3 <-  glmer.nb(cells ~ treat*poly(day,2) + (1|day), data=df4)

AIC(nb.model1,nb.model2,nb.model3)
anova(nb.model1,nb.model2) #model 1 is the best of these models

anova(nb.model1,null.nb.model1) #model 1 is better than the null model

summary(nb.model1) #Shows that the effect of treatment and time on cell count is highly significant

plot(allEffects(nb.model1))
summ(nb.model1,scale=TRUE,digits=3, confint=T)



