#DavidHCollins
#Created: 1 July 2020



#This script covers figure 1c and the effect of treatment on average egg count over time throughout the course of the experiment. Censoring queens should not matter for these analyses, because these analyses are about using the data for every individual up until the point it died/was censored)


#First you must always reset R
rm(list=ls())

#make these packages (and associated functions) available to use in the script

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
#packageVersion("xxxx") tells you the version of the package xxxx you are currently using or sessionInfo() for all packages

#individual level data for fecundity and cell counts

setwd("Documents/UEA Work documents/Post Doc_LHT/Objective 1 Experiment 1/Analysis")
cells <- read_csv("cells.csv")
View(cells)
df1<-cells


#do you have the right data?
glimpse(df1)

#Turn treat into a factor and reorder so R comes first
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

#Before modeling, rescale the variables in case the models are unidentifiable


df4 <- df4 %>% 
  mutate(treat.logic = ifelse(treat=="Treatment",T,F),
         treat.scale = scale(treat.logic),
         day.scale= scale(day),
         eggs.scale=scale(eggs))


# Find the peak egg laying for each treatment
se <- function(x, ...) sqrt(var(x, ...)/length(x))


mean_se_df <- df4 %>% 
  group_by(treat,day) %>% 
  summarise(mean_day = mean(eggs),
            SD_day = sd(eggs),
            SE_day = se(eggs),
            n_day = n())%>%
  mutate(lower.ci = mean_day - qt(1 - (0.05 / 2), n_day - 1) * SE_day,
         upper.ci = mean_day + qt(1 - (0.05 / 2), n_day - 1) * SE_day)


table2<-dat %>%
  group_by(type) %>%
  summarise(N.num=n(),
            mean.num = mean(dat$num),
            sd.num = sd(dat$num),
            "Percent"=n_perc(num > 0)) %>%
  mutate(se.num = sd.num / sqrt(N.num),
         lower.ci = 100*(mean.num - qt(1 - (0.05 / 2), N.num - 1) * se.num),
         upper.ci = 100*(mean.num + qt(1 - (0.05 / 2), N.num - 1) * se.num))


# First make the figure, this will be figure 1b (note counts more than 68 are filtered because after that count there was only 1 colony left, at which point box plots become meaningless

df2$day.factor<-as.factor(df2$day)
is.factor(df2$day.factor)
is.factor(df2)

(plot_bp.time.eggs<-df2 %>% 
  filter(queenright == T) %>% 
  group_by(day.factor, treat) %>% 
  filter(n() > 1)%>% 
  ungroup()%>% 
  ggplot(aes(x = day.factor, y = eggs)) +
  geom_boxplot(aes(fill=treat), outlier.shape=NA)+
  theme(plot.title.position = "plot", 
        text = element_text(size=11),
        axis.text.x = element_text(),
        axis.line= element_line(),
        legend.position = "none",
        panel.border = element_blank(), 
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.spacing.x = unit(0.3,"line"))+ 
  scale_fill_manual(values = c("orange", "gray68"))+
  labs(x ="Day", y = "Colony fertility (eggs produced)", fill="Treatment")+ 
  scale_x_discrete(breaks=seq(5,155,5))+
  scale_y_continuous(breaks=seq(0,250,50))+
  geom_vline(xintercept=2.5, size = 0.5, linetype='dashed', col = 'red')+
  geom_vline(xintercept=13.5, size = 0.5, linetype='dashed', col = 'black'))
  
ggsave("Treatment_time_boxplot.pdf", plot_bp.time.eggs, width = 16, height = 10, units = "cm")
ggsave("Treatment_time_boxplot.png", plot_bp.time.eggs, width = 16, height = 10, units = "cm")
ggsave("Treatment_time_boxplot.svg", plot_bp.time.eggs, width = 16, height = 6, units = "cm")



# The model for these data should be a glm, with a quadratic term (to deal with the hump)

#Re-order the factors so that R is contrasted against C (rather than the other way around)

df4$treat <- factor(df4$treat, levels = c("C", "R"))


# As the data are count data we shall start with glmm from the poisson family. As I was curious, I included treatment as a random effect in these models, but as there are only two levels to treatment it makes much more sense to only include it as a fixed effect. Therefore the only models included in the manuscript include treatment as a fixed effect.

#Make null models

null.p.model1 <- glmer( eggs ~ 1 + (1|ID), family=poisson, data=df4)
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

overdisp_fun(null.p.model1) #we want the ratio of the residual deviance to the residual degrees of freedom to be 1 (or less than 1), here the ratio is much higher than 1 so the data is overdispersed  - this is a model written by Bolker to test for overdiespersion in glmer models. You can use it but cite Bolker, 2016 (see Liz thesis for full ref and when/how to cite him) - therefore need to retest this model when correcting for overdispersion

#Make a new column in the data frame with a random factor called 'obs'
df4$obs <- factor(1:nrow(df4))
df4$obs
df4

#Make null model with observation as a random factor

null.p.model2 <- glmer( eggs ~ 1 + (1|ID)+(1|obs), family=poisson, data=df4)

overdisp_fun(null.p.model2) #data no longer overdispersed

#Also want to check if over random effects make for a better fit to the data
null.p.model3 <- glmer( eggs ~ 1 + (1|ID)+(1|obs)+(1|day), family=poisson, data=df4)
null.p.model4 <- glmer( eggs ~ 1 + (1|ID)+(1|obs)+(1|day)+(1|treat), family=poisson, data=df4)

anova(null.p.model2, null.p.model3, null.p.model4) #nullmodel3 provides a better fit, and nullmodel4 provides an even better fit

#Make a null model where the effect of treatment on egg count is dependent on ID. I.e. the effect of treatment on egg count is blocked with respect to ID.

null.p.model5 <- glmer(eggs ~ 1 +(1|ID)+(1|obs)+(1|day)+(treat|ID), family=poisson, data=df4)
summary(null.p.model5) #Note the convergence warning, try rescaling variables

null.p.model5 <- glmer(eggs ~ 1 +(1|ID)+(1|obs)+(1|day.scale)+(treat.scale|ID), family=poisson, data=df4)#convergence warning still there


anova(null.p.model4, null.p.model5) #nullmodel4 provides a better fit to the data, problematic warning for nullmodel5, and the random terms treat|ID explain a tiny proportion of the variance


## Full model of fixed effects, days as a fixed factor and id as random factors. Continuous random effects go in brackets (continuous variable 1|categorical variable 1), categorical random factors in individual brackets (i.e. . (1|category1)+(1|category2)., fixed factor in the main equation e.g. response variable1 ~ predictor variable1*fixedfactor1) 

p.model1<- glmer(eggs ~ treat*day + (1|ID)+(1|obs)+(1|day)+(1|treat), family=poisson, data=df4)
#convergence warning, recalculate with scaled day variable
p.model1<- glmer(eggs ~ treat*day.scale + (1|ID)+(1|obs)+(1|day.scale)+(1|treat), family=poisson, data=df4)

summary(p.model1) #Note. Convergence warnings, none of the variance is explained by having treatment as a random factor
overdisp_fun(p.model1) 
#The ratio is now less than 1 so the data is underdispersed. Fine for this class of models

#Is the full model better than the null model?
anova(null.p.model4,p.model1) #much better fit for the full model

#Handy function in the MuMIn package which allows you to generate R2 values so you can work out how well your model fits your data (should do model simplification or model averaging before this step however). This function uses the method of Nakagawa, S., Johnson, P. C. D., and Schielzeth, H. 2012. The coefficient of determination R2 and intra-class correlation coefficient from generalized linear mixed-effects models revisited and expanded. Journal of the Royal Society Interface. 14(134): 20170213 to calculate r2 from glmms (and should be cited)
r.squaredGLMM(p.model1)



#Now let's try coding it as a quadratic model (as the relationship between egg count and time is curved), note. poly(day,2) is the same effect as including an effect of (day^2) to the model

p.model2 <- glmer(eggs ~ treat*poly(day,2) + (1|ID)+(1|obs)+(1|day)+(1|treat), family=poisson, data=df4)#the model (warning message 'failure to converge' might need to simplify)

##rescale treatment and day
p.model2 <- glmer(eggs ~ treat.scale*poly(day.scale,2) + (1|ID)+(1|obs)+(1|day.scale)+(1|treat.scale), family=poisson, data=df4)#the model (warning message 'failure to converge' might need to simplify)
overdisp_fun(p.model2)#not overdispersed 
summary(p.model2) #model is a significant fit to the data, and shows that Treatment and time have a strongly significant effect on egg count, and there is a significant interaction between treatment and time on eggcount. Note that treatment as a random effect explains very little.
r.squaredGLMM(p.model2)



#model simplification steps
anova(p.model1,p.model2)#model 2 is the best explanation for the data

p.model3 <- glmer(eggs ~ treat*poly(day,2)+(1|ID)+(1|obs)+(1|day), family=poisson, data=df4)#the model with treat removed as a random effect, this model doesn't have the convergence warnings
p.model4 <- glmer(eggs ~ treat*poly(day,2)+(1|ID)+(1|obs)+(1|treat), family=poisson, data=df4)#the model with day removed as a random effec, this model doesn't have the convergence warnings
p.model5 <- glmer(eggs ~ treat*poly(day,2)+(1|ID)+(1|obs), family=poisson, data=df4)#the model with both treat and day removed as a random effect, this model doesn't have the convergence warnings
anova(p.model2,p.model3,p.model4,p.model5) # model 2 and 3 are performing best and are significantly different from the other models
anova(p.model2,p.model3) # No significant difference between model 2 and 3, but as model 3 has lowest AIC, and as it's a simpler model, and does not have the convergence warnings, I will use model 3

r.squaredGLMM(p.model3) #marginal R2 is 0.404 and conditional r2 is 0.988
r.squaredGLMM(p.model4)
r.squaredGLMM(p.model5) #models 4 and 5 explain more of the variation but have lower liklihood scores, R2 often increases as the number of variables increase, but that doesn't make them better models (i.e. they might be overfit)


#test if data is explained by treatment varying by allowing the slopes for day and treatment to vary for a given ID
p.model6 <- glmer(eggs ~ treat*poly(day,2) + (1|ID) + (1|obs) + (1|day) + (day+treat|ID), family=poisson, data=df4)
r.squaredGLMM(p.model6)

p.model7 <- glmer(eggs ~ treat*poly(day,2) + (1|ID) + (1|obs) + (1|day) + (day|ID), family=poisson, data=df1)
r.squaredGLMM(p.model7)

p.model8 <- glmer(eggs ~ treat*poly(day,2) + (1|ID) + (1|obs) + (1|day) + (treat|ID), family=poisson, data=df1)

anova(p.model3,p.model6,p.model7,p.model8) #Model 7 is better, but I don't like the warning messages
anova(p.model3,p.model7) #Model 7 has a lower AIC, however it explains less of the variation (lower R2 values) and the warning message might be due to not enough data. These warnings are not due to scaling issues. So I will use model 3 instead

#Ultimately model 3 is the best performing model
#Chosen model parameters
summary(p.model3)
r.squaredGLMM(p.model3) #marginal R2 is 0.404 and conditional r2 is 0.988
plot(p.model3)#Looks horribly funnel shaped

plot(p.model3,residuals(.) ~log(fitted(.)))#Not sure if logging makes it look any better

ggplot(fortify(p.model3),
       aes(x=.fitted,y=sqrt(abs(.scresid))))+geom_point()+
  geom_smooth(colour="red",alpha=0.3)

#Check the residual distributions for individuals
plot(p.model3,ID~resid(.,type="pearson"))#E5 is quite an outlier, see what happens if removed
df4.or <- df4 %>% 
  filter(ID!="E5")

p.model3.or <- glmer(eggs ~ treat*poly(day,2)+(1|ID)+(1|obs)+(1|day), family=poisson, data=df4.or)#model3 with E5 removed
summary(p.model3.or) #very little difference to the results
r.squaredGLMM(p.model3.or) #marginal R2 is 0.402 and conditional r2 is 0.988 - again very little difference
plot(p.model3.or)

qqnorm(resid(p.model3)) #qqplot doesn't look normal either, so this is probably wrong, it's the random effects that are supposed to show a normal distribution


####Negative Binomial Models ####

#Will now try modeling the data using negative binomial rather than Poisson, we shall see if this provides a better fit to the data - note, these models are computationally intensive so they are hidden in hashtags. Only run the ones that you think are necessary to recreate the analysis. As model 2 is the 'best' model that is left without hashtags

#Start with null model

null.nb.model1 <- glmer.nb(eggs ~ 1 +(1|ID)+(1|obs)+(1|day), family=poisson, data=df4)

#Also look at partial model which does not account for treatment
#null.nb.model2 <- glmer.nb(eggs ~ poly(day,2) +(1|ID)+(1|obs)+(1|day), family=poisson, data=df4)

#nb.model1 <- glmer.nb(eggs ~ treat*poly(day,2) + (1|ID) + (1|obs) + (1|day), data=df4)#this model takes a long time to run and comes with the following Warning messages:
#Warning messages:
# 1: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#                  Model failed to converge with max|grad| = 0.00462676 (tol = 0.002, component 1)
#               2: In theta.ml(Y, mu, weights = object@resp$weights, limit = limit,  :
#                               iteration limit reached           iteration limit reached

#anova(null.nb.model1,nb.model1) #full model is better than the null model
#anova(null.nb.model2,nb.model2) #model is better when treatment is included as a fixed effect

#(p.model3,nb.model1)# the poisson and negative binomial models are not significantly different from each other

#lrtest(p.model3,nb.model1)# models aren't significantly different, but as negative binomial seems to perform better it is probably the better model

#Check if the observation term in the negative binomial can be dropped, and it remain under-dispersed (should do because negative binomial should not be overdispersed)
nb.model2 <-  glmer.nb(eggs ~ treat*poly(day,2) + (1|ID) + (1|day), data=df4)
summary(nb.model2)
r.squaredGLMM(nb.model2)
overdisp_fun(nb.model2) #model is not overdispersed (as a negative binomial model this doesn't really matter anyway)
getME(nb.model2, "glmer.nb.theta") #theta is 1.208918
plot(nb.model2) #plot looks better for negative binomial than poisson
plot(nb.model2,residuals(.) ~log(fitted(.)))#logging makes it look better
qqnorm(resid(nb.model2)) #note. using poisson errors to model count data, assumes a certain amount of heteroscedasicity. since the poisson ditribution has just one parameter it assumes that the variation increases with the linear predictor (i.e. this is heteroskedasticity), BUT for a given predicted value the variation that is expected around it (dispersion) is fixed. Therefore overdispersion matters, but other assumptions that are relevant to other models don't really

#nb.model3 <-  glmer.nb(eggs ~ treat*poly(day,2) + (1|ID) + (1|day) + (treat|ID), data=df1)

#nb.model4 <-  glmer.nb(eggs ~ treat*poly(day,2) + (1|ID), data=df1)

#anova(nb.model2, nb.model3, nb.model4)#no significant difference between models 2 and 3, however model 2 has lower AIC and therefore performs better than the other 2 models


#nb.model5 <- glmer.nb(eggs ~ treat*poly(day,2) + (1|ID) + (1|day), data=df1, control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=2e5))) - no significant difference between this and model2

#Compare best negative binomial model with its null model
null.best <- glmer.nb(eggs ~ 1 + (1|ID) + (1|day), data=df4)
anova(nb.model2, null.best) 



#Of the best of the negative binomial models, and the poisson models, the negative binomial has a much better fit to the data. In addition, it is not overdispersed even without the observation level effect, and seems like a more sensible model given my data. However, the two models are not significantly different to each other
pchisq(2 * (logLik(nb.model2) - logLik(p.model3)), df = 1, lower.tail = FALSE) #Small value again implies non significant difference


#Get the model outputs as confidence intervals
#confint(nb.model2) #Note. takes a very long time to run, confidence intervals for the fixed effects written below:
#                                    2.5 %      97.5 %
#.sig01                         0.09515878   0.2057957
#.sig02                         0.22636768   0.3592742
#(Intercept)                    2.75281420   2.9774944
#treatTreatment                 0.99115933   1.3023273
#poly(day, 2)1                -24.81239103 -18.4731615
#poly(day, 2)2                 -4.37631019   1.6046228
#treatTreatment:poly(day, 2)1  10.73627161  19.5500906
#treatTreatment:poly(day, 2)2 -11.68949875  -3.5338499

#Use allEffects function to more easily visualise the output of the model
plot(allEffects(nb.model2))

#Use summ to get the confidence intervals (works much quicker than the above function), the model outputs, the ICC values of the grouping variables
summ(nb.model2,scale=TRUE,digits=3, confint=T)

#visualise the summs output, and compare between the best models
plot_summs(nb.model2,p.model3)

#ICC values for the entire models
performance::icc(nb.model2)

#Get the summ values for the null model
summ(null.best, confint = T, digits=3)

#summ doesnt work if the confit argument is true with a null model so use the confint function
confint(null.best) 

#Export to word
#export_summs(nb.model2, scale = TRUE, error_format = "[{conf.low}, {conf.high}]")
#export_summs(nb.model2, scale = TRUE, to.file = "docx", file.name = "test.docx")

#As the 'best' model worked well for the egg count data, we shall try the same model with the cell count data (which is a similar datatype with a similar structure)
nb.cell.null1<-  glmer.nb(cells ~ 1 + (1|ID) + (1|day), data=df4)
nb.cell.null2<-  glmer.nb(cells ~ poly(day,2) + (1|ID) + (1|day), data=df4)
nb.cell.model2 <-  glmer.nb(cells ~ treat*poly(day,2) + (1|ID) + (1|day), data=df4)
plot(nb.cell.model2)
summ(nb.cell.model2,scale=TRUE,digits=3, confint=T)
anova(nb.cell.model2,nb.cell.null1)
anova(nb.cell.model2,nb.cell.null2)
#The full model explains more of the variation than the two best null models.
