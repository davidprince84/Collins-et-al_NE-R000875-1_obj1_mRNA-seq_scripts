#David Collins
#22/08/2020
#This script is to analyse data on the queen activity levels from NER008751 Objective 1 Experiment 1.


#packages required
library(tidyverse)
library(effects)
library(jtools)



#First you must always reset R
rm(list=ls())


#set file directory
setwd("C:/Users/dhcollins500/OneDrive - University of East Anglia/Documents/Post Doc_LHT/Objective 1 Experiment 1/Analysis/")

aggression_data <- read_csv("aggression.csv")
View(aggression_data)
df1<-aggression_data


#Make treat a factor and reorder it so R is the first level
df1$treat <- as.factor(df1$treat)
is.factor(df1$treat)
df1$treat <- factor(df1$treat, levels = c("R", "C"))


#Make a new version of df2 with the missing values omitted
df2<-na.omit(df1)

distinct(df1,active)
distinct(df2,active)


######ACTIVITY########

##Activity Plot:

# Make a new dataframe that summarises the queen activity levels on each day
df3 <- df2 %>%
  group_by(day,day_abc,treat_activity) %>%
  summarise (
    total_individuals = n(),
    total_activity = sum(active == TRUE),
    total_inactivity = sum (active == FALSE),
    total_response = sum(response== TRUE),
    total_non_response = sum(response==FALSE))

#Mutate so that there is a seperate column for just treatment, for the purposes of this script then the removal (R) queens are referred to as 'A' queens, and the control(C) queens are referred to as 'C' queens. The reason for this is that the levels argument for the barplot shown below is complicated and it always defaults to alphabetical even when levels is specified as R coming before C. The only solution I have found is to change the data on treatment so the R queens come alphabetically before the C queens. If necessary R queens true name can then be specifically specified in the labelling arguments
df3 <- df3 %>%
  mutate(treat = case_when(treat_activity == 'C_Active' ~ 'C',
                               treat_activity == 'C_Inactive' ~ 'C',
                               treat_activity == 'R_Active' ~ 'A',
                               treat_activity == 'R_Inactive' ~ 'A'))


#Make treat a factor and reorder it so A (referred to as R in the paper) is the first level
df3$treat <- as.factor(df3$treat)
is.factor(df3$treat)
df3$treat <- factor(df3$treat, levels = c("A", "C"))


df3$treat_activity <- as.factor(df3$treat_activity)
df3$treat_activity <- factor(df3$treat_activity, levels = c("R_Active", "R_Inactive","C_Active", "C_Inactive"))


# Make figure for activity
(active.bc<-ggplot(data = df3,
                  aes(x=paste0(day_abc, "/",treat), #Use paste0 to split up a axis by treatment
                      y=total_individuals, fill=treat_activity)) +
  geom_bar(stat="identity", position="stack")+
  scale_x_discrete (labels=c("2", "", "6", "","10", "", "14", "","18", "", "22", "","26", "", "30", "","34", "", "38", "","42", "", "46", "","50", "", "54", "","58", "", "62", "","66", "", "70", "","74", "", "78", "","82", "", "86", "","90", "", "94", "","98", "",  "102", "", "106", "","110", "", "114", "","118", "", "122", "","126","", "130", "","134", "", "138", "","142", "", "150",""))+
  theme(plot.title.position = "plot", 
        text = element_text(size=10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line= element_line(),
        panel.border = element_blank(), 
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.spacing.x = unit(0.3,"line"),
        legend.position = c(0.89,0.8))+
        labs(x="Day", y="Number of individuals")+
  scale_fill_manual(values = c("orange", "orangered","gray68", "gray2"),
                    name="Treatment & Activity",
                    breaks=c("R_Active","R_Inactive",
                             "C_Active", "C_Inactive"),
                    labels=c("R queen, Active","R queen, Inactive",
                             "C queen, Active", "C queen, Inactive")))

#Might be useful in future to have the name and treatment of every box:
# (labels=c("2 (C)", "2 (T)", "6 (C)", "6 (T)","10 (C)", "10 (T)", "14 (C)", "14 (T)","18 (C)", "18 (T)", "22 (C)", "22 (T)","26 (C)", "26 (T)", "30 (C)", "30 (T)","34 (C)", "34 (T)", "38 (C)", "38 (T)","42 (C)", "42 (T)", "46 (C)", "46 (T)","50 (C)", "50 (T)", "54 (C)", "54 (T)","58 (C)", "58 (T)", "62 (C)", "62 (T)","66 (C)", "66 (T)", "70 (C)", "70 (T)","74 (C)", "74 (T)", "78 (C)", "78 (T)","82 (C)", "82 (T)", "86(C)", "86 (T)","90 (C)", "90 (T)", "94(C)", "94 (T)","98 (C)", "98 (T)",  "102 (C)", "102 (T)", "106 (C)", "106 (T)","110 (C)", "110 (T)", "114 (C)", "114 (T)","118 (C)", "118 (T)", "122 (C)", "122 (T)","126 (C)", "126 (T)", "130 (C)", "130 (T)","134 (C)", "134 (T)", "138(C)", "138 (T)","142 (C)", "146 (C)", "150 (C)","154 (C)")

ggsave("Treatment_obs_activity.png", active.bc, width = 16, height = 10, units = "cm")

##Activity Model:

is.factor(df2$active)

df2$active <- as.factor(df2$active)


null.model.a1 <- glmer( active ~ 1 + (1|ID) + (1|day), family=binomial(link = "logit"), data=df2)
null.model.a2 <- glmer( active ~ 1 + (1|day), family=binomial(link = "logit"), data=df2)
null.model.a3 <- glmer( active ~ 1 + (1|ID), family=binomial(link = "logit"), data=df2)

anova (null.model.a1, null.model.a2, null.model.a3) # null model 1 has the best fit to the data


#full model with fixed effects

model.a1 <- glmer( active ~ treat + (1|ID) + (1|day), family=binomial(link = "logit"), data=df2)
model.a2 <- glmer( active ~ day + (1|ID) + (1|day), family=binomial(link = "logit"), data=df2)

#day needs to be rescaled, can attempt this in two different ways (log makes more sense to me)

df2 <- mutate(df2, day.1 = day/10,
              logday = log(day))

#redo null models with rescaled data
null.model.a1 <- glmer( active ~ 1 + (1|ID) + (1|logday), family=binomial(link = "logit"), data=df2)
null.model.a2 <- glmer( active ~ 1 + (1|logday), family=binomial(link = "logit"), data=df2)
null.model.a3 <- glmer( active ~ 1 + (1|ID), family=binomial(link = "logit"), data=df2)

#full models with rescaled data
model.a1 <- glmer( active ~ treat + (1|ID) + (1|logday), family=binomial(link = "logit"), data=df2)
model.a2 <- glmer( active ~ logday + (1|ID) + (1|logday), family=binomial(link = "logit"), data=df2)
model.a3 <- glmer( active ~ treat*logday + (1|ID) + (1|logday), family=binomial(link = "logit"), data=df2)
model.a4 <- glmer( active ~ treat*logday + (1|ID), family=binomial(link = "logit"), data=df2)


anova(null.model.a1,model.a1) #no different to the null model

anova(null.model.a1,model.a2) #no better than the null model

anova(null.model.a1,model.a3) #no better than the null model (and not significant)

anova(null.model.a1,model.a4) #no different to the null model

anova(model.a2,model.a3,model.a4) #model a3 does better than a4, but overall model a2 provides the best fit to the data. This may indicate that treatment has no effect on activity level, where as day has a slight negative effect. I'm most interested in model 3 simply because I want to establish properly that treatment has *no* effect on activity. In addition, it is useful to see that there was no activity:day interaction

#Check what happens when you remove the interaction term
model.a5 <- glmer( active ~ treat + logday + (1|ID) + (1|logday), family=binomial(link = "logit"), data=df2)
summary(model.a5) #No indication of anything significant here.
#
anova(null.model.a1,model.a3,model.a5) #The null model does the best, but model a5 does slightly better than model a3. None of these models have very different AIC values from each other.


#Note: I would check for overdispersion (variance greater than mean) here, but because the outcome is binary, the data cannot be overdispersed



summary(model.a2)
summary(model.a3)
summ(model.a3)

plot(model.a2)
plot(model.a3)

plot(allEffects(model.a2))
plot(allEffects(model.a3)) #again shows graphically how activity falls over time, but is not affected by treatment.


qqnorm(resid(model.a3)) #qqplot doesn't look normal either, so this is probably wrong, it's the random effects that are supposed to show a normal distribution

#Note. I also tried model a3 as a polynomial just to check, and found the linear model performed better. No need to chase it up as no reason to suspect that treatment/day are having non-linear effeccts on activity levels

#Use prediction functions to estimate how well predicted values correlate with the observed values
Pred <- predict(model.a3, type = "response")
Pred <- if_else(Pred > 0.5, 1, 0)
ConfusionMatrix <- table(Pred, pull(df2, active)) #`pull` results in a vector
#correct classification rate
sum(diag(ConfusionMatrix))/sum(ConfusionMatrix) #This shows that the model correctly predicts whether an individual will be active or not roughly 66% of the time
ConfusionMatrix # A matrix of predictions, and whether they are correct, 1027 observations were correctly predicted to be inactive, 518 were incorrectly classified as active, 14 were incorrectly classified as inactive, and 24 were correctly classified as active


# Compute AUC for predicting activity with the model
Prob <- predict(model.a3, type="response")
Pred <- prediction(Prob, as.vector(pull(df2, active)))
AUC <- performance(Pred, measure = "auc")
AUC <- AUC@y.values[[1]]
AUC # if AUC = 0.5 means the model performs no better than random chance, the value of 0.696 given shows the model is okay (though not great)



### QUEEN RESPONSE ###


# Make a new dataframe that summarises the queen response levels on each day
df4 <- df2 %>%
  group_by(day_abc,day,treat_response) %>%
  summarise (
    total_individuals = n(),
    total_response = sum(response== TRUE),
    total_non_response = sum(response==FALSE))

#Mutate so that there is a seperate column for just treatment
df4 <- df4 %>%
  mutate(treat = case_when(treat_response == 'C_Response' ~ 'C',
                               treat_response == 'C_NoResponse' ~ 'C',
                               treat_response == 'R_Response' ~ 'A',
                               treat_response == 'R_NoResponse' ~ 'A'))

# Make figure for response


(response.bc<-ggplot(data = df4,
                    aes(x=paste0(day_abc, "/",treat), #Use paste0 to split up a axis by treatment
                        y=total_individuals, fill=treat_response)) +
  geom_bar(stat="identity", position="stack")+
  scale_x_discrete (labels=c("2", "", "6", "","10", "", "14", "","18", "", "22", "","26", "", "30", "","34", "", "38", "","42", "", "46", "","50", "", "54", "","58", "", "62", "","66", "", "70", "","74", "", "78", "","82", "", "86", "","90", "", "94", "","98", "",  "102", "", "106", "","110", "", "114", "","118", "", "122", "","126","", "130", "","134", "", "138", "","142", "", "150",""))+
  theme(plot.title.position = "plot", 
        text = element_text(size=10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line= element_line(),
        panel.border = element_blank(), 
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.spacing.x = unit(0.3,"line"),
        legend.position = c(0.87,0.8))+
  labs(x="Day", y="Number of individuals")+
  scale_fill_manual(values = c("orange", "orangered","gray68", "gray2"),
                    name=
"Treatment & 
Response level"
,
                    breaks=c("R_Response","R_NoResponse",
                             "C_Response", "C_NoResponse"),
                    labels=c("R queen, Response","R queen, Non-response",
                              "C queen, Response", "C queen, Non-response")))

ggsave("Treatment_obs_response.png", response.bc, width = 16, height = 10, units = "cm")


##Response Model:

is.factor(df2$response)

df2$response <- as.factor(df2$response)


null.model.r1 <- glmer( response ~ 1 + (1|ID) + (1|logday), family=binomial(link = "logit"), data=df2)
null.model.r2 <- glmer( response ~ 1 + (1|logday), family=binomial(link = "logit"), data=df2)
null.model.r3 <- glmer( response ~ 1 + (1|ID), family=binomial(link = "logit"), data=df2)

anova (null.model.r1, null.model.r2, null.model.r3) # null model 1 has the best fit to the data


#full model with fixed effects

model.r1 <- glmer( response ~ treat + (1|ID) + (1|logday), family=binomial(link = "logit"), data=df2)
model.r2 <- glmer( response ~ logday + (1|ID) + (1|logday), family=binomial(link = "logit"), data=df2)
model.r3 <- glmer( response ~ treat*logday + (1|ID) + (1|logday), family=binomial(link = "logit"), data=df2)
model.r4 <- glmer( response ~ treat*logday + (1|ID), family=binomial(link = "logit"), data=df2)

anova(model.r1,model.r2,model.r3,model.r4) #model 2 and model 3 are the best, retain model 3 (as we are most interested in the effect of treatment)

anova(null.model.r1,model.r3) #model 3 is significantly better than the null model
anova(null.model.r1,model.r4) #model 4 is significantly better than the null model

anova(model.r3, model.r4) #model r3 is the best

#Test whether it is worth dropping the interaction term:
model.r5 <- glmer( response ~ treat + logday + (1|ID) + (1|logday), family=binomial(link = "logit"), data=df2)

anova(model.r3, model.r5) #model r5 performs no better, so stick to the original model 

summary(model.r3) #This shows no significant effect of treatment on queen's response levels. As time went on both treatments became significantly less responsive to disturbance.

plot(model.r3)
plot(allEffects(model.r3))
qqnorm(resid(model.r3)) #qqplot doesn't look normal either, so this is probably wrong, it's the random effects that are supposed to show a normal distribution

Pred <- predict(model.r3, type = "response")
Pred <- if_else(Pred > 0.5, 1, 0)
ConfusionMatrix <- table(Pred, pull(df2, response)) #`pull` results in a vector
#correct classification rate
sum(diag(ConfusionMatrix))/sum(ConfusionMatrix) #This shows that the model correctly predicts whether an individual will be active or not roughly 90.8% of the time... pretty good
ConfusionMatrix # A matrix of predictions, and whether they are correct, 48 observations were correctly predicted to be non-responsive, 12 were incorrectly classified as responsive, 133 were incorrectly classified as non-responsive, and 1390 were correctly identified as responsive

# Compute AUC for predicting response with the model
Prob <- predict(model.r3, type="response")
Pred <- prediction(Prob, as.vector(pull(df2, response)))
AUC <- performance(Pred, measure = "auc")
AUC <- AUC@y.values[[1]]
AUC # if AUC = 0.5 means the model performs no better than random chance, the value of 0.895 given shows the predictions of the model are pretty good

