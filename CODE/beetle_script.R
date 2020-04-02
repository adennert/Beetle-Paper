####R Script for Beetle Paper
####Originally written in fall 2018 by NR; updated in spring 2020 by NR and AD
####For submission to Ecology & Evolution
#TEST merge, 26 March 2020 3:46 pm

#### Libraries ####

#install.packages("datasets")
library(datasets)
#install.packages("stats")
library(stats)
#install.packages("doBy")
library(doBy)
#install.packages("lme4")
library(lme4)
#install.packages("Matrix")
library(Matrix)
help(Matrix)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("car")
library(car)
#install.packages("ggpubr")
library(ggpubr)
#install.packages("Hmisc")
library(Hmisc)
#install.packages("DHARMa")
library(DHARMa)

#### 0. REPEATABILITY ####
# dataset
repeat.data <- read.csv("DATA/(0) repeat.csv", header = TRUE, ",", strip.white = TRUE)
View(repeat.data)
str(repeat.data)
summary(repeat.data)
names(repeat.data)

# (0) repeatability model
repeatmodel <- lmer(length ~ 1 + (1|insectid), data = repeat.data)
summary(repeatmodel)

# check residuals
ggplot(repeat.data, aes(x = fitted(repeatmodel), y = resid(repeatmodel))) +
  geom_point() +
  theme_classic() +
  geom_line(y=0, colour="red") +
  labs(x="Fitted values", y= "Residuals")
qqnorm(as.vector(resid(repeatmodel)))
qqline(as.vector(resid(repeatmodel)), col = "blue")

# calculate repeatability factor (use VarCorr to extract variance components)
VarCorr(repeatmodel)
# R  = (2.099381)^2 / ((2.099381^2) + (0.051741^2))
#    = 0.99939295
#    Repeatability factor (R) is 0.99 so very high repeatability
#    (very low within group variation compared to b/w group variation)

#### 1. SOIL NUTRIENTS ####

# dataset
soiln.data <- read.csv("DATA/(1) soiln.csv", header = TRUE, ",", strip.white = TRUE)
View(soiln.data)
str(soiln.data) #check that distance is an integer or numeric
levels(soiln.data$transect) #check levels of categorical variable
#soiln.data$distance <- as.numeric(soiln.data$distance)
summary(soiln.data)
names(soiln.data)

# visualize raw soil data
library(ggplot2)
ggplot(soiln.data, aes(distance, soild15N)) +
  geom_point() + theme_classic() + 
  theme(text = element_text(size = 16)) +
  theme(axis.text.x = element_text(size = 16, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.text.y = element_text(size = 16, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  coord_cartesian(ylim = c(0, 12)) +
  labs(x="Distance from the Kunsoot River (m)", 
       y = expression(paste("Soil δ"^"15"*"N"*" (‰)")))

# (1) soil model
# cannot include random effect of transect due to singularity issues
# must standardize these variables as they are both continuous with differing
# variance, range, and units
soiln.data$distance.std <- c(scale(soiln.data$distance, center = TRUE, scale = TRUE))
soiln.data$moisture.std <- c(scale(soiln.data$moisture, center = TRUE, scale = TRUE))

# create the standardized soil model
soilnmodel <- lm(soild15N ~ distance.std * moisture.std, data = soiln.data)
summary(soilnmodel)

# check for correlation between distance and moisture
cor(soiln.data$distance, soiln.data$moisture, method = c("pearson"))
# Pearson's correlation coefficient: r = -0.006, where 0 means there is no association.

# check residuals
ggplot(soiln.data, aes(x = fitted(soilnmodel), y = resid(soilnmodel))) +
  geom_point() +
  theme_classic() +
  geom_line(y=0, colour="red") +
  labs(x="Fitted values", y= "Residuals")
qqnorm(as.vector(resid(soilnmodel)))
qqline(as.vector(resid(soilnmodel)), col = "blue")

##USING DHARMa PACKAGE TO INTERPRET RESIDUALS (March 2020):
# set simulations constant 
set.seed(1)
# calculate scaled residuals
library(DHARMa)
sim <- simulateResiduals(fittedModel = soilnmodel, n = 500) 
# the calculated residuals are stored in sim$scaledResiduals
# plot residuals against the other predictors
plotResiduals(soiln.data$distance, sim$scaledResiduals)
plotResiduals(soiln.data$moisture, sim$scaledResiduals)
# shows QQ plot, dispersion, outliers in 1 plot
testResiduals(sim)

# soil model coefficient plot with parameters with strong effects in black, 
# and parameters with weak/no effect in grey
library(sjPlot)
plot_model(soilnmodel, type = "est", title = "", group.terms = c(1,2,2),
           order.terms = c(1,2,3), colors = c("black", "grey"), axis.labels = 
             c("Distance * Moisture", "Moisture", "Distance")) +
  theme_classic() + geom_hline(yintercept = 0, lty = 2, colour = "gray") 
  
# create a new data frame to make model predictions
newdat.soil <- expand.grid(distance.std = with(soiln.data, 
                                             seq(min(distance.std), max(distance.std), 
                                                 length.out = 30)),
                           moisture.std = mean(soiln.data$moisture.std))
modelexp <- predict(soilnmodel, newdat.soil, type = "response", 
                    se.fit = TRUE, re.form = NA, full = T)
newdat.soil$fit <- modelexp$fit
newdat.soil$lwr <- modelexp$fit - 1.96 * modelexp$se.fit
newdat.soil$upr <- modelexp$fit + 1.96 * modelexp$se.fit

# Code to display the x-axis unstandardized
#unstandardize <- function() {
  #function(x) format(x*sd(soiln.data$distance) + mean(soiln.data$distance), digits = 0) 
#}

# plot the raw data points with model on top
# note that moisture is displayed held at its mean
ggplot() +
  geom_point(data = soiln.data, aes(x = distance.std, y = soild15N), alpha = 0.6, 
             pch = 16, size = 2, colour = "black") + 
  theme_classic() + 
  geom_line(data = newdat.soil, aes(x = distance.std, y = fit), size = 1, 
            colour = "black") +
  geom_ribbon(data = newdat.soil, aes(x = distance.std, ymin = lwr, ymax = upr), fill = "black", 
              alpha = .25) +
  theme(text = element_text(size = 16)) +
  theme(axis.text.x = element_text(size = 16, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.text.y = element_text(size = 16, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  coord_cartesian(ylim = c(0, 12)) +
  labs(x="Distance from the Kunsoot River (m)", 
       y = expression(paste("Soil δ"^"15"*"N"*" (‰)"))) +
  scale_x_continuous(breaks = c(-1.2247, 0, 1.2247), label = c("0", "25", "50"))
  #scale_x_continuous(labels = unstandardize())
                   

#### 2. BODY N (SIA) ####

# dataset
sia.data <- read.csv("DATA/(2) sia.csv", header = TRUE, ",", strip.white = TRUE)
View(sia.data)
str(sia.data) 
summary(sia.data)
names(sia.data)

# subsets
library(dplyr)
weevil.subset <- sia.data %>% filter(trophic == "Curculionidae")
View(weevil.subset)
carabid.subset <- sia.data %>% filter(trophic == "Carabidae")
View(carabid.subset)

# visualize raw data 
library(ggplot2)
sia.data$distance <- as.factor(sia.data$distance) #make distance a factor for the boxplot 
str(sia.data)
ggplot(sia.data, aes(distance, bodyd15N, fill = species)) +
  stat_boxplot() + theme_classic() +
  theme(text = element_text(size = 16)) +
  theme(axis.text.x = element_text(size = 16, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.text.y = element_text(size = 16, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  coord_cartesian(ylim = c(-2, 12)) +
  theme(legend.position="top") + 
  theme(legend.title = element_text(size = 14, color = "black", face = "bold")) +
  theme(legend.text = element_text(size = 14, face = "italic")) +
  labs(x="Distance from the Kunsoot River (m)", 
       y = expression(paste("Body δ"^"15"*"N"*" (‰)")))

# (2a) weevil body N model - using weevil.subset
hist(weevil.subset$bodyd15N) #check distribution of response
levels(weevil.subset$trophic) #check levels of categorical variables
weevil.subset$distance <- as.integer(weevil.subset$distance) #make distance an integer 
str(weevil.subset)
weevilbodynmodel <- lmer(bodyd15N ~ distance + (1|transect), data = weevil.subset)

# check residuals
ggplot(weevil.subset, aes(x = fitted(weevilbodynmodel), y = resid(weevilbodynmodel))) +
  geom_point() +
  theme_classic() +
  geom_line(y=0, colour="red") +
  labs(x="Fitted values", y= "Residuals")
qqnorm(as.vector(resid(weevilbodynmodel)))
qqline(as.vector(resid(weevilbodynmodel)), col = "blue")

## using DHARMa to interpret residuals
# set simulations constant 
set.seed(1)
# calculate scaled residuals
library(DHARMa)
sim <- simulateResiduals(fittedModel = weevilbodynmodel, n = 500) # the calculated residuals are stored in sim$scaledResiduals
# plot the scaled residuals (Observed vs Expected)
plot(sim)
# plot residuals against the other predictors
plotResiduals(weevil.subset$distance, sim$scaledResiduals) 
# test outliers
testOutliers(sim)
# test dispersion 
testDispersion(sim)
# shows QQ plot, dispersion, outliers in 1 plot
testResiduals(sim)

# (2a) coefficient plot with parameters with strong effects in black, 
# and parameters with weak/no effect in grey
library(sjPlot)
 plot_model(weevilbodynmodel, type = "std2", title = "", 
           group.terms = c(1), order.terms = c(1), colors = c("black"), 
           axis.labels = c("Distance"), axis.lim = c(-1,1)) + 
  theme_classic() + geom_hline(yintercept = 0, lty = 2, colour = "gray") 

# (2b) carabid body N model - using carabid.subset
hist(carabid.subset$bodyd15N) #check distribution of response
levels(carabid.subset$trophic) #check levels of categorical variables
carabid.subset$distance <- as.integer(carabid.subset$distance) #make distance an integer 
str(carabid.subset)
carabidbodynmodel <- lmer(bodyd15N ~ distance + (1|transect), data = carabid.subset)
summary(carabidbodynmodel)

# check residuals
ggplot(carabid.subset, aes(x = fitted(carabidbodynmodel), y = resid(carabidbodynmodel))) +
  geom_point() +
  theme_classic() +
  geom_line(y=0, colour="red") +
  labs(x="Fitted values", y= "Residuals")
qqnorm(as.vector(resid(carabidbodynmodel)))
qqline(as.vector(resid(carabidbodynmodel)), col = "blue")

## using DHARMa to interpret residuals 
# set simulations constant 
set.seed(1)
# calculate scaled residuals
library(DHARMa)
sim <- simulateResiduals(fittedModel = carabidbodynmodel, n = 500) # the calculated residuals are stored in sim$scaledResiduals
# plot the scaled residuals (Observed vs Expected)
plot(sim)
# plot residuals against the other predictors
plotResiduals(carabid.subset$distance, sim$scaledResiduals) 
# test outliers
testOutliers(sim)
# test dispersion 
testDispersion(sim)
# shows QQ plot, dispersion, outliers in 1 plot
testResiduals(sim)

# (2b) coefficient plot with parameters with strong effects in black, 
# and parameters with weak/no effect in grey
library(sjPlot)
plot_model(carabidbodynmodel, type = "std2", title = "", 
           group.terms = c(1), order.terms = c(1), colors = c("grey"), 
           axis.labels = c("Distance"), axis.lim = c(-1,1)) +
  theme_classic() + geom_hline(yintercept = 0, lty = 2, colour = "gray") 

 
#### 3. BODY SIZE (SIA) ####

# dataset
sia.data <- read.csv("DATA/(2) sia.csv", header = TRUE, ",", strip.white = TRUE)
View(sia.data)
str(sia.data) 
summary(sia.data)
names(sia.data)

# subsets
library(dplyr)
weevil.subset <- sia.data %>% filter(trophic == "Curculionidae")
View(weevil.subset)
carabid.subset <- sia.data %>% filter(trophic == "Carabidae")
View(carabid.subset)

#visualize all raw data
library(ggplot2)
ggplot(sia.data, aes(bodyd15N, median, color = species)) +
  stat_smooth(method = "lm") + geom_point(size = 0.7) + theme_classic() +
  theme(text = element_text(size = 16)) +
  theme(axis.text.x = element_text(size = 16, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.text.y = element_text(size = 16, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  coord_cartesian(ylim = c(NULL)) +
  theme(legend.position="top") + 
  theme(legend.title = element_text(size = 14, color = "black", face = "bold")) +
  theme(legend.text = element_text(size = 14, face = "italic")) +
  labs(x = expression(paste("Body δ"^"15"*"N"*" (‰)")), 
       y = "Elytron length (mm)")

# (3a) weevil body size SIA subset model - using weevil.subset
hist(weevil.subset$median) #check distribution of response
levels(weevil.subset$trophic) #check levels of categorical variables
weevilbodysizemodel <- lmer(median ~ bodyd15N + (1|transect), data = weevil.subset) 
summary(weevilbodysizemodel)

# check residuals
ggplot(weevil.subset, aes(x = fitted(weevilbodysizemodel), y = resid(weevilbodysizemodel))) +
  geom_point() +
  theme_classic() +
  geom_line(y=0, colour="red") +
  labs(x="Fitted values", y= "Residuals")
qqnorm(as.vector(resid(weevilbodysizemodel)))
qqline(as.vector(resid(weevilbodysizemodel)), col = "blue")

## using DHARMa to interpret residuals
# set simulations constant 
set.seed(1)
# calculate scaled residuals
library(DHARMa)
sim <- simulateResiduals(fittedModel = weevilbodysizemodel, n = 500) # the calculated residuals are stored in sim$scaledResiduals
# plot the scaled residuals (Observed vs Expected)
plot(sim)
# plot residuals against the other predictors
plotResiduals(weevil.subset$bodyd15N, sim$scaledResiduals) 
# test outliers
testOutliers(sim)
# test dispersion 
testDispersion(sim)
# shows QQ plot, dispersion, outliers in 1 plot
testResiduals(sim)

# (3a) coefficient plot with parameters with strong effects in black, 
# and parameters with weak/no effect in grey
library(sjPlot)
plot_model(weevilbodysizemodel, type = "std2", title = "", 
           group.terms = c(1), order.terms = c(1), colors = c("grey"), 
           axis.labels = c(paste("Body δ15N (‰)")), axis.lim = c(-1,1)) +
  theme_classic() + geom_hline(yintercept = 0, lty = 2, colour = "gray") 


# (3b) carabid body size SIA subset model - using carabid.subset
hist(carabid.subset$median) #check distribution of response
levels(carabid.subset$trophic) #check levels of categorical variables
carabidbodysizemodel <- lmer(median ~ bodyd15N + (1|transect), data = carabid.subset) 
summary(carabidbodysizemodel)

# check residuals
ggplot(carabid.subset, aes(x = fitted(carabidbodysizemodel), y = resid(carabidbodysizemodel))) +
  geom_point() +
  theme_classic() +
  geom_line(y=0, colour="red") +
  labs(x="Fitted values", y= "Residuals")
qqnorm(as.vector(resid(carabidbodysizemodel)))
qqline(as.vector(resid(carabidbodysizemodel)), col = "blue")

## using DHARMa to interpret residuals
# set simulations constant 
set.seed(1)
# calculate scaled residuals
library(DHARMa)
sim <- simulateResiduals(fittedModel = carabidbodysizemodel, n = 500) # the calculated residuals are stored in sim$scaledResiduals
# plot the scaled residuals (Observed vs Expected)
plot(sim)
# plot residuals against the other predictors
plotResiduals(carabid.subset$bodyd15N, sim$scaledResiduals) 
# test outliers
testOutliers(sim)
# test dispersion 
testDispersion(sim)
# shows QQ plot, dispersion, outliers in 1 plot
testResiduals(sim)

# (3b) coefficient plot with parameters with strong effects in black, 
# and parameters with weak/no effect in grey
library(sjPlot)
plot_model(carabidbodysizemodel, type = "std2", title = "", 
           group.terms = c(1), order.terms = c(1), colors = c("grey"), 
           axis.labels = c(paste("Body δ15N (‰)")), axis.lim = c(-1,1)) +
  theme_classic() + geom_hline(yintercept = 0, lty = 2, colour = "gray") 


#### 4. BODY SIZE (FULL MODEL) ####

# dataset 
bodysize.data <- read.csv("DATA/(3) bodysize.csv", header = TRUE, ",", strip.white = TRUE)
View(bodysize.data)
str(bodysize.data) 
bodysize.data$round <- ordered(bodysize.data$round, levels = 1:4) #round as ordinal 
summary(bodysize.data)
names(bodysize.data)

# subset
# create a subset that eliminates N. incomptus (n=1) and L. ferruginosus (n=2)
# because insuffient sample size to estimate model parameters for these species
library(dplyr)
fullbodysize.subset <- filter(bodysize.data, species %in% c("Steremnius tuberosus","Cychrus tuberculatus", 
                                                            "Pterostichus amethystinus","Pterostichus crenicollis", 
                                                            "Scaphinotus angusticollis","Steremnius carinatus", 
                                                            "Zacotus matthewsii"))
View(fullbodysize.subset)

# (4) full bodysize model
hist(fullbodysize.subset$median) #check distribution of response
levels(fullbodysize.subset$trophic) #check levels of categorical variables
levels(fullbodysize.subset$sex)
levels(fullbodysize.subset$species)
levels(fullbodysize.subset[,"species"]) # check S. tuberosus first (intercept)
fullbodysize.subset$species = factor(fullbodysize.subset$species, levels(fullbodysize.subset$species)[c(2,3,4,1,5,6,7,8,9)])
levels(fullbodysize.subset[, "sex"]) # check that F first (intercept)
fullbodysize.subset$sex = factor(fullbodysize.subset$sex, levels(fullbodysize.subset$sex)[c(2,3,1)])
str(fullbodysize.subset) #check 'round' is ordinal

fullbodysize.subset$round <- ordered(fullbodysize.subset$round, levels = 1:4, labels=c("1", "2", "3", "4"))

fullbodysizemodel <- lmer(median ~ distance*species + trophic + sex + round + 
                            (1|transect), data = fullbodysize.subset)
summary(fullbodysizemodel)

# check residuals
ggplot(fullbodysize.subset, aes(x = fitted(fullbodysizemodel), y = resid(fullbodysizemodel))) +
  geom_point() +
  theme_classic() +
  geom_line(y=0, colour="red") +
  labs(x="Fitted values", y= "Residuals")
qqnorm(as.vector(resid(fullbodysizemodel)))
qqline(as.vector(resid(fullbodysizemodel)), col = "blue")

## using DHARMa to interpret residuals
# set simulations constant 
set.seed(1)
# calculate scaled residuals
library(DHARMa)
sim <- simulateResiduals(fittedModel = fullbodysizemodel, n = 500) # the calculated residuals are stored in sim$scaledResiduals
# plot the scaled residuals (Observed vs Expected)
plot(sim)
# plot residuals against the other predictors
plotResiduals(fullbodysize.subset$distance, sim$scaledResiduals) 
plotResiduals(fullbodysize.subset$species, sim$scaledResiduals) #2 spp. have been eliminated
plotResiduals(fullbodysize.subset$trophic, sim$scaledResiduals)
plotResiduals(fullbodysize.subset$sex, sim$scaledResiduals)
plotResiduals(fullbodysize.subset$round, sim$scaledResiduals)
# test outliers
testOutliers(sim)
# test dispersion 
testDispersion(sim)
# shows QQ plot, dispersion, outliers in 1 plot
testResiduals(sim)

# (4) coefficient plot with parameters with strong effects in black, 
# and parameters with weak/no effect in grey
library(sjPlot)
plot_model(fullbodysizemodel, type = "est", title = "", #choose correct type =
           group.terms = c(1,2,1,2,2,2,2,2,1,1,1,2,1,1,1,1,1,1), 
           order.terms = c(1,13,14,15,16,17,18,12,11,10,8,9,2,3,4,5,6,7), 
           colors = c("grey", "black"), 
           axis.labels = c("")) +
  theme_classic() + geom_hline(yintercept = 0, lty = 2, colour = "gray") 


#testing out labels:
# (4) coefficient plot with parameters with strong effects in black, 
# and parameters with weak/no effect in grey
library(sjPlot)
plot_model(fullbodysizemodel, type = "est", title = "", #choose correct type =
           group.terms = c(1,2,1,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1), 
           order.terms = c(1,13,14,15,16,17,18,12,11,10,8,9,2,3,4,5,6,7), 
           colors = c("grey", "black"), 
           axis.labels = c("Z. matthewsii", "S. carinatus", "S. angusticollis", 
                           "P. crenicollis", "P. amethystinus", "C. tuberculatus", 
                           "Sex NA", "Male", "Sampling Round 4", "Sampling Round 3", 
                           "Sampling Round 2", "Distance * Z. matthewsii", "Distance * S. carinatus", 
                           "Distance * S. angusticollis", "Distance * P. crenicollis", 
                           "Distance * P. amethystinus", "Distance * C. tuberculatus", 
                           "Distance"), axis.lim = c(-3,9)) +
  theme_classic() + geom_hline(yintercept = 0, lty = 2, colour = "gray") 

#### 5. POST HOC (SIA) ####

# dataset
sia.data <- read.csv("DATA/(2) sia.csv", header = TRUE, ",", strip.white = TRUE)
View(sia.data)
str(sia.data) 
summary(sia.data)
names(sia.data)

# subsets
library(dplyr)
weevil.subset <- sia.data %>% filter(trophic == "Curculionidae")
View(weevil.subset)
carabid.subset <- sia.data %>% filter(trophic == "Carabidae")
View(carabid.subset)

# visualize raw data
library(ggplot2)
ggplot(sia.data, aes(bodyd15N, bodyncombust, color = species)) +
  stat_smooth(method = "lm") + geom_point(size = 0.7) + theme_classic() +
  theme(text = element_text(size = 16)) +
  theme(axis.text.x = element_text(size = 16, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.text.y = element_text(size = 16, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  coord_cartesian(ylim = c(NULL)) +
  theme(legend.position="top") + 
  theme(legend.title = element_text(size = 14, color = "black", face = "bold")) +
  theme(legend.text = element_text(size = 14, face = "italic")) +
  labs(x = expression(paste("Body δ"^"15"*"N"*" (‰)")), 
       y = "Body N (%)")

# (5a) weevil post hoc model
hist(weevil.subset$bodyncombust) #check distribution of response
levels(weevil.subset$trophic) #check levels of categorical variables
weevilposthocmodel <- lmer(bodyncombust ~ bodyd15N + (1|transect), data = weevil.subset)

# check residuals
ggplot(weevil.subset, aes(x = fitted(weevilposthocmodel), y = resid(weevilposthocmodel))) +
  geom_point() +
  theme_classic() +
  geom_line(y=0, colour="red") +
  labs(x="Fitted values", y= "Residuals")
qqnorm(as.vector(resid(weevilposthocmodel)))
qqline(as.vector(resid(weevilposthocmodel)), col = "blue")

## using DHARMa to interpret residuals
# set simulations constant 
set.seed(1)
# calculate scaled residuals
library(DHARMa)
sim <- simulateResiduals(fittedModel = weevilposthocmodel, n = 500) # the calculated residuals are stored in sim$scaledResiduals
# plot the scaled residuals (Observed vs Expected)
plot(sim)
# plot residuals against the other predictors
plotResiduals(weevil.subset$bodyd15N, sim$scaledResiduals) 
# test outliers
testOutliers(sim)
# test dispersion 
testDispersion(sim)
# shows QQ plot, dispersion, outliers in 1 plot
testResiduals(sim)

# (5a) coefficient plot with parameters with strong effects in black, 
# and parameters with weak/no effect in grey
library(sjPlot)
plot_model(weevilposthocmodel, type = "std2", title = "", 
           group.terms = c(1), order.terms = c(1), colors = c("black"), 
           axis.labels = c(paste("Body δ15N (‰)")), axis.lim = c(-1,1)) +
  theme_classic() + geom_hline(yintercept = 0, lty = 2, colour = "gray") 

# (5b) carabid post hoc model
hist(carabid.subset$bodyncombust) #check distribution of response
levels(carabid.subset$trophic) #check levels of categorical variables
carabidposthocmodel <- lmer(bodyncombust ~ bodyd15N + (1|transect), data = carabid.subset)
summary(carabidposthocmodel)

# check residuals
ggplot(carabid.subset, aes(x = fitted(carabidposthocmodel), y = resid(carabidposthocmodel))) +
  geom_point() +
  theme_classic() +
  geom_line(y=0, colour="red") +
  labs(x="Fitted values", y= "Residuals")
qqnorm(as.vector(resid(carabidposthocmodel)))
qqline(as.vector(resid(carabidposthocmodel)), col = "blue")

## using DHARMa to interpret residuals
# set simulations constant 
set.seed(1)
# calculate scaled residuals
library(DHARMa)
sim <- simulateResiduals(fittedModel = carabidposthocmodel, n = 500) # the calculated residuals are stored in sim$scaledResiduals
# plot the scaled residuals (Observed vs Expected)
plot(sim)
# plot residuals against the other predictors
plotResiduals(carabid.subset$bodyd15N, sim$scaledResiduals) 
# test outliers
testOutliers(sim)
# test dispersion 
testDispersion(sim)
# shows QQ plot, dispersion, outliers in 1 plot
testResiduals(sim)

# (5b) coefficient plot with parameters with strong effects in black, 
# and parameters with weak/no effect in grey
library(sjPlot)
plot_model(carabidposthocmodel, type = "std2", title = "", 
           group.terms = c(1), order.terms = c(1), colors = c("grey"), 
           axis.labels = c(paste("Body δ15N (‰)")), axis.lim = c(-1,1)) +
  theme_classic() + geom_hline(yintercept = 0, lty = 2, colour = "gray") 
