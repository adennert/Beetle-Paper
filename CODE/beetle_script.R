####R Script for Beetle Paper
####Originally written in fall 2018 by NR; updated in spring 2020 by NR and AD
####For submission to Ecology & Evolution

#### Library Loading & Set-up ####

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
#install.packages("patchwork")
library(patchwork)
#install.packages("png")
library(png)
#install.packages("grid")
library(grid)
#install.packages("dplyr")
library(dplyr)

## load figure silohouettes of weevils and carabids
## images downloaded under Creative Commons from phylopic.org
wee.pic <- readPNG("FIGURES/weevil.png") 
car.pic <- readPNG("FIGURES/carabid.png") 

w <- rasterGrob(wee.pic, interpolate = TRUE, 
                width=unit(1,'cm'),
                x = unit(1,"npc"), y = unit(1,"npc"),
                hjust = 1, vjust = 1)

c <- rasterGrob(car.pic, interpolate = TRUE, 
                width=unit(1.5,'cm'),
                x = unit(1,"npc"), y = unit(1,"npc"),
                hjust = 1, vjust = 1)

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
str(soiln.data) #check that distance is an integer
levels(soiln.data$transect) #check levels of categorical variable
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
plotResiduals(soiln.data$distance.std, sim$scaledResiduals)
plotResiduals(soiln.data$moisture.std, sim$scaledResiduals)
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

# soil model coefficient plot with parameters with strong effects in black, 
# and parameters with weak/no effect in grey
library(sjPlot)
p1 <- plot_model(soilnmodel, type = "est", title = "", group.terms = c(1,2,2),
           order.terms = c(1,2,3), colors = c("black", "grey"), axis.labels = 
             c("Distance *\nMoisture", "Moisture", "Distance"), dot.size = 3) +
  theme_classic(30) + geom_hline(yintercept = 0, lty = 2, colour = "gray") +
  labs(tag = "B")
  
# create a plot showing the model predictions and raw data
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


# plot the raw data points with model on top
# note that moisture is displayed held at its mean
p2 <- ggplot() +
  geom_point(data = soiln.data, aes(x = distance.std, y = soild15N), alpha = 0.6, 
             pch = 16, size = 3, colour = "black") + 
  theme_classic(30) + 
  geom_line(data = newdat.soil, aes(x = distance.std, y = fit), size = 1, 
            colour = "black") +
  geom_ribbon(data = newdat.soil, aes(x = distance.std, ymin = lwr, ymax = upr), fill = "black", 
              alpha = .25) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  coord_cartesian(ylim = c(0, 12)) +
  labs(x="Distance from river (m)", 
       y = expression(paste("Soil δ"^"15"*"N"*" (‰)"))) +
  scale_x_continuous(breaks = c(-1.2247, 0, 1.2247), label = c("0", "25", "50")) +
  #scale_x_continuous(labels = unstandardize())
  labs(tag = "A")

#use patchwork package to merge the coefficient plot with the model/raw data plot
p2 + p1
ggsave("FIGURES/fig2.png",  height=6, width=14, dpi = "retina")


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

str(sia.data)
ggplot(sia.data, aes(as.factor(distance), bodyd15N, fill = species)) +
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
weevilbodynmodel <- lm(bodyd15N ~ distance, data = weevil.subset)
summary(weevilbodynmodel)

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

# (2b) carabid body N model - using carabid.subset
hist(carabid.subset$bodyd15N) #check distribution of response
levels(carabid.subset$trophic) #check levels of categorical variables
carabid.subset$distance <- as.integer(carabid.subset$distance) #make distance an integer 
str(carabid.subset)
carabidbodynmodel <- lm(bodyd15N ~ distance, data = carabid.subset)
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

# create plots showing the model predictions and raw data
# create a new data frame to make model predictions
newdat.bodyn <- expand.grid(distance = with(sia.data, 
                                               seq(min(distance), max(distance), 
                                                   length.out = 58)))

# make model predictions with upr/lower confidence intervals for carabid model
carmodelexp <- predict(carabidbodynmodel, newdat.bodyn, type = "response", 
                    se.fit = TRUE, re.form = NA, full = T)
newdat.bodyn$car.fit <- carmodelexp$fit
newdat.bodyn$car.lwr <- carmodelexp$fit - 1.96 * carmodelexp$se.fit
newdat.bodyn$car.upr <- carmodelexp$fit + 1.96 * carmodelexp$se.fit

# make model predictions with upr/lower confidence intervals for weevil model
weemodelexp <- predict(weevilbodynmodel, newdat.bodyn, type = "response", 
                       se.fit = TRUE, re.form = NA, full = T)
newdat.bodyn$wee.fit <- weemodelexp$fit
newdat.bodyn$wee.lwr <- weemodelexp$fit - 1.96 * weemodelexp$se.fit
newdat.bodyn$wee.upr <- weemodelexp$fit + 1.96 * weemodelexp$se.fit

# plot the weevil raw data points with weevil model on top
p1 <- ggplot() +
  #add raw data points
  geom_point(data = weevil.subset, aes(x = distance, y = bodyd15N), alpha = 0.6, 
             pch = 16, size = 3, colour = "black") + 
  theme_classic(30) + 
  #add weevil model fit
  geom_line(data = newdat.bodyn, aes(x = distance, y = wee.fit), size = 1, 
            colour = "black") +
  geom_ribbon(data = newdat.bodyn, aes(x = distance, ymin = wee.lwr, ymax = 
                                         wee.upr), fill = "black", alpha = .25) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  coord_cartesian(ylim = c(0, 12)) +
  labs(x="Distance from river (m)", 
       y = expression(paste("Body δ"^"15"*"N"*" (‰)"))) +
  scale_x_continuous(breaks = c(0, 25, 50), label = c("0", "25", "50")) +
  scale_y_continuous(breaks = c(0,4,8,12), label = c("0","4","8","12")) +
  labs(tag = "A") +
  annotation_custom(grob = w) 

# plot the carabid raw data points with carabid model on top
p2 <- ggplot() +
  #add raw data points
  geom_point(data = carabid.subset, aes(x = distance, y = bodyd15N), alpha = 0.6, 
             pch = 16, size = 3, colour = "black") + 
  theme_classic(30) + 
  #add carabid model fit
  geom_line(data = newdat.bodyn, aes(x = distance, y = car.fit), size = 1, 
            colour = "black") +
  geom_ribbon(data = newdat.bodyn, aes(x = distance, ymin = car.lwr, ymax = 
                                         car.upr), fill = "black", alpha = .25) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  coord_cartesian(ylim = c(0, 12)) +
  labs(x="Distance from river (m)", 
       y = expression(paste("Body δ"^"15"*"N"*" (‰)"))) +
  scale_x_continuous(breaks = c(0, 25, 50), label = c("0", "25", "50")) +
  scale_y_continuous(breaks = c(0,4,8,12), label = c("0","4","8","12")) +
  labs(tag = "B") +
  annotation_custom(grob = c)

p1 + p2
ggsave("FIGURES/fig3.png",  height=6, width=14, dpi = "retina")

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
weevilbodysizemodel <- lm(median ~ bodyd15N, data = weevil.subset) 
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
sim <- simulateResiduals(fittedModel = weevilbodysizemodel, n = 500) 
# the calculated residuals are stored in sim$scaledResiduals
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

# (3b) carabid body size SIA subset model - using carabid.subset
hist(carabid.subset$median) #check distribution of response
levels(carabid.subset$trophic) #check levels of categorical variables
carabidbodysizemodel <- lm(median ~ bodyd15N, data = carabid.subset) 
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
sim <- simulateResiduals(fittedModel = carabidbodysizemodel, n = 500) 
# the calculated residuals are stored in sim$scaledResiduals
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


# create plots showing the model predictions and raw data
# create a new data frame to make model predictions
newdat.bodysize.wee <- expand.grid(bodyd15N = with(weevil.subset, 
                                            seq(min(bodyd15N), max(bodyd15N), 
                                                length.out = 28)))
newdat.bodysize.car <- expand.grid(bodyd15N = with(carabid.subset, 
                                                   seq(min(bodyd15N), max(bodyd15N), 
                                                       length.out = 30)))

# make model predictions with upr/lower confidence intervals for carabid model
carsizemodelexp <- predict(carabidbodysizemodel, newdat.bodysize.car, type = "response", 
                       se.fit = TRUE, re.form = NA, full = T)
newdat.bodysize.car$car.fit <- carsizemodelexp$fit
newdat.bodysize.car$car.lwr <- carsizemodelexp$fit - 1.96 * carsizemodelexp$se.fit
newdat.bodysize.car$car.upr <- carsizemodelexp$fit + 1.96 * carsizemodelexp$se.fit

# make model predictions with upr/lower confidence intervals for weevil model
weesizemodelexp <- predict(weevilbodysizemodel, newdat.bodysize.wee, type = "response", 
                       se.fit = TRUE, re.form = NA, full = T)
newdat.bodysize.wee$wee.fit <- weesizemodelexp$fit
newdat.bodysize.wee$wee.lwr <- weesizemodelexp$fit - 1.96 * weesizemodelexp$se.fit
newdat.bodysize.wee$wee.upr <- weesizemodelexp$fit + 1.96 * weesizemodelexp$se.fit

# plot the weevil raw data points with weevil model on top
p1 <- ggplot() +
  #add raw data points
  geom_point(data = weevil.subset, aes(x = bodyd15N, y = median), alpha = 0.6, 
             pch = 16, size = 3, colour = "black") + 
  theme_classic(30) + 
  #add weevil model fit
  geom_line(data = newdat.bodysize.wee, aes(x = bodyd15N, y = wee.fit), size = 1, 
            colour = "black") +
  geom_ribbon(data = newdat.bodysize.wee, aes(x = bodyd15N, ymin = wee.lwr, ymax = 
                                         wee.upr), fill = "black", alpha = .25) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  coord_cartesian(ylim = c(4, 6)) +
  labs(y="Elytron length (mm)", 
       x = expression(paste("Body δ"^"15"*"N"*" (‰)"))) +
  scale_y_continuous(breaks = c(4,5,6), label = c("4","5","6")) +
  labs(tag = "A") +
  annotation_custom(grob = w) 


# plot the carabid raw data points with carabid model on top
p2 <- ggplot() +
  #add raw data points
  geom_point(data = carabid.subset, aes(x = bodyd15N, y = median), alpha = 0.6, 
             pch = 16, size = 3, colour = "black") + 
  theme_classic(30) + 
  #add carabid model fit
  geom_line(data = newdat.bodysize.car, aes(x = bodyd15N, y = car.fit), size = 1, 
            colour = "black") +
  geom_ribbon(data = newdat.bodysize.car, aes(x = bodyd15N, ymin = car.lwr, ymax = 
                                         car.upr), fill = "black", alpha = .25) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  coord_cartesian(ylim = c(8, 11)) +
  labs(y="Elytron length (mm)", 
       x = expression(paste("Body δ"^"15"*"N"*" (‰)"))) +
  scale_y_continuous(breaks = c(8,9,10,11), label = c("8","9","10","11")) +
  labs(tag = "B") +
  annotation_custom(grob = c) 

p1 + p2
ggsave("FIGURES/fig4.png",  height=6, width=14, dpi = "retina")


#### 4. BODY SIZE (FULL MODEL) ####

# dataset 
bodysize.data <- read.csv("DATA/(3) bodysize.csv", header = TRUE, ",", strip.white = TRUE)
View(bodysize.data)
str(bodysize.data) 
summary(bodysize.data)
names(bodysize.data)

# must standardize sampling round and distance variables as they are both 
#continuous with differing variance, range, and units
bodysize.data$distance.std <- c(scale(bodysize.data$distance, center = TRUE, 
                                      scale = TRUE))
bodysize.data$round.std <- c(scale(bodysize.data$round, center = TRUE, 
                                      scale = TRUE))
str(bodysize.data)

# subsets
library(dplyr)
#weevil subset:
#create a subset that includes all weevils EXCEPT N. incomptus (n=1) 
# because insuffient sample size to estimate model parameters for this species
fullbodysizeweevil.subset <- filter(bodysize.data, species %in% c("Steremnius tuberosus", 
                                              "Steremnius carinatus")) %>% 
  droplevels()
View(fullbodysizeweevil.subset) #n=1010

#carabid subset:
#create a subset that includes all SEXED carabids 
# (EXCEPT L. ferruginosus (n=2) and C. tuberculatus (n=3))
# because insuffient sample size to estimate model parameters for this species
fullbodysizecarabid.subset <- filter(bodysize.data, species %in% c("Pterostichus amethystinus",
                                                                   "Pterostichus crenicollis",
                                                                   "Scaphinotus angusticollis",
                                                                   "Zacotus matthewsii")) %>% 
  droplevels()
fullbodysizesexedcarabid.subset <- filter(fullbodysizecarabid.subset, sex %in% c("F", "M")) %>% 
  droplevels()
View(fullbodysizesexedcarabid.subset) #n=290

# (4a) weevil full bodysize model
hist(fullbodysizeweevil.subset$median) #check distribution of response 
levels(fullbodysizeweevil.subset$species) #check levels of categorical variables
levels(fullbodysizeweevil.subset[,"species"]) #make S. tuberosus first (intercept)
fullbodysizeweevil.subset$species <- relevel(fullbodysizeweevil.subset$species,
                                             "Steremnius tuberosus")

#make 'round' an integer - but note not important as a factor either
fullbodysizeweevil.subset$round <- as.integer(fullbodysizeweevil.subset$round) 
str(fullbodysizeweevil.subset) 

fullbodysizeweevilmodel <- lmer(median ~ distance.std*species + round.std + 
                          (1|transect), data = fullbodysizeweevil.subset)
summary(fullbodysizeweevilmodel)

# (4a) check residuals
ggplot(fullbodysizeweevil.subset, aes(x = fitted(fullbodysizeweevilmodel), 
                                      y = resid(fullbodysizeweevilmodel))) +
  geom_point() +
  theme_classic() +
  geom_line(y=0, colour="red") +
  labs(x="Fitted values", y= "Residuals")
qqnorm(as.vector(resid(fullbodysizeweevilmodel)))
qqline(as.vector(resid(fullbodysizeweevilmodel)), col = "blue")

## (4a) using DHARMa to interpret residuals
# set simulations constant 
set.seed(1)
# calculate scaled residuals
library(DHARMa)
sim <- simulateResiduals(fittedModel = fullbodysizeweevilmodel, n = 500) 
# the calculated residuals are stored in sim$scaledResiduals
# plot the scaled residuals (Observed vs Expected)
plot(sim)
# plot residuals against the other predictors
plotResiduals(fullbodysizeweevil.subset$distance.std, sim$scaledResiduals) 
plotResiduals(fullbodysizeweevil.subset$species, sim$scaledResiduals) 
plotResiduals(fullbodysizeweevil.subset$round.std, sim$scaledResiduals)
# test outliers
testOutliers(sim)
# test dispersion 
testDispersion(sim)
# shows QQ plot, dispersion, outliers in 1 plot
testResiduals(sim)

# (4a) full coefficient plot for Supplemental 
#parameters with strong effects in black, parameters with weak/no effect in grey
library(sjPlot)
p3 <- plot_model(fullbodysizeweevilmodel, type = "est", title = "",
           group.terms = c(1,2,1,1), 
           order.terms = c(1,4,3,2), 
           colors = c("grey","black"), 
           axis.labels = c("S. carinatus", "Sampling \nWeek", 
                           "Distance * \nS. carinatus", "Distance")) +
  theme_classic(15) + geom_hline(yintercept = 0, lty = 2, colour = "gray") +
  labs(tag = "A") +
  annotation_custom(grob = w) 

# (4a) select coefficient plot for Body Text
# parameters with strong effects in black, parameters with weak/no effect in grey
# 'sampling round' and 'species' removed
library(sjPlot)
p1 <- plot_model(fullbodysizeweevilmodel, type = "est", title = "",
           group.terms = c(1,1,1), 
           order.terms = c(1,3,2), 
           colors = c("grey"), 
           rm.terms = c("speciesSteremnius carinatus"),
           axis.labels = c("Sampling \nWeek", "Distance * \nS. carinatus", 
                           "Distance")) +
  theme_classic(20) + geom_hline(yintercept = 0, lty = 2, colour = "gray") +
  labs(tag = "A") +
  annotation_custom(grob = w) 

# (4b) carabid full bodysize model 
hist(fullbodysizesexedcarabid.subset$median) #check distribution of response 
levels(fullbodysizesexedcarabid.subset$sex) #check levels of categorical variables
levels(fullbodysizesexedcarabid.subset$species) #check levels of categorical variables

levels(fullbodysizesexedcarabid.subset[,"species"]) #make P. amethystinus first (intercept) 
fullbodysizesexedcarabid.subset$species <- 
  relevel(fullbodysizesexedcarabid.subset$species, "Pterostichus amethystinus")

levels(fullbodysizesexedcarabid.subset[, "sex"]) # make F first (intercept)
fullbodysizesexedcarabid.subset$sex = factor(fullbodysizesexedcarabid.subset$sex, 
                                             levels(fullbodysizesexedcarabid.subset$sex)[c(1,2)])
str(fullbodysizesexedcarabid.subset) 
#make 'round' an integer - but note not important as a factor either
fullbodysizecarabid.subset$round <- as.integer(fullbodysizecarabid.subset$round)  

fullbodysizecarabidmodel <- lmer(median ~ distance.std*species + sex + round.std 
                                 + (1|transect), 
                                 data = fullbodysizesexedcarabid.subset)
summary(fullbodysizecarabidmodel)

# (4b) check residuals
ggplot(fullbodysizesexedcarabid.subset, aes(x = fitted(fullbodysizecarabidmodel), 
                                            y = resid(fullbodysizecarabidmodel))) +
  geom_point() +
  theme_classic() +
  geom_line(y=0, colour="red") +
  labs(x="Fitted values", y= "Residuals")
qqnorm(as.vector(resid(fullbodysizecarabidmodel)))
qqline(as.vector(resid(fullbodysizecarabidmodel)), col = "blue")

## (4b) using DHARMa to interpret residuals
# set simulations constant 
set.seed(1)
# calculate scaled residuals
library(DHARMa)
sim <- simulateResiduals(fittedModel = fullbodysizecarabidmodel, n = 500) 
# the calculated residuals are stored in sim$scaledResiduals
# plot the scaled residuals (Observed vs Expected)
plot(sim)
# plot residuals against the other predictors
plotResiduals(fullbodysizesexedcarabid.subset$distance.std, sim$scaledResiduals) 
plotResiduals(fullbodysizesexedcarabid.subset$species, sim$scaledResiduals) 
plotResiduals(fullbodysizesexedcarabid.subset$sex, sim$scaledResiduals)
plotResiduals(fullbodysizesexedcarabid.subset$round.std, sim$scaledResiduals)
# test outliers
testOutliers(sim)
# test dispersion 
testDispersion(sim)
# shows QQ plot, dispersion, outliers in 1 plot
testResiduals(sim)

# (4b) full coefficient plot for Supplemental 
#parameters with strong effects in black, parameters with weak/no effect in grey
library(sjPlot)
p4 <- plot_model(fullbodysizecarabidmodel, type = "est", title = "",
           group.terms = c(1,2,2,2,2,1,1,1,1), 
           order.terms = c(1,7,8,9,6,5,2,3,4), 
           colors = c("grey", "black"), 
           wrap.labels = 50,
           rm.terms = c(),
           axis.labels = c("Z. matthewsii", 
                           "S. angusticollis", 
                           "P. crenicollis", 
                           "Sex: Male", 
                           "Sampling \n Week", 
                           "Distance * \n Z. matthewsii", 
                           "Distance * \n S. angusticollis", 
                           "Distance * \n P. crenicollis", 
                           "Distance")) +
  theme_classic(15) + geom_hline(yintercept = 0, lty = 2, colour = "gray") +
  labs(tag = "B") +
  annotation_custom(grob = c) 

#(4b) select coefficient plot for Body Text 
#parameters with strong effects in black, parameters with weak/no effect in grey
library(sjPlot)
p2 <- plot_model(fullbodysizecarabidmodel, type = "est", title = "",
           group.terms = c(1,2,1,1,1,1), 
           order.terms = c(1,4,5,6,3,2), 
           colors = c("grey", "black"), 
           wrap.labels = 50,
           rm.terms = c("speciesPterostichus crenicollis", 
                        "speciesScaphinotus angusticollis",
                        "speciesZacotus matthewsii"),
           axis.labels = c("Sex: Male", 
                           "Sampling \n Week", 
                           "Distance * \n Z. matthewsii", 
                           "Distance * \n S. angusticollis", 
                           "Distance * \n P. crenicollis", 
                           "Distance")) +
  theme_classic(20) + geom_hline(yintercept = 0, lty = 2, colour = "gray") +
  labs(tag = "B") +
  annotation_custom(grob = c) 

p1 + p2
ggsave("FIGURES/fig5_simplified.png",  height=6, width=14, dpi = "retina")

p3 + p4
ggsave("FIGURES/fig5_full.png",  height=6, width=14, dpi = "retina")

### create plots showing the model predictions and raw data
### create a new data frame to make model predictions
#newdat.4a.wee <- expand.grid(distance.std = with(fullbodysizeweevil.subset, 
#                                                   seq(min(distance.std), max(distance.std), 
#                                                       length.out = 1010)),
#                                       species = c("Steremnius carinatus"),
#                                                   #"Steremnius tuberosus"), 
#                                       round.std = mean(fullbodysizeweevil.subset$round.std))
#
#newdat.4b.car <- expand.grid(distance.std = with(fullbodysizesexedcarabid.subset, 
#                                                   seq(min(distance.std), max(distance.std), 
#                                                       length.out = 290)),
#                                      species = c("Zacotus matthewsii"),
#                                       round.std = 
#                                        mean(fullbodysizesexedcarabid.subset$round.std),
#                                       sex = c("F"))
#                                     
### make model predictions with upr/lower confidence intervals for weevil model
#modelexp <- predict(fullbodysizeweevilmodel, 
#                         newdata = newdat.4a.wee,
#                           type = "response", re.form = NA)
#newdat.4a.wee$wee.fit <- modelexp
#
#boot <- bootMer(fullbodysizeweevilmodel, predict, nsim = 100, re.form = NA,
#               type = "parametric", seed = 101)
#std.err <- apply(boot$t, 2, sd)
#newdat.4a.wee$wee.lwr <- newdat.4a.wee$wee.fit - 1.96 * std.err
#newdat.4a.wee$wee.upr <- newdat.4a.wee$wee.fit + 1.96 * std.err
#
### make model predictions with upr/lower confidence intervals for carabid model
#modelexp <- predict(fullbodysizecarabidmodel, 
#                    newdata = newdat.4b.car,
#                    type = "response", re.form = NA)
#newdat.4b.car$car.fit <- modelexp
#
#boot <- bootMer(fullbodysizecarabidmodel, predict, nsim = 10000, re.form = NA,
#                type = "parametric")
#std.err <- apply(boot$t, 2, sd)
#newdat.4b.car$car.lwr <- newdat.4b.car$car.fit - 1.96 * std.err
#newdat.4b.car$car.upr <- newdat.4b.car$car.fit + 1.96 * std.err
#
#

#
## plot the weevil raw data points with weevil model on top
#ggplot() +
#  #add raw data points
#  geom_jitter(data = fullbodysizeweevil.subset, aes(x = distance.std, y = median), 
#             alpha = 0.6, pch = 16, size = 3, colour = "black", width = 0.1) + 
#  theme_classic(30) + 
#  #add weevil model fit
#  geom_line(data = newdat.4a.wee, aes(x = distance.std, y = wee.fit), 
#            size = 1, colour = "black") +
#  geom_ribbon(data = newdat.4a.wee, aes(x = distance.std, ymin = wee.lwr, 
#                                                    ymax = wee.upr), 
#             fill = "black", alpha = 0.25) +
#  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
#  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
#  #coord_cartesian(ylim = c(4, 6)) +
#  labs(y="Elytron length (mm)", 
#       x = expression(paste("Distance from river (m)"))) +
#  scale_x_continuous(breaks = c(-1.2247, 0, 1.2247), label = c("0", "25", "50")) +
#  #scale_y_continuous(breaks = c(4,5,6), label = c("4","5","6")) +
#  labs(tag = "A") +
#  annotation_custom(grob = w) 
#
## plot the weevil raw data points against distance
#data_summary <- function(x) {
#  m <- mean(x)
#  ymin <- m-(1.96*plotrix::std.error(x))
#  ymax <- m+(1.96*plotrix::std.error(x))
#  return(c(y=m,ymin=ymin,ymax=ymax)) }
#
#ggplot() +
#  #add raw data points
#  geom_jitter(data = fullbodysizeweevil.subset, aes(x = distance, y = median), 
#              alpha = 0.6, pch = 16, size = 3, colour = "black", width = 1.5) + 
#  theme_classic(30) + 
#  #add weevil model fit
#  #geom_line(data = newdat.4a.wee, aes(x = distance.std, y = wee.fit), 
#            #size = 1, colour = "black") +
#  #geom_ribbon(data = newdat.4a.wee, aes(x = distance.std, ymin = wee.lwr, 
#                                        #ymax = wee.upr), 
#              #fill = "black", alpha = 0.25) +
#  #stat_smooth(data = fullbodysizeweevil.subset, aes(x = distance, y = median),
#              #method = lm) +
#  #geom_errorbar(data = fullbodysizeweevil.subset, aes(x = distance, y = median) +
#  stat_summary(data = fullbodysizeweevil.subset, aes(x = distance, y = median),
#    fun.data=data_summary, color="red", size = 0.5) +
#  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
#  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
#  #coord_cartesian(ylim = c(4, 6)) +
#  labs(y="Elytron length (mm)", 
#       x = expression(paste("Distance from river (m)"))) +
#  scale_x_continuous(breaks = c(0, 25, 50), label = c("0", "25", "50")) +
#  #scale_y_continuous(breaks = c(4,5,6), label = c("4","5","6")) +
#  labs(tag = "A") +
#  annotation_custom(grob = w) 
#
## plot the carabid raw data points with carabid model on top
#ggplot() +
#  #add raw data points
#  geom_jitter(data = fullbodysizesexedcarabid.subset, aes(x = distance.std, 
#                                                         y = median), 
#             alpha = 0.6, pch = 16, size = 3, colour = "black", width = 0.1) + 
#  theme_classic(30) + 
#  #add carabid model fit
#  geom_line(data = newdat.4b.car, aes(x = distance.std, y = car.fit), size = 1, 
#            colour = "black") +
#  geom_ribbon(data = newdat.4b.car, aes(x = distance.std, ymin = car.lwr, ymax = 
#                                                car.upr), fill = "black", alpha = .25) +
#  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
#  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
#  #coord_cartesian(ylim = c(4, 6)) +
#  labs(y="Elytron length (mm)", 
#       x = expression(paste("Distance from river (m)"))) +
#  #scale_y_continuous(breaks = c(4,5,6), label = c("4","5","6")) +
#  scale_x_continuous(breaks = c(-1.2247, 0, 1.2247), label = c("0", "25", "50")) +
#  labs(tag = "B") +
#  annotation_custom(grob = c) 

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
weevilposthocmodel <- lm(bodyncombust ~ bodyd15N, data = weevil.subset)
summary(weevilposthocmodel)

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


# (5b) carabid post hoc model
hist(carabid.subset$bodyncombust) #check distribution of response
levels(carabid.subset$trophic) #check levels of categorical variables
carabidposthocmodel <- lm(bodyncombust ~ bodyd15N, data = carabid.subset)
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

# create plots showing the model predictions and raw data
# create a new data frame to make model predictions
newdat.post.wee <- expand.grid(bodyd15N = with(weevil.subset, 
                                                   seq(min(bodyd15N), max(bodyd15N), 
                                                       length.out = 28)))
newdat.post.car <- expand.grid(bodyd15N = with(carabid.subset, 
                                                   seq(min(bodyd15N), max(bodyd15N), 
                                                       length.out = 30)))


# make model predictions with upr/lower confidence intervals for carabid model
carpostmodelexp <- predict(carabidposthocmodel, newdat.post.car, type = "response", 
                           se.fit = TRUE, re.form = NA, full = T)
newdat.post.car$car.fit <- carpostmodelexp$fit
newdat.post.car$car.lwr <- carpostmodelexp$fit - 1.96 * carpostmodelexp$se.fit
newdat.post.car$car.upr <- carpostmodelexp$fit + 1.96 * carpostmodelexp$se.fit

# make model predictions with upr/lower confidence intervals for weevil model
weepostmodelexp <- predict(weevilposthocmodel, newdat.post.wee, type = "response", 
                           se.fit = TRUE, re.form = NA, full = T)
newdat.post.wee$wee.fit <- weepostmodelexp$fit
newdat.post.wee$wee.lwr <- weepostmodelexp$fit - 1.96 * weepostmodelexp$se.fit
newdat.post.wee$wee.upr <- weepostmodelexp$fit + 1.96 * weepostmodelexp$se.fit

# plot the weevil raw data points with weevil model on top
p1 <- ggplot() +
  #add raw data points
  geom_point(data = weevil.subset, aes(x = bodyd15N, y = bodyncombust), alpha = 0.6, 
             pch = 16, size = 3, colour = "black") + 
  theme_classic(30) + 
  #add weevil model fit
  geom_line(data = newdat.post.wee, aes(x = bodyd15N, y = wee.fit), size = 1, 
            colour = "black") +
  geom_ribbon(data = newdat.post.wee, aes(x = bodyd15N, ymin = wee.lwr, ymax = 
                                                wee.upr), fill = "black", alpha = .25) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  #coord_cartesian(ylim = c(4, 6)) +
  labs(y="Body N (%)", 
       x = expression(paste("Body δ"^"15"*"N"*" (‰)"))) +
  scale_y_continuous(breaks = c(8,8.5,9, 9.5), label = c("8","8.5","9", "9.5")) +
  labs(tag = "A") +
  annotation_custom(grob = w) 


# plot the carabid raw data points with carabid model on top
p2 <- ggplot() +
  #add raw data points
  geom_point(data = carabid.subset, aes(x = bodyd15N, y = bodyncombust), alpha = 0.6, 
             pch = 16, size = 3, colour = "black") + 
  theme_classic(30) + 
  #add carabid model fit
  geom_line(data = newdat.post.car, aes(x = bodyd15N, y = car.fit), size = 1, 
            colour = "black") +
  geom_ribbon(data = newdat.post.car, aes(x = bodyd15N, ymin = car.lwr, ymax = 
                                                car.upr), fill = "black", alpha = .25) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  #coord_cartesian(ylim = c(8, 11)) +
  labs(y="Body N (%)", 
       x = expression(paste("Body δ"^"15"*"N"*" (‰)"))) +
  scale_y_continuous(breaks = c(10,11,12), label = c("10","11","12")) +
  labs(tag = "B") +
  annotation_custom(grob = c) 

p1 + p2
ggsave("FIGURES/fig6.png",  height=6, width=14, dpi = "retina")

