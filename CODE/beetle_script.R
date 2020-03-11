####R Script for Beetle Paper
####Originally written in fall 2018 by NR; updated in spring 2020 by NR and AD
####For submission to Ecology & Evolution -AD

#### BELLA BELLA BEETLE ANALYSIS 2018 ####

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

# for bootstrapped confidence intervals
library(visreg)
library(devtools)
devtools::install_github("remkoduursma/bootpredictlme4")
library(bootpredictlme4)


#### DECEMBER 2018 NEW ANALYSIS ####
#### 0. REPEATABILITY ####

# data - repeatability (n = 1338 * 3 = 4014)
repeat.data <- read.csv(file.choose(), header = TRUE, ",")
View(repeat.data)
str(repeat.data)
summary(repeat.data)
names(repeat.data)

# (0) repeatability model
repeatmodel <- lmer(length ~ 1 + (1|insectid), data = repeat.data)
summary(repeatmodel)
df.residual(repeatmodel, type = c("lmer")) #this needs to be corrected, df can't be 4011

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

# data - soil (n = 30)
soiln.data <- read.csv(file.choose(), header = TRUE, ",")
View(soiln.data)
str(soiln.data) #distance as an integer for the model 
summary(soiln.data)
names(soiln.data)

soilnfig.data <- read.csv(file.choose(), header = TRUE, ",")
soilnfig.data$distance <- as.factor(soilnfig.data$distance)
View(soilnfig.data)
str(soilnfig.data) #distance as a factor for the boxplot figure
summary(soilnfig.data)
names(soilnfig.data)

# visualize raw soil data
library(ggplot2)
ggplot(soilnfig.data, aes(distance, soild15N)) +
  stat_boxplot(geom = "boxplot") + theme_classic() + 
  theme(text = element_text(size = 16)) +
  theme(axis.text.x = element_text(size = 16, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.text.y = element_text(size = 16, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  coord_cartesian(ylim = c(0, 12)) +
  labs(x="Distance from the Kunsoot River (m)", 
       y = expression(paste("Soil δ"^"15"*"N"*" (‰)")))

# (1) soil model
soilnmodel <- lmer(soild15N ~ distance*moisture + (1|transect), data = soiln.data)
summary(soilnmodel)
df.residual(soilnmodel, type = c("lmer")) 

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
sim <- simulateResiduals(fittedModel = soilnmodel, n = 500) # the calculated residuals are stored in sim$scaledResiduals
# plot the scaled residuals (Observed vs Expected)
plot(sim)
# plot residuals against the other predictors
plotResiduals(soiln.data$distance, sim$scaledResiduals)
plotResiduals(soiln.data$moisture, sim$scaledResiduals)
# test outliers
testOutliers(sim)
# test dispersion 
testDispersion(sim)
# shows QQ plot, dispersion, outliers in 1 plot
testResiduals(sim)

# soil model coefficient plot
library(sjPlot)
plot_model(soilnmodel, type = "std2", colors = "bw", title = "soil δ15N") +
  theme_classic() + geom_hline(yintercept = 0, lty = 2, colour = "gray")

# get the standardized coefficients: "std" = forest-plot of standardized beta values
library(sjPlot)
get_model_data(soilnmodel, type = c("std"), show.df = TRUE)


#### 2. BODY NUTRIENTS ####

# data - weevil body nutrients (n = 28)
weevilbodyn.data <- read.csv(file.choose(), header = TRUE, ",")
View(weevilbodyn.data)
str(weevilbodyn.data) #distance as an integer for model
summary(weevilbodyn.data)
names(weevilbodyn.data)

weevilbodynfig.data <- read.csv(file.choose(), header = TRUE, ",")
weevilbodynfig.data$distance <- as.factor(weevilbodynfig.data$distance)
View(weevilbodynfig.data)
str(weevilbodynfig.data) #distance as a factor for figure
summary(weevilbodynfig.data)
names(weevilbodynfig.data)

# data - carabid body nutrients (n = 30)
carabidbodyn.data <- read.csv(file.choose(), header = TRUE, ",")
View(carabidbodyn.data)
str(carabidbodyn.data) #distance as an integer for model
summary(carabidbodyn.data)
names(carabidbodyn.data)

carabidbodynfig.data <- read.csv(file.choose(), header = TRUE, ",")
carabidbodynfig.data$distance <- as.factor(carabidbodynfig.data$distance)
View(carabidbodynfig.data)
str(carabidbodynfig.data) #distance as a factor for figure
summary(carabidbodynfig.data)
names(carabidbodynfig.data)

# data - combined for figure
bothbodynfig.data <- read.csv(file.choose(), header = TRUE, ",")
bothbodynfig.data$distance <- as.factor(bothbodynfig.data$distance)
View(bothbodynfig.data)
str(bothbodynfig.data) #distance as factor for figure (once run ^ as.factor)
names(bothbodynfig.data)

# visualize raw weevil data 
library(ggplot2)
ggplot(weevilbodynfig.data, aes(distance, bodyd15N)) +
  stat_boxplot() + theme_classic() +
  theme(text = element_text(size = 16)) +
  theme(axis.text.x = element_text(size = 16, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.text.y = element_text(size = 16, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  coord_cartesian(ylim = c(NULL)) +
  labs(x="Distance from the Kunsoot River (m)", 
       y = expression(paste("Curculionidae body δ"^"15"*"N"*" (‰)")))

# (2a) weevil body nutrients model
weevilbodynmodel <- lmer(bodyd15N ~ distance + (1|transect), data = weevilbodyn.data) # singular fit message because only 2 T5s and 2T7s (all others have 3 each)
summary(weevilbodynmodel)
df.residual(weevilbodynmodel)

# check residuals
ggplot(weevilbodyn.data, aes(x = fitted(weevilbodynmodel), y = resid(weevilbodynmodel))) +
  geom_point() +
  theme_classic() +
  geom_line(y=0, colour="red") +
  labs(x="Fitted values", y= "Residuals")
qqnorm(as.vector(resid(weevilbodynmodel)))
qqline(as.vector(resid(weevilbodynmodel)), col = "blue")

##USING DHARMa PACKAGE TO INTERPRET RESIDUALS (March 2020):
# set simulations constant 
set.seed(1)
# calculate scaled residuals
library(DHARMa)
sim <- simulateResiduals(fittedModel = weevilbodynmodel, n = 500) # the calculated residuals are stored in sim$scaledResiduals
# plot the scaled residuals (Observed vs Expected)
plot(sim)
# plot residuals against the other predictors
plotResiduals(weevilbodyn.data$distance, sim$scaledResiduals) 
# test outliers
testOutliers(sim)
# test dispersion 
testDispersion(sim)
# shows QQ plot, dispersion, outliers in 1 plot
testResiduals(sim)

# coefficient plot - distance affects weevil body nutrients
library(sjPlot)
plot_model(weevilbodynmodel, type = "std2", colors = "bw", title = "body d15N") +
  theme_classic() +
  geom_hline(yintercept = 0, lty = 2, colour = "gray")

# get standardized coefficients: "std" = forest-plot of standardized beta values
library(sjPlot)
get_model_data(weevilbodynmodel, type = c("std"), show.df = TRUE)

# visualize raw carabid data 
library(ggplot2)
ggplot(carabidbodynfig.data, aes(distance, bodyd15N)) +
  stat_boxplot() + theme_classic() +
  theme(text = element_text(size = 16)) +
  theme(axis.text.x = element_text(size = 16, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.text.y = element_text(size = 16, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  coord_cartesian(ylim = c(0, 12)) +
  labs(x="Distance from the Kunsoot River (m)", 
       y = expression(paste("Carabidae body δ"^"15"*"N"*" (‰)")))

# (2b) carabid body nutrients model
carabidbodynmodel <- lmer(bodyd15N ~ distance + (1|transect), data = carabidbodyn.data)
summary(carabidbodynmodel)
df.residual(carabidbodynmodel)

# check residuals
ggplot(carabidbodyn.data, aes(x = fitted(carabidbodynmodel), y = resid(carabidbodynmodel))) +
  geom_point() +
  theme_classic() +
  geom_line(y=0, colour="red") +
  labs(x="Fitted values", y= "Residuals")
qqnorm(as.vector(resid(carabidbodynmodel)))
qqline(as.vector(resid(carabidbodynmodel)), col = "blue")

##USING DHARMa PACKAGE TO INTERPRET RESIDUALS (March 2020):
# set simulations constant 
set.seed(1)
# calculate scaled residuals
library(DHARMa)
sim <- simulateResiduals(fittedModel = carabidbodynmodel, n = 500) # the calculated residuals are stored in sim$scaledResiduals
# plot the scaled residuals (Observed vs Expected)
plot(sim)
# plot residuals against the other predictors
plotResiduals(carabidbodyn.data$distance, sim$scaledResiduals) 
# test outliers
testOutliers(sim)
# test dispersion 
testDispersion(sim)
# shows QQ plot, dispersion, outliers in 1 plot
testResiduals(sim)

# coefficient plot - no distance effect
library(sjPlot)
plot_model(carabidbodynmodel, type = "std2", colors = "bw", title = "body d15N") +
  theme_classic() + geom_hline(yintercept = 0, lty = 2, colour = "gray")

# get standardized coefficients: "std" = forest-plot of standardized beta values
library(sjPlot)
get_model_data(carabidbodynmodel, type = c("std"), show.df = TRUE)

# see both weevils and carabids together for plot
library(ggplot2)
ggplot(bothbodynfig.data, aes(distance, bodyd15N, fill = species)) +
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


#### 3. SIA body subset ####
# data - weevil SIA body size (average of 3 individuals) (n = 28)
weevilsia.data <- read.csv(file.choose(), header = TRUE, ",")
View(weevilsia.data)
str(weevilsia.data) 
summary(weevilsia.data)
names(weevilsia.data)

# data - carabid SIA body size (1 indiv. per sample) (n = 30)
carabidsia.data <- read.csv(file.choose(), header = TRUE, ",")
View(carabidsia.data)
str(carabidsia.data) 
summary(carabidsia.data)
names(carabidsia.data)

# data - both combined for figure
bothsia.data <- read.csv(file.choose(), header = TRUE, ",")
View(bothsia.data)
str(bothsia.data) 
summary(bothsia.data)
names(bothsia.data)

# visualize raw weevil data 
library(ggplot2)
ggplot(weevilsia.data, aes(bodyd15N, avglength)) +
  stat_smooth() + theme_classic() +
  theme(text = element_text(size = 16)) +
  theme(axis.text.x = element_text(size = 16, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.text.y = element_text(size = 16, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  coord_cartesian(ylim = c(NULL)) +
  labs(x = expression(paste("Curculionidae body δ"^"15"*"N"*" (‰)")), 
       y = "Elytron length (mm)")

# (3a) weevil sia model
weevilsiamodel <- lmer(avglength ~ bodyd15N + (1|transect), data = weevilsia.data)
summary(weevilsiamodel)
df.residual(weevilsiamodel)

# check residuals
ggplot(weevilsia.data, aes(x = fitted(weevilsiamodel), y = resid(weevilsiamodel))) +
  geom_point() +
  theme_classic() +
  geom_line(y=0, colour="red") +
  labs(x="Fitted values", y= "Residuals")
qqnorm(as.vector(resid(weevilsiamodel)))
qqline(as.vector(resid(weevilsiamodel)), col = "blue")

##USING DHARMa PACKAGE TO INTERPRET RESIDUALS (March 2020):
# set simulations constant 
set.seed(1)
# calculate scaled residuals
library(DHARMa)
sim <- simulateResiduals(fittedModel = weevilsiamodel, n = 500) # the calculated residuals are stored in sim$scaledResiduals
# plot the scaled residuals (Observed vs Expected)
plot(sim)
# plot residuals against the other predictors
plotResiduals(weevilsia.data$bodyd15N, sim$scaledResiduals) 
# test outliers
testOutliers(sim)
# test dispersion 
testDispersion(sim)
# shows QQ plot, dispersion, outliers in 1 plot
testResiduals(sim)

# coefficient plot - distance affects weevil body nutrients
library(sjPlot)
plot_model(weevilsiamodel, type = "std2", colors = "bw", title = "elytron length") +
  theme_classic() + geom_hline(yintercept = 0, lty = 2, colour = "gray")

# get the standardized coefficients: "std" = forest-plot of standardized beta values
library(sjPlot)
get_model_data(weevilsiamodel, type = c("std"), show.df = TRUE)         

# visualize raw carabid data 
library(ggplot2)
ggplot(carabidsia.data, aes(bodyd15N, medianelytronlength)) +
  stat_smooth() + theme_classic() +
  theme(text = element_text(size = 16)) +
  theme(axis.text.x = element_text(size = 16, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.text.y = element_text(size = 16, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  coord_cartesian(ylim = c(NULL)) +
  labs(x = expression(paste("Carabidae body δ"^"15"*"N"*" (‰)")), 
       y = "Elytron length (mm)")

# (3b) carabid sia model
carabidsiamodel <- lmer(medianelytronlength ~ bodyd15N + (1|Transect), data = carabidsia.data)
summary(carabidsiamodel)
df.residual(carabidsiamodel)

# check residuals
ggplot(carabidsia.data, aes(x = fitted(carabidsiamodel), y = resid(carabidsiamodel))) +
  geom_point() +
  theme_classic() +
  geom_line(y=0, colour="red") +
  labs(x="Fitted values", y= "Residuals")
qqnorm(as.vector(resid(carabidsiamodel)))
qqline(as.vector(resid(carabidsiamodel)), col = "blue")

##USING DHARMa PACKAGE TO INTERPRET RESIDUALS (March 2020):
# set simulations constant 
set.seed(1)
# calculate scaled residuals
library(DHARMa)
sim <- simulateResiduals(fittedModel = carabidsiamodel, n = 500) # the calculated residuals are stored in sim$scaledResiduals
# plot the scaled residuals (Observed vs Expected)
plot(sim)
# plot residuals against the other predictors
plotResiduals(carabidsia.data$bodyd15N, sim$scaledResiduals) 
# test outliers
testOutliers(sim)
# test dispersion 
testDispersion(sim)
# shows QQ plot, dispersion, outliers in 1 plot
testResiduals(sim)

# coefficient plot
library(sjPlot)
plot_model(carabidsiamodel, type = "std2", colors = "bw", title = "elytron length") +
  theme_classic() + geom_hline(yintercept = 0, lty = 2, colour = "gray")

# get standardized coefficients: "std" = forest-plot of standardized beta values
library(sjPlot)
get_model_data(carabidsiamodel, type = c("std"), show.df = TRUE)

# see together for plot
library(ggplot2)
ggplot(bothsia.data, aes(bodyd15N, median, color = species)) +
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


#### 4. BODY SIZE ####

# data - weevils (n = 1011)
weevil.data <- read.csv(file.choose(), header = TRUE, ",")
View(weevil.data)
str(weevil.data) 
summary(weevil.data)
names(weevil.data)
levels(weevil.data[,"species"]) # make sure order is S. tuberosus first (intercept)
weevil.data$species = factor(weevil.data$species, levels(weevil.data$species)[c(3,1,2)])

# data - sexed carabids (n = 295)
sexedcarabid.data <- read.csv(file.choose(), header = TRUE, ",")
View(sexedcarabid.data)
str(sexedcarabid.data)
summary(sexedcarabid.data)
names(sexedcarabid.data)
levels(sexedcarabid.data[,"species"]) # order is Z. matthewsii first (intercept)
sexedcarabid.data$species = factor(sexedcarabid.data$species, levels(sexedcarabid.data$species)[c(6,1,2,3,4,5)])

# (4a) weevil body size model (March 2020, include sampling round)
weevilmodel <- lmer(median ~ distance*species + round + (1|transect), data = weevil.data)
summary(weevilmodel)
df.residual(weevilmodel)

# OLD (4a) weevil body size model (December 2018)
weevilmodel <- lmer(median ~ distance*species + (1|transect), data = weevil.data)
summary(weevilmodel)
df.residual(weevilmodel)

# check residuals
ggplot(weevil.data, aes(x = fitted(weevilmodel), y = resid(weevilmodel))) +
  geom_point() +
  theme_classic() +
  geom_line(y=0, colour="red") +
  labs(x="Fitted values", y= "Residuals")
qqnorm(as.vector(resid(weevilmodel)))
qqline(as.vector(resid(weevilmodel)), col = "blue")

##USING DHARMa PACKAGE TO INTERPRET RESIDUALS (March 2020):
# set simulations constant 
set.seed(1)
# calculate scaled residuals
library(DHARMa)
sim <- simulateResiduals(fittedModel = weevilmodel, n = 500) # the calculated residuals are stored in sim$scaledResiduals
# plot the scaled residuals (Observed vs Expected)
plot(sim)
# plot residuals against the other predictors
plotResiduals(weevil.data$distance, sim$scaledResiduals) 
plotResiduals(weevil.data$species, sim$scaledResiduals)
plotResiduals(weevil.data$round, sim$scaledResiduals) 
# test outliers
testOutliers(sim)
# test dispersion 
testDispersion(sim)
# shows QQ plot, dispersion, outliers in 1 plot
testResiduals(sim)

# coefficient plot for weevil model
library(sjPlot)
plot_model(weevilmodel, type = "std2", colors = "bw") +
  theme_classic() + geom_hline(yintercept = 0, lty = 2, colour = "gray")

# standardized coefficients 
get_model_data(weevilmodel, type = c("std"))

# (4b) carabid body size model (March 2020, including sampling round and sex as additive effects)
carabidmodel <- lmer(median ~ distance*species + sex + round + (1|transect), data = sexedcarabid.data)
summary(carabidmodel)
df.residual(carabidmodel)

# OLD (4b) carabid body size model (December 2018, using sexed carabids n = 295 for a 3-way interaction)
carabidmodel <- lmer(median ~ distance*sex*species + (1|transect), data = sexedcarabid.data)
summary(carabidmodel)
df.residual(carabidmodel)

# check residuals
ggplot(sexedcarabid.data, aes(x = fitted(carabidmodel), y = resid(carabidmodel))) +
  geom_point() +
  theme_classic() +
  geom_line(y=0, colour="red") +
  labs(x="Fitted values", y= "Residuals")
qqnorm(as.vector(resid(carabidmodel)))
qqline(as.vector(resid(carabidmodel)), col = "blue")

##USING DHARMa PACKAGE TO INTERPRET RESIDUALS (March 2020):
# set simulations constant 
set.seed(1)
# calculate scaled residuals
library(DHARMa)
sim <- simulateResiduals(fittedModel = carabidmodel, n = 500) # the calculated residuals are stored in sim$scaledResiduals
# plot the scaled residuals (Observed vs Expected)
plot(sim)
# plot residuals against the other predictors
plotResiduals(sexedcarabid.data$distance, sim$scaledResiduals) 
plotResiduals(sexedcarabid.data$species, sim$scaledResiduals) 
plotResiduals(sexedcarabid.data$sex, sim$scaledResiduals) 
plotResiduals(sexedcarabid.data$round, sim$scaledResiduals) 
# test outliers
testOutliers(sim)
# test dispersion 
testDispersion(sim)
# shows QQ plot, dispersion, outliers in 1 plot
testResiduals(sim)

# coefficient plot 
library(sjPlot)
plot_model(carabidmodel, type = "std2", colors = "bw") +
  theme_classic() + geom_hline(yintercept = 0, lty = 2, colour = "gray")

# standardized coefficients 
get_model_data(carabidmodel, type = c("std"))

#### 5. POST HOC (DISCUSSION) ####

# data - weevil body nutrients (n = 28) - note same dataframe as 2. BODY NUTRIENTS
weevilbodyn.data <- read.csv(file.choose(), header = TRUE, ",")
View(weevilbodyn.data)
str(weevilbodyn.data)
summary(weevilbodyn.data)
names(weevilbodyn.data)

# data - carabid body nutrients (n = 30) - note same dataframe as 2. BODY NUTRIENTS
carabidbodyn.data <- read.csv(file.choose(), header = TRUE, ",")
View(carabidbodyn.data)
str(carabidbodyn.data)
summary(carabidbodyn.data)
names(carabidbodyn.data)

# data - combined for figure - note same dataframe as 2. BODY NUTRIENTS except don't want distance as a factor
bothbodynfig.data <- read.csv(file.choose(), header = TRUE, ",")
View(bothbodynfig.data)
str(bothbodynfig.data)
summary(bothbodynfig.data)
names(bothbodynfig.data)

# visualize raw weevil data 
library(ggplot2)
ggplot(weevilbodyn.data, aes(bodyd15N, bodyncombust)) +
  stat_smooth() + theme_classic() +
  theme(text = element_text(size = 16)) +
  theme(axis.text.x = element_text(size = 16, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.text.y = element_text(size = 16, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  coord_cartesian(ylim = c(NULL)) +
  labs(x = expression(paste("Curculionidae body δ"^"15"*"N"*" (‰)")), 
       y = "Body %N")

# (5a) weevil post hoc
weevilposthocmodel <- lmer(bodyncombust ~ bodyd15N + (1|transect), data = weevilbodyn.data)
summary(weevilposthocmodel)
df.residual(weeviladhocmodel)

# check residuals
ggplot(weevilbodyn.data, aes(x = fitted(weevilposthocmodel), y = resid(weevilposthocmodel))) +
  geom_point() +
  theme_classic() +
  geom_line(y=0, colour="red") +
  labs(x="Fitted values", y= "Residuals")
qqnorm(as.vector(resid(weevilposthocmodel)))
qqline(as.vector(resid(weevilposthocmodel)), col = "blue")

##USING DHARMa PACKAGE TO INTERPRET RESIDUALS (March 2020):
# set simulations constant 
set.seed(1)
# calculate scaled residuals
library(DHARMa)
sim <- simulateResiduals(fittedModel = weevilposthocmodel, n = 500) # the calculated residuals are stored in sim$scaledResiduals
# plot the scaled residuals (Observed vs Expected)
plot(sim)
# plot residuals against the other predictors
plotResiduals(weevilbodyn.data$bodyd15N, sim$scaledResiduals) 
# test outliers
testOutliers(sim)
# test dispersion 
testDispersion(sim)
# shows QQ plot, dispersion, outliers in 1 plot
testResiduals(sim)

# coefficient plot 
library(sjPlot)
plot_model(weevilposthocmodel, type = "std2", colors = "bw", title = "body %N") +
  theme_classic() + geom_hline(yintercept = 0, lty = 2, colour = "gray")

# get the standardized coefficients: "std" = forest-plot of standardized beta values
library(sjPlot)
get_model_data(weevilposthocmodel, type = c("std"), show.df = TRUE)         

# (5b) carabid post hoc
carabidposthocmodel <- lmer(bodyncombust ~ bodyd15N + (1|transect), data = carabidbodyn.data)
summary(carabidposthocmodel)
df.residual(carabidposthocmodel)

# check residuals
ggplot(carabidbodyn.data, aes(x = fitted(carabidposthocmodel), y = resid(carabidposthocmodel))) +
  geom_point() +
  theme_classic() +
  geom_line(y=0, colour="red") +
  labs(x="Fitted values", y= "Residuals")
qqnorm(as.vector(resid(carabidposthocmodel)))
qqline(as.vector(resid(carabidposthocmodel)), col = "blue")

##USING DHARMa PACKAGE TO INTERPRET RESIDUALS (March 2020):
# set simulations constant 
set.seed(1)
# calculate scaled residuals
library(DHARMa)
sim <- simulateResiduals(fittedModel = carabidposthocmodel, n = 500) # the calculated residuals are stored in sim$scaledResiduals
# plot the scaled residuals (Observed vs Expected)
plot(sim)
# plot residuals against the other predictors
plotResiduals(carabidbodyn.data$bodyd15N, sim$scaledResiduals) 
# test outliers
testOutliers(sim)
# test dispersion 
testDispersion(sim)
# shows QQ plot, dispersion, outliers in 1 plot
testResiduals(sim)

# coefficient plot 
library(sjPlot)
plot_model(carabidposthocmodel, type = "std2", colors = "bw", title = "body %N") +
  theme_classic() + geom_hline(yintercept = 0, lty = 2, colour = "gray")

# get the standardized coefficients: "std" = forest-plot of standardized beta values
library(sjPlot)
get_model_data(carabidposthocmodel, type = c("std"), show.df = TRUE)         

# see together for plot
library(ggplot2)
ggplot(bothbodynfig.data, aes(bodyd15N, bodyncombust, color = species)) +
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
