####R Script for Beetle Paper
####Originally written in fall 2018 by NR; updated in spring 2020 by NR and AD
####For submission to Ecology & Evolution

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
names(repeat.data)

# (0) repeatability model
repeatmodel <- lmer(length ~ 1 + (1|insectid), data = repeat.data)
df.residual(repeatmodel, type = c("lmer")) #this needs to be corrected, df can't be 4011
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

# data - soil (n = 30)
soiln.data <- read.csv(file.choose(), header = TRUE, ",")
View(soiln.data)
str(soiln.data)
summary(soiln.data)

soilnfig.data <- read.csv(file.choose(), header = TRUE, ",")
soilnfig.data$distance <- as.factor(soilnfig.data$distance)
str(soilnfig.data)

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
df.residual(soilnmodel, type = c("lmer")) 
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
summary(weevilbodyn.data)

weevilbodynfig.data <- read.csv(file.choose(), header = TRUE, ",")
weevilbodynfig.data$distance <- as.factor(weevilbodynfig.data$distance)
str(weevilbodynfig.data)

# data - carabid body nutrients (n = 30)
carabidbodyn.data <- read.csv(file.choose(), header = TRUE, ",")
View(carabidbodyn.data)
summary(carabidbodyn.data)

carabidbodynfig.data <- read.csv(file.choose(), header = TRUE, ",")
carabidbodynfig.data$distance <- as.factor(carabidbodynfig.data$distance)
str(carabidbodynfig.data)

# data - combined for figure
bothbodynfig.data <- read.csv(file.choose(), header = TRUE, ",")
bothbodynfig.data$distance <- as.factor(bothbodynfig.data$distance)
str(bothbodynfig.data)

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
weevilbodynmodel <- lmer(bodyd15N ~ distance + (1|transect), data = weevilbodyn.data)
df.residual(weevilbodynmodel)

# check residuals
ggplot(weevilbodyn.data, aes(x = fitted(weevilbodynmodel), y = resid(weevilbodynmodel))) +
  geom_point() +
  theme_classic() +
  geom_line(y=0, colour="red") +
  labs(x="Fitted values", y= "Residuals")
qqnorm(as.vector(resid(weevilbodynmodel)))
qqline(as.vector(resid(weevilbodynmodel)), col = "blue")

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
ggplot(carabidbodynfig.data, aes(group = distance, y = bodyd15N)) +
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
df.residual(carabidbodynmodel)

# check residuals
ggplot(carabidbodyn.data, aes(x = fitted(carabidbodynmodel), y = resid(carabidbodynmodel))) +
  geom_point() +
  theme_classic() +
  geom_line(y=0, colour="red") +
  labs(x="Fitted values", y= "Residuals")
qqnorm(as.vector(resid(carabidbodynmodel)))
qqline(as.vector(resid(carabidbodynmodel)), col = "blue")

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
summary(weevilsia.data)

# data - carabid SIA body size (1 indiv. per sample) (n = 30)
carabidsia.data <- read.csv(file.choose(), header = TRUE, ",")
View(carabidsia.data)
summary(carabidsia.data)

# data - both combined for figure
bothsia.data <- read.csv(file.choose(), header = TRUE, ",")
View(bothsia.data)
summary(bothsia.data)

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
df.residual(weevilsiamodel)

# check residuals
ggplot(weevilsia.data, aes(x = fitted(weevilsiamodel), y = resid(weevilsiamodel))) +
  geom_point() +
  theme_classic() +
  geom_line(y=0, colour="red") +
  labs(x="Fitted values", y= "Residuals")
qqnorm(as.vector(resid(weevilsiamodel)))
qqline(as.vector(resid(weevilsiamodel)), col = "blue")

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
df.residual(carabidsiamodel)

# check residuals
ggplot(carabidsia.data, aes(x = fitted(carabidsiamodel), y = resid(carabidsiamodel))) +
  geom_point() +
  theme_classic() +
  geom_line(y=0, colour="red") +
  labs(x="Fitted values", y= "Residuals")
qqnorm(as.vector(resid(carabidsiamodel)))
qqline(as.vector(resid(carabidsiamodel)), col = "blue")

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
levels(weevil.data[,"species"]) # make sure order is S. tuberosus first (intercept)
weevil.data$species = factor(weevil.data$species, levels(weevil.data$species)[c(3,1,2)])

# data - sexed carabids (n = 295)
sexedcarabid.data <- read.csv(file.choose(), header = TRUE, ",")
View(sexedcarabid.data)
str(sexedcarabid.data)
levels(sexedcarabid.data[,"species"]) # order is Z. matthewsii first (intercept)
sexedcarabid.data$species = factor(sexedcarabid.data$species, levels(sexedcarabid.data$species)[c(6,1,2,3,4,5)])

# (4a) weevil body size model (March 2020, include sampling round)
weevilmodel <- lmer(median ~ distance*species + round + (1|transect), data = weevil.data)
df.residual(weevilmodel)

# OLD (4a) weevil body size model (December 2018)
weevilmodel <- lmer(median ~ distance*species + (1|transect), data = weevil.data)
df.residual(weevilmodel)

# check residuals
ggplot(weevil.data, aes(x = fitted(weevilmodel), y = resid(weevilmodel))) +
  geom_point() +
  theme_classic() +
  geom_line(y=0, colour="red") +
  labs(x="Fitted values", y= "Residuals")
qqnorm(as.vector(resid(weevilmodel)))
qqline(as.vector(resid(weevilmodel)), col = "blue")

# coefficient plot for weevil model
library(sjPlot)
plot_model(weevilmodel, type = "std2", colors = "bw") +
  theme_classic() + geom_hline(yintercept = 0, lty = 2, colour = "gray")

# standardized coefficients 
get_model_data(weevilmodel, type = c("std"))

# (4b) carabid body size model (March 2020, including sampling round and sex as additive effects)
carabidmodel <- lmer(median ~ distance*species + sex + round + (1|transect), data = sexedcarabid.data)
df.residual(carabidmodel)

# OLD (4b) carabid body size model (December 2018, using sexed carabids n = 295 for a 3-way interaction)
carabidmodel <- lmer(median ~ distance*sex*species + (1|transect), data = sexedcarabid.data)
df.residual(carabidmodel)

# check residuals
ggplot(sexedcarabid.data, aes(x = fitted(carabidmodel), y = resid(carabidmodel))) +
  geom_point() +
  theme_classic() +
  geom_line(y=0, colour="red") +
  labs(x="Fitted values", y= "Residuals")
qqnorm(as.vector(resid(carabidmodel)))
qqline(as.vector(resid(carabidmodel)), col = "blue")

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
names(weevilbodyn.data)

# data - carabid body nutrients (n = 30) - note same dataframe as 2. BODY NUTRIENTS
carabidbodyn.data <- read.csv(file.choose(), header = TRUE, ",")
View(carabidbodyn.data)
names(carabidbodyn.data)

# data - combined for figure - note same dataframe as 2. BODY NUTRIENTS except don't want distance as a factor
bothbodynfig.data <- read.csv(file.choose(), header = TRUE, ",")
View(bothbodyn.data)
str(bothbodyn.data)

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
weeviladhocmodel <- lmer(bodyncombust ~ bodyd15N + (1|transect), data = weevilbodyn.data)
df.residual(weeviladhocmodel)

# check residuals
ggplot(weevilbodyn.data, aes(x = fitted(weeviladhocmodel), y = resid(weeviladhocmodel))) +
  geom_point() +
  theme_classic() +
  geom_line(y=0, colour="red") +
  labs(x="Fitted values", y= "Residuals")
qqnorm(as.vector(resid(weeviladhocmodel)))
qqline(as.vector(resid(weeviladhocmodel)), col = "blue")

# coefficient plot 
library(sjPlot)
plot_model(weeviladhocmodel, type = "std2", colors = "bw", title = "body %N") +
  theme_classic() + geom_hline(yintercept = 0, lty = 2, colour = "gray")

# get the standardized coefficients: "std" = forest-plot of standardized beta values
library(sjPlot)
get_model_data(weeviladhocmodel, type = c("std"), show.df = TRUE)         

# (5b) carabid post hoc
carabidadhocmodel <- lmer(bodyncombust ~ bodyd15N + (1|transect), data = carabidbodyn.data)
df.residual(carabidadhocmodel)

# check residuals
ggplot(carabidbodyn.data, aes(x = fitted(carabidadhocmodel), y = resid(carabidadhocmodel))) +
  geom_point() +
  theme_classic() +
  geom_line(y=0, colour="red") +
  labs(x="Fitted values", y= "Residuals")
qqnorm(as.vector(resid(carabidadhocmodel)))
qqline(as.vector(resid(carabidadhocmodel)), col = "blue")

# coefficient plot 
library(sjPlot)
plot_model(carabidadhocmodel, type = "std2", colors = "bw", title = "body %N") +
  theme_classic() + geom_hline(yintercept = 0, lty = 2, colour = "gray")

# get the standardized coefficients: "std" = forest-plot of standardized beta values
library(sjPlot)
get_model_data(carabidadhocmodel, type = c("std"), show.df = TRUE)         

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



