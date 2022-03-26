# set wd, load data, install ggplot2
rm(list=ls())
install.packages("readr")
install.packages("ggplot2")
install.packages("plyr")
install.packages("dplyr")
install.packages("tidyr")
install.packages("car")
install.packages("MuMIn")
install.packages("ggeffects")
library(readr)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(lme4)
library(car)
library(MuMIn)
library(tidyverse)
library(ggeffects)

# Set working directory
setwd("/Users/lucy/Desktop/Imperial/Winter/Valewood Preparation/Valewood 2/Data/")

# summary function
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, !!measurevar := mean)
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


### Tea bags
# load data
tea_k <- read.csv(file = "teabags_k_values.csv")
tea_k$Site <- as.factor(tea_k$Site)

# split into green and rooibos
tea_green <- filter(tea_k, Type == "Green")
tea_r <- filter(tea_k, Type == "Rooibos")

# simplify df and pivot
tea_green <- dplyr::select(tea_green, Site, Location, Type, k_R1, k_R2, k_R3, k_R4) %>% 
  dplyr::rename(R1 = k_R1,
                R2 = k_R2,
                R3 = k_R3,
                R4 = k_R4,
  ) %>% 
  pivot_longer(
    cols = starts_with("R"), 
    names_to = "Replicate", 
    values_to = "k",
    values_drop_na = TRUE
  )

tea_r <- dplyr::select(tea_r, Site, Location, Type, k_R1, k_R2, k_R3, k_R4) %>% 
  dplyr::rename(
    R1 = k_R1,
    R2 = k_R2,
    R3 = k_R3,
    R4 = k_R4,
  ) %>% 
  pivot_longer(
    cols = starts_with("R"), 
    names_to = "Replicate", 
    values_to = "k",
    values_drop_na = TRUE
  )

## Control & Enclosure analysis
# create separate dataframe for 'control & enclosure' analysis
green_CE <- dplyr::filter(tea_green, Location == "Control" | Location == "Enclosure")
rooibos_CE <- dplyr::filter(tea_r, Location == "Control" | Location == "Enclosure")

# LMMs
# Green
hist(green_CE$k)
lmer_g_CE <- lmer(k ~ Location + (1|Site) + (1|Replicate), data = green_CE)
summary(lmer_g_CE)
Anova(lmer_g_CE, test.statistic = "F")
# check model assumptions
plot(lmer_g_CE)

# random-slope and random-intercept model
lmer_g_CE_Slope <- lmer(k ~ Location + (1 + Location|Site) + (1 + Location|Replicate), data = green_CE)
summary(lmer_g_CE_Slope)
anova(lmer_g_CE, lmer_g_CE_Slope)
# no significant difference between models
AIC(lmer_g_CE, lmer_g_CE_Slope, k = T)
# intercept-only model has lower AIC

# pairwise
lmer_g_sites <- lmer(k ~ Site + (1|Replicate), data = green_CE)
em <- emmeans(lmer_g_sites, c("Site"))
contrast(em, method = "pairwise")

# Rooibos
hist(rooibos_CE$k)
lmer_r_CE <- lmer(k ~ Location + (1|Site) + (1|Replicate), data = rooibos_CE)
summary(lmer_r_CE)
Anova(lmer_r_CE, test.statistic = "F")

# check model assumptions
plot(lmer_r_CE)

# random-slope and random-intercept model
lmer_r_CE_Slope <- lmer(k ~ Location + (1 + Location|Site) + (1 + Location|Replicate), data = rooibos_CE)
summary(lmer_r_CE_Slope)
anova(lmer_r_CE, lmer_r_CE_Slope)
# no significant difference between models
AIC(lmer_r_CE, lmer_r_CE_Slope, k = 2)
# intercept-only model has lower AI

# pairwise
lmer_r_sites <- lmer(k ~ Site + (1|Replicate), data = rooibos_CE)
em <- emmeans(lmer_r_sites, c("Site"))
contrast(em, method = "pairwise")

## enclosure & downstream
# load data and split into Green and Rooibos for  separate analysis
tea_green_D <- filter(tea_green, Location == "Downstream")
tea_r_D <- filter(tea_r, Location == "Downstream")
tea_green_E <- filter(tea_green, Location == "Enclosure")
tea_r_E <- filter(tea_r, Location == "Enclosure")

# add distance
dist <- c(3.2, 3.2, 3.2, 3.2, 1.6, 1.6, 1.6, 1.6, 0.8, 0.8, 0.8, 0.8, 0.4, 0.4, 0.4, 0.4, 0.2, 0.2, 0.2, 0.2)
tea_green_D$Distance <- dist
tea_r_D$Distance <- dist
tea_green_E$Distance <- 0
tea_r_E$Distance <- 0

# combine E and D
tea_green_ED <- bind_rows(tea_green_D, tea_green_E)
tea_r_ED <- bind_rows(tea_r_D, tea_r_E)

# LMMs
# Green
lm_green_ED <- lmer(k ~ Distance + (1|Replicate), data = tea_green_ED)
summary(lm_green_ED)
Anova(lm_green_ED, test.statistic = "F")
# check model assumptions
plot(lm_green_ED)

# random-slope and random-intercept model
lm_green_ED_Slope <- lmer(k ~ Distance + (1 + Distance|Replicate), data = tea_green_ED)
summary(lm_green_ED_Slope)
anova(lm_green_ED, lm_green_ED_Slope)
# no significant difference between models
AIC(lm_green_ED, lm_green_ED_Slope, k = 2)
# intercept-only model has lower AIC

# Rooibos
lm_r_ED <- lmer(k ~ Distance + (1|Replicate), data = tea_r_ED)
summary(lm_r_ED)
Anova(lm_r_ED, test.statistic = "F")
# check model assumptions
plot(lm_r_ED)

# random-slope and random-intercept model
lm_r_ED_Slope <- lmer(k ~ Distance + (1 + Distance|Replicate), data = tea_r_ED)
summary(lm_r_ED_Slope)
anova(lm_r_ED, lm_r_ED_Slope)
# no significant difference between models
AIC(lm_r_ED, lm_r_ED_Slope, k = 2)
# intercept-only model has lower AIC

### Leaf Litter Analysis ###
# load data
leaves <- read.csv(file = "all_leaves.csv")
leaves$Site <- as.factor(leaves$Site)
leaves$Location <- as.factor(leaves$Location)
leaves$Leaf <- as.factor(leaves$Leaf)
leaves$Mesh <- as.factor(leaves$Mesh)

# simplify df and pivot
leaves <- dplyr::select(leaves, Site, Location, Leaf, Mesh, k_R1, k_R2, k_R3, k_R4) %>% 
  dplyr::rename(R1 = k_R1,
                R2 = k_R2,
                R3 = k_R3,
                R4 = k_R4,
  ) %>% 
  pivot_longer(
    cols = starts_with("R"), 
    names_to = "Replicate", 
    values_to = "k",
    values_drop_na = TRUE
  )

# create separate dataframes for alder and oak, and control & enclosure and enclosure & downstream for separate analysis
leaves_CE <- dplyr::filter(leaves, Location == "Control" | Location == "Enclosure")
alder_CE <- dplyr::filter(leaves_CE, Leaf == "Alder")
oak_CE <- dplyr::filter(leaves_CE, Leaf == "Oak")

# LMEs for 'control and enclosure'
# Alder
lmer_alder_CE <- lmer(k ~ Location*Mesh + (1|Site/Replicate), data = alder_CE)
summary(lmer_alder_CE)
Anova(lmer_alder_CE, test.statistic = "F")
# check model assumptions
plot(lmer_alder_CE)

# model for alder fine only
alder_CE_f <- filter(alder_CE, Mesh == "Fine")
lmer_alder_CE_f <- lmer(k ~ Location + (1|Replicate), data = alder_CE_f)
summary(lmer_alder_CE_f)
Anova(lmer_alder_CE_f, test.statistic = "F")

# model for alder coarse only
alder_CE_c <- filter(alder_CE, Mesh == "Coarse")
lmer_alder_CE_c <- lmer(k ~ Location + (1|Replicate), data = alder_CE_c)
summary(lmer_alder_CE_c)
Anova(lmer_alder_CE_c, test.statistic = "F")

# Oak
lmer_oak_CE <- lmer(k ~ Location*Mesh + (1|Site/Replicate), data = oak_CE)
summary(lmer_oak_CE)
Anova(lmer_oak_CE, test.statistic = "F")
# check model assumptions
plot(lmer_oak_CE)

# model for oak fine only
oak_CE_f <- filter(oak_CE, Mesh == "Fine")
lmer_oak_CE_f <- lmer(k ~ Location + (1|Replicate), data = oak_CE_f)
summary(lmer_oak_CE_f)
Anova(lmer_oak_CE_f, test.statistic = "F")

# model for oak coarse only
oak_CE_c <- filter(oak_CE, Mesh == "Coarse")
lmer_oak_CE_c <- lmer(k ~ Location + (1|Replicate), data = oak_CE_c)
summary(lmer_oak_CE_c)
Anova(lmer_oak_CE_c, test.statistic = "F")

# enclosure & downstream
oak_D <- filter(leaves, Location == "Downstream", Leaf == "Oak")
alder_D <- filter(leaves, Location == "Downstream", Leaf == "Alder")
oak_E <- filter(leaves, Location == "Enclosure", Leaf == "Oak")
alder_E <- filter(leaves, Location == "Enclosure", Leaf == "Alder")

# add distance
dist <- c(3.2, 3.2, 3.2, 3.2, 1.6, 1.6, 1.6, 1.6, 0.8, 0.8, 0.8, 0.8, 0.4, 0.4, 0.4, 0.4, 0.2, 0.2, 0.2, 0.2,3.2, 3.2, 3.2, 3.2, 1.6, 1.6, 1.6, 1.6, 0.8, 0.8, 0.8, 0.8, 0.4, 0.4, 0.4, 0.4, 0.2, 0.2, 0.2, 0.2)
oak_D$Distance <- dist
alder_D$Distance <- dist
oak_E$Distance <- 0
alder_E$Distance <- 0

# combine E and D
oak_ED <- bind_rows(oak_D, oak_E)
alder_ED <- bind_rows(alder_D, alder_E)

# Stats
# Oak
lm_oak_ED <- lmer(k ~ Distance*Mesh + (1|Replicate), data = oak_ED)
summary(lm_oak_ED)
Anova(lm_oak_ED, test.statistic = "F")
# check model assumptions
plot(lm_oak_ED)

# random-slope and random-intercept model
lm_oak_ED_Slope <- lmer(k ~ Distance*Mesh + (1 + Distance|Replicate), data = oak_ED)
summary(lm_oak_ED_Slope)
anova(lm_oak_ED, lm_oak_ED_Slope)
# intercept-only model has lower AIC

# Alder
lm_alder_ED <- lmer(k ~ Distance*Mesh + (1|Replicate), data = alder_ED)
summary(lm_alder_ED)
Anova(lm_alder_ED, test.statistic = "F")
# check model assumptions
plot(lm_alder_ED)
# qq plot
qqnorm(resid(lm_alder_ED))
qqline(resid(lm_alder_ED))
# skewed right residuals

# transform k values
alder_ED$Kt <- log10(alder_ED$k)
# refit LMM
lm_alder_ED <- lmer(Kt ~ Distance*Mesh + (1|Replicate), data = alder_ED, na.action = na.omit)
summary(lm_alder_ED)
Anova(lm_alder_ED, test.statistic = "F")
# check model assumptions
plot(lm_alder_ED)

# random-slope and random-intercept model
lm_alder_ED_Slope <- lmer(Kt ~ Distance*Mesh + (1 + Distance|Replicate), data = alder_ED, na.action = na.omit)
summary(lm_alder_ED_Slope)
anova(lm_alder_ED, lm_alder_ED_Slope)
# models not significantly different, intercept-only model has lower AIC


### invertebrates ###
# load data and format
inverts1 <- read.csv(file = "leaf_k_invertebrates.csv")
inverts1$Site <- as.factor(inverts1$Site)
inverts1$all_caddis <- inverts1$cased_caddis + inverts1$caseless_caddis
inverts1 <- dplyr::select(inverts1, Site, Location, Leaf, Mesh, Replicate, k, Gammarus, cased_caddis, caseless_caddis, large_bodied, Xt)
inverts1$all_caddis <- inverts1$cased_caddis + inverts1$caseless_caddis

# new column with correction to give abundance per unit mass of leaf litter remaining in bag
inverts1$abundance_gammarus <- inverts1$Gammarus/inverts1$Xt
inverts1$abundance_all_caddis <- inverts1$all_caddis/inverts1$Xt
inverts1$abundance_large_bodied <- inverts1$large_bodied/inverts1$Xt

# remove 0s
inverts1[inverts1==0] <- NA

# average for each site
inverts_avg <- inverts1 %>%
  group_by(Site, Leaf, Mesh) %>%
  summarise_at(vars(-Location,-Replicate), funs(mean(., na.rm=TRUE)))

# C and E
inverts_leaves <- filter(inverts1, Mesh == "Coarse", Location == "Control" | Location == "Enclosure")

lmer_gammarus_leaves <- lmer(k~abundance_gammarus*Leaf*Location + (1|Site/Replicate), data = inverts_leaves, na.action = na.omit)
summary(lmer_gammarus_leaves)
Anova(lmer_gammarus_leaves)
# check model assumptions
plot(lmer_gammarus_leaves)




######### Valewood 1 (V1) analysis - for Supplementary Materials
# load data
v1_tea_k <- read.csv(file = "v1_teabags.csv")
v1_tea_k <- v1_tea_k %>% drop_na()
v1_tea_k$Site <- as.factor(v1_tea_k$Site)

# split into green and rooibos
v1_tea_green <- filter(v1_tea_k, Type == "Green")
v1_tea_r <- filter(v1_tea_k, Type == "Rooibos")

# simplify df and pivot
v1_tea_green <- dplyr::select(v1_tea_green, Site, Location, Type, k)
v1_tea_r <- dplyr::select(v1_tea_r, Site, Location, Type, k)

# combine tea_r and tea_green
v1_tea_all <- rbind(v1_tea_green, v1_tea_r)

# LMMs for control and enclosure
# Green
# histogram of response variable
v1_green_CE <- filter(v1_tea_green, Location !="Downstream")
v1_lm_g_CE <- lm(k ~ Location, data = v1_green_CE)
summary(v1_lm_g_CE)
# check model assumptions
plot(v1_lm_g_CE)

# Rooibos
v1_rooibos_CE <- filter(v1_tea_r, Location !="Downstream")
v1_lm_r_CE <- lm(k ~ Location, data = v1_rooibos_CE)
summary(v1_lm_r_CE)
# check model assumptions
plot(v1_lm_r_CE)

## enclosure & downstream
# load data and split into Green and Rooibos for  separate analysis
v1_tea_green_D <- filter(v1_tea_green, Location == "Downstream")
v1_tea_r_D <- filter(v1_tea_r, Location == "Downstream")
v1_tea_green_E <- filter(v1_tea_green, Location == "Enclosure")
v1_tea_r_E <- filter(v1_tea_r, Location == "Enclosure")

# add distance
dist <- c(1.6, 0.8, 0.4, 0.2)
v1_tea_green_D$Distance <- dist
v1_tea_r_D$Distance <- dist
v1_tea_green_E$Distance <- 0
v1_tea_r_E$Distance <- 0

# combine E and D
v1_tea_green_ED <- bind_rows(v1_tea_green_D, v1_tea_green_E)
v1_tea_r_ED <- bind_rows(v1_tea_r_D, v1_tea_r_E)

# LMMs
# Green
lm_v1_green_ED <- lm(k ~ Distance, data = v1_tea_green_ED)
summary(lm_v1_green_ED)

# Rooibos
lm_v1_r_ED <- lm(k ~ Distance, data = v1_tea_r_ED)
summary(lm_v1_r_ED)

########## leaf analysis

# load data
v1_leaves <- read.csv(file = "valewood_1_leaves.csv")
v1_leaves$Site <- as.factor(v1_leaves$Site)
v1_leaves$Location <- as.factor(v1_leaves$Location)
v1_leaves$Leaf <- as.factor(v1_leaves$Leaf)
v1_leaves$Mesh <- as.factor(v1_leaves$Mesh)

# simplify df and pivot
v1_leaves <- dplyr::select(v1_leaves, Site, Location, Leaf, Mesh, k_R1, k_R2, k_R3, k_R4) %>% 
  dplyr::rename(R1 = k_R1,
                R2 = k_R2,
                R3 = k_R3,
                R4 = k_R4,
  ) %>% 
  pivot_longer(
    cols = starts_with("R"), 
    names_to = "Replicate", 
    values_to = "k",
    values_drop_na = TRUE
  )

# create separate dataframes for alder and oak, and control & enclosure and enclosure & downstream for separate analysis
v1_leaves_CE <- dplyr::filter(v1_leaves, Location == "Control" | Location == "Enclosure")
v1_alder_CE <- dplyr::filter(v1_leaves_CE, Leaf == "Alder")
v1_oak_CE <- dplyr::filter(v1_leaves_CE, Leaf == "Oak")

# LMEs for 'control and enclosure'
# Alder
v1_alder_CE <- as.tibble(v1_alder_CE)
v1_lmer_alder_CE <- lmer(k ~ Location*Mesh + (1|Replicate), data = v1_alder_CE)
summary(v1_lmer_alder_CE)
Anova(v1_lmer_alder_CE, test.statistic = "F")

# check model assumptions
plot(v1_lmer_alder_CE)
# qq plot
qqnorm(resid(v1_lmer_alder_CE))
qqline(resid(v1_lmer_alder_CE))
# normal residuals

# Oak
v1_lmer_oak_CE <- lmer(k ~ Location*Mesh + (1|Site/Replicate), data = v1_oak_CE)
summary(v1_lmer_oak_CE)
Anova(v1_lmer_oak_CE, test.statistic = "F")

# check model assumptions
plot(v1_lmer_oak_CE)
# qq plot
qqnorm(resid(v1_lmer_oak_CE))
qqline(resid(v1_lmer_oak_CE))
# normal residuals

# enclosure & downstream
v1_oak_D <- filter(v1_leaves, Location == "Downstream", Leaf == "Oak")
v1_alder_D <- filter(v1_leaves, Location == "Downstream", Leaf == "Alder")
v1_oak_E <- filter(v1_leaves, Location == "Enclosure", Leaf == "Oak")
v1_alder_E <- filter(v1_leaves, Location == "Enclosure", Leaf == "Alder")

# add distance
dist <- c(1.6, 1.6, 1.6, 1.6, 0.8, 0.8, 0.8, 0.8, 0.4, 0.4, 0.4, 0.4, 0.2, 0.2, 0.2, 0.2, 1.6, 1.6, 1.6, 1.6, 0.8, 0.8, 0.8, 0.8, 0.4, 0.4, 0.4, 0.4, 0.2, 0.2, 0.2, 0.2)
v1_oak_D$Distance <- dist
v1_alder_D$Distance <- dist
v1_oak_E$Distance <- 0
v1_alder_E$Distance <- 0

# combine E and D
v1_oak_ED <- bind_rows(v1_oak_D, v1_oak_E)
v1_alder_ED <- bind_rows(v1_alder_D, v1_alder_E)

# Stats
# Oak
lm_v1_oak_ED <- lmer(k ~ Distance*Mesh + (1|Replicate), data = v1_oak_ED)
summary(lm_v1_oak_ED)
Anova(lm_v1_oak_ED, test.statistic = "F")
# check model assumptions
plot(lm_v1_oak_ED)

# random-slope and random-intercept model
lm_v1_oak_ED_Slope <- lmer(k ~ Distance*Mesh + (1 + Distance|Replicate), data = v1_oak_ED)
summary(lm_v1_oak_ED_Slope)
anova(lm_v1_oak_ED, lm_v1_oak_ED_Slope)
# intercept-only model has lower AIC

# Alder
lm_v1_alder_ED <- lmer(k ~ Distance*Mesh + (1|Replicate), data = v1_alder_ED)
summary(lm_v1_alder_ED)
Anova(lm_v1_alder_ED, test.statistic = "F")
# check model assumptions
plot(lm_v1_alder_ED)

### invertebrates ###
# load data and format
v1_inverts <- read.csv(file = "v1_iverts.csv")
v1_inverts$Site <- as.factor(v1_inverts$Site)
v1_inverts$all_caddis <- v1_inverts$cased_caddis + v1_inverts$caseless_caddis

# new column with correction to give abundance pr unit mass of leaf litter remaining in bag
v1_inverts$abundance_gammarus <- v1_inverts$Gammarus/v1_inverts$Xt
v1_inverts$abundance_all_caddis <- v1_inverts$all_caddis/v1_inverts$Xt
v1_inverts$abundance_large_bodied <- v1_inverts$large_bodied/v1_inverts$Xt

# remove 0s
v1_inverts[v1_inverts==0] <- NA
v1_inverts_leaves <- filter(v1_inverts, Mesh == "Coarse", Location == "Control" | Location == "Enclosure")

lmer_gammarus_leaves <- lmer(k~abundance_gammarus*Leaf*Location + (1|Site/Replicate), data = v1_inverts_leaves, na.action = na.omit)
summary(lmer_gammarus_leaves)
Anova(lmer_gammarus_leaves)

