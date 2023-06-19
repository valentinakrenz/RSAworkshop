# R script for analyzing trial-wise, roi-based RSA data
# Valentina Krenz 2023



# set paths ####

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Set the working directory to the script location

output_dir <- "./results"
input_dir <- "./results"
file_name <- "ERS"

# read in data ####
if(!require(readxl)){
  install.packages("readxl")
}
library(readxl) #prepare data #needs to be loaded after Rmisc
if(!require(tidyverse)){
  install.packages("tidyverse")
}
library(tidyverse) #prepare data #needs to be loaded after Rmisc

df <- read_excel(file.path(input_dir, paste0(file_name, ".xlsx"))) 

# run LMM ####
if(!require(lme4)){
  install.packages("lme4")
}
library(lme4) #for lmer / glmer -> mixed effects model

if(!require(lmerTest)){
  install.packages("lmerTest")
}
library(lmerTest) #show p_values in mixed effetcs model

sub_df <- subset(df, ROI == "AmyBil" | ROI == "HCBil") %>%
  mutate(ROI = factor(ROI, levels=c("HCBil", "AmyBil"), labels=c("hippocampus","amygdala")),
         emotion = factor(emotion, levels=c("neutral","negative"))
  )

lmm <- lmer(corr ~ emotion * ROI + (1 + emotion + ROI | sj) + (1 | item), data = sub_df)
summary(lmm)
anova(lmm)

step_model <- lmerTest::step(lmm, reduce.fixed=FALSE, reduce.random=TRUE)
step_model

lmm <- lmer(corr ~ emotion + ROI + (1 | item) + (1 | sj) + emotion:ROI, data = sub_df)
summary(lmm)

# post hoc tests
if(!require(emmeans)){
  install.packages("emmeans")
}
library(emmeans) #post-hoc tests for ANOVAs and  (generalized) linear mixed models
model.contrast <- emmeans(lmm, pairwise ~ ROI | emotion, lmer.df = "satterthwaite")
summary(model.contrast, adjust="FDR")

# plot without emotion factor
if(!require(sjPlot)){
  install.packages("sjPlot")
}
library(sjPlot)#plotting (generalized) linear mixed models 

p <- plot_model(lmm, type = "pred", terms = c("emotion","ROI"),
                show.data = FALSE, value.offset = TRUE, jitter = TRUE, 
                dot.size = 4, grid = FALSE, line.size = 2, 
                axis.title = c("stimulus emotionality", "Fisher-transformed r"),
                color=c("#0072B2","firebrick")) + theme_minimal()
p

# Save the plot in png format 
png_file_name <- paste0(file_name, "_allROIs.png")
png_file_path <- file.path(output_dir, png_file_name)
png(png_file_path)
print(p)
dev.off()

# run ANOVA #####

if(!require(afex)){
  install.packages("afex")
}
library(afex) #prepare data #needs to be loaded after Rmisc

agg_df <- aggregate(corr ~ emotion + sj + ROI, data = sub_df, FUN = mean) # aggregate over items

# run ANOVA for each region
model <- aov_ez(
  "sj" # subject indicator
  ,"corr"#define dependent variable
  ,agg_df
  ,within=c("emotion","ROI") # define within-factors
  #,between="group",
  ,anova_table="pes") # gives partial eta square

summary(model)

model.contrast <- emmeans(model, pairwise ~ ROI | emotion)
summary(model.contrast, adjust="FDR")

model.contrast <- emmeans(model, pairwise ~ emotion * ROI)
summary(model.contrast, adjust="bonferroni")



