## ----setup, include=FALSE---------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE)
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(dplyr)


## ----Load data--------------------------------------------------------------------------
alldata <- read.csv("groupedvirulence7days.csv", sep=',', header=T)


## ----Subset for including stock and excluding unfamiliar and 1XP:-----------------------
famtest <- subset(alldata, Familiarity=="FAM" & Treatment !="STOCK")

famtest <- subset(famtest, Treatment !="1XP")


## ----Leveling - relevel to test pooling on each host------------------------------------
famtest$Host <- relevel(factor(famtest$Host), ref = "BALBF")
famtest$Treatment <- relevel(factor(famtest$Treatment), ref = "1X")
#& Host != "BL6F"

#Scale and center
famtest$Day <- scale(famtest$Day, center = TRUE, scale = TRUE)


## ----Run model--------------------------------------------------------------------------

# Levels to iterate over
hosts <- c("BALBF", "BALBM", "BL6F")

# Open the file for writing all summaries
sink("LMM-Pooling_summary_treatment_plus_Host.txt")


# Iterate over each level
for (level in hosts) {
  
  # Relevel HostFam
  famtest$Host <- relevel(factor(famtest$Host), ref = level)
  
  # Run model
  model.pooling <- lmer(Pct ~ Day * Treatment + Host + (Day|Line/ID), data = famtest, REML = TRUE)
  
  # Write the title indicating the level
  cat("##############################################\n")
  cat("            Leveled on", level, "\n")
  cat("##############################################\n\n")
  
  # Output summaries to the file
  print(summary(model.pooling))
}

# Save and close the file
sink()



## ----1X 2X 4X in each host separately---------------------------------------------------
# Levels to iterate over
hosts <- c("BALBF", "BALBM", "BL6F")

# Open the file for writing all summaries
sink("LMM-Pooling_summary_treatment_by_Host.txt")


# Iterate over each level
for (level in hosts) {
  
  # Relevel HostFam
  famtest2 <- subset(famtest, Host == level)
  
  # Run model
  model.pooling <- lmer(Pct ~ Day * Treatment + (Day|Line/ID), data = famtest2, REML = TRUE)
  
  # Write the title indicating the level
  cat("##############################################\n")
  cat("            Pooling differences for", level, "\n")
  cat("##############################################\n\n")
  
  # Output summaries to the file
  print(summary(model.pooling))
}

# Save and close the file
sink()



## ----1XP vs 1X tests, all strains, familiar.--------------------------------------------
fam1xp <- subset(alldata, Familiarity=="FAM" & Treatment != "2X" & Treatment != "4X" & Treatment !="STOCK")

fam1xp$Host <- relevel(factor(fam1xp$Host), ref = "BALBF")
fam1xp$Treatment <- relevel(factor(fam1xp$Treatment), ref = "1X")
#& Host != "BL6F"

#Scale and center
fam1xp$Day <- scale(fam1xp$Day, center = TRUE, scale = TRUE)


# Levels to iterate over
hosts <- c("BALBF", "BALBM", "BL6F")

# Open the file for writing all summaries
sink("LMM-Pooling_summary_1X_vs_1XP_plus_host.txt")


# Iterate over each level
for (level in hosts) {
  
  # Relevel HostFam
  fam1xp$Host <- relevel(factor(fam1xp$Host), ref = level)
  
  # Run model
  model.1xp <- lmer(Pct ~ Day * Treatment + Host + (Day|Line/ID), data = fam1xp, REML = TRUE)
  
  # Write the title indicating the level
  cat("##############################################\n")
  cat("            Leveled on", level, "\n")
  cat("##############################################\n\n")
  
  # Output summaries to the file
  print(summary(model.1xp))
}

# Save and close the file
sink()


## ----By host----------------------------------------------------------------------------
fam1xp <- subset(alldata, Familiarity=="FAM" & Treatment != "2X" & Treatment != "4X" & Treatment !="STOCK")

fam1xp$Host <- relevel(factor(fam1xp$Host), ref = "BALBF")
fam1xp$Treatment <- relevel(factor(fam1xp$Treatment), ref = "1X")
#& Host != "BL6F"

#Scale and center
fam1xp$Day <- scale(fam1xp$Day, center = TRUE, scale = TRUE)


# Levels to iterate over
hosts <- c("BALBF", "BALBM", "BL6F")

# Open the file for writing all summaries
sink("LMM-Pooling_summary_1X_vs_1XP_by_host.txt")


# Iterate over each level
for (level in hosts) {
  
  # Relevel HostFam
   fam1xp2 <- subset(fam1xp, Host == level)
  
  # Run model
  model.1xp <- lmer(Pct ~ Day * Treatment + (Day|Line/ID), data = fam1xp2, REML = TRUE)
  
  # Write the title indicating the level
  cat("##############################################\n")
  cat("            Leveled on", level, "\n")
  cat("##############################################\n\n")
  
  # Output summaries to the file
  print(summary(model.1xp))
}

# Save and close the file
sink()

