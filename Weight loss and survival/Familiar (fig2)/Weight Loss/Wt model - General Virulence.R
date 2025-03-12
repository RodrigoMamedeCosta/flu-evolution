## ----setup, include=FALSE---------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE)
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(dplyr)


## ---------------------------------------------------------------------------------------
alldata <- read.csv("groupedvirulence7days.csv", sep=',', header=T)


## ---------------------------------------------------------------------------------------
# Subset data
famtest <- subset(alldata, Familiarity == "FAM" & Treatment != "1XP" )

# Scale and center Day
famtest$Day <- scale(famtest$Day, center = TRUE, scale = TRUE)

# Hosts to iterate over
hosts <- c("BALBF", "BALBM", "BL6F")

# Open a connection to the output file
output_file <- "LMM_Results_Host_passaged_vs_unpassaged.txt"
sink(output_file)

cat("LMM results for ancestral vs evolved virus in each host")
# Iterate through each host as the reference level
for (host in hosts) {
  # Relevel Host
  famtest$Host <- relevel(factor(famtest$Host), ref = host)
  
  # Fit the model
  model <- lmer(Pct ~ Day * Status * Host + (Day | Line/ID), data = famtest, REML = TRUE)
  
  # Print header for clarity
  cat("\n\nResults with Host as Reference Level:", host, "\n")
  cat("---------------------------------------------------\n")
  
  # Output summary
  print(summary(model))
}

# Close the connection
sink()



## ---------------------------------------------------------------------------------------
# Subset data
famtest <- subset(alldata, Familiarity == "FAM" & Treatment != "1XP" & Status == "Passaged")

# Scale and center Day
famtest$Day <- scale(famtest$Day, center = TRUE, scale = TRUE)

# Hosts to iterate over
hosts <- c("BALBF", "BALBM", "BL6F")

# Open a connection to the output file
output_file <- "LMM_Results_Host_Differences.txt"
sink(output_file)

cat("LMM results for weight loss comparisons between hosts [ Pct ~ Day * Host  + (Day | Line/ID)]")

# Iterate through each host as the reference level
for (host in hosts) {
  # Relevel Host
  famtest$Host <- relevel(factor(famtest$Host), ref = host)
  
  # Fit the model
  model <- lmer(Pct ~ Day * Host + (Day | Line/ID), data = famtest, REML = TRUE)
  
  # Print header for clarity
  cat("\n\nResults with Host as Reference Level:", host, "\n")
  cat("---------------------------------------------------\n")
  
  # Output summary
  print(summary(model))
}

# Close the connection
sink()


