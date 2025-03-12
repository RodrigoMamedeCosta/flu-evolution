## ----setup, include=FALSE---------------------------------------------------------------
knitr::opts_chunk$set(warning = FALSE, include = FALSE)
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(MuMIn)


## ---------------------------------------------------------------------------------------
alldata <- read.csv("groupedvirulence7days.csv", sep=',', header=T)


## ----BALB_FAM vs BALB_UNFAM vs BL6FAM vs BL6UNFAM---------------------------------------
# Remove data that is not relevant
tradeoffs <- subset(alldata, Treatment != "STOCK" & Treatment != "1XP" & Sex == "F")

# Scale and center
tradeoffs$Day <- scale(tradeoffs$Day, center = TRUE, scale = TRUE)

# Create a new column with Host and familiarity
tradeoffs$HostFam <- paste0(tradeoffs$Host,"_",tradeoffs$Familiarity)
tradeoffs$HostFam <- as.factor(tradeoffs$HostFam)

# Levels to iterate over
host_fams <- c("BALBF_FAM", "BALBF_UNF", "BL6F_FAM", "BL6F_UNF")

# Open the file for writing all summaries
sink("LMM-Tradeoffs_summary.txt")

# Iterate over each level
for (level in host_fams) {
  
  # Relevel HostFam
  tradeoffs$HostFam <- relevel(factor(tradeoffs$HostFam), ref = level)
  
  # Run model
  model.tradeoffs <- lmer(Pct ~ Day * HostFam + (Day|Line/ID), data = tradeoffs, REML = TRUE)
  
  # Write the title indicating the level
  cat("\n\n##############################################\n")
  cat("            Leveled on", level, "\n")
  cat("##############################################\n\n")
  
  # Output summaries to the file
  print(summary(model.tradeoffs))
}

# Save and close the file
sink()



## ----BALB-adapted vs BL6-adapted--------------------------------------------------------
#Remove data that is not relevant
overallvir <- subset(alldata, Treatment !="STOCK" & Sex == "F" & Treatment !="1XP" )

#Level of 1X and familiar

overallvir$Host <- relevel(factor(overallvir$Host), ref = "BALBF")
overallvir$Familiarity <- relevel(factor(overallvir$Familiarity), ref = "FAM")

#Scale and center.
overallvir$Day <- scale(overallvir$Day, center = TRUE, scale = TRUE)

#Run model
model.overallvir <- lmer(Pct ~ Day * Host + (Day|Line/ID), data = overallvir, REML = TRUE)

summary(model.overallvir)

# Open the file for writing
sink("LMM-Overall_virulence_BALBF_vs_BL6F.txt")

# Output summaries to the file
summary(model.overallvir)
# Save
sink()



## ----STRAIN MORE SUSCEPTIBLE TO FAM - DELETE--------------------------------------------
#Remove data that is not relevant
overallvir <- subset(alldata, Treatment !="STOCK" & Treatment != "2X" & Treatment != "1XP" & Sex == "F")

#Scale and center.
overallvir$Day <- scale(overallvir$Day, center = TRUE, scale = TRUE)

#overallvir <- subset(overallvir, Familiarity == "FAM")
overallvir <- subset(overallvir, Familiarity == "UNF")


#Run model
model.overallsusc <- lmer(Pct ~ Day * HostStrain + (Day|Line/ID), data = overallvir, REML = TRUE)

summary(model.overallsusc)

# Open the file for writing
sink("LMM-1X&4X_susceptibility_BALBC_vs_BL6_Unfamiliar.txt")

# Output summaries to the file
summary(model.overallsusc)
# Save
sink()



## ----POOLING 1X vs 4X TESTS-------------------------------------------------------------
#Remove data that is not relevant
treattest <- subset(alldata, Treatment !="STOCK" & Treatment != "2X" & Sex == "F" & Treatment != "1XP")

# Subset to specify comparisons

 #treattest <- subset(treattest, Host == "BL6F")
 #treattest <- subset(treattest, HostStrain == "C57BL/6")
 treattest <- subset(treattest, Familiarity == "UNF")

#Level on 1X

treattest$Treatment <- relevel(factor(treattest$Treatment), ref = "1X")

#Scale and center.
treattest$Day <- scale(treattest$Day, center = TRUE, scale = TRUE)

#Run model
model.treat <- lmer(Pct ~ Day * Treatment + (Day|Line/ID), data = treattest, REML = TRUE)

# Open the file for writing
sink("LMM-Tradeoff-1X-vs-4X.txt")

# Output summaries to the file
summary(model.treat)
# Save
sink()



## ----POOLING 1X vs 1XP tests------------------------------------------------------------
#Remove data that is not relevant
poolingtest <- subset(alldata, Treatment !="STOCK" & Treatment != "2X" & Sex == "F" & Treatment != "4X")

# Subset to specify comparisons

 #poolingtest <- subset(poolingtest, Host == "BL6F")
 #poolingtest <- subset(poolingtest, HostStrain == "BALB/C")
 poolingtest <- subset(poolingtest, Familiarity == "UNF")

#Level on 1X

poolingtest$Treatment <- relevel(factor(poolingtest$Treatment), ref = "1X")

#Scale and center.
poolingtest$Day <- scale(poolingtest$Day, center = TRUE, scale = TRUE)

#Run model
model.pooling <- lmer(Pct ~ Day * Treatment + Host + (Day|Line/ID), data = poolingtest, REML = TRUE)

# Open the file for writing
sink("LMM-Tradeoff-1X-vs-1XP_plus_host.txt")

# Output summaries to the file
summary(model.pooling)
# Save
sink()



## ----1XP vs combined 1X 2X 4X-----------------------------------------------------------
#Remove data that is not relevant
poolingtest <- subset(alldata, Treatment !="STOCK" & Sex == "F")

# Subset to specify comparisons

 #poolingtest <- subset(poolingtest, Host == "BL6F")
 #poolingtest <- subset(poolingtest, HostStrain == "BALB/C")
 poolingtest <- subset(poolingtest, Familiarity == "UNF")

#Level on 1X

poolingtest$Treatment <- ifelse(poolingtest$Treatment %in% c("1X", "2X", "4X"), "1X2X4X", "1XP")

poolingtest$Treatment <- relevel(factor(poolingtest$Treatment), ref = "1X2X4X")

#Scale and center.
poolingtest$Day <- scale(poolingtest$Day, center = TRUE, scale = TRUE)

#Run model
model.pooling <- lmer(Pct ~ Day * Treatment + Host + (Day|Line/ID), data = poolingtest, REML = TRUE)

# Open the file for writing
sink("LMM-Tradeoff-1X2X4X-vs-1XP_plus_host.txt")

# Output summaries to the file
summary(model.pooling)
# Save
sink()



## ----OVERALL VIRULENCE DIFFERENCES------------------------------------------------------
#Remove data that is not relevant
poolingtest <- subset(alldata, Treatment !="STOCK" & Sex =="F")

# Subset to specify comparisons

 #poolingtest <- subset(poolingtest, Host == "BL6F")
 #poolingtest <- subset(poolingtest, HostStrain == "BALB/C")
 #poolingtest <- subset(poolingtest, Familiarity == "UNF")

#Level on 1X

poolingtest$Treatment <- relevel(factor(poolingtest$Treatment), ref = "1X")

#Scale and center.
poolingtest$Day <- scale(poolingtest$Day, center = TRUE, scale = TRUE)

#Run model
model.pooling <- lmer(Pct ~ Day * Treatment + Host + (Day|Line/ID), data = poolingtest, REML = TRUE)

# Open the file for writing
sink("LMM-ALL_Treatments_in_allHosts_Females_only.txt")

# Output summaries to the file
summary(model.pooling)
# Save
sink()


