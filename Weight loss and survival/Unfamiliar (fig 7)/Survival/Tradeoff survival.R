## ----setup, include=FALSE---------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(lme4)
library(lmerTest)
library(survival)
library(coxphf)


## ---------------------------------------------------------------------------------------
pvalue <- function(x, ...) UseMethod("pvalue")

pvalue.survdiff <- function (x, ...) {
  if (length(x$n) == 1) {
    df <- 1
    pval <- as.numeric(pchisq(x$chisq, 1, lower.tail = FALSE))
  } else {
    if (is.matrix(x$obs)) {
      otmp <- rowSums(x$obs)
      etmp <- rowSums(x$exp)
    } else {
      otmp <- x$obs
      etmp <- x$exp
    }
    df <- sum(etmp > 0) - 1
    pval <- as.numeric(pchisq(x$chisq, df, lower.tail = FALSE))
  }
  # Ensure the list is explicitly returned
  result <- list(chisq = x$chisq, p.value = pval, df = df)
  return(result)  # Explicit return
}



## ---------------------------------------------------------------------------------------
survdata<- read.csv("Groupedvirulence.master.csv", sep=',', header=T)


## ---------------------------------------------------------------------------------------
survdata <- subset(survdata, Treatment!="STOCK" )
survdata$Host <- relevel(factor(survdata$Host), ref = "BALBF")


## ---------------------------------------------------------------------------------------
tradeoffs <- subset(survdata, Treatment != "1XP" & Sex =="F")

# Create a new column with Host and familiarity
tradeoffs$HostFam <- paste0(tradeoffs$Host,"_",tradeoffs$Familiarity)
tradeoffs$HostFam <- as.factor(tradeoffs$HostFam)

# Levels to iterate over, excluding BL6F_UNF
host_fams <- c("BALBF_FAM", "BALBF_UNF", "BL6F_FAM", "BL6F_UNF")

# Open the file for writing all summaries
sink("Survival-Tradeoffs_summary.txt")
#sink("Survival-Tradeoffs_NOBL6UNF_summary.txt")

# Iterate over each level, excluding BL6F_UNF from analysis
for (level in host_fams) {
  
  # Subset to exclude BL6F_UNF
  #filtered_data <- subset(tradeoffs, HostFam != "BL6F_UNF")
  
  # Relevel HostFam
  tradeoffs$HostFam <- relevel(factor(tradeoffs$HostFam), ref = level)
  
  # Run analyses
  cox <- coxph(Surv(tradeoffs$TTD, tradeoffs$Event) ~ factor(HostFam), data = tradeoffs)
  
  # Write the title indicating the level
  cat("\n\n##############################################\n")
  cat("            Leveled on", level, "\n")
  cat("##############################################\n\n")
  
  # Output summaries to the file
  print(summary(cox))
}

# Save and close the file
sink()



## ---------------------------------------------------------------------------------------
# Open the file for writing all summaries
sink("Survival-Tradeoffs_LOGRANK.txt")

  # Create the survival object
  surv_obj <- Surv(tradeoffs$TTD, tradeoffs$Event)
  
  # Perform the log-rank test
  logrank <- survdiff(surv_obj ~ HostFam, data = tradeoffs)
  
  # Write the results to the file
  
  print(logrank)

# Close the file
sink()



## ---------------------------------------------------------------------------------------

# Open the file for writing all summaries
sink("Survival-Tradeoffs_FAMvsUNF_LOGRANK.txt")

# Define the groups to compare
comparisons <- list(
  BALBF = c("BALBF_FAM", "BALBF_UNF"),
  BL6F = c("BL6F_FAM", "BL6F_UNF")
)

# Iterate through each comparison
for (comparison_name in names(comparisons)) {
  
  # Subset the data to include only the relevant groups
  relevant_groups <- comparisons[[comparison_name]]
  tradeoffs_subset <- subset(tradeoffs, HostFam %in% relevant_groups)
  
  # Relevel the HostFam factor so the first group is the reference
  tradeoffs_subset$HostFam <- relevel(factor(tradeoffs_subset$HostFam), ref = relevant_groups[1])
  
  # Create the survival object
  surv_obj <- Surv(tradeoffs_subset$TTD, tradeoffs_subset$Event)
  
  # Perform the log-rank test
  logrank <- survdiff(surv_obj ~ HostFam, data = tradeoffs_subset)
  
    # Calculate p-value using the custom function
  pval_info <- pvalue(logrank)
  
  # Write the results to the file
  cat("\n##############################################\n")
  cat("            ", relevant_groups[1], "vs", relevant_groups[2], "\n")
  cat("##############################################\n")
  print(logrank)
  cat("Chi-squared:", pval_info$chisq, "\n")
  cat("P-value:", signif(pval_info$p.value, 4), "\n") # Extracted p.value for formatting
  cat("Degrees of Freedom:", pval_info$df, "\n")
}

# Close the file
sink()



## ---------------------------------------------------------------------------------------
#Remove data that is not relevant
overallvir <- subset(survdata, Treatment !="1XP" & Sex =="F")

#Level of 1X and familiar

overallvir$Host <- relevel(factor(overallvir$Host), ref = "BALBF")
overallvir$Familiarity <- relevel(factor(overallvir$Familiarity), ref = "FAM")

#Run model
model.overallvir <- coxph(Surv(tradeoffs$TTD, tradeoffs$Event) ~ factor(Familiarity) + Host, data = overallvir)

summary(model.overallvir)

# Open the file for writing
sink("Survival-Overall_virulence_FAM_vs_UNFAM.txt")

# Output summaries to the file
summary(model.overallvir)
# Save
sink()


## ---------------------------------------------------------------------------------------
#Remove data that is not relevant
overallvir <- subset(survdata, Treatment !="1XP" & Sex =="F")

#Level of 1X and familiar

overallvir$Host <- relevel(factor(overallvir$Host), ref = "BALBF")
overallvir$Familiarity <- relevel(factor(overallvir$Familiarity), ref = "FAM")

#Run model
model.overallvir <- coxph(Surv(tradeoffs$TTD, tradeoffs$Event) ~ factor(Host), data = overallvir)

summary(model.overallvir)

# Open the file for writing
sink("Survival-Overall_virulence_BALBF_vs_BL6F.txt")

# Output summaries to the file
summary(model.overallvir)
# Save
sink()


## ---------------------------------------------------------------------------------------

overallmort <- subset(survdata, Familiarity =="UNF" & Sex =="F")


overallmort<- subset (overallmort, Treatment !="1XP")


#overallmort$Treatment <- relevel(factor(overallmort$Treatment), ref="1XP")


# Open the file for writing all summaries
sink("Survival-Unfamiliar_1X_2X_4X_plus_host.txt")

#Run analyses
  overallmort$TTD
  table(overallmort$Event,overallmort$TTD)
  Surv(overallmort$TTD,overallmort$Event)
  
  cox <- coxph(Surv(overallmort$TTD,overallmort$Event) ~ Treatment+Host, data=overallmort)

  # Output summaries to the file
  print(summary(cox))

# Save and close the file
sink()



## ---------------------------------------------------------------------------------------

overallmort <- subset(survdata, Familiarity =="UNF")


overallmort<- subset (overallmort, Treatment !="2X" & Treatment != "4X")


#overallmort$Treatment <- relevel(factor(overallmort$Treatment), ref="1XP")


# Open the file for writing all summaries
sink("Survival-Unfamiliar_1X_vs_1XP_plus_host.txt")

#Run analyses
  overallmort$TTD
  table(overallmort$Event,overallmort$TTD)
  Surv(overallmort$TTD,overallmort$Event)
  
  cox <- coxph(Surv(overallmort$TTD,overallmort$Event) ~ Treatment+Host, data=overallmort)

  # Output summaries to the file
  print(summary(cox))

# Save and close the file
sink()



## ---------------------------------------------------------------------------------------
# Filter for "UNF" familiarity
overallmort <- subset(survdata, Familiarity == "FAM")

# Create a new treatment category: Combine 1X, 2X, and 4X into "1X2X4X"
overallmort$Treatment <- ifelse(overallmort$Treatment %in% c("1X", "2X", "4X"), "1X2X4X", "1XP")

# Ensure the Treatment column is a factor and relevel with 1XP as the reference
overallmort$Treatment <- factor(overallmort$Treatment, levels = c("1XP", "1X2X4X"))

# Open the file for writing all summaries
sink("Survival-Familiar_1X2X4X_vs_1XP_plus_host.txt")

# Run analyses
print("Time to Death (TTD):")
print(overallmort$TTD)

print("Event Table (Event vs TTD):")
print(table(overallmort$Event, overallmort$TTD))

print("Survival Object:")
print(Surv(overallmort$TTD, overallmort$Event))

cox <- coxph(Surv(overallmort$TTD, overallmort$Event) ~ Treatment + Host, data = overallmort)

# Output summaries to the file
print(summary(cox))

# Save and close the file
sink()



## ---------------------------------------------------------------------------------------

hostmort <- subset(survdata, Treatment != "1XP" & Sex =="F")


hostmort$HostStrain <- relevel(factor(hostmort$HostStrain), ref = "BALB/C")

# Open the file for writing all summaries
sink("Survival-host-susceptibility_summary.txt")

#Run analyses
  hostmort$TTD
  table(hostmort$Event,hostmort$TTD)
  Surv(hostmort$TTD,hostmort$Event)
  
  cox <- coxph(Surv(hostmort$TTD,hostmort$Event) ~ factor(HostStrain), data=hostmort)

  # Output summaries to the file
  print(summary(cox))

# Save and close the file
sink()




## ---------------------------------------------------------------------------------------
#Remove data that is not relevant
poolingtest <- subset(survdata, Treatment !="STOCK" & Sex =="F")

#Level on 1X

poolingtest$Treatment <- relevel(factor(poolingtest$Treatment), ref = "1X")

#Run analyses
  poolingtest$TTD
  table(poolingtest$Event,poolingtest$TTD)
  Surv(poolingtest$TTD,poolingtest$Event)
  
  cox.pool <- coxph(Surv(poolingtest$TTD,poolingtest$Event) ~ factor(Treatment), data=poolingtest)

# Open the file for writing
sink("Survival-ALL_Treatments_in_allHosts.txt")

# Output summaries to the file
summary(cox.pool)
# Save
sink()


