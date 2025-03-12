## ----setup, include=FALSE---------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(lme4)
library(lmerTest)
library(survival)


## ---------------------------------------------------------------------------------------
pvalue <- function(x, ...) UseMethod("pvalue")
pvalue.survdiff <- function (x, ...) {
  if (length(x$n) == 1) {
    df <- 1
    pval <- pchisq(x$chisq, 1, lower.tail = FALSE)
  } else {
    if (is.matrix(x$obs)) {
      otmp <- rowSums(x$obs)
      etmp <- rowSums(x$exp)
    } else {
      otmp <- x$obs
      etmp <- x$exp
    }
    df <- sum(etmp > 0) - 1
    pval <- pchisq(x$chisq, df, lower.tail = FALSE)
  }
  list(chisq = x$chisq, p.value = pval, df = df)
  cat("Chi-squared:", x$chisq, "\nP-value:", pval, "\nDegrees of Freedom:", df, "\n")
}



## ---------------------------------------------------------------------------------------
survdata<- read.csv("Groupedvirulence.master.csv", sep=',', header=T)


## ---------------------------------------------------------------------------------------
#All the data, including unfamiliar
survdata2 <- subset(survdata, Treatment!="1XP" & Treatment !="STOCK") 
#Only familiar infections
survdatafam <- subset(survdata, Treatment!="1XP" & Treatment!="STOCK" & Familiarity=="FAM")

#Set base levels
survdatafam$Host <- relevel(factor(survdatafam$Host), ref = "BALBF")
survdatafam$Treatment <- relevel(factor(survdatafam$Treatment), ref = "1X")



## ---------------------------------------------------------------------------------------
# Subset your data
survdata3 <- subset(survdata, Treatment!="1XP" & Familiarity=="FAM")
# Define the host levels and the output file name
host_levels <- c("BALBF", "BL6F", "BALBM")
output_file <- "logrank_results_combined.txt"

# Open the output file for writing
sink(output_file)

# Iterate over the host levels
for (host in host_levels) {
  
  # Create a survival object
  surv_obj <- Surv(survdata3$TTD, survdata3$Event)
  
  # Run the log-rank test
  logrank_result <- survdiff(surv_obj ~ factor(Status), data = survdata3)
  
  # Print the result to the file
  cat("Log-rank test results for Host = ", host, ":\n", sep = "")
  print(logrank_result)
  cat("\n")
  
  # Print detailed p-value
  cat(pvalue(logrank_result), "\n\n", sep = "")
}

# Close the file
sink()



## ---------------------------------------------------------------------------------------
# Define the host levels and the output file name
host_levels <- c("BALBF", "BL6F", "BALBM")
output_file <- "CoxPH_Host_results_combined.txt"

# Open the output file for writing
sink(output_file)

# Iterate over the host levels
for (host in host_levels) {
  
  # Relevel Host variable to set the current host as the reference level
  survdatafam$Host <- relevel(factor(survdatafam$Host), ref = host)
  
  # Fit the CoxPH model for the current host level
  cox_model <- coxph(Surv(TTD, Event) ~ Host, data = survdatafam)
  
  # Print the result to the file
  cat("\n\nCox Proportional Hazards Model for Host = ", host, " (reference = ", host, "):\n", sep = "")
  print(cox_model)
}

# Close the file
sink()


## ---------------------------------------------------------------------------------------
# Fit a CoxPH model for Treatment only
survdatafam$Host <- relevel(factor(survdatafam$Host), ref = "BALBF")
survdatafam$Treatment <- relevel(factor(survdatafam$Treatment), ref = "1X")
cox_model <- coxph(Surv(TTD, Event) ~ Treatment + Host, data = survdatafam)

# Print the summary of the model to see the results
summary(cox_model)

# Save to txt (no csv possible)
sink("CoxPH-Treatment.txt")
print(cox_model)
sink()



## ---------------------------------------------------------------------------------------
# Fit the CoxPH model for the current host level
survdatatradeoff <- subset(survdata, Sex=="F" & Treatment != "STOCK" & Treatment !="1XP")

survdatatradeoff$Host <- relevel(factor(survdatatradeoff$Host), ref = "BL6F")

cox_model <- coxph(Surv(TTD, Event) ~ Host, data = survdatatradeoff)
  
sink("CoxPH_Familiar+Unfamiliar_results_combined.txt")
  
# Print the result to the file
cat("Cox Proportional Hazards Model")
print(cox_model)


# Close the file
sink()


## ---------------------------------------------------------------------------------------
# Fit the CoxPH model for the current host level
survdatatradeoff <- subset(survdata, Sex=="F" & Treatment != "STOCK" & Treatment !="1XP")

survdatatradeoff$Host <- relevel(factor(survdatatradeoff$HostStrain), ref = "C57BL/6")

cox_model <- coxph(Surv(TTD, Event) ~ HostStrain, data = survdatatradeoff)
  
sink("CoxPH_General_Host_Susceptibility_results.txt")
  
# Print the result to the file
cat("Cox Proportional Hazards Model")
print(cox_model)


# Close the file
sink()


## ---------------------------------------------------------------------------------------
# Remove data that is not relevant
tradeoffs <- subset(survdata, Treatment != "STOCK" & Treatment != "2X" & Treatment != "1XP" & Sex == "F")

# Create a new column with Status and Host -> Type
tradeoffs$HostFam <- paste0(tradeoffs$Host,"_",tradeoffs$Familiarity)

# Levels to iterate over
host_fams <- c("BALBF_FAM", "BALBF_UNF", "BL6F_FAM", "BL6F_UNF")

# Open the file for writing all summaries
sink("CoxPH-Tradeoffs_1X_4X_only_summary.txt")

# Iterate over each level
for (level in host_fams) {
  
  # Relevel HostFam
  
  tradeoffs$HostFam <- relevel(factor(tradeoffs$HostFam), ref = level)
  
  # Run model

  model.tradeoffs <- coxph(Surv(TTD, Event) ~ HostFam, data = tradeoffs)
  
  # Write the title indicating the level
  cat("\n\n##############################################\n")
  cat("            Leveled on", level, "\n")
  cat("##############################################\n\n")
  
  # Output summaries to the file
  print(summary(model.tradeoffs))
}

# Save and close the file
sink()



## ---------------------------------------------------------------------------------------
survdatafem <- subset(survdata, Treatment=="1X" | Treatment=="1XP") #Only 1X and 1XP
survdatafem <- subset(survdatafem, Sex=="F")
survdatafemfam <- subset(survdatafem, Familiarity == "FAM")


## ---------------------------------------------------------------------------------------
# Convert Test to a factor
survdatafem$Test <- factor(survdatafem$Test)

#Relevel 
survdatafem$Test <- relevel(survdatafem$Test, ref="BL6F.FAM")

#Run analyses
survdatafem$TTD
table(survdatafem$Event,survdatafem$TTD)
Surv(survdatafem$TTD,survdatafem$Event)



cox <- coxph(Surv(survdatafem$TTD,survdatafem$Event) ~ factor(Treatment), data=survdatafem)
summary(cox)

#Save
sink("coxph-1XP-female-only.txt")
print(cox)
sink()


## ---------------------------------------------------------------------------------------
# Convert Test to a factor
survdatafem$Test <- factor(survdatafem$Test)

#Relevel 
survdatafem$Test <- relevel(survdatafem$Test, ref="BL6F.FAM")

#Run analyses
survdatafem$TTD
table(survdatafem$Event,survdatafem$TTD)
Surv(survdatafem$TTD,survdatafem$Event)



cox <- coxph(Surv(survdatafem$TTD,survdatafem$Event) ~ factor(Treatment) * factor(Familiarity), data=survdatafem)
summary(cox)

sink("coxph-1XP-treatment-familiarity.txt")
print(cox)
sink()


## ---------------------------------------------------------------------------------------
# Convert Test to a factor
survdatafem$Test <- factor(survdatafem$Test)
survdatafem$Treatment <- factor(survdatafem$Treatment)
#Relevel 
survdatafem$Test <- relevel(survdatafem$Test, ref="BL6F.FAM")
survdatafem$Treatment <- relevel(survdatafem$Treatment, ref="1X")

#Run analyses
survdatafem$TTD
table(survdatafem$Event,survdatafem$TTD)
Surv(survdatafem$TTD,survdatafem$Event)



cox <- coxph(Surv(survdatafem$TTD,survdatafem$Event) ~ factor(Test) * factor(Treatment), data=survdatafem)
summary(cox)

sink("coxph-1XP-test.txt")
print(cox)
sink()

