## ----setup, include=FALSE---------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(lme4)
library(lmerTest)
library(survival)
library(survminer)


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
survdatafam <- subset(survdata, Treatment!="STOCK" & Familiarity=="FAM")

#Set base levels
survdatafam$Host <- relevel(factor(survdatafam$Host), ref = "BALBF")
survdatafam$Treatment <- relevel(factor(survdatafam$Treatment), ref = "1X")



## ----FOR 1X 2X 4X (Optionally can include 1XP) BY HOST----------------------------------

survdatafam2 <- subset(survdata, Treatment!="STOCK" & Treatment!="1XP" & Familiarity == "FAM")

# Define the host levels and the output file name
host_levels <- c("BALBF", "BL6F", "BALBM")
output_file <- "Survival_Treatment_Results_By_Host_1X_2X_4X.txt"

# Open the output file for writing
sink(output_file)

# Iterate over the host levels
for (host in host_levels) {
  
  cat("\n\n##############################################\n")
  cat("Analysis for Host = ", host, "\n")
  cat("##############################################\n")
  
  #subset host
  host_subset <- subset(survdatafam2, Host == host)
  
  # Fit the CoxPH model
  cox_model <- tryCatch({
    coxph(Surv(TTD, Event) ~ Treatment, data = host_subset)
  }, warning = function(w) {
    cat("\nWARNING: CoxPH model for Host =", host, "issued the following warning:\n", w$message, "\n")
    cat("REFER TO LOGRANK\n")
    NULL
  }, error = function(e) {
    cat("\nERROR: CoxPH model for Host =", host, "failed with the following error:\n", e$message, "\n")
    NULL
  })
  
  if (!is.null(cox_model)) {
    # Print CoxPH model results if it converged
    cat("\nCox Proportional Hazards Model Results:\n")
    print(summary(cox_model))
  }
  
  # Perform the log-rank test (survdiff)
  surv_obj <- Surv(host_subset$TTD, host_subset$Event)
  logrank <- tryCatch({
    survdiff(surv_obj ~ Treatment, data = host_subset)
  }, error = function(e) {
    cat("\nERROR: Log-rank test for Host =", host, "failed with the following error:\n", e$message, "\n")
    NULL
  })
  
  if (!is.null(logrank)) {
    # Print log-rank test results
    cat("\nLog-Rank Test Results:\n")
    print(logrank)
    
    # Perform pairwise comparisons using pairwise_survdiff
    pairwise_results <- tryCatch({
      pairwise_survdiff(Surv(TTD, Event) ~ Treatment, data = host_subset)
    }, error = function(e) {
      cat("\nERROR: Pairwise comparisons for Host =", host, "failed with the following error:\n", e$message, "\n")
      NULL
    })
    
    if (!is.null(pairwise_results)) {
      cat("\nPairwise Log-Rank Test Results:\n")
      print(pairwise_results)
    }
  }
}

# Close the file
sink()



## ----FOR 1X vs 1XP BY HOST--------------------------------------------------------------

survdatafam3 <- subset(survdata, Treatment!="STOCK" & Treatment!="2X" & Treatment!="4X" & Familiarity == "FAM")

# Define the host levels and the output file name
host_levels <- c("BALBF", "BL6F", "BALBM")
output_file <- "Survival_Treatment_Results_By_Host_1X_vs_1XP.txt"

# Open the output file for writing
sink(output_file)

# Iterate over the host levels
for (host in host_levels) {
  
  cat("\n\n##############################################\n")
  cat("Analysis for Host = ", host, "\n")
  cat("##############################################\n")
  
  #subset host
  host_subset <- subset(survdatafam3, Host == host)
  
  # Fit the CoxPH model
  cox_model <- tryCatch({
    coxph(Surv(TTD, Event) ~ Treatment, data = host_subset)
  }, warning = function(w) {
    cat("\nWARNING: CoxPH model for Host =", host, "issued the following warning:\n", w$message, "\n")
    cat("REFER TO LOGRANK\n")
    NULL
  }, error = function(e) {
    cat("\nERROR: CoxPH model for Host =", host, "failed with the following error:\n", e$message, "\n")
    NULL
  })
  
  if (!is.null(cox_model)) {
    # Print CoxPH model results if it converged
    cat("\nCox Proportional Hazards Model Results:\n")
    print(summary(cox_model))
  }
  
  # Perform the log-rank test (survdiff)
  surv_obj <- Surv(host_subset$TTD, host_subset$Event)
  logrank <- tryCatch({
    survdiff(surv_obj ~ Treatment, data = host_subset)
  }, error = function(e) {
    cat("\nERROR: Log-rank test for Host =", host, "failed with the following error:\n", e$message, "\n")
    NULL
  })
  
  if (!is.null(logrank)) {
    # Print log-rank test results
    cat("\nLog-Rank Test Results:\n")
    print(logrank)
    
    # Perform pairwise comparisons using pairwise_survdiff
    pairwise_results <- tryCatch({
      pairwise_survdiff(Surv(TTD, Event) ~ Treatment, data = host_subset)
    }, error = function(e) {
      cat("\nERROR: Pairwise comparisons for Host =", host, "failed with the following error:\n", e$message, "\n")
      NULL
    })
    
    if (!is.null(pairwise_results)) {
      cat("\nPairwise Log-Rank Test Results:\n")
      print(pairwise_results)
    }
  }
}

# Close the file
sink()


