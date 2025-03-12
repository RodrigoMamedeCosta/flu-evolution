# Load necessary libraries
library(dplyr)
library(dunn.test)
library(Hmisc)



titers <- read.csv("GroupedtiterR10.csv", sep=',', header=T)

r10titers <- subset(titers, Round != "STOCK")

#OR INCLUDE STOCK:
# r10titers <- titers


### TITER STATS###

# Prepare the output file for capturing the analysis results
sink("TCID50_kruskal_dunn_results.txt")

# Kruskal-Wallis test
cat("Kruskal-Wallis Test:\n")
kruskal_result <- kruskal.test(TCID50 ~ Host, data = r10titers)
print(kruskal_result)

# Check if Kruskal-Wallis test is significant before proceeding to Dunn's test
if (kruskal_result$p.value < 0.05) {
  cat("Since p-value is less than 0.05, proceeding with Dunn's test for pairwise comparisons.\n")
  
  # Dunn's test for pairwise comparisons
  dunn_result <- dunn.test(r10titers$TCID50, r10titers$Host, method="bonferroni")
  print(dunn_result)
} else {
  cat("Kruskal-Wallis test p-value is not less than 0.05. No significant difference detected among groups.\n")
}

cat("\n") # Add extra line for readability

# Redirect output back to the console
sink()



### COPY NUMBER STATS###

# Prepare the output file for capturing the analysis results
sink("COPYN_kruskal_dunn_results.txt")

# Kruskal-Wallis test
cat("Kruskal-Wallis Test:\n")
kruskal_result <- kruskal.test(Copy.Number ~ Host, data = r10titers)
print(kruskal_result)

# Check if Kruskal-Wallis test is significant before proceeding to Dunn's test
if (kruskal_result$p.value < 0.05) {
  cat("Since p-value is less than 0.05, proceeding with Dunn's test for pairwise comparisons.\n")
  
  # Dunn's test for pairwise comparisons
  dunn_result <- dunn.test(r10titers$Copy.Number, r10titers$Host, method="bonferroni")
  print(dunn_result)
} else {
  cat("Kruskal-Wallis test p-value is not less than 0.05. No significant difference detected among groups.\n")
}

cat("\n") # Add extra line for readability

# Redirect output back to the console
sink()
