library(tidyverse)
library(stats)
library(dunn.test)
library(Hmisc)


#Load data

data <- read.csv("Input Data/ancdercounts-unfiltered.csv")

#Subset for round 10 samples with homogenate balb female data (not purified).
#Also remove lower coverage 4X biological replicates
data<- data %>% filter(Round== "10" & SampleID != 10318 & SampleID != 9770)

# Extract segment names dynamically from column names
segment_column_list <- names(data)[10:length(names(data))] %>%
  gsub("^([^.]+)\\..*$", "\\1", .) %>%
  unique()

####Set parameters for analysis####

#Set parameters
frequency_cutoff <- 0.03 #Minimum frequency of indels
length_cutoff <- 3 #Minimum length of indels
coverage_cutoff <- 100 #Average coverage at PA (eliminates samples with poor coverage)

####Wrangle data####
detect_stretches <- function(df) {
  df <- df %>%
    group_by(SampleID, HostStrain, Sex, Host, Treatment, Line, Round, Virulence, Medium, Segment) %>%
    mutate(Frequency = Der / (Anc + Der)) %>%
    filter(Frequency > 0.005) %>%
    group_by(SampleID, HostStrain, Sex, Host, Treatment, Line, Round, Virulence, Medium, Segment, add = TRUE) %>%
    mutate(Group = cumsum(c(TRUE, diff(as.numeric(Position)) != 1))) %>%
    group_by(SampleID, HostStrain, Sex, Host, Treatment, Line, Round, Virulence, Medium, Segment, Group, add = TRUE) %>%
    filter(n() > 1) %>%
    summarise(
      Start = if (length(Position) > 0) min(as.numeric(Position)) else NA_integer_,
      End = if (length(Position) > 0) max(as.numeric(Position)) else NA_integer_,
      Length = n(),
      Indel_Frequency = mean(Frequency),
      .groups = 'drop'
    )
  
  return(df)
}

results <- list()
for (Segment in segment_column_list) {
  sample_info_columns <- data %>%
    select(SampleID, HostStrain, Sex, Host, Treatment, Line, Round, Virulence, Medium)  # Ensure these are the columns you need
  
  anc_columns <- data %>%
    select(matches(paste0("^", Segment, "\\.\\d+\\.Anc$")))
  der_columns <- data %>%
    select(matches(paste0("^", Segment, "\\.\\d+\\.Der$")))
  
  if (ncol(anc_columns) > 0 && ncol(der_columns) > 0 && nrow(anc_columns) == nrow(der_columns)) {
    names(anc_columns) <- gsub(paste0(Segment, "\\.(\\d+)\\.(.*)$"), "\\1_Anc", names(anc_columns))
    names(der_columns) <- gsub(paste0(Segment, "\\.(\\d+)\\.(.*)$"), "\\1_Der", names(der_columns))
    
    combined_data <- cbind(sample_info_columns, Segment = Segment, data.frame(anc_columns, der_columns, check.names = FALSE))
    
    long_data <- pivot_longer(
      combined_data,
      cols = matches("_Anc$|_Der$"),
      names_to = c("Position", ".value"),
      names_sep = "_"
    )
    
    results[[Segment]] <- detect_stretches(long_data)
  } else {
    cat(sprintf("Skipping Segment %s due to insufficient data\n", Segment))
  }
}


temp_results <- bind_rows(results)

head(temp_results)

####Determine sample mean coverage#####

# Calculate average coverage for PA segment
data$Total_Coverage_PA <- rowSums(data[grep("^PA\\.\\d+\\.Anc$", names(data))] + data[grep("^PA\\.\\d+\\.Der$", names(data))], na.rm = TRUE) / length(grep("^PA\\.\\d+\\.Anc$", names(data)))

# Create a dataframe with mean coverage for each SampleID
PA_coverage <- data %>%
  group_by(SampleID) %>%
  summarise(PAcov = mean(Total_Coverage_PA, na.rm = TRUE)) %>%
  ungroup()

# Filter out SampleIDs with PA coverage less than 100
valid_sample_ids <- PA_coverage %>%
  filter(PAcov >= coverage_cutoff) %>%
  pull(SampleID)

# Filter to include and exclude certain rows based on conditions
filtered_results_above_pct <- temp_results %>%
  filter(Length >= length_cutoff, Indel_Frequency >= frequency_cutoff, SampleID %in% valid_sample_ids)

excluded_results <- temp_results %>%
  filter(Length < length_cutoff | Indel_Frequency < frequency_cutoff | !(SampleID %in% valid_sample_ids))

# Calculate grouped by SampleID
sample_calculations <- filtered_results_above_pct %>%
  group_by(SampleID) %>%
  summarise(
    Indel_Count_Sample = n(),
    Average_Length_Sample = mean(Length),
    Average_Frequency_Sample = mean(Indel_Frequency)
  ) %>%
  ungroup()

#Save grouped by sampleID to get per sample Stats.
# Convert frequency cutoff to a percentage string for the filename
frequency_pct <- paste0(format(frequency_cutoff * 100, nsmall = 0), "%")
# Write the final combined results to CSV with dynamic filename based on frequency cutoff
filename <- paste("Per_sample_stats-", frequency_pct, ".csv", sep = "")
write.csv(sample_calculations, filename, row.names = FALSE)

# Join sample calculations back to the filtered data
filtered_with_samples <- filtered_results_above_pct %>%
  left_join(sample_calculations, by = "SampleID")


# Calculate grouped by Host
host_calculations <- filtered_with_samples %>%
  group_by(Host) %>%
  summarise(
    Indel_Total_Host = n(),
    Average_Count_Host = mean(Indel_Count_Sample, na.rm = TRUE), # Using mean from sample_calculations
    Average_Length_Host = mean(Length),
    Average_Frequency_Host = mean(Indel_Frequency)
  ) %>%
  ungroup()

# Join host calculations back to the main filtered data
final_results <- filtered_with_samples %>%
  left_join(host_calculations, by = "Host")

# Convert numeric columns to character in the filtered results for compatibility in binding
final_results <- filtered_results_above_pct %>%
  left_join(sample_calculations, by = "SampleID") %>%
  left_join(host_calculations, by = "Host") %>%
  mutate(across(c(Indel_Count_Sample, Average_Length_Sample, Average_Frequency_Sample,
                  Indel_Total_Host, Average_Count_Host, Average_Length_Host, Average_Frequency_Host), as.character))

# Combine the above 1pct results with the excluded entries
final_combined_results <- bind_rows(
  final_results,
  excluded_results %>%
    mutate(
      Indel_Count_Sample = "Excluded",
      Average_Length_Sample = "Excluded",
      Average_Frequency_Sample = "Excluded",
      Indel_Total_Host = "Excluded",
      Average_Count_Host = "Excluded",
      Average_Length_Host = "Excluded",
      Average_Frequency_Host = "Excluded"
    )
)

# Add a column to indicate whether the sample was included
final_combined_results <- final_combined_results %>%
  mutate(Sample_included = if_else(SampleID %in% valid_sample_ids, "Yes", "No"))

#save
# Convert frequency cutoff to a percentage string for the filename
frequency_pct <- paste0(format(frequency_cutoff * 100, nsmall = 0), "%")

# Write the final combined results to CSV with dynamic filename based on frequency cutoff
filename <- paste("Nucleotide-indels-R10-homogenates-", frequency_pct, ".csv", sep = "")
write.csv(final_combined_results, filename, row.names = FALSE)

#Convert everything important to numeric before doing stats

# Ensure columns are numeric (convert non-numeric entries that may cause issues)
final_combined_results2 <- final_combined_results %>%
  mutate(
    Average_Count_Host = as.numeric(Average_Count_Host),
    Average_Length_Host = as.numeric(Average_Length_Host),
    Average_Frequency_Host = as.numeric(Average_Frequency_Host)
  )

# Select unique host summaries and round values
unique_host_summaries <- final_combined_results2 %>%
  select(Host, Average_Count_Host, Average_Length_Host, Average_Frequency_Host) %>%
  distinct(Host, .keep_all = TRUE) %>%
  mutate(
    Average_Count_Host = round(Average_Count_Host, 1),
    Average_Length_Host = round(Average_Length_Host, 1),
    Average_Frequency_Host = round(Average_Frequency_Host, 3)
  )

# Continue with your sink() and output formatting
output_file <- paste("Results/host_indel_summaries_", frequency_pct, ".txt", sep = "")

sink(output_file)
cat("Cutoff Parameters:\n")
cat("Indel Frequency Cutoff = ", format(frequency_cutoff * 100, digits = 2), "%\n")
cat("Indel Length Cutoff = ", length_cutoff, "\n")
cat("Polymerase Mean Coverage Cutoff = ", coverage_cutoff, "\n\n")

cat("Average number of indels:\n")
for (i in 1:nrow(unique_host_summaries)) {
  cat(unique_host_summaries$Host[i], "-", unique_host_summaries$Average_Count_Host[i], "\n")
}

cat("\nAverage length of indels:\n")
for (i in 1:nrow(unique_host_summaries)) {
  cat(unique_host_summaries$Host[i], "-", unique_host_summaries$Average_Length_Host[i], "\n")
}

cat("\nAverage frequency of indels:\n")
for (i in 1:nrow(unique_host_summaries)) {
  cat(unique_host_summaries$Host[i], "-", unique_host_summaries$Average_Frequency_Host[i], "\n")
}

sink()  # Close the sink function to stop redirection

####stats###
#Convert to numeric, if not done yet.
final_results$Indel_Count_Sample <- as.numeric(final_results$Indel_Count_Sample)
final_results$Average_Length_Sample <- as.numeric(final_results$Average_Length_Sample)
final_results$Average_Frequency_Sample <- as.numeric(final_results$Average_Frequency_Sample)

#Load data with set cutoff.
data<- final_results %>% filter(Indel_Frequency >= frequency_cutoff)

#Order data
desired_order <- c("BALBF", "BALBM", "BL6F", "BL6M")

# Set the factor levels for 'Host' according to the desired order
data <- data %>%
  mutate(Host = factor(Host, levels = desired_order))

#Calculate means for each variable, then prepare to print
# Calculate means and counts by Host
means_by_host <- data %>%
  group_by(Host) %>%
  summarise(
    `Mean number of deletions` = n() / n_distinct(SampleID),
    Length = mean(Length, na.rm = TRUE),
    Deletion_Frequency = mean(Indel_Frequency, na.rm = TRUE)
  )

# Calculate means and counts by HostStrain
means_by_hoststrain <- data %>%
  group_by(HostStrain) %>%
  summarise(
    `Mean number of deletions` = n() / n_distinct(SampleID),
    Length = mean(Length, na.rm = TRUE),
    Deletion_Frequency = mean(Indel_Frequency, na.rm = TRUE)
  )

# Calculate means and counts by Sex
means_by_sex <- data %>%
  group_by(Sex) %>%
  summarise(
    `Mean number of deletions` = n() / n_distinct(SampleID),
    Length = mean(Length, na.rm = TRUE),
    Deletion_Frequency = mean(Indel_Frequency, na.rm = TRUE)
  )
# Function to print means to sink
print_means <- function(means, group) {
  group_column <- sym(group)
  
  cat(paste("Mean number of deletions above",frequency_pct, "frequency by", group, ":\n"))
  means %>% select(!!group_column, `Mean number of deletions`) %>%
    mutate(`Mean number of deletions` = round(`Mean number of deletions`, 2)) %>%
    rowwise() %>%
    mutate(output = paste(!!group_column, "-", `Mean number of deletions`)) %>%
    pull(output) %>%
    cat(sep = "\n")
  cat("\n")
  
  cat(paste("Mean length of deletions above",frequency_pct, "frequency by", group, ":\n"))
  means %>% select(!!group_column, Length) %>%
    mutate(Length = round(Length, 2)) %>%
    rowwise() %>%
    mutate(output = paste(!!group_column, "-", Length)) %>%
    pull(output) %>%
    cat(sep = "\n")
  cat("\n")
  
  cat(paste("Mean frequency of deletions above",frequency_pct, "frequency by", group, ":\n"))
  means %>% select(!!group_column, Deletion_Frequency) %>%
    mutate(Deletion_Frequency = round(Deletion_Frequency, 3)) %>%
    rowwise() %>%
    mutate(output = paste(!!group_column, "-", Deletion_Frequency)) %>%
    pull(output) %>%
    cat(sep = "\n")
  cat("\n")
}

# Create a function to perform Kruskal-Wallis and Dunn's tests
perform_kruskal_dunn <- function(data, variable, group) {
  kruskal_result <- kruskal.test(as.formula(paste(variable, "~", group)), data = data)
  print(kruskal_result)
  
  if(kruskal_result$p.value < 0.05) {
    dunn_result <- dunn.test(data[[variable]], data[[group]], method="bonferroni")
    print(dunn_result)
  } else {
    cat(paste("Kruskal-Wallis test p-value for", variable, "by", group, "is not less than 0.05. No significant difference detected among groups.\n"))
  }
}

# Create a function to perform Mann-Whitney U test
perform_mann_whitney <- function(data, variable, group) {
  mann_whitney_result <- wilcox.test(as.formula(paste(variable, "~", group)), data = data)
  print(mann_whitney_result)
  
  if(mann_whitney_result$p.value < 0.05) {
    cat(paste("Mann-Whitney test p-value for", variable, "by", group, "is less than 0.05. Significant difference detected among groups.\n"))
  } else {
    cat(paste("Mann-Whitney test p-value for", variable, "by", group, "is not less than 0.05. No significant difference detected among groups.\n"))
  }
}

# Print the results to console and sink to a file
KWDunnsresults <- paste("Results/kruskal_dunn_results_", frequency_pct, "_", coverage_cutoff, "X", ".txt", sep = "")
sink(KWDunnsresults)

cat("\n-------------------------------------------------------")
cat("\nParameters for statistical analysis of genome deletions\n")
cat("-------------------------------------------------------\n\n")

cat("Indel Frequency Cutoff =", format(frequency_cutoff * 100, digits = 2), "%\n")
cat("Indel Length Cutoff =", length_cutoff, "nucleotides\n")
cat("Polymerase Mean Coverage Cutoff =", coverage_cutoff, "X\n")

cat("\n-------------------------------------------------------")
cat("\n Comparisons for Host (BALBF vs BALBM vs BL6F vs BL6M)\n")
cat("-------------------------------------------------------\n")
# Length and Indel Frequency by Host
cat("\n--- Number of Indels by Host ---\n")
perform_kruskal_dunn(data, "Indel_Count_Sample", "Host")
cat("\n--- Length by Host ---\n")
perform_kruskal_dunn(data, "Length", "Host")
cat("\n--- Frequency by Host ---\n")
perform_kruskal_dunn(data, "Indel_Frequency", "Host")
cat("\n--- Means for each Host ---\n\n")
print_means(means_by_host, "Host")

cat("\n-------------------------------------------------------")
cat("\n    Comparisons for HostStrain (BALB/c vs C57BL/6)\n")
cat("-------------------------------------------------------\n")

# Length and Indel Frequency by HostStrain
cat("\n--- Number of Indels by HostStrain ---\n")
perform_mann_whitney(data, "Indel_Count_Sample","HostStrain")
cat("\n--- Length by HostStrain ---\n")
perform_mann_whitney(data, "Length", "HostStrain")
cat("\n--- Frequency by HostStrain ---\n")
perform_mann_whitney(data, "Indel_Frequency", "HostStrain")
cat("\n--- Means for each HostStrain ---\n\n")
print_means(means_by_hoststrain, "HostStrain")

cat("\n-------------------------------------------------------")
cat("\n           Comparisons for Sex (M vs F)\n")
cat("-------------------------------------------------------\n")
# Length and Indel Frequency by Sex
cat("\n--- Number of Indels by Sex ---\n")
perform_mann_whitney(data, "Indel_Count_Sample","Sex")
cat("\n--- Length by Sex ---\n")
perform_mann_whitney(data, "Length", "Sex")
cat("\n--- Frequency by Sex ---\n")
perform_mann_whitney(data, "Indel_Frequency", "Sex")
cat("\n--- Means for each Sex ---\n\n")
print_means(means_by_sex, "Sex")

sink()

cat("Stats complete")


####Correlations between deletions measurements.

# Calculate correlation and p-values for each group
calculate_correlation <- function(data, group_var) {
  data %>%
    group_by_at(group_var) %>%
    summarise(correlation = cor(Length, Indel_Frequency, method = "spearman", use = "complete.obs"),
              p_value = cor.test(Length, Indel_Frequency, method = "spearman")$p.value)
}

# Correlation for each group
correlation_by_host <- calculate_correlation(data, "Host")
correlation_by_hoststrain <- calculate_correlation(data, "HostStrain")
correlation_by_sex <- calculate_correlation(data, "Sex")

# Print correlations with p-values
cat("Correlation between Length and Frequency by Host:\n")
print(correlation_by_host)
cat("\nCorrelation between Length and Frequency by HostStrain:\n")
print(correlation_by_hoststrain)
cat("\nCorrelation between Length and Frequency by Sex:\n")
print(correlation_by_sex)

# Scatter plots
plot_relationship <- function(data, group_var) {
  ggplot(data, aes(x = Length, y = Indel_Frequency, color = !!sym(group_var))) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    labs(title = paste("Relationship between Length and Frequency by", group_var),
         x = "Deletion Length",
         y = "Deletion Frequency") +
    theme_minimal()
}

# Plot for each group
plot_relationship(data, "Host")
plot_relationship(data, "HostStrain")
plot_relationship(data, "Sex")


####Make stats summaries####

# Deletions above 3%
data <- final_results %>% filter(Indel_Frequency >= 0.03)

# Order data
desired_order <- c("BALBF", "BALBM", "BL6F", "BL6M")

# Set the factor levels for 'Host' according to the desired order
data <- data %>%
  mutate(Host = factor(Host, levels = desired_order))

# Calculate means and counts by Host
means_by_host <- data %>%
  group_by(Host) %>%
  summarise(
    `Mean number of deletions` = n() / n_distinct(SampleID),
    Length = mean(Length, na.rm = TRUE),
    Deletion_Frequency = mean(Indel_Frequency, na.rm = TRUE)
  )

# Create a function to perform Kruskal-Wallis and capture results
perform_kruskal_dunn_and_capture <- function(data, variable, group) {
  kruskal_result <- kruskal.test(as.formula(paste(variable, "~", group)), data = data)
  kw_statistic <- kruskal_result$statistic
  p_value <- kruskal_result$p.value
  
  # Capture means by group
  means <- data %>%
    group_by(!!sym(group)) %>%
    summarise(
      Mean = mean(!!sym(variable), na.rm = TRUE)
    ) %>%
    mutate(Analysis = variable, KW_statistic = kw_statistic, p_value = p_value)
  
  return(means)
}

# Perform Kruskal-Wallis and capture results for each variable
mean_number_of_deletions <- perform_kruskal_dunn_and_capture(data, "Indel_Count_Sample", "Host")
mean_length_of_deletions <- perform_kruskal_dunn_and_capture(data, "Length", "Host")
mean_frequency_of_deletions <- perform_kruskal_dunn_and_capture(data, "Indel_Frequency", "Host")

# Combine results into a single table
summary_table <- bind_rows(
  mean_number_of_deletions %>% mutate(Measure = "Mean number of deletions"),
  mean_length_of_deletions %>% mutate(Measure = "Mean length of deletions"),
  mean_frequency_of_deletions %>% mutate(Measure = "Mean frequency of deletions")
) %>%
  select(Measure, Host, Mean, KW_statistic, p_value)

# Save to CSV
write_csv(summary_table, "Results/summary_table_above_3_percent.csv",row.names = FALSE)

# Print table to console
print(summary_table)



# Save to CSV
write_csv(data, "Deletion_Table_3_percent_100cov.csv")


####Does pooling treatment affect deletions?####

# Function to perform Kruskal-Wallis and Dunn's tests
perform_kruskal_dunn_tests <- function(data, variable, hoststrain) {
  cat("Analysis for HostStrain:", hoststrain, "and variable:", variable, "\n")
  
  # Kruskal-Wallis test for Treatment
  cat("Kruskal-Wallis Test for Treatment:\n")
  kw_treatment <- kruskal.test(as.formula(paste(variable, "~ Treatment")), data = data)
  print(kw_treatment)
  
  # Dunn's test for Treatment if Kruskal-Wallis is significant
  if (kw_treatment$p.value < 0.05) {
    cat("Since p-value is less than 0.05, proceeding with Dunn's test for Treatment pairwise comparisons.\n")
    dunn_result_treatment <- dunn.test(data[[variable]], data$Treatment, method = "bonferroni")
    print_dunn_test_results(dunn_result_treatment)
  } else {
    cat("Kruskal-Wallis test p-value for Treatment is not less than 0.05. No significant difference detected among groups.\n")
  }
  
  cat("\n") # Add extra line for readability
}

# Function to print Dunn's test results
print_dunn_test_results <- function(dunn_result) {
  cat("Dunn's test results:\n")
  print(dunn_result$Z)  # Z scores
  print(dunn_result$P.adjusted)  # Adjusted p-values
}

# Split the data by HostStrain
data_BALBC <- data %>% filter(HostStrain == "BALB/C")
data_C57BL6 <- data %>% filter(HostStrain == "C57BL/6")

# Prepare the output file for capturing the analysis results
sink("kruskal_dunn__treatment_by_hoststrain_results.txt")

# Variables to analyze
variables <- c("Indel_Count_Sample", "Length", "Indel_Frequency")

# Perform the analysis for each variable and HostStrain
for (variable in variables) {
  perform_kruskal_dunn_tests(data_BALBC, variable, "BALB/C")
  perform_kruskal_dunn_tests(data_C57BL6, variable, "C57BL/6")
}

# Redirect output back to the console
sink()
