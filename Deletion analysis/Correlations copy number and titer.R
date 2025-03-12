# Load necessary libraries
library(dplyr)
library(ggplot2)
library(scales)
library(MASS)


data <- read.csv("Input Data/GroupedtiterR10.csv", sep=',', header=T)

data <- data %>%
  filter(!is.na(Indel.count.3pct))

# Ensure TCID50 and Copy Number are numeric
data <- data %>%
  mutate(TCID50 = as.numeric(TCID50),
         Copy.Number = as.numeric(`Copy.Number`),
         Host = as.factor(Host))

data <- data %>%
  mutate(log.TCID50 = log10(TCID50)) %>%
  mutate(log.Copy.Number = log10(Copy.Number))

# Calculate Spearman's correlation for the entire dataset (no grouping by Host)
spearman_results_combined <- data %>%
  summarise(
    spearman_rho_TCID50_Indel = cor(log.TCID50, Indel.count.3pct, method = "spearman"),
    p_value_TCID50_Indel = suppressWarnings(cor.test(log.TCID50, Indel.count.3pct, method = "spearman")$p.value),
    
    spearman_rho_TCID50_Length = cor(log.TCID50, Average_Length_Sample.3pct, method = "spearman"),
    p_value_TCID50_Length = suppressWarnings(cor.test(log.TCID50, Average_Length_Sample.3pct, method = "spearman")$p.value),
    
    spearman_rho_TCID50_Frequency = cor(log.TCID50, Average_Frequency_Sample.3pct, method = "spearman"),
    p_value_TCID50_Frequency = suppressWarnings(cor.test(log.TCID50, Average_Frequency_Sample.3pct, method = "spearman")$p.value),
    
    spearman_rho_CopyNumber_Indel = cor(log.Copy.Number, Indel.count.3pct, method = "spearman"),
    p_value_CopyNumber_Indel = suppressWarnings(cor.test(log.Copy.Number, Indel.count.3pct, method = "spearman")$p.value),
    
    spearman_rho_CopyNumber_Length = cor(log.Copy.Number, Average_Length_Sample.3pct, method = "spearman"),
    p_value_CopyNumber_Length = suppressWarnings(cor.test(log.Copy.Number, Average_Length_Sample.3pct, method = "spearman")$p.value),
    
    spearman_rho_CopyNumber_Frequency = cor(log.Copy.Number, Average_Frequency_Sample.3pct, method = "spearman"),
    p_value_CopyNumber_Frequency = suppressWarnings(cor.test(log.Copy.Number, Average_Frequency_Sample.3pct, method = "spearman")$p.value)
  )

# Specify the output file
output_file <- "TCID50_CopyNumber_vs_Deletions_combined_spearman_results.txt"

# Sink the results to a text file
sink(output_file)
cat("Spearman's Correlation Results (Combined Analysis)\n")
cat("-----------------------------------------------\n")
cat(paste0(
  "TCID50 vs Indel count - Spearman's rho: ", round(spearman_results_combined$spearman_rho_TCID50_Indel, 3), 
  ", p-value: ", format.pval(spearman_results_combined$p_value_TCID50_Indel, digits = 3), "\n",
  "TCID50 vs Average Length - Spearman's rho: ", round(spearman_results_combined$spearman_rho_TCID50_Length, 3), 
  ", p-value: ", format.pval(spearman_results_combined$p_value_TCID50_Length, digits = 3), "\n",
  "TCID50 vs Average Frequency - Spearman's rho: ", round(spearman_results_combined$spearman_rho_TCID50_Frequency, 3), 
  ", p-value: ", format.pval(spearman_results_combined$p_value_TCID50_Frequency, digits = 3), "\n",
  "Copy Number vs Indel count - Spearman's rho: ", round(spearman_results_combined$spearman_rho_CopyNumber_Indel, 3), 
  ", p-value: ", format.pval(spearman_results_combined$p_value_CopyNumber_Indel, digits = 3), "\n",
  "Copy Number vs Average Length - Spearman's rho: ", round(spearman_results_combined$spearman_rho_CopyNumber_Length, 3), 
  ", p-value: ", format.pval(spearman_results_combined$p_value_CopyNumber_Length, digits = 3), "\n",
  "Copy Number vs Average Frequency - Spearman's rho: ", round(spearman_results_combined$spearman_rho_CopyNumber_Frequency, 3), 
  ", p-value: ", format.pval(spearman_results_combined$p_value_CopyNumber_Frequency, digits = 3), "\n"
))

sink()  # End sink

# Confirmation message
cat(paste("Results saved to", output_file, "\n"))




##### PLOTS #####



######TCID50######

# List of TCID50 comparisons
tcid50_comparisons <- list(
  list(x = "Indel.count.3pct", y = "TCID50", x_title = "Deletion Number", y_title = expression(bold(TCID[50]))),
  list(x = "Average_Length_Sample.3pct", y = "TCID50", x_title = "Deletion Length", y_title = expression(bold(TCID[50]))),
  list(x = "Average_Frequency_Sample.3pct", y = "TCID50", x_title = "Deletion Frequency", y_title = expression(bold(TCID[50])))
)

# Loop through each TCID50 comparison to generate the plots
for (comp in tcid50_comparisons) {
  
  # Create the plot for each comparison with LM trendline
  p <- ggplot(data, aes_string(x = comp$x, y = comp$y)) +
    geom_point(aes(color = Host, shape = Host), size = 4, stroke = 1.2) + # Scatter plot
    geom_smooth(method = "rlm", se = FALSE, color = "black", size = 0.5) + # LM trendline
    scale_y_log10(limits = c(500, 1e7),
                  breaks = c(1e3, 1e4, 1e5, 1e6, 1e7),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +  # Apply log scale to y-axis
    scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
    geom_hline(yintercept = 1000, linetype = "dotted", color = "black", size = 1) + # Dotted line at y = 1000
    labs(x = comp$x_title, y = comp$y_title) +  # Use custom titles for x and y axes
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 24, face = "bold", color = "black", 
                                  margin = margin(t = 15)), # Increase distance for x-axis title
      axis.title.y = element_text(size = 24, face = "bold", color = "black", 
                                  margin = margin(r = 15)),  # Increase distance for y-axis title
      axis.text = element_text(size = 18, color = "black"),
      axis.line = element_line(size = 1, color = "black"),
      axis.ticks = element_line(size = 1),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 18),
      legend.title = element_blank()
    ) +
    scale_color_manual(values = c("#f94040", "#0000FF", "#C5944E", "#0F99B2")) +
    scale_shape_manual(values = c(1, 22, 16, 15)) +
    guides(color = guide_legend(override.aes = list(linetype = "blank")),
           shape = guide_legend(override.aes = list(linetype = "blank")))
  
  # Print the plot
  print(p)
  
  # Save as SVG and PNG in the "Images" folder with simplified filenames
  filename <- paste0("Images/TCID50_vs_", gsub(" ", "_", comp$x_title))
  
  svg(paste0(filename, ".svg"), width = 10, height = 6)
  print(p) 
  dev.off()
  
  ggsave(filename = paste0(filename, ".png"),
         bg = 'white', width = 8, height = 5, device = 'png', dpi = 300)
}

cat("TCID50 plots saved in the 'Images' folder as SVG and PNG files.\n")


######Copy Number#####

# List of Copy Number comparisons
copy_number_comparisons <- list(
  list(x = "Indel.count.3pct", y = "Copy.Number", x_title = "Deletion Number", y_title = "Copy Number"),
  list(x = "Average_Length_Sample.3pct", y = "Copy.Number", x_title = "Deletion Length", y_title = "Copy Number"),
  list(x = "Average_Frequency_Sample.3pct", y = "Copy.Number", x_title = "Deletion Frequency", y_title = "Copy Number")
)

# Loop through each Copy Number comparison to generate the plots
for (comp in copy_number_comparisons) {
  
  # Create the plot for each comparison with LM trendline
  p <- ggplot(data, aes_string(x = comp$x, y = comp$y)) +
    geom_point(aes(color = Host, shape = Host), size = 4, stroke = 1.2) + # Scatter plot
    geom_smooth(method = "rlm", se = FALSE, color = "black", size = 0.5) + # LM trendline
    scale_y_log10(limits = c(5e6, 1.2e9),
                  breaks = c(1e7, 1e8, 1e9),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +  # Apply log scale to y-axis
    scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
    labs(x = comp$x_title, y = comp$y_title) +  # Use custom titles for x and y axes
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 24, face = "bold", color = "black", 
                                  margin = margin(t = 15)), # Increase distance for x-axis title
      axis.title.y = element_text(size = 24, face = "bold", color = "black", 
                                  margin = margin(r = 15)),  # Increase distance for y-axis title
      axis.text = element_text(size = 18, color = "black"),
      axis.line = element_line(size = 1, color = "black"),
      axis.ticks = element_line(size = 1),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 18),
      legend.title = element_blank()
    ) +
    scale_color_manual(values = c("#f94040", "#0000FF", "#C5944E", "#0F99B2")) +
    scale_shape_manual(values = c(1, 22, 16, 15)) +
    guides(color = guide_legend(override.aes = list(linetype = "blank")),
           shape = guide_legend(override.aes = list(linetype = "blank")))
  
  # Print the plot
  print(p)
  
  # Save as SVG and PNG in the "Images" folder with simplified filenames
  filename <- paste0("Images/Copy_Number_vs_", gsub(" ", "_", comp$x_title))
  
  svg(paste0(filename, ".svg"), width = 10, height = 6)
  print(p) 
  dev.off()
  
  ggsave(filename = paste0(filename, ".png"),
         bg = 'white', width = 8, height = 5, device = 'png', dpi = 300)
}

cat("Copy Number plots saved in the 'Images' folder as SVG and PNG files.\n")
