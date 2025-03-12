# Load necessary libraries
library(dplyr)
library(ggplot2)
library(scales)
library(MASS)


## ---------------------------------------------------------------------------------------
data <- read.csv("GroupedtiterR10.csv", sep=',', header=T)

data <- data %>%
  filter(!is.na(Copy.Number)) %>%
   filter(!ID %in% c(9770, 10318))

# Ensure TCID50 and Copy Number are numeric
data <- data %>%
  mutate(TCID50 = as.numeric(TCID50),
         Copy.Number = as.numeric(`Copy.Number`),
         Host = as.factor(Host))

data <- data %>%
  mutate(log.TCID50 = log10(TCID50 + 1)) %>%
  mutate(log.Copy.Number = log10(Copy.Number))



## ---------------------------------------------------------------------------------------
# Calculate Spearman's correlation for each Host
spearman_results <- data %>%
  group_by(Host) %>%
  summarise(
    spearman_rho = cor(log.TCID50, log.Copy.Number, method = "spearman"),
    p_value = suppressWarnings(cor.test(log.TCID50, log.Copy.Number, method = "spearman")$p.value) # Suppress warnings
  )

# Specify the output file
output_file <- "TCID50_vs_CopyNumber_spearman_results.txt"

# Sink the results to a text file
sink(output_file)
cat("Spearman's Correlation Results by Host\n")
cat("--------------------------------------\n")
for (i in 1:nrow(spearman_results)) {
  cat(paste0(
    "Host: ", spearman_results$Host[i], 
    ", Spearman's rho: ", round(spearman_results$spearman_rho[i], 3), 
    ", p-value: ", format.pval(spearman_results$p_value[i], digits = 3), "\n"
  ))
}
sink()  # End sink

# Confirmation message
cat(paste("Results saved to", output_file, "\n"))


###PLOTS
## ---------------------------------------------------------------------------------------
p <- ggplot(data, aes(x = Copy.Number, y = TCID50, color = Host, shape = Host, linetype = Host)) +
  geom_point(size = 4.5, stroke = 1.2) +
  geom_smooth(data = subset(data, Host %in% c("BALBF", "BALBM")), method = "rlm", se = FALSE, size = 0.7) +  # Thicker dashed lines for BALBF & BALBM
  geom_smooth(data = subset(data, Host %in% c("BL6F", "BL6M")), method = "rlm", se = FALSE, size = 0.5) +  # Thinner solid lines for BL6F & BL6M
  scale_x_log10(
    limits = c(9e5, 1.2e9),
    breaks = c(1e6, 1e7, 1e8, 1e9),
    labels = scales::trans_format("log10", scales::math_format(10^.x))  # Use proper scientific formatting
  ) +
  scale_y_log10(
    limits = c(500, 1e7),
    breaks = c(1e3, 1e4, 1e5, 1e6, 1e7),
    labels = scales::trans_format("log10", scales::math_format(10^.x))  # Use proper scientific formatting
  ) +
  geom_hline(yintercept = 1e3, linetype = "dotted", size = 1) +
  labs(x = "Copy number", y = expression(bold(TCID[50]))) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 24, face = "bold", color = "black", 
                                margin = margin(t = 15)),
    axis.title.y = element_text(size = 24, face = "bold", color = "black", 
                                margin = margin(r = 15)),
    axis.text = element_text(size = 20, face = "bold", color = "black"),
    axis.line = element_line(size = 1.3, color = "black"),
    axis.ticks = element_line(size = 1.2),
    axis.ticks.length = unit(8, "pt"),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 18),
    legend.title = element_blank()
  ) +
  scale_color_manual(values = c("BALBF" = "#f94040", "BALBM" = "#0000FF", "BL6F" = "#C5944E", "BL6M" = "#0F99B2")) +
  scale_shape_manual(values = c("BALBF" = 1, "BALBM" = 22, "BL6F" = 16, "BL6M" = 15)) +
  scale_linetype_manual(values = c("BALBF" = "dashed", "BALBM" = "dashed", "BL6F" = "solid", "BL6M" = "solid")) +
  guides(color = guide_legend(override.aes = list(linetype = "blank")),
         shape = guide_legend(override.aes = list(linetype = "blank")))

print(p)

svg("TCID50vsCopyNumber.svg", width = 10, height = 6)
print(p) # Replace `p` with your ggplot object
dev.off()

ggsave(filename = "TCID50vsCopyNumber.png", bg = 'white', width = 8, height = 5, device='png', dpi=300)


