getwd()
setwd("C:/Documents/MiCoDa/1.MiCoDa_V2_EVECluster/2.DetectPrimer_FelipeScript/0.FinalTables_Joao/Results/13_Jun/Github/Figure_S1")

################################################################################
################################################################################
# FIGURE S1: A) Variation in gene coverage across sequences, B) Number of samples in which primer assignment matched the metadata.
################################################################################
################################################################################
# Packages
library(ggplot2)
library(cowplot)
library(dplyr)
#install.packages("readr")
library(readr)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("cowplot")
library(cowplot)
#install.packages("ggpmisc")
library(ggpmisc)
library(RColorBrewer)

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################### PANEL A ####################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# BOXPLOTS
################################################################################
# Table
combined_data <- read.csv("MiCoDaV2_FinalTableComplete_20250509_HvRegLocator.csv")
names(combined_data)

combined_data <-  combined_data  %>% dplyr::filter(Threshold == 0.6)

################################################################################
# Alignment end
df_plot <- combined_data %>%
  select(Predicted_HV_region_end, Coverage_HV_region_end)

hv_levels <- c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9")

v1_dummy <- data.frame(
  Predicted_HV_region_end = factor("V1", levels = hv_levels),
  Coverage_HV_region_end = NA
)

# Combine with original data
combined_data_aug <- rbind(df_plot, v1_dummy)

################################################################################
# Plot
end <- ggplot(combined_data_aug, aes(x = Predicted_HV_region_end, y = Coverage_HV_region_end, fill = Predicted_HV_region_end)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FFFFB3",
                               "#BEBADA",
                               "#FB8072",
                               "#80B1D3",
                               "#FDB462",
                               "#B3DE69",
                               "#FCCDE5",
                               "#D9D9D9")) +
  #scale_fill_brewer(palette = "Set3") +
  labs(
    title = "",
    x = "Predicted HV region end",
    y = "Coverage HV region"
  ) +
  theme_light() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1, face = "bold"),
    legend.position = "none"
  )
end

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################### PANEL B ####################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# BARPLOT
################################################################################
# Table
combined_data <- read.csv("MiCoDaV2_FinalTableComplete_20250509_HvRegLocator.csv")
names(combined_data)

################################################################################
# Summarize counts per HV_region_end and Correct_Assignment
sample_count <- combined_data %>%
  group_by(Predicted_HV_region_end, Threshold, True_end) %>%
  summarise(Count = n(), .groups = "drop")  %>%
  filter(Threshold == 0.6)
names(sample_count)

################################################################################
# Create the new row as a data frame
new_row <- data.frame(
  Predicted_HV_region_end = "V1",
  Threshold = 0.6,
  True_end = TRUE,#factor("TRUE", levels = levels(sample_count$True_end)),
  Count = 0
)

# Add the row to the existing sample_count table
sample_count <- bind_rows(sample_count, new_row)

sample_count$True_end <- factor(sample_count$True_end, levels = c("TRUE", "FALSE"))

################################################################################
# Simple bar plot with colors by Correct_Assignment_end
ggplot(sample_count, aes(x = Predicted_HV_region_end, y = Count, fill = True_end)) +
  geom_bar(stat = "identity", position = "dodge") +   
  scale_fill_manual(values = c("TRUE" = "#a1d76a", "FALSE" = "#e9a3c9")) +  # Custom colors
  #facet_wrap(~Threshold) +  # Facet by Threshold
  labs(title = " ",
       x = "HV region end",
       y = "Number of samples",
       fill = "Correct assignment") +  # Fix legend label
  theme_light()

################################################################################
# No Legend
p2 <- ggplot(sample_count, aes(x = Predicted_HV_region_end, y = Count, fill = True_end)) +
  geom_bar(stat = "identity", position = "dodge") +   
  scale_fill_manual(values = c("TRUE" = "#a1d76a", "FALSE" = "#e9a3c9")) +  # Custom colors
  #facet_wrap(~Threshold) +  # Facet by Threshold
  labs(title = " ",
       x = "HV region end",
       y = "",
       fill = "Correct assignment") +  # Fix legend label
  theme_light() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)
  )  # Remove legend
p2

################################################################################
################################################################################