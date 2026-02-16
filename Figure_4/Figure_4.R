getwd()
setwd("C:/nextcloud/Documents/MiCoDa/1.MiCoDa_V2_EVECluster/2.DetectPrimer_FelipeScript/0.FinalTables_Joao/Results/13_Jun/Github/Figure_4/")

################################################################################
################################################################################
# FIGURE 4: Application of HVRLocator for the selection of V4 16S rRNA gene metabarcoding samples from MiCoDa V2
################################################################################
################################################################################
# Packages
library(ggplot2)
library(cowplot)
library(dplyr)
#install.packages("readr")
library(readr)
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
################################################################################
# Table
combined_data <- read.csv("MiCoDaV2_FinalTableComplete_20250509_HvRegLocator.csv")
names(combined_data)

################################################################################
combined_data_long <- tidyr::pivot_longer(combined_data, cols = c(Media_Length), names_to = "Variable", values_to = "Value") %>%
  filter(Threshold == 0.6)
names(combined_data_long)

# DESCENDING: V9 (start ~1500) to V1 (start ~0)
combined_data_long$mean_alignment <- rowMeans(combined_data_long[, c("Median_Alignment_start", "Median_Alignment_end")])
combined_data_long_sorted <- combined_data_long[order(-combined_data_long$mean_alignment), ]
combined_data_long_sorted$seq_index <- seq_len(nrow(combined_data_long_sorted))

unique(combined_data_long_sorted$Platform_2)

my_colors <- c("Illumina" = "#8da0cb", "Ion Torrent" = "#a6d854", "454 pyrosequencing" = "#e7298a", "PacBio" = "#fc8d62", "BGISEQ" = "#e6ab02")

top_y <- max(combined_data_long_sorted$seq_index) + 50

################################################################################
# Plot
ggplot(combined_data_long_sorted, aes(y = seq_index, color = Platform_2)) +
  geom_segment(aes(x = Median_Alignment_start, xend = Median_Alignment_end,
                   yend = seq_index),
               linewidth = 1.2, alpha = 0.5) + # set alpha here
  scale_color_manual(name = "Platform", values = my_colors) + # colors
  
  # Regiones V1–V9 como segmentos negros en distintas alturas
  #annotate("segment", x = 0,    xend = 96,   y = 7500, yend = 7500, color = "black", size = 2.0) +
  annotate("segment", x = 66,    xend = 99,   y = 42317, yend = 42317, color = "black", size = 2.0) +
  #annotate("text",    x = 48,   y = 8000, label = "V1", size = 5, fontface = "bold") +
  
  #annotate("segment", x = 97,   xend = 299,  y = 10000, yend = 10000, color = "black", size = 2.0) +
  annotate("segment", x = 137,   xend = 242,  y = 42317, yend = 42317, color = "black", size = 2.0) +
  #annotate("text",    x = 198,  y = 8000, label = "V2", size = 5, fontface = "bold") +
  
  #annotate("segment", x = 300,  xend = 479,  y = 7500, yend = 7500, color = "black", size = 2.0) +
  annotate("segment", x = 433,  xend = 497,  y = 42317, yend = 42317, color = "black", size = 2.0) +
  #annotate("text",    x = 389.5,y = 8000, label = "V3", size = 5, fontface = "bold") +
  
  #annotate("segment", x = 480,  xend = 811,  y =10000, yend = 10000, color = "black", size = 2.0) +
  annotate("segment", x = 576,  xend = 682,  y =42317, yend = 42317, color = "black", size = 2.0) +
  #annotate("text",    x = 630,y =8000, label = "V4", size = 5, fontface = "bold") +
  
  #annotate("segment", x = 812,  xend = 885,  y = 7500, yend = 7500, color = "black", size = 2.0) +
  annotate("segment", x = 822,  xend = 879,  y = 42317, yend = 42317, color = "black", size = 2.0) +
  #annotate("text",    x = 848.5,y = 8000, label = "V5", size = 5, fontface = "bold") +
  
  #annotate("segment", x = 886,  xend = 1065, y = 10000, yend = 10000, color = "black", size = 2.0) +
  annotate("segment", x = 986,  xend = 1043, y = 42317, yend = 42317, color = "black", size = 2.0) +
  #annotate("text",    x = 975.5,y = 8000, label = "V6", size = 5, fontface = "bold") +
  
  #annotate("segment", x = 1066, xend = 1180, y = 7500, yend = 7500, color = "black", size = 2.0) +
  annotate("segment", x = 1117, xend = 1173, y = 42317, yend = 42317, color = "black", size = 2.0) +
  #annotate("text",    x = 1123,y = 8000, label = "V7", size = 5, fontface = "bold") +
  
  #annotate("segment", x = 1181, xend = 1372, y = 10000, yend = 10000, color = "black", size = 2.0) +
  annotate("segment", x = 1243, xend = 1294, y = 42317, yend = 42317, color = "black", size = 2.0) +
  #annotate("text",    x = 1276.5,y = 8000, label = "V8", size = 5, fontface = "bold") +
  
  #annotate("segment", x = 1373, xend = 8000, y = 7500, yend = 7500, color = "black", size = 2.0) +
  annotate("segment", x = 1435, xend = 1465, y = 42317, yend = 42317, color = "black", size = 2.0) +
  #annotate("text",    x = 1436.5,y =8000, label = "V9", size = 5, fontface = "bold") +
  
  geom_vline(xintercept = c(96, 299, 479, 811, 885, 1065,1180,1372,1500), linetype = "dashed", color = "black", linewidth = 0.6) +
  labs(
    x = "Alignment Position (bp)",
    y = "Sequence ordered by alignment position start",
    title = " "
  ) +
  theme_light(base_size = 14) +
  scale_x_continuous(breaks = seq(0, 1600, 100), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0.00, 0.00))) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )



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
# PIE CHART
################################################################################
# Table
combined_data <- read.csv("MiCoDaV2_FinalTableComplete_20250509_HvRegLocator.csv")
names(combined_data)

################################################################################
# Summarize counts per HV_region_start and Correct_Assignment
sample_count <- combined_data %>%
  group_by(Predicted_HV_region_start, Threshold) %>%
  summarise(Count = n(), .groups = "drop") %>%
  filter(Threshold == 0.6) %>%
  mutate(Percentage = round(100 * Count / sum(Count), 1),
         Label = paste0(Predicted_HV_region_start, " (", Percentage, "%) ", Count))
names(sample_count)

################################################################################
# Plot 
ggplot(sample_count, aes(x = "", y = Count, fill = Label)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Set3") +
  labs(fill = "Predicted HV Region Start",
       title = "") +
  theme_void()

p3 <- ggplot(sample_count, aes(x = "", y = Count, fill = Label)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Set3") +
  labs(fill = "Predicted HV Region Start",
       title = "") +
  theme_void()


plot_grid(p3, nrow = 1, labels = c("B"))



################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################### PANEL C ####################################
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
# Alignment start
start <- ggplot(combined_data, aes(x = Predicted_HV_region_start, y = Coverage_HV_region_start, fill = as.factor(Predicted_HV_region_start))) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set3") +
  labs(
    title = "",
    x = "Predicted HV region start",
    y = "Coverage HV region"
  ) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1,face ="bold" ),
        legend.position = "none")
start




################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################### PANEL D ####################################
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
# Summarize counts per HV_region_start and Correct_Assignment
sample_count <- combined_data %>%
  group_by(Predicted_HV_region_start, Threshold, True_start) %>%
  summarise(Count = n(), .groups = "drop") %>%
  filter(Threshold == 0.6)
names(sample_count)

sample_count$True_start <- factor(sample_count$True_start, levels = c("TRUE", "FALSE"))

################################################################################
# Legend
ggplot(sample_count, aes(x = Predicted_HV_region_start, y = Count, fill = True_start)) +
  geom_bar(stat = "identity", position = "dodge") +   
  scale_fill_manual(values = c("TRUE" = "#a1d76a", "FALSE" = "#e9a3c9")) +  # colors
  labs(title = " ",
       x = "HV region start",
       y = "Number of samples",
       fill = "Correct assignment") +  # Fix legend label
  theme_light()

################################################################################
# No Legend
p1 <- ggplot(sample_count, aes(x = Predicted_HV_region_start, y = Count, fill = True_start)) +
  geom_bar(stat = "identity", position = "dodge") +   
  scale_fill_manual(values = c("TRUE" = "#a1d76a", "FALSE" = "#e9a3c9")) +  # colors
  labs(title = " ",
       x = "HV region start",
       y = "Number of samples",
       fill = "Correct assignment") +  # Fix legend label
  theme_light() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)
  )  # Remove legend
p1

################################################################################
################################################################################