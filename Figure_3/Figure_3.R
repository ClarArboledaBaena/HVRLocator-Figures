################################################################################
################################################################################
# FIGURE 3: Differences in 16S rRNA gene coverage using the same primer set (Primer 515R-806R for V4 region) but different sequencing setups.
################################################################################
################################################################################
# Packages
library(ggplot2)
library(cowplot)

################################################################################
# Table
combined_data <- read.csv("1.MergeCombineData&Metadata_Set2_Figure3_EliminateNoV4.csv")
names(combined_data)

################################################################################
# Colors
my_colors <- c("V4" = "#7fcdbb", "V5" = "#41b6c4", "V6" = "#1d91c0", "V7" = "#225ea8", "V8" = "#253494", "V9" = "#081d58")

################################################################################
start <- ggplot(combined_data, aes(x = Predicted_HV_region_start, y = Coverage_HV_region_start, fill = as.factor(Predicted_HV_region_start))) +
  geom_boxplot() +
  #geom_jitter(width = 0.2) + # Points over the boxplot
  #scale_fill_brewer(palette = "Set3") +
  scale_fill_manual(values = my_colors) +
  facet_wrap(~ Threshold) +
  labs(
    title = "",
    x = "",
    y = "Coverage HV region"
  ) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1,face ="bold" ),
        legend.position = "none",
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(color = "black", face = "bold")
  )
start

################################################################################
end <- ggplot(combined_data, aes(x = Predicted_HV_region_end, y = Coverage_HV_region_end, fill = Predicted_HV_region_end)) +
  geom_boxplot() +
  scale_fill_manual(values = my_colors) +
  #scale_fill_brewer(palette = "Set3") +
  facet_wrap(~ Threshold) +
  labs(
    title = "",
    x = "",
    y = "Coverage HV region"
  ) +
  theme_light() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1, face = "bold"),
    legend.position = "none",
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(color = "black", face = "bold")
    
  )
end

################################################################################
# Merge the graphs
plot_grid(start, end, ncol = 1, nrow = 2, labels = c("A", "B"))

################################################################################
################################################################################
