getwd()
setwd("C:/nextcloud/Documents/MiCoDa/1.MiCoDa_V2_EVECluster/2.DetectPrimer_FelipeScript/0.FinalTables_Joao/Results/13_Jun/Github/Figure_2/")

################################################################################
################################################################################
# FIGURE 2: Alignment positions across 16S rRNA regions and sequencing setups
################################################################################
################################################################################
# Packages
# Load required packages
#install.packages("ggplot2")
library(ggplot2)
library(dplyr)
#install.packages("readr")
library(readr)
library(tidyverse)
#install.packages("cowplot")
library(cowplot)
#install.packages("ggpmisc")
library(ggpmisc)

################################################################################
################################################################################
################################################################################
################################################################################
################################  SET 1 FIGURE 2A  #############################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# Upload new table
combined_data <- read.csv("Dataset_A_1.MergeCombineData&Metadata_Set1_EMP_Final.csv")
names(combined_data)

combined_data_long <-  combined_data  %>% dplyr::filter(Threshold == 0.6)

################################################################################
################################################################################
################################################################################
################################################################################
# Final graph for the paper :) EMP
# DESCENDING: V9 (start ~1500) to V1 (start ~0)
combined_data_long$mean_alignment <- rowMeans(combined_data_long[, c("Median_Alignment_start", "Median_Alignment_end")])
combined_data_long_sorted <- combined_data_long[order(-combined_data_long$mean_alignment), ]
combined_data_long_sorted$seq_index <- seq_len(nrow(combined_data_long_sorted))

# Updated plot
unique(combined_data_long_sorted$Platform_2)

my_colors <- c("Illumina" = "#8da0cb", "Ion Torrent" = "#d4ec3f", "454 pyrosequencing" = "#e7298a", "PacBio" = "#fc8d62", "BGISEQ" = "#1b9e77")

top_y <- max(combined_data_long_sorted$seq_index) + 50

################################################################################
################################################################################
# OPTION WITHOUT BLACK LINES & DASHES
Figure2A <- ggplot(combined_data_long_sorted, aes(y = seq_index, color = Platform_2)) +
  geom_segment(aes(x = Median_Alignment_start, xend = Median_Alignment_end,
                   yend = seq_index),
               linewidth = 1.2, alpha = 0.5) +
  scale_color_manual(name = "Platform", values = my_colors) +
  
  # The nine hypervariable regions spanned nucleotides, respectively [numbering based on the E. coli system of nomenclature (Brosius et al., 1978)]
  
  # v1 69–99,
  # V2 137–242,
  # V3 433–497,
  # v4 576–682,
  # V5 822–879,
  # V6 986–1043,
  # V7 1117–1173,
  # V8 1243–1294,
  # V9 1435–1465
  
  # V1–V9 regions as black segments at different heights
  annotate("segment", x = 66,    xend = 99,   y = 16060, yend = 16060, color = "black", size = 2.0) +
  annotate("segment", x = 137,   xend = 242,  y = 16060, yend = 16060, color = "black", size = 2.0) +
  annotate("segment", x = 433,   xend = 497,  y = 16060, yend = 16060, color = "black", size = 2.0) +
  annotate("segment", x = 576,   xend = 682,  y = 16060, yend = 16060, color = "black", size = 2.0) +
  annotate("segment", x = 822,   xend = 879,  y = 16060, yend = 16060, color = "black", size = 2.0) +
  annotate("segment", x = 986,   xend = 1043, y = 16060, yend = 16060, color = "black", size = 2.0) +
  annotate("segment", x = 1117,  xend = 1173, y = 16060, yend = 16060, color = "black", size = 2.0) +
  annotate("segment", x = 1243,  xend = 1294, y = 16060, yend = 16060, color = "black", size = 2.0) +
  annotate("segment", x = 1435,  xend = 1465, y = 16060, yend = 16060, color = "black", size = 2.0) +
  
  geom_vline(xintercept = c(96, 299, 479, 811, 885, 1065,1180,1372,1500),
             linetype = "dashed", color = "black", linewidth = 0.6) +
  
  # General settings
  labs(x = "", y = "", title = " ") +
  theme_light(base_size = 14) +
  scale_x_continuous(limits = c(0, 1600), breaks = seq(0, 1600, 100), expand = c(0,0)) +
  scale_y_continuous(expand = expansion(mult = c(0.00, 0.00))) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    plot.margin = margin(t = 5, r = 20, b = 5, l = 5)
  )
Figure2A





################################################################################
################################################################################
################################################################################
################################################################################
################################  SET 4 FIGURE 2B  #############################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# Upload new table
combined_data <- read.csv("Dataset_B_1.Combine_data_Set5_Set6_FINAL.csv")
names(combined_data)

combined_data_long <-  combined_data  %>% dplyr::filter(Threshold == 0.6)

################################################################################
################################################################################
################################################################################
################################################################################
# Final graph for the paper :) Lebre & Wasimuddi
# DESCENDING: V9 (start ~1500) to V1 (start ~0)
combined_data_long$mean_alignment <- rowMeans(
  combined_data_long[, c("Median_Alignment_start", "Median_Alignment_end")]
)
combined_data_long_sorted <- combined_data_long[order(-combined_data_long$mean_alignment), ]
combined_data_long_sorted$seq_index <- seq_len(nrow(combined_data_long_sorted))

# Updated plot
unique(combined_data_long_sorted$Platform_2)

my_colors <- c(
  "Illumina" = "#8da0cb",
  "Ion Torrent" = "#d4ec3f",
  "454 pyrosequencing" = "#e7298a",
  "PacBio" = "#fc8d62",
  "BGISEQ" = "#1b9e77"
)

top_y <- max(combined_data_long_sorted$seq_index) + 50

################################################################################
################################################################################
# OPTION WITHOUT BLACK LINES & DASHES
Figure2B <- ggplot(combined_data_long_sorted, aes(y = seq_index, color = Platform_2)) +
  geom_segment(
    aes(x = Median_Alignment_start, xend = Median_Alignment_end, yend = seq_index),
    linewidth = 1.2, alpha = 0.5
  ) +
  scale_color_manual(name = "Platform", values = my_colors) +
  
  # The nine hypervariable regions spanned nucleotides, respectively
  # [numbering based on the E. coli system of nomenclature (Brosius et al., 1978)]
  
  # V1–V9 regions as black segments at different heights
  annotate("segment", x = 66,   xend = 99,   y = 240, yend = 240, color = "black", size = 1.5) +
  annotate("segment", x = 137,  xend = 242,  y = 240, yend = 240, color = "black", size = 1.5) +
  annotate("segment", x = 433,  xend = 497,  y = 240, yend = 240, color = "black", size = 1.5) +
  annotate("segment", x = 576,  xend = 682,  y = 240, yend = 240, color = "black", size = 1.5) +
  annotate("segment", x = 822,  xend = 879,  y = 240, yend = 240, color = "black", size = 1.5) +
  annotate("segment", x = 986,  xend = 1043, y = 240, yend = 240, color = "black", size = 1.5) +
  annotate("segment", x = 1117, xend = 1173, y = 240, yend = 240, color = "black", size = 1.5) +
  annotate("segment", x = 1243, xend = 1294, y = 240, yend = 240, color = "black", size = 1.5) +
  annotate("segment", x = 1435, xend = 1465, y = 240, yend = 240, color = "black", size = 1.5) +
  
  geom_vline(
    xintercept = c(96, 299, 479, 811, 885, 1065, 1180, 1372, 1500),
    linetype = "dashed",
    color = "black",
    linewidth = 0.6
  ) +
  
  # General settings
  labs(x = "", y = "", title = " ") +
  theme_light(base_size = 14) +
  scale_x_continuous(
    limits = c(0, 1600),
    breaks = seq(0, 1600, 100),
    expand = c(0, 0)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.00, 0.00))) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    plot.margin = margin(t = 5, r = 20, b = 5, l = 5)
  )
Figure2B





################################################################################
################################################################################
################################################################################
################################################################################
################################  SET 2  FIGURE 2C #############################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# Upload new table
#combined_data <- read.csv("Set_2_MiCoDaV1/1.MergeCombineData&Metadata_Set2_Figure3_EliminateNoV4.csv") # OLD, before checking
combined_data <- readxl::read_xlsx("Dataset_C_1.MergeCombineData&Metadata_Set2_Figure3_EliminateNoV4_CheckSequences.xlsx", sheet = 1)
names(combined_data)

combined_data_long <-  combined_data  %>% dplyr::filter(Threshold == 0.6)

################################################################################
################################################################################
################################################################################
################################################################################
# Final graph for the paper :) MICODA V1
# DESCENDING: V9 (start ~1500) to V1 (start ~0)
combined_data_long$mean_alignment <- rowMeans(combined_data_long[, c("Median_alignment_start", "Median_alignment_end")])
combined_data_long_sorted <- combined_data_long[order(-combined_data_long$mean_alignment), ]
combined_data_long_sorted$seq_index <- seq_len(nrow(combined_data_long_sorted))
names(combined_data_long_sorted)

# Updated plot
unique(combined_data_long_sorted$Platform_2)

my_colors <- c("Illumina" = "#8da0cb", "Ion Torrent" = "#d4ec3f", "454 pyrosequencing" = "#e7298a", "PacBio" = "#fc8d62", "BGISEQ" = "#1b9e77")

top_y <- max(combined_data_long_sorted$seq_index) + 50

################################################################################
################################################################################
# OPTION WITHOUT BLACK LINES & DASHES
Figure2C <- ggplot(combined_data_long_sorted, aes(y = seq_index, color = Platform_2)) +
  geom_segment(aes(x = Median_alignment_start, xend = Median_alignment_end,
                   yend = seq_index),
               linewidth = 1.2, alpha = 0.5) +
  scale_color_manual(name = "Platform", values = my_colors) +
  
  # The nine hypervariable regions spanned nucleotides, respectively [numbering based on the E. coli system of nomenclature (Brosius et al., 1978)]
  
  # v1 69–99,
  # V2 137–242,
  # V3 433–497,
  # v4 576–682,
  # V5 822–879,
  # V6 986–1043,
  # V7 1117–1173,
  # V8 1243–1294,
  # V9 1435–1465
  
  # V1–V9 regions as black segments at different heights
  annotate("segment", x = 66,    xend = 99,   y = 15050, yend = 15050, color = "black", size = 2.0) +
  annotate("segment", x = 137,   xend = 242,  y = 15050, yend = 15050, color = "black", size = 2.0) +
  annotate("segment", x = 433,   xend = 497,  y = 15050, yend = 15050, color = "black", size = 2.0) +
  annotate("segment", x = 576,   xend = 682,  y = 15050, yend = 15050, color = "black", size = 2.0) +
  annotate("segment", x = 822,   xend = 879,  y = 15050, yend = 15050, color = "black", size = 2.0) +
  annotate("segment", x = 986,   xend = 1043, y = 15050, yend = 15050, color = "black", size = 2.0) +
  annotate("segment", x = 1117,  xend = 1173, y = 15050, yend = 15050, color = "black", size = 2.0) +
  annotate("segment", x = 1243,  xend = 1294, y = 15050, yend = 15050, color = "black", size = 2.0) +
  annotate("segment", x = 1435,  xend = 1465, y = 15050, yend = 15050, color = "black", size = 2.0) +
  
  geom_vline(xintercept = c(96, 299, 479, 811, 885, 1065,1180,1372,1500),
             linetype = "dashed", color = "black", linewidth = 0.6) +
  
  # General settings
  labs(
    x = "",
    y = "",
    title = " "
  ) +
  theme_light(base_size = 14) +
  scale_x_continuous(limits = c(0, 1600), breaks = seq(0, 1600, 100), expand = c(0,0)) +
  scale_y_continuous(expand = expansion(mult = c(0.00, 0.00))) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none", # or "none" when I want to remove the legend
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    plot.margin = margin(t = 5, r = 20, b = 5, l = 5)
  )
Figure2C





################################################################################
################################################################################
################################################################################
################################################################################
################################  SET 3  FIGURE 2D #############################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# Upload new table
#combined_data <- read.csv("Set_3_DatathonMiCoDaV2/1.MergeCombineData&Metadata_Set3_Final.csv") # OLD, before checking
combined_data <- readxl::read_xlsx("Dataset_D_1.CombineTotalSamples_AnalysesResultsCaseStudy.xlsx", sheet = 1)
names(combined_data)

combined_data_long <-  combined_data  %>% dplyr::filter(Threshold == 0.6)

################################################################################
################################################################################
################################################################################
################################################################################
# Final graph for the paper :) Datathon
# DESCENDING: V9 (start ~1500) to V1 (start ~0)
combined_data_long$mean_alignment <- rowMeans(combined_data_long[, c("Median_Alignment_start", "Median_Alignment_end")])
combined_data_long_sorted <- combined_data_long[order(-combined_data_long$mean_alignment), ]
combined_data_long_sorted$seq_index <- seq_len(nrow(combined_data_long_sorted))

# Updated plot
unique(combined_data_long_sorted$Platform_2)

my_colors <- c(
  "Illumina" = "#8da0cb",
  "Ion Torrent" = "#d4ec3f",
  "454 pyrosequencing" = "#e7298a",
  "PacBio" = "#fc8d62",
  "BGISEQ" = "#1b9e77"
)

top_y <- max(combined_data_long_sorted$seq_index) + 50

################################################################################
################################################################################
# OPTION WITHOUT BLACK LINES & DASHES
Figure2D <- ggplot(combined_data_long_sorted, aes(y = seq_index, color = Platform_2)) +
  geom_segment(
    aes(x = Median_Alignment_start, xend = Median_Alignment_end, yend = seq_index),
    linewidth = 1.2, alpha = 0.5
  ) +
  scale_color_manual(name = "Platform", values = my_colors) +
  
  # The nine hypervariable regions spanned nucleotides, respectively
  # [numbering based on the E. coli system of nomenclature (Brosius et al., 1978)]
  
  # V1–V9 regions as black segments at different heights
  annotate("segment", x = 66,   xend = 99,   y = 5114, yend = 5114, color = "black", size = 2.0) +
  annotate("segment", x = 137,  xend = 242,  y = 5114, yend = 5114, color = "black", size = 2.0) +
  annotate("segment", x = 433,  xend = 497,  y = 5114, yend = 5114, color = "black", size = 2.0) +
  annotate("segment", x = 576,  xend = 682,  y = 5114, yend = 5114, color = "black", size = 2.0) +
  annotate("segment", x = 822,  xend = 879,  y = 5114, yend = 5114, color = "black", size = 2.0) +
  annotate("segment", x = 986,  xend = 1043, y = 5114, yend = 5114, color = "black", size = 2.0) +
  annotate("segment", x = 1117, xend = 1173, y = 5114, yend = 5114, color = "black", size = 2.0) +
  annotate("segment", x = 1243, xend = 1294, y = 5114, yend = 5114, color = "black", size = 2.0) +
  annotate("segment", x = 1435, xend = 1465, y = 5114, yend = 5114, color = "black", size = 2.0) +
  
  geom_vline(
    xintercept = c(96, 299, 479, 811, 885, 1065, 1180, 1372, 1500),
    linetype = "dashed",
    color = "black",
    linewidth = 0.6
  ) +
  
  # General settings
  labs(x = "", y = "", title = " ") +
  theme_light(base_size = 14) +
  scale_x_continuous(limits = c(0, 1600), breaks = seq(0, 1600, 100), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0.00, 0.00))) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none", # or "bottom" or "none" when I want to remove the legend
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    plot.margin = margin(t = 5, r = 20, b = 5, l = 5)
  )
Figure2D


################################################################################
################################################################################
################################################################################
################################################################################
################################# FINAL GRAPH ##################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# Merge the graphs
# Package cowplot
library(cowplot)
plot_grid(Figure2A, Figure2B, Figure2C, Figure2D, ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"))

################################################################################
################################################################################
################################################################################