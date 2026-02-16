################################################################################
################################################################################
# FIGURE 1: Correlation between number of samples and runtime
################################################################################
################################################################################
# Packages
library(ggplot2)

################################################################################
# Open data frame
df <- readxl::read_xlsx("Processing_time_all_datasets.xlsx", sheet = 7)
names(df)
str(df)
df$Threshold <- as.character(df$Threshold)

################################################################################
# Final plot
ggplot(df, aes(x = samples, y = runtime_minutes, color = Threshold, group = Threshold)) +
  geom_point(size = 4) +
  geom_line(linewidth = 2) +
  scale_x_log10(breaks = df$samples) +
  scale_color_grey(start = 0.8, end = 0.2) +  # You can adjust the gray range here
  labs(
    title = " ",
    x = "Number of samples",
    y = "Run time (minutes)"
  ) +
  theme_light(base_size = 14)
################################################################################
################################################################################
