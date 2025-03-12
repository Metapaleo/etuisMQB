library(data.table)
library(ggplot)
library(cowplot)

## PCA EBV Type 1

setwd("${WORKING}/resultsEBVT1/plink")

pca_data <- fread("EBV.evec", header = TRUE, sep = "\t", select = c(1:14))

# Replace the names of the PC columns
colnames(pca_data) <- c("Sample", paste0("PC", 1:(ncol(pca_data)-4)), "Geography", "EBV_Type", "Color")

plot1 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Geography, shape = EBV_Type)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_color_manual(values = setNames(pca_data$Color, pca_data$Geography)) +
    scale_shape_manual(values = c("EBV Type 1" = 16, "EBV Type 2" = 17, "Recombinant" = 15)) +
    geom_text(data = subset(pca_data, Sample == "1998_5_4"),  # Add label to sample `1998_5_4`
              aes(label = Sample),
              hjust = -0.2, vjust = -0.5, size = 4, color = "black") +
    labs(x = "PC1 (11.023%)", y = "PC2 (8.315%)") +
    theme_minimal() +
    theme(legend.position = "right", legend.text = element_text(size = 11),
          # Panel border
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          # Remove gridline
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # Increase label size
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          # Increase tick labels
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12)
    )

plot2 <- ggplot(pca_data, aes(x = PC2, y = PC3, color = Geography, shape = EBV_Type)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_color_manual(values = setNames(pca_data$Color, pca_data$Geography)) +
    scale_shape_manual(values = c("EBV Type 1" = 16, "EBV Type 2" = 17, "Recombinant" = 15)) +
    geom_text(data = subset(pca_data, Sample == "1998_5_4"),  # Add label to sample `1998_5_4`
              aes(label = Sample),
              hjust = -0.2, vjust = -0.5, size = 4, color = "black") +
    labs(x = "PC2 (8.315%)", y = "PC3 (6.523%)") +
    theme_minimal() +
    theme(legend.position = "none",
          # Panel border
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          # Remove gridline
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # Increase label size
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          # Increase tick labels
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12)
    )

# Extract legend from the plot
legend <- get_legend(plot1)

# Remove legend from plot1
plot1 <- plot1 + theme(legend.position = "none")
# Combine both plots and the legend adjusting the ratio
combined_plot <- plot_grid(plot1, plot2, legend, ncol = 3, rel_widths = c(1, 1, 0.5))
# Remove the legend for manuscript figure
#combined_plot <- plot_grid(plot1, plot2, ncol = 2, rel_widths = c(1, 1))
# Title of the general plot
title <- ggdraw() +
  draw_label("PCA Plot EBV Type 1 Reference", fontface = "bold", x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))

final_plot <- plot_grid(title, combined_plot, ncol = 1, rel_heights = c(0.1, 1))

pdf("../PCA_EBV_Type1_Recombinant.pdf", height = 7, width = 10)
  print(final_plot)
dev.off()


## PCA EBV Type 2

setwd("${WORKINGDIR}/EBV_T2/resultsEBVT2/plink")

# Read the .evec file
pca_data <- fread("EBV.evec", header = TRUE, sep = "\t", select = c(1:14))
colnames(pca_data) <- c("Sample", paste0("PC", 1:(ncol(pca_data)-4)), "Geography", "EBV_Type", "Color")

# Create two plots, PC1 vs PC2 and PC2 vs PC3
plot1 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Geography, shape = EBV_Type)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = setNames(pca_data$Color, pca_data$Geography)) +
  scale_shape_manual(values = c("EBV Type 1" = 16, "EBV Type 2" = 17, "Recombinant" = 15)) +
  geom_text(data = subset(pca_data, Sample == "1998_5_4"),  # Add label to sample `1998_5_4`
              aes(label = Sample),
              hjust = -0.2, vjust = -0.5, size = 4, color = "black") +
  labs(x = "PC1 (10.762%)", y = "PC2 (8.137%)") +
  theme_minimal() +
  theme(legend.position = "right", legend.text = element_text(size = 11),
        # Panel border
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        # Remove gridline
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Increase label size
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        # Increase tick labels
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)
  )

plot2 <- ggplot(pca_data, aes(x = PC2, y = PC3, color = Geography, shape = EBV_Type)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = setNames(pca_data$Color, pca_data$Geography)) +
  scale_shape_manual(values = c("EBV Type 1" = 16, "EBV Type 2" = 17, "Recombinant" = 15)) +
  geom_text(data = subset(pca_data, Sample == "1998_5_4"),  # Add label to sample `1998_5_4`
              aes(label = Sample),
              hjust = -0.2, vjust = -0.5, size = 3, color = "black") +
  labs(x = "PC2 (8.137%)", y = "PC3 (6.827%)") +
  theme_minimal() +
  theme(legend.position = "none", legend.text = element_text(size = 11),
        # Panel border
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        # Remove gridline
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Increase label size
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        # Increase tick labels
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)
  )

# Extract legend from the plot
legend <- get_legend(plot1)
# Remove legend from plot1
plot1 <- plot1 + theme(legend.position = "none")
# Combine both plots and the legend adjusting the ratio
combined_plot <- plot_grid(plot1, plot2, legend, ncol = 3, rel_widths = c(1, 1, 0.5))
# Title of the general plot
title <- ggdraw() +
  draw_label("PCA Plot EBV Type 2 Reference", fontface = "bold", x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))

final_plot <- plot_grid(title, combined_plot, ncol = 1, rel_heights = c(0.1, 1))

pdf("../PCA_EBV_Type2_Recombinant.pdf", height = 7, width = 12)
  print(final_plot)
dev.off()
