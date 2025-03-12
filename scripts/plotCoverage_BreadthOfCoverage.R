library(ggplot2)
library(data.table)
library(dplyr)
library(cowplot)

## Plot EBV Type 1
setwd("${WORKINGDIR}/resultsEBVT1/bamBaseCoverage")

# Selected genomes to plot
genomes <- c("1998_5_4", "KC207814", "MG298911")

# Define regions for EBV Type 1
regions_T1 <- data.frame(
  Region = c("EBNA1", "EBNA2", "EBNA3A", "EBNA3B", "EBNA3C"),
  Start = c(55189, 36098, 79950, 83065, 86083),
  End = c(55361, 37739, 82960, 89502, 89135)
)

# Read and process selected genome files
df_combined_T1 <- bind_rows(lapply(genomes, function(strain) {
  # Read data
  df <- fread(paste0(strain, "_depth.rmdup.tsv"), header = FALSE, sep = "\t", col.names = c("Strain", "Position", "Coverage"))

  # Expand dataframe to handle overlapping regions
  df_expanded <- df %>%
    left_join(regions_T1, by = character()) %>%  # Cross join to add all regions
    filter(Position >= Start & Position <= End) %>%  # Keep only positions in regions
    mutate(Strain = strain)  # Add strain column (since it's missing from the TSV)

  return(df_expanded)
}))

# Plot with independent X and Y axes per facet
p1 <- ggplot(df_combined_T1, aes(x = Position, y = Coverage)) +
  geom_line() +
  facet_grid(Strain ~ Region, scales = "free") +  # Allow each region and strain to have independent axes
  labs(title = "EBNA genes Mapping Coverage, EBV Type 1 reference",
       x = "Genomic Position",
       y = "Coverage") +
  theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7))



## Plot EBV Type 2

setwd("${WORKINGDIR}/EBV_T2/resultsEBVT2/bamBaseCoverage")

# Selected genomes to plot
genomes <- c("1998_5_4", "KC207814", "MG298911")

# Define regions for EBV Type 2
regions_T2 <- data.frame(
  Region = c("EBNA1", "EBNA2", "EBNA3A", "EBNA3B", "EBNA3C"),
  Start = c(96492, 36201, 80026, 83074, 86654),
  End = c(98417, 37565, 82888, 86532, 89937)
)

df_combined_T2 <- bind_rows(lapply(genomes, function(strain) {
  # Read data
  df <- fread(paste0(strain, "_depth.rmdup.tsv"), header = FALSE, sep = "\t", col.names = c("Strain", "Position", "Coverage"))

  # Expand dataframe to handle overlapping regions
  df_expanded <- df %>%
    left_join(regions_T2, by = character()) %>%  # Cross join to add all regions
    filter(Position >= Start & Position <= End) %>%  # Keep only positions in regions
    mutate(Strain = strain)  # Add strain column (since it's missing from the TSV)

  return(df_expanded)
}))

# Plot with independent X and Y axes per facet
p2 <- ggplot(df_combined_T2, aes(x = Position, y = Coverage)) +
  geom_line() +
  facet_grid(Strain ~ Region, scales = "free") +  # Allow each region and strain to have independent axes
  labs(title = "EBNA genes Mapping Coverage, EBV Type 2 reference",
       x = "Genomic Position",
       y = "Coverage") +
  theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7))

## Plot breadth of coverage
setwd("${WORKINGDIR}")

## Read file with breadth of coverage information, 4 columns as in example:
# Strain  EBV_Type1_Reference_(B95-8) EBV_Type2_Reference_(AG876) EBV_Type
# 1998_5_4  80.1301 75.8798 EBV_Type1
# NC_009334 80.1557 85.3285 EBV_Type2

breadthOfCov <- fread("EBV_Type_breadhtOfCov.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

p3 <- ggplot(breadthOfCov, aes(x = EBV_Type, y = EBV_Type1_Reference, fill = EBV_Type)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  # Semi-transparent boxplot
  geom_point(position = position_jitter(width = 0.25, height = 0),
             shape = 21, size = 1, color = "black", fill = "black", alpha = 0.5) +
  labs(title = "Breadth of Coverage mapping to EBV Type 1",
       x = "EBV Type",
       y = "Breadth of Coverage") +
  theme()

p4 <- ggplot(breadthOfCov, aes(x = EBV_Type, y = EBV_Type2_Reference, fill = EBV_Type)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  # Semi-transparent boxplot
  geom_point(position = position_jitter(width = 0.25, height = 0),
             shape = 21, size = 1, color = "black", fill = "black", alpha = 0.5) +
  labs(title = "Breadth of Coverage mapping to EBV Type 2",
       x = "EBV Type",
       y = "Breadth of Coverage") +
  theme()

# Save mapping coverage for T1 and T2, as well as BoC
pdf("EBV_CoverageMapping_BoC.pdf", height = 7, width = 11)
  plot_grid(p1, p3, p2, p4, labels = c("A", "B", "C", "D"), ncol = 2, rel_widths = c(1,0.6))
dev.off()
