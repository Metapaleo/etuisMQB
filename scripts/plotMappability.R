library(data.table)
library(ggplot2)

# Read mappability against reference EBV Type 1
setwd("${WORKINGDIR}/resultsEBVT1")

EBVmappabilityT1 <- fread("Mappability_EBVT1.tsv", header = FALSE, sep = "\t", col.names = c("Strain", "Mappability"))
EBVmappabilityT1$Type <- rep("EBV Type 1 reference vs EBV Type 1", nrow(EBVmappabilityT1))

Type2_Strains <- scan("${WORKINGDIR}/listEbvType2.txt", what = "character")

coords <- which(EBVmappabilityT1$Strain %in% Type2_Strains)
EBVmappabilityT1$Type[coords] <- c("EBV Type 1 reference vs EBV Type 2")

# Mappability of reads mapping to T1 and T2

EBVmappabilityT1T2 <- fread("Mappability_EBVT1_T2.tsv", header = FALSE, sep = "\t", col.names = c("Strain", "Mappability"))
EBVmappabilityT1T2$Type <- rep("z_EBV Type 1", nrow(EBVmappabilityT1T2))

coords <- which(EBVmappabilityT1T2$Strain %in% Type2_Strains)
EBVmappabilityT1T2$Type[coords] <- c("z_EBV Type 2")

# Read mappability against reference EBV Type 2
setwd("${WORKINGDIR}/EBV_T2/resultsEBVT2")

EBVmappabilityT2 <- fread("Mappability_EBVT2.tsv", header = FALSE, sep = "\t", col.names = c("Strain", "Mappability"))
EBVmappabilityT2$Type <- rep("EBV Type 2 reference vs EBV Type 1", nrow(EBVmappabilityT2))

coords <- which(EBVmappabilityT2$Strain %in% Type2_Strains)
EBVmappabilityT2$Type[coords] <- c("EBV Type 2 reference vs EBV Type 2")

# Merge the two data sets and plot them together
EBVmappability <- rbind(EBVmappabilityT1, EBVmappabilityT2, EBVmappabilityT1T2)

pdf("${WORKINGDIR}/Mappability_EBVT1_EBVT2.pdf", height = 7, width = 11)
  ggplot(EBVmappability, aes(x=Type, y=Mappability, fill=Type)) +
    geom_violin() +  
    labs(title = "Mappability against EBV type 1 and type 2", 
         x = "", 
         y = "Mappability") +
    theme(axis.text.x = element_text(size = 7, angle = 0, hjust = 0.5), 
          legend.position = "bottom",
          legend.box = "vertical")
dev.off()
