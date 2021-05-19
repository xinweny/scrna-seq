#### Packages ####
library(monocle3)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(glue)

#### Functions ####

#### Config ####
# Set working directory
setwd("~/mrc/project/scrna-seq")

# Parameters
inhibitors <- c("Luminespib", "Alvespimycin", "Tanespimycin",
                "Trichostatin A")

#### Load data ####
col.data <- read.csv("./GSE139944/data/GSM4150378_sciPlex3_pData.txt", 
                     sep=" ", quote='"')
protein.coding <- read.csv("./GSE139944/data/gencode.v27.transcripts.bed", 
                           sep="\t", quote='', header=FALSE) %>% 
  filter(V9 == "protein_coding") %>%
  .$V7 %>% unique()

# Format sample names
col.data$cell_product_dose <- paste0(col.data$cell_type, "_",
                                     col.data$product_name, "_", 
                                     col.data$dose)

#### Filtering ####
# Keep valid cells
col.data <- col.data[scan("./GSE139944/data/sciPlex3_valid_cells.tsv", character(), quote=""), ]

# Output frequency table of no. of cells per sample
cell.freq <- table(col.data$cell_product_dose)
col.data$cell_freq <- cell.freq[match(col.data$cell_product_dose, rownames(cell.freq))]

# Filter cell metadata for selected inhibitors
filt.col.data <- rbind(col.data[col.data$vehicle, ], 
                       col.data[grep(paste(inhibitors, sep="|", collapse="|"), col.data$product_name), ])

# Filter out mouse genes and non-HSP90 targets
filt.counts <- counts(readRDS("./GSE139944/data/GSM4150378_sciPlex3_cds_all_cells.RDS"))[protein.coding, rownames(filt.col.data)]

filt.gene.data <- fData(readRDS("./GSE139944/data/GSM4150378_sciPlex3_cds_all_cells.RDS"))[protein.coding, ]

# Set rownames for matrix and gene metadata
rownames(filt.counts) <- filt.gene.data$gene_short_name
rownames(filt.gene.data) <- filt.gene.data$gene_short_name

# Keep cells with ≥ 500 different transcripts expressed
filt.counts <- filt.counts[, colSums(filt.counts > 0) >= 500]
filt.col.data <- filt.col.data[colnames(filt.counts), ]

#### Mean expression and CV ####
# Parameters
cell.type <- "MCF7"
product <- "Alvespimycin"
d <- c(10, 100, 1000, 10000)

# Calculate gene expression summary for vehicle
vehicle.cells <- filt.col.data[filt.col.data$vehicle, ] %>%
  filter(cell_type == cell.type, dose == 0) %>%
  rownames()

vehicle.counts <- filt.counts[, vehicle.cells]
vehicle.rpm <- vehicle.counts / (colSums(vehicle.counts) / 10^6)
means.vehicle <- rowMeans(vehicle.rpm)
cv.vehicle <- rowSds(vehicle.rpm) / means.vehicle

cv.df <- data.frame(mean_rpm=means.vehicle,
                    cv=cv.vehicle,
                    gene_short_name=names(means.vehicle),
                    sample=rep("Vehicle_0", length(means.vehicle)))

# Calculate gene expression summary for product across different doses
for (i in 1:length(d)) {
  sample.cells <- filt.col.data[grep(product, filt.col.data$product_name), ] %>%
    filter(cell_type == cell.type, dose == d[i]) %>% 
    rownames()
  
  sample.counts <- filt.counts[, sample.cells]
  sample.rpm <- sample.counts / (colSums(sample.counts) / 10^6)
  means.sample <- rowMeans(sample.rpm)
  cv.sample <- rowSds(sample.rpm) / means.sample
  
  cv.df <- rbind(cv.df,
                 data.frame(mean_rpm=means.sample,
                            cv=cv.sample,
                            gene_short_name=names(means.sample),
                            sample=rep(paste0(product, "_", d[i]), length(means.sample))))
}

samples <- unique(cv.df$sample)

# Get cell number information
ncells <- c(stack(cell.freq)[grep(paste0(cell.type, "_Vehicle_0"), stack(cell.freq)$ind), "values"])

for (i in 1:length(d)) {
  ncells[i + 1] <- stack(cell.freq)[grep(glue("{cell.type}_{product}.+_{d[i]}$"), stack(cell.freq)$ind), "values"]
}

names(ncells) <- paste0(cell.type, "_", samples)

# CV-mean RPM plot
for (i in 1:length(samples)) {
  cv.rpm.plot <- ggplot(cv.df %>% filter(sample == samples[i]),
                        aes(x=log10(mean_rpm + 1), y=cv)) +
    geom_point(size=0.2) +
    xlab("log10(Mean RPM + 1)") +
    ylab("Coefficient of variation (CV)") +
    labs(title=glue("CV-mean plot of N={cv.df %>% filter(sample == samples[i]) %>% nrow()} genes
                    for {cell.type} {samples[i]} cells (nCells={ncells[[paste0(cell.type, '_', samples[i])]]})"))
  ggsave(cv.rpm.plot,
         file=glue("./GSE139944/transcriptional_noise/sciPlex3_CVvsRPM_{cell.type}_{samples[i]}.png"),
         width=30, height=20, units="cm")
}

# Filter out genes not expressed in any sample
keep.genes <- cv.df %>%
  group_by(gene_short_name) %>%
  summarize(mean=mean(mean_rpm)) %>%
  filter(mean > 0) %>%
  .$gene_short_name %>% unique()

cv.df <- cv.df[cv.df$gene_short_name %in% keep.genes, ]

# Set order of samples
cv.df$sample <- factor(cv.df$sample, levels=samples)

#### Visualisation ####
# Plot mean gene expression
png(glue("./GSE139944/transcriptional_noise/sciPlex3_meanRPM_{cell.type}_{product}.png"),
    width=3000, height=2000, res=300)
ggplot(cv.df, aes(x=sample, y=log10(mean_rpm + 1))) +
  geom_violin() +
  stat_summary(fun=mean, geom="errorbar",
               aes(ymax = ..y.., ymin = ..y.., group = factor(sample)),
               width=0.75, colour="red") +
  xlab("") +
  ylab("log10(Mean RPM + 1)") +
  scale_x_discrete(labels=paste0(samples, " (", ncells, ")")) +
  ggtitle(glue("Mean RPM of N={length(keep.genes)} genes across {product} treatments
               in {cell.type} cells"))
dev.off()

# Plot coefficient of variation
png(glue("./GSE139944/transcriptional_noise/sciPlex3_CV_{cell.type}_{product}.png"),
    width=3000, height=2000, res=300)
ggplot(cv.df, aes(x=sample, y=cv)) +
  geom_boxplot() +
  ylim(1, 5) +
  xlab("") +
  ylab("CV") +
  scale_x_discrete(labels=paste0(samples, " (", ncells, ")")) +
  ggtitle(glue("Coefficient of variation of genes across {product} treatments
               in {cell.type} cells"))
dev.off()

# Plot CV vs. CV of 2 conditions
dose <- 100

cv.cv <- data.frame(cv.df[cv.df$sample == "Vehicle_0", "cv"],
                    cv.df[cv.df$sample == paste0(product, "_", dose), "cv"],
                    cv.df[cv.df$sample == "Vehicle_0", "gene_short_name"])
names(cv.cv) <- c("vehicle_cv", "treatment_cv", "gene_short_name")

cv.cv$direction <- ifelse(cv.cv$treatment_cv > cv.cv$vehicle_cv, "UP", "DOWN")

png(glue("./GSE139944/transcriptional_noise/sciPlex3_CVCV_{cell.type}_{product}_{dose}_vs_Vehicle_0.png"),
    width=3000, height=2000, res=300)
cv.cv %>% drop_na(vehicle_cv, treatment_cv) %>%
ggplot(aes(x=log2(vehicle_cv), y=log2(treatment_cv))) +
  geom_abline(intercept=0, 
              slope=1,
              linetype="dashed", size=0.5, color="gray") +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(aes(fill=direction),
             colour="black", pch=21) +
  xlab("log2(CV (Vehicle_0))") + 
  ylab(glue("log2(CV ({product}_{dose}))")) +
  labs(title=glue("Coefficients of variation of N={cv.cv %>% drop_na(vehicle_cv, treatment_cv) %>% nrow()} genes in {cell.type} cells"),
       subtitle=glue("{product}_{dose} (nCells={ncells[[paste0(cell.type, '_', product, '_', dose)]]}) vs. Vehicle_0 (nCells={ncells[[paste0(cell.type, '_', product, '_', dose)]]})
                     nUP={cv.cv %>% drop_na(vehicle_cv, treatment_cv) %>% filter(direction == 'UP') %>% nrow()}, nDOWN={cv.cv %>% drop_na(vehicle_cv, treatment_cv) %>% filter(direction == 'DOWN') %>% nrow()}"),
       fill="Direction (Treatment vs. Vehicle)")
dev.off()
  

