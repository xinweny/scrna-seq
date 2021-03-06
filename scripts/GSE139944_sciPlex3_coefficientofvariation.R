#### Packages ####
library(monocle3)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(ggrepel)
library(tidyverse)
library(glue)

#### Functions ####

#### Config ####
# Set working directory
setwd("~/mrc/project/scrna-seq")

set.seed(1)

# Parameters
inhibitors <- c("Luminespib", "Alvespimycin", "Tanespimycin",
                "Trichostatin A")

#### Load data ####
col.data <- read.csv("./data/GSE139944/data/GSM4150378_sciPlex3_pData.txt", 
                     sep=" ", quote='"')
protein.coding <- read.csv("./data/GSE139944/data/gencode.v27.transcripts.bed", 
                           sep="\t", quote='', header=FALSE) %>% 
  filter(V9 == "protein_coding") %>%
  .$V7 %>% unique()

# Format sample names
col.data$cell_product_dose <- paste0(col.data$cell_type, "_",
                                     col.data$product_name, "_", 
                                     col.data$dose)

#### Filtering ####
# Keep valid cells
col.data <- col.data[scan("./data/GSE139944/data/sciPlex3_valid_cells.tsv", character(), quote=""), ]

# Filter cell metadata for selected inhibitors
filt.col.data <- rbind(col.data[col.data$vehicle, ], 
                       col.data[grep(paste(inhibitors, sep="|", collapse="|"), col.data$product_name), ])

# Filter out mouse genes and non-HSP90 targets
filt.counts <- counts(readRDS("./data/GSE139944/data/GSM4150378_sciPlex3_cds_all_cells.RDS"))[protein.coding, rownames(filt.col.data)]

filt.gene.data <- fData(readRDS("./data/GSE139944/data/GSM4150378_sciPlex3_cds_all_cells.RDS"))[protein.coding, ]

# Set rownames for matrix and gene metadata
rownames(filt.counts) <- filt.gene.data$gene_short_name
rownames(filt.gene.data) <- filt.gene.data$gene_short_name

# Output frequency table of no. of cells per sample
cell.freq <- table(filt.col.data$cell_product_dose)
filt.col.data$cell_freq <- cell.freq[match(filt.col.data$cell_product_dose, rownames(cell.freq))]

#### Coefficient of variation ####
# Parameters
cell.type <- "K562"
product <- "Alvespimycin"
d <- c(10, 100, 1000, 10000)
ncells.vehicle <- 250

samples <- c("Vehicle_0", paste0(product, "_", d))

# Initialise data frame
cv.df <- data.frame(mean_rpm=numeric(),
                    cv=numeric(),
                    gene_short_name=character(),
                    sample=character())

# Calculate mean RPM and CV for product across different doses
for (i in 1:length(samples)) {
  if (samples[i] == "Vehicle_0") {
    # Randomly sample vehicle cells
    sample.cells <- filt.col.data[filt.col.data$vehicle, ] %>%
      filter(cell_type == cell.type, dose == 0) %>%
      sample_n(ncells.vehicle) %>%
      rownames()
  } else {
    sample.cells <- filt.col.data[grep(product, filt.col.data$product_name), ] %>%
      filter(cell_type == cell.type, dose == d[i - 1]) %>% 
      rownames() 
  }
  
  sample.counts <- filt.counts[, sample.cells]
  sample.rpm <- sample.counts / (colSums(sample.counts) / 10^6)
  means.sample <- rowMeans(sample.rpm)
  cv.sample <- rowSds(sample.rpm) / means.sample
  
  cv.df <- rbind(cv.df,
                 data.frame(mean_rpm=means.sample,
                            cv=cv.sample,
                            gene_short_name=names(means.sample),
                            sample=rep(samples[i], length(means.sample))))
}

# Get cell number information
# ncells <- c(stack(cell.freq)[grep(glue("{cell.type}_Vehicle_0"), stack(cell.freq)$ind), "values"])
ncells <- c(ncells.vehicle)

for (i in 1:length(d)) {
  ncells[i + 1] <- stack(cell.freq)[grep(glue("{cell.type}_{product}.+_{d[i]}$"), stack(cell.freq)$ind), "values"]
}

names(ncells) <- paste0(cell.type, "_", samples)

# CV-mean RPM plot
cv.rpm.plot.list <- list()
for (i in 1:length(samples)) {
  cv.rpm.plot <- ggplot(cv.df %>% filter(sample == samples[i]),
                        aes(x=log10(mean_rpm + 1), y=cv)) +
    geom_point(size=0.2) +
    xlab("log10(Mean RPM + 1)") +
    ylab("Coefficient of variation (CV)") +
    labs(title=glue("{samples[i]} (nCells={ncells[[paste0(cell.type, '_', samples[i])]]}), N={cv.df %>% drop_na(cv) %>% filter(sample == samples[i]) %>% nrow()}"))
  
  cv.rpm.plot.list[[i]] <- cv.rpm.plot
}
ggsave(glue("./data/GSE139944/transcriptional_noise/CVvsRPM/sciPlex3_CVvsRPM_{cell.type}_{product}.png"),
       grid.arrange(grobs=cv.rpm.plot.list, 
                    ncol=2, nrow=3,
                    top=textGrob(glue("Coefficient of variation (CV) vs. mean RPM of protein coding genes
                                      in {cell.type} cells for {product} across different doses and Vehicle"),
                                 gp=gpar(fontsize=20, font=2))),
       width=40, height=60, units="cm")

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
# Plot mean RPM distribution across different doses
png(glue("./data/GSE139944/transcriptional_noise/meanRPM/sciPlex3_meanRPM_{cell.type}_{product}.png"),
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

# Plot coefficient of variation distribution across different doses
png(glue("./data/GSE139944/transcriptional_noise/CV/sciPlex3_CV_{cell.type}_{product}.png"),
    width=3000, height=2000, res=300)
ggplot(cv.df, aes(x=sample, y=cv)) +
  geom_boxplot() +
  xlab("") +
  ylab("CV") +
  scale_x_discrete(labels=paste0(samples, " (", ncells, ")")) +
  ggtitle(glue("Coefficient of variation of genes across {product} treatments
               in {cell.type} cells"))
dev.off()

# Plot CV vs. CV of product_dose vs. vehicle
cv.cv.plot.list <- list()
for (i in 1:length(d)) {
  cv.cv <- data.frame(cv.df[cv.df$sample == "Vehicle_0", "cv"],
                      cv.df[cv.df$sample == paste0(product, "_", d[i]), "cv"],
                      cv.df[cv.df$sample == "Vehicle_0", "gene_short_name"])
  names(cv.cv) <- c("vehicle_cv", "treatment_cv", "gene_short_name")
  
  cv.cv$direction <- ifelse(cv.cv$treatment_cv > cv.cv$vehicle_cv, "UP", "DOWN")
  
  cv.cv.plot <- cv.cv %>% drop_na(vehicle_cv, treatment_cv) %>%
    ggplot(aes(x=log2(vehicle_cv), y=log2(treatment_cv))) +
    geom_abline(intercept=0, 
                slope=1,
                linetype="dashed", size=0.5, color="gray") +
    geom_hline(yintercept=0) +
    geom_vline(xintercept=0) +
    geom_point(aes(fill=direction),
               colour="black", pch=21) +
    xlab("log2(CV (Vehicle_0))") + 
    ylab(glue("log2(CV ({product}_{d[i]}))")) +
    labs(title=glue("{product}_{d[i]} (nCells={ncells[[paste0(cell.type, '_', product, '_', d[i])]]}) vs. Vehicle_0 (nCells={ncells[[paste0(cell.type, '_Vehicle_0')]]})"),
         subtitle=glue("N={cv.cv %>% drop_na(vehicle_cv, treatment_cv) %>% nrow()}, nUP={cv.cv %>% drop_na(vehicle_cv, treatment_cv) %>% filter(direction == 'UP') %>% nrow()}, nDOWN={cv.cv %>% drop_na(vehicle_cv, treatment_cv) %>% filter(direction == 'DOWN') %>% nrow()}"),
         fill="Direction")
  
  cv.cv.plot.list[[i]] <- cv.cv.plot
}
ggsave(glue("./data/GSE139944/transcriptional_noise/CVvsCV/sciPlex3_CVvsCV_{cell.type}_{product}.png"),
       grid.arrange(grobs=cv.cv.plot.list, 
                    ncol=2, nrow=2,
                    top=textGrob(glue("Coefficient of variation (CV) of protein coding genes
                                      in {cell.type} cells for {product} vs. Vehicle across different doses"),
                                 gp=gpar(fontsize=20, font=2))),
       width=60, height=40, units="cm")

