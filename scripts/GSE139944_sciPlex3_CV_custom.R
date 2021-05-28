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

proteo.list <- read.csv("./data/GSE139944/data/proteostasis_gene_list_16_03_21_NON_CORE0_CORE1.csv",
                        sep="\t")
proteo.genes <- proteo.list[proteo.list$CORE == "CORE", c("Human_gene_ID")]

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

filt.gene.data$id <- gsub("\\.[0-9_A-Z]+$", "", filt.gene.data$id)

# Output frequency table of no. of cells per sample
cell.freq <- table(filt.col.data$cell_product_dose)
filt.col.data$cell_freq <- cell.freq[match(filt.col.data$cell_product_dose, rownames(cell.freq))]

# RPM normalisation
filt.rpms <- filt.counts / (colSums(filt.counts) / 10^6)

#### Parameters ####
cell.type <- "K562"
product <- "Alvespimycin"
d <- c(10, 100, 1000, 10000)

ncells.vehicle <- 250
nrep.vehicle <- 3

vehicles <- paste0("Vehicle_0_rep", 1:nrep.vehicle)
samples <- c(vehicles, paste0(product, "_", d))

#### Coefficient of variation ####
# Get list of cells for each sample
sample.cells.list <- list()

for (i in 1:length(samples)) {
  if (grepl("Vehicle_0", samples[i])) {
    # Randomly sample vehicle cells
    sample.cells.list[[samples[i]]] <- filt.col.data[filt.col.data$vehicle, ] %>%
      filter(cell_type == cell.type, dose == 0) %>%
      sample_n(ncells.vehicle) %>%
      rownames()
  } else {
    sample.cells.list[[samples[i]]] <- filt.col.data[grep(product, filt.col.data$product_name), ] %>%
      filter(cell_type == cell.type, dose == d[i - nrep.vehicle]) %>% 
      rownames() 
  }
}

# Get cell number information
ncells <- c(rep(ncells.vehicle, nrep.vehicle),
            stack(cell.freq)[grep(glue("{cell.type}_{product}.+_10+$"), 
                                  stack(cell.freq)$ind), 
                             "values"])
names(ncells) <- paste0(cell.type, "_", samples)

# Calculate CV and mean RPM for each gene in each sample
cv.df <- data.frame(mean_rpm=numeric(),
                    cv=numeric(),
                    gene_short_name=character(),
                    sample=character())
for (i in 1:length(samples)) {
  sample.rpm <- filt.rpms[, sample.cells.list[[samples[i]]]]
  means.sample <- rowMeans(sample.rpm)
  cv.sample <- rowSds(sample.rpm) / means.sample
  
  cv.df <- rbind(cv.df,
                 data.frame(mean_rpm=means.sample,
                            cv=cv.sample,
                            gene_short_name=names(means.sample),
                            sample=rep(samples[i], length(means.sample))))
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
# Plot mean RPM distribution across different samples
png(glue("./data/GSE139944/transcriptional_noise/meanRPM/sciPlex3_meanRPM_custom_{cell.type}_{product}.png"),
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

# Plot coefficient of variation distribution across different samples
png(glue("./data/GSE139944/transcriptional_noise/CV/sciPlex3_CV_custom_{cell.type}_{product}.png"),
    width=3000, height=2000, res=300)
ggplot(cv.df, aes(x=sample, y=cv)) +
  geom_violin() +
  xlab("") +
  ylab("CV") +
  scale_x_discrete(labels=paste0(samples, " (", ncells, ")")) +
  ggtitle(glue("Coefficient of variation of genes across {product} treatments
               in {cell.type} cells"))
dev.off()

# Plot CV vs. CV of product_dose vs. vehicle
cv.cv.plot.list <- list()
cv.cv <- data.frame(vehicle_cv=numeric(),
                    treatment_cv=numeric(),
                    vehicle_rpm=numeric(),
                    treatment_rpm=numeric(),
                    sample=character(),
                    gene_short_name=character())

cv.diff.thresh <- 1.5
rpm.diff.thresh <- 0.2

for (i in 1:length(vehicles)) {
  cv.cv <- rbind(cv.cv,
                 data.frame(vehicle_cv=cv.df[cv.df$sample == vehicles[i], "cv"],
                            treatment_cv=cv.df[cv.df$sample == paste0(product, "_", d), "cv"],
                            vehicle_rpm=cv.df[cv.df$sample == vehicles[i], "mean_rpm"],
                            treatment_rpm=cv.df[cv.df$sample == paste0(product, "_", d), "mean_rpm"],
                            sample=cv.df[cv.df$sample == vehicles[i], "sample"],
                            gene_short_name=cv.df[cv.df$sample == vehicles[i], "gene_short_name"]))
}

# Remove genes with max CV
# max.thresh <- 0.1
# cv.cv <- cv.cv %>% filter(vehicle_cv < max(vehicle_cv, na.rm=TRUE) - max.thresh &
#                           treatment_cv < max(treatment_cv, na.rm=TRUE) - max.thresh)

cv.cv$cv_diff <- cv.cv$treatment_cv / cv.cv$vehicle_cv
cv.cv$cv_direction <- ifelse(cv.cv$treatment_cv > cv.cv$vehicle_cv, "UP", "DOWN")
cv.cv$rpm_diff <- (cv.cv$treatment_rpm - cv.cv$vehicle_rpm) / cv.cv$vehicle_rpm

cv.cv <- cv.cv %>% 
  mutate(label=ifelse((cv_diff > cv.diff.thresh | cv_diff < 1 / cv.diff.thresh) & abs(rpm_diff) < rpm.diff.thresh, 
                      gene_short_name, ""))
cv.cv <- drop_na(cv.cv, vehicle_cv, treatment_cv)

count <- 0

for (i in 1:length(vehicles)) {
  for (j in 1:length(d)) {
    cv.cv.plot <- cv.cv %>% filter(sample == vehicles[i]) %>%
      ggplot(aes(x=log2(vehicle_cv), y=log2(treatment_cv))) +
      geom_abline(intercept=0, 
                  slope=1,
                  linetype="dashed", size=0.5, color="gray") +
      geom_hline(yintercept=0) +
      geom_vline(xintercept=0) +
      geom_point(aes(fill=cv_direction),
                 colour="black", pch=21) +
      geom_label_repel(data=cv.cv %>% filter(sample == vehicles[i]),
                       aes(label=label),
                       size=3, segment.alpha=0.5, fill="white",
                       box.padding=0.25,
                       force=5,
                       max.overlaps=Inf) +
      xlab(glue("log2(CV ({vehicles[i]}))")) + 
      ylab(glue("log2(CV ({product}_{d}))")) +
      labs(title=glue("{product}_{d} (nCells={ncells[[paste0(cell.type, '_', product, '_', d)]]}) vs. {vehicles[i]} (nCells={ncells[[paste0(cell.type, '_', vehicles[i])]]})"),
           subtitle=glue("N={cv.cv %>% filter(sample == vehicles[i]) %>% drop_na(vehicle_cv, treatment_cv) %>% nrow()}, nUP={cv.cv %>% filter(sample == vehicles[i]) %>% drop_na(vehicle_cv, treatment_cv) %>% filter(cv_direction == 'UP') %>% nrow()}, nDOWN={cv.cv %>% filter(sample == vehicles[i]) %>% drop_na(vehicle_cv, treatment_cv) %>% filter(cv_direction == 'DOWN') %>% nrow()}
                       {cv.cv %>% filter(sample == vehicles[i]) %>% .$label %>% .[!is.na(.) & . != ''] %>% length()} genes with CV fold-change > {cv.diff.thresh} and RPM % difference < {rpm.diff.thresh * 100}%"),
           fill="CV direction")
    
    cv.cv.plot.list[[count]] <- cv.cv.plot
    count <- count + 1
  }
}
ggsave(glue("./data/GSE139944/transcriptional_noise/CVvsCV/sciPlex3_CVvsCV_custom_{cell.type}_{product}_{d}.png"),
       grid.arrange(grobs=cv.cv.plot.list, 
                    ncol=3, nrow=4,
                    top=textGrob(glue("Coefficient of variation (CV) of protein coding genes
                                      in {cell.type} cells for {product} vs. Vehicle (Nrep={nrep.vehicle})"),
                                 gp=gpar(fontsize=20, font=2))),
       width=30, height=90, units="cm")

# Output overlapping genes between vehicles for each direction
gene.overl.list <- list()

for (d in c("UP", "DOWN")) {
  gene.overl.list[[d]] <- sort(table(cv.cv %>% filter(cv_direction == d) %>% .$label),
                              decreasing=TRUE)[-1]
}


##### Distribution of RPMs ####
select.genes <- c( "TNIK", "RPL37A", "YES1", "APEX2", "PRNP")
density.plot.list <- list()

rpms.df <- data.frame(rpm=numeric(),
                      sample=character(),
                      gene_short_name=character())

for (i in 1:length(select.genes)) {
  for (j in 1:length(samples)) {
    select.rpms <- filt.rpms[select.genes[i], sample.cells.list[[samples[j]]]]
    rpms.df <- rbind(rpms.df,
                     data.frame(rpm=select.rpms,
                                sample=rep(samples[j], length(select.rpms)),
                                gene_short_name=rep(select.genes[i], length(select.rpms))))
  }
  
  density.plot <- rpms.df %>% filter(gene_short_name == select.genes[i]) %>% 
    ggplot(aes(x=log10(rpm + 1), color=factor(sample, levels=samples))) +
    geom_density() +
    xlab("log10(RPM + 1)") +
    ylab("Density") +
    scale_color_discrete(name="Samples",
                         labels=paste0(samples, " (", ncells, ")")) +
    labs(title=select.genes[i])
  density.plot.list[[i]] <- density.plot
}
ggsave(glue("./data/GSE139944/transcriptional_noise/RPMdist/sciPlex3_RPMdist_{cell.type}_{product}_{d}.png"),
       grid.arrange(grobs=density.plot.list, 
                    ncol=1, nrow=length(select.genes),
                    top=textGrob(glue("Distribution of RPM values for N={length(select.genes)} genes
                                      in {cell.type} cells for {product}_{d} vs. Vehicle (Nrep={nrep.vehicle})"),
                                 gp=gpar(fontsize=20, font=2))),
       width=30, 
       height=15 * length(select.genes), 
       units="cm")

#### UMAP ####
deseq.output <- read.table(glue("./data/GSE139944/deseq/DEGtable/sciPlex3_DESeq_{cell.type}_{product}_{d}_vs_Vehicle_0.txt"), 
                           header=TRUE, sep='\t', row.names=1, check.names=FALSE)
rownames(deseq.output) <- gsub("\\.[0-9_A-Z]+$", "", rownames(deseq.output))
proteo.deseq <- deseq.output[rownames(deseq.output) %in% proteo.genes &
                             deseq.output$log2FoldChange > 0, ]

umap.genes <- c("HSPA1A", "HSP90AB1", "HSP90AA1")
# c("HSPA1A", "HSP90AB1", "HSP90AA1")
# c("RPL37A", "RYBP", "TNIK", "VWDE")

umap.plot.list <- list()

for (i in 1:length(samples)) {
  umap.cds <- new_cell_data_set(expression_data=filt.counts[rownames(filt.gene.data), sample.cells.list[[samples[i]]]],
                                cell_metadata=filt.col.data[sample.cells.list[[samples[i]]], ],
                                gene_metadata=filt.gene.data) 
  
  umap.cds <- preprocess_cds(umap.cds,
                             method="PCA",
                             norm_method="log",
                             pseudo_count=1)
  
  # Reduce dimensions using UMAP
  umap.cds <- reduce_dimension(umap.cds,
                               reduction_method="UMAP",
                               preprocess_method="PCA")
  
  # Cluster cells
  umap.cds <- cluster_cells(umap.cds,
                            reduction_method="UMAP",
                            num_iter=100)
  
  umap.plot <- plot_cells(umap.cds,
                          genes=umap.genes,
                          group_cells_by="cluster",
                          show_trajectory_graph=FALSE,
                          label_cell_groups=FALSE,
                          cell_size=1) +
    ggtitle(glue("{samples[i]} (Ncells={ncells[i]})"))
  
  
  umap.plot.list[[i]] <- umap.plot
}
ggsave(glue("./data/GSE139944/UMAP/sciPlex3_UMAP_custom_{cell.type}_{product}_{d}_{paste(umap.genes, sep='-', collapse='-')}.png"),
       grid.arrange(grobs=umap.plot.list, 
                    ncol=2, nrow=2,
                    top=textGrob(glue("UMAP of protein-coding genes (N={nrow(umap.cds)})
                                      in {cell.type} cells for {product}_{d} vs. Vehicle (Nrep={length(vehicles)})
                                      labelled by {paste(umap.genes, sep=', ', collapse=', ')} expressions"),
                                 gp=gpar(fontsize=20, font=2))),
       width=40, 
       height=20, 
       units="cm")

#### RPM vs. RPM ####
x.genes <- c("HSP90AA1", "HSP90AA1", "HSP90AB1")
y.genes <- c("HSP90AB1", "HSPA1A", "HSPA1A")

rpm.rpm.plot.list <- list()

for (j in 1:length(x.genes)) {
  rpm.rpm.df <- data.frame(x.gene=numeric(),
                           y.gene=numeric(),
                           sample=character(),
                           row.names=character())
  
  for (i in 1:length(samples)) {
    cells <- sample.cells.list[[samples[i]]]
    
    rpm.rpm.df <- rbind(rpm.rpm.df,
                        data.frame(x.gene=filt.rpms[x.genes[j], cells],
                                   y.gene=filt.rpms[y.genes[j], cells],
                                   sample=rep(samples[i], length(cells)),
                                   row.names=cells))
  }
  
  rpm.rpm.plot <- rpm.rpm.df %>%
    ggplot(aes(x=log10(x.gene + 1), 
               y=log10(y.gene + 1), 
               color=factor(sample, levels=samples))) +
    geom_point(alpha=1) +
    scale_color_manual("Samples",
                       labels=paste0(samples, " (", ncells, ")"),
                       breaks=samples,
                       values=c(rep("gray", nrep.vehicle), "yellow", "orange", "red", "darkred")) +
    xlab(glue("{x.genes[j]} log10(RPM + 1)")) +
    ylab(glue("{y.genes[j]} log10(RPM + 1)")) +
    labs(title=glue("{y.genes[j]} vs. {x.genes[j]}"))
  
  rpm.rpm.plot.list[[j]] <- rpm.rpm.plot
}
ggsave(glue("./data/GSE139944/transcriptional_noise/RPMvsRPM/sciPlex3_RPMvsRPM_custom_{cell.type}_{product}.png"),
       grid.arrange(grobs=rpm.rpm.plot.list, 
                    ncol=1, nrow=length(x.genes),
                    top=textGrob(glue("RPM vs. RPM of N={length(x.genes)} pairwise gene combinations
                                      in {cell.type} cells for {product} vs. Vehicle (Nrep={length(vehicles)})"),
                                 gp=gpar(fontsize=20, font=2))),
       width=25, 
       height=20 * length(x.genes), 
       units="cm")

