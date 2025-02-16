## R for Metagenomics: ZINB-WaVE and DESeq2 Analysis

# Set working directory
setwd("~/Microbiome Analysis/Metagenomics")

# Load necessary libraries
install.packages("beanplot")
install.packages("RColorBrewer")
install.packages("Seurat")
install.packages("zinbwave")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")  # Install BiocManager if not already installed
BiocManager::install("DESeq2")       # Install DESeq2

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")  # Install BiocManager if not already installed
BiocManager::install("zinbwave")    # Install zinbwave

library("beanplot")
library("RColorBrewer")
library("Seurat")
library("DESeq2")
library("zinbwave")

# Load and Preprocess Count Data
# Load count data
counts_raw <- read.table('kaiju_species.tsv', header=TRUE, sep='\t', quote="")

# Clean up sample names
counts_raw[,1] <- gsub('_kaiju.*', '', gsub('.*/', '', counts_raw[,1]))

# Transform data from long to wide format
counts <- reshape(counts_raw, idvar='taxon_id', timevar='file', v.names='reads', direction='wide', drop=c('percent', 'taxon_name'))

# Replace NAs with 0
counts[is.na(counts)] <- 0

# Set row names and remove the taxon_id column
rownames(counts) <- paste0('taxon_id_', counts[,1])
counts <- counts[,-1]

# Clean up column names
colnames(counts) <- sub('reads.', '', colnames(counts))

# Convert to matrix
counts <- as.matrix(counts)

# Filter out extremely rare taxa
counts <- counts[apply(counts, 1, function(x) sum(x > 5) > 5),]

# Load and Preprocess meta Data
# Load metadata
meta <- read.table('meta_data_sra_acc.tsv', header=TRUE, sep='\t')

# Set row names
rownames(meta) <- meta$Sample

# Ensure metadata is in the same order as counts
meta <- meta[colnames(counts),]

# Convert variables to appropriate types
meta$Gender <- as.factor(meta$Gender)
meta$Group <- relevel(as.factor(meta$Group), 'Term')
meta$Individual <- as.factor(meta$Individual)
meta$DeliveryMode <- as.factor(meta$DeliveryMode)
meta$EnteralFeeds_2mo <- as.factor(meta$EnteralFeeds_2mo)

# Clean up levels of 'EnteralFeeds_2mo'
levels(meta$EnteralFeeds_2mo)[gsub(' ', '', levels(meta$EnteralFeeds_2mo)) == 'Breastmilk'] <- 'Breastmilk'

# ZINB-WaVE Analysis
# Create a SummarizedExperiment object
zi <- SummarizedExperiment(assays = SimpleList(counts = counts), colData = meta)

# Run ZINB-WaVE
zi <- zinbwave(zi, K = 5, X = ~ Gender + Group + Individual, observationalWeights = TRUE, epsilon=1e12)

# Extract latent variables
W <- reducedDim(zi)

# Plot latent variables
plot(W[,1:2], col=c('red','blue')[meta$DeliveryMode])

# Plot with more complex coloring
factorcolors <- RColorBrewer::brewer.pal(5, "Set1")
plot(W, col = factorcolors[meta$EnteralFeeds_2mo])

# Clustering with Seurat
# Convert to Seurat object
seu <- as.Seurat(x=zi, counts = 'counts', data = 'counts')

# Run UMAP
seu <- RunUMAP(seu, dims = 1:ncol(W), reduction = 'zinbwave', n.components = 2)

# Find neighbors and clusters
seu <- FindNeighbors(seu, reduction = "umap", dims = 1:2)
seu <- FindClusters(object = seu)

# Extract clusters
clusts <- as.factor(seu$seurat_clusters)

# Plot with cluster coloring
factorcolors <- RColorBrewer::brewer.pal(9, "Set1")
plot(W[,1:2], col=factorcolors[clusts])

# Add clusters to metadata
colData(zi) <- cbind(colData(zi), clusts)

# Differential Abundance Analysis with DESeq2
# Create DESeqDataSet object
deseq_object <- DESeqDataSet(zi, design = ~ Gender + Group + clusts)

# Fit DESeq2 models
res_group <- DESeq(deseq_object, test='LRT', sfType='poscounts', reduced = ~ Gender + clusts, useT=TRUE, minmu=1e-6, minReplicatesForReplace=Inf)

# Extract results
ps_group <- results(res_group)
ps_group <- ps_group[order(ps_group[,'padj']),]

# Plot top differential taxon
topdiff <- rownames(ps_group)[[1]]
relabund <- apply(counts, 2, function(x) x/sum(x))
merged_rel <- cbind(colData(zi), t(relabund))

# Boxplot
boxplot(as.formula(paste0(topdiff, ' ~ Group')), merged_rel)

# Log-transformed boxplot
relabund_0na <- relabund
relabund_0na[relabund_0na == 0] <- NA
merged_rel_0na <- cbind(meta, t(relabund_0na))
boxplot(as.formula(paste0(topdiff, ' ~ Group')), merged_rel_0na, log='y')

# Beanplot
beanplot::beanplot(as.formula(paste0(topdiff, ' ~ Group')), merged_rel_0na, log='y', las=2, cex.axis=0.5)

# Interaction plot
beanplot::beanplot(as.formula(paste0(topdiff, ' ~ Gender + Group')), merged_rel_0na, log='y', las=2, cex.axis=0.5)

# Identify the top differential taxon
counts_raw$taxon_name[counts_raw$taxon_id == sub('taxon_id_', '', topdiff)]

# Control for Individual Effects
# Re-code Individual as a factor nested within Group
ind_group <- rep('', nrow(colData(zi)))
for (g in unique(colData(zi)$Group)) {
  i <- 1
  for (h in unique(colData(zi)$Individual[colData(zi)$Group == g])) {
    ind_group[colData(zi)$Group == g & colData(zi)$Individual == h] <- i
    i <- i + 1
  }
}
ind_group <- as.factor(ind_group)

# Add to metadata
colData(zi) <- cbind(colData(zi), ind_group)

# Create custom model matrix
mm <- model.matrix(~ Group + Group:ind_group, data=colData(zi))
mm <- cbind(1, mm[, apply(mm, 2, function(x) sd(x) > 0)])

# Create DESeq object
deseq_object <- DESeqDataSet(zi, design = mm)

# Fit DESeq model
res_group_cont <- DESeq(deseq_object, sfType='poscounts', useT=TRUE, minmu=1e-6, minReplicatesForReplace=Inf)

# Shrink estimates
ps_group_cont <- lfcShrink(res_group_cont, coef="GroupEarly.only.abx", type = "normal")

# Order results
ps <- ps_group_cont[order(ps_group_cont[,'padj']),]

# Plot top differential taxon
topdiff <- rownames(ps)[[1]]
beanplot::beanplot(as.formula(paste0(topdiff, ' ~ Group + Gender')), merged_rel_0na, log='y', las=2, cex.axis=0.5)

# Identify the top differential taxon
counts_raw$taxon_name[counts_raw$taxon_id == sub('taxon_id_', '', topdiff)]

# # Re-code Individual as a factor nested within Group
ind_group <- rep('', nrow(colData(zi)))
for (g in unique(colData(zi)$Group)) {
  i <- 1
  for (h in unique(colData(zi)$Individual[colData(zi)$Group == g])) {
    ind_group[colData(zi)$Group == g & colData(zi)$Individual == h] <- i
    i <- i + 1
  }
}
ind_group <- as.factor(ind_group)

# Add to metadata
colData(zi) <- cbind(colData(zi), ind_group)

# Create custom model matrix
mm <- model.matrix(~ Group + Group:ind_group, data=colData(zi))
mm <- cbind(1, mm[, apply(mm, 2, function(x) sd(x) > 0)])

# Create DESeq object
library(DESeq2)
deseq_object <- DESeqDataSet(zi, design = mm)

# Fit DESeq model
res_group_cont <- DESeq(deseq_object, sfType='poscounts', useT=TRUE, minmu=1e-6, minReplicatesForReplace=Inf)

# Shrink estimates
ps_group_cont <- lfcShrink(res_group_cont, coef="GroupEarly.only.abx", type = "normal")

# Order results
ps <- ps_group_cont[order(ps_group_cont[,'padj']),]

# Plot top differential taxon
topdiff <- rownames(ps)[[1]]
beanplot::beanplot(as.formula(paste0(topdiff, ' ~ Group + Gender')), merged_rel_0na, log='y', las=2, cex.axis=0.5)

# Identify the top differential taxon
counts_raw$taxon_name[counts_raw$taxon_id == sub('taxon_id_', '', topdiff)]

# Save Session
# Create the 'outputs' directory if it doesn't exist
if (!dir.exists("outputs")) {
  dir.create("outputs")
}

# Save the session
save.image(file = 'outputs/res_deseq_tax_20210726.RData')
