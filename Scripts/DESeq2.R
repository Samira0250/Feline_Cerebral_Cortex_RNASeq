# === Load required libraries ===
library(DESeq2)
library(ggplot2)
library(ggrepel)  # For labeling top genes nicely

# === Set working directory ===
setwd("C:/Users/.../Downloads/DESeq2/final/WT vs KO analyses") #Customize to your own directory

# === Load count matrix ===
countData <- read.csv("gene_count_matrix.csv", row.names = 1, check.names = FALSE)
countData <- as.matrix(countData)
mode(countData) <- "numeric"

# === Define sample metadata ===
sample_ids <- c(
  "SRR8092570", "SRR8092571", "SRR8092574", "SRR8092575",  # WT
  "SRR8092587", "SRR8092589", "SRR8092590", "SRR8092591", "SRR8092592",  # CR
  "SRR8092566", "SRR8092567", "SRR8092569", "SRR8092572", "SRR8092573", "SRR8092588"   # KO
)
genotypes <- c(rep("WT", 4), rep("CR", 5), rep("KO", 6))

# Save metadata
ordered_metadata <- data.frame(sampleID = sample_ids, genotypes = genotypes)
write.csv(ordered_metadata, "ordered_sample_metadata.csv", row.names = FALSE)

# === Load and align metadata ===
colData <- read.csv("ordered_sample_metadata.csv", row.names = 1)
countData <- countData[, colnames(countData) != "SRR8092586"]
colData <- colData[colnames(countData), , drop = FALSE]
colnames(colData)[1] <- "genotypes"
colData$genotypes <- factor(colData$genotypes, levels = c("WT", "CR", "KO"))
stopifnot(all(colnames(countData) == rownames(colData)))

# === Create DESeq2 dataset and run DESeq ===
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ genotypes)
dds <- dds[rowSums(counts(dds)) > 100, ]
dds <- DESeq(dds)
res <- results(dds)

# === Save DE results ===
resOrdered <- res[order(res$padj), ]
sigGenes_0.05 <- subset(resOrdered, padj < 0.05)
write.csv(sigGenes_0.05, "DEG_results_padj_0.05.csv")

# === Volcano Plot with Top 10 Genes Labeled ===
volcano_df <- as.data.frame(resOrdered)
volcano_df$significant <- "Not Significant"
volcano_df$significant[volcano_df$padj < 0.05 & volcano_df$log2FoldChange > 0] <- "Upregulated"
volcano_df$significant[volcano_df$padj < 0.05 & volcano_df$log2FoldChange < 0] <- "Downregulated"

# Label top 10 genes by smallest padj
top10_genes <- rownames(volcano_df)[1:10]

ggplot(volcano_df, aes(x = log2FoldChange, y = -log10(pvalue), color = significant,
                       label = ifelse(rownames(volcano_df) %in% top10_genes, rownames(volcano_df), ""))) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(title = "Volcano Plot: WT vs KO (Top Genes Labeled)",
       x = "Log2 Fold Change",
       y = "-Log10(p-value)") +
  theme(legend.position = "top")

ggsave("Volcano_Plot_WT_vs_KO.png", width = 10, height = 8, dpi = 300)

# === Top Gene Expression - WT vs KO Only ===
# Assuming topGene_df is already created from:
# topGene_df <- plotCounts(dds, gene = "PEA15", intgroup = "genotypes", returnData = TRUE)

# First, add sample IDs if not done yet
topGene_df$sampleID <- rownames(topGene_df)

# Now filter only WT and KO samples (exclude CR)
topGene_WT_KO <- subset(topGene_df, genotypes %in% c("WT", "KO"))

# Now re-plot
library(ggplot2)

ggplot(topGene_WT_KO, aes(x = genotypes, y = count, label = sampleID, color = genotypes)) +
  geom_point(size = 3) +
  geom_text(vjust = -1.2, size = 3) +
  scale_color_manual(values = c(WT = "green", KO = "red")) +  # No CR color anymore
  scale_y_continuous() +
  labs(
    title = "Normalized Expression of PEA15 by Genotype (WT vs KO)",
    x = "Genotype",
    y = "Normalized Count"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Save the plot
ggsave("PEA15_expression_WT_vs_KO_only.png", width = 10, height = 6, dpi = 300)

