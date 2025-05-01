# === Load libraries ===
library(dplyr)
library(ggplot2)

# === Set working directory ===
setwd("C:/Users/szs0394/Downloads/DESeq2/final/WT vs KO analyses")

# === Read the GSEA result files ===
gsea_pos <- read.delim("gsea_report_for_na_pos_1745614753040.tsv", header = TRUE)
gsea_neg <- read.delim("gsea_report_for_na_neg_1745614753040.tsv", header = TRUE)

# === Combine positive and negative results ===
gsea_combined <- bind_rows(gsea_pos, gsea_neg)

# === Create 'Direction' column (Upregulated vs Downregulated) ===
gsea_combined$Direction <- ifelse(gsea_combined$NES > 0, "Upregulated", "Downregulated")

# === Clean the pathway names (remove KEGG_ or HALLMARK_ prefixes) ===
gsea_combined$CleanName <- gsub("^HALLMARK_|^KEGG_", "", gsea_combined$NAME)

# === Sort by FDR q-value (most significant first) ===
gsea_combined <- gsea_combined %>% arrange(FDR.q.val)

# === Select Top 20 most significant pathways ===
gsea_top20 <- gsea_combined %>% head(20)

# === Plot ===
gsea_plot <- ggplot(gsea_top20, aes(x = reorder(CleanName, NES), y = NES, fill = Direction)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("Downregulated" = "blue", "Upregulated" = "red")) +
  labs(
    x = NULL,
    y = "Normalized Enrichment Score (NES)",
    title = "Top 20 Enriched Pathways (WT vs KO)"
  ) +
  theme_minimal(base_size = 14) +  # Increase overall font size
  theme(
    legend.position = "top",
    axis.text.y = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# === Save the plot ===
ggsave("GSEA_Top20_Barplot_WT_vs_KO.png", plot = gsea_plot, width = 10, height = 8, dpi = 300)

# === Print it to RStudio Viewer ===
print(gsea_plot)
