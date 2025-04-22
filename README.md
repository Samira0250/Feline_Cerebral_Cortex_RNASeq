# Feline_Cerebral_Cortex_RNASeq

🧠 Feline_Cerebral_Cortex_RNASeq
This repository contains the full pipeline for analyzing RNA-seq data from the feline cerebral cortex, with a focus on:

Differential Gene Expression (DGE) using DESeq2

Differential Exon Usage (DEU) using DEXSeq

The analysis is performed on mapped and counted RNA-seq reads, and includes downstream steps like normalization, visualization, gene ranking for GSEA, and exon-level quantification.

📁 Repository Structure
graphql
Copy
Edit
Feline_Cerebral_Cortex_RNASeq/
├── DESeq2_analysis.R           # Main DESeq2 script for DGE
├── DEXSeq_analysis.R           # DEXSeq pipeline for differential exon usage
├── submit_deseq2.pbs           # PBS job script for running DESeq2 on ASC
├── data/                       # Count matrices and metadata
├── results/                    # Output files: CSVs, plots, heatmaps
├── gsea/                       # Ranked list and files for GSEA PreRanked analysis
└── README.md                   # This file
🛠️ Requirements
The pipeline is designed to run on a high-performance computing environment (Alabama Supercomputer). Make sure the following R packages are installed:

DESeq2

DEXSeq

pheatmap

RColorBrewer

ggplot2

biomaRt (optional: for cross-species gene mapping)


📊 Outputs
Normalized gene expression tables

DE results with adjusted p-values and log2 fold changes

MA plots, PCA plots, and heatmaps of variable genes

Ranked gene list for GSEA

Differential exon usage results from DEXSeq

🧬 Project Goal
The goal of this project is to characterize transcriptional changes in the feline cerebral cortex, particularly in the context of neurodevelopmental disorders and genotype effects (PEA15 knockout and heterozygote).
