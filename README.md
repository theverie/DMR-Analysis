# DMR-Analysis
DMR analysis of Alzheimer Related Genes (APOE, MAPT, APP, PSEN1, PSEN2, TREM2)


# Loading required libraries for DMR analysis
library(dplyr)
library(ggplot2)
library(data.table)
library(ChIPseeker)
library(GenomicRanges)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)
library(fgsea)  # Added for fgseaMultilevel

# Install and load KEGGREST if not already installed
> if (!requireNamespace("KEGGREST", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install("KEGGREST")
library(KEGGREST)

if (!requireNamespace("pathview", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  cat("Installing pathview package...\n")
  BiocManager::install("pathview", update = TRUE, ask = FALSE)
  if (!requireNamespace("pathview", quietly = TRUE)) {
    stop("Failed to install pathview package. Please install manually with BiocManager::install('pathview') and restart R.")
  }
}
library(pathview)

# Reading the CSV data
dmr_data <- fread("/Users/admin/Downloads/GSE244352_DifferentialMethylation.csv")

# Cleaning and preprocessing the data
dmr_data <- dmr_data %>%
  filter(!is.na(pvalue) & !is.na(qvalue) & !is.na(meth.diff)) %>%
  mutate(
    chr = as.factor(chr),
    strand = as.factor(strand),
    significant = qvalue < 0.05 & abs(meth.diff) > 15
  )
# Summarizing significant DMRs
summary_stats <- dmr_data %>%
  summarise(
    total_regions = n(),
    significant_regions = sum(significant),
    hypermethylated = sum(significant & meth.diff > 0),
    hypomethylated = sum(significant & meth.diff < 0),
    mean_meth_diff = mean(abs(meth.diff[significant]), na.rm = TRUE),
    median_pvalue = median(pvalue[significant], na.rm = TRUE)
  )

# Printing summary statistics
cat("Summary of DMR Analysis:\n")
print(summary_stats)


# Creating a Manhattan plot for visualization
dmr_data$significance <- ifelse(dmr_data$qvalue < 0.05 & abs(dmr_data$meth.diff) > 15,
                                "Significant", "Not Significant") #additional code to set meth.diff to 15%
manhattan_plot <- ggplot(dmr_data, aes(x = start, y = -log10(pvalue), color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  facet_wrap(~ chr, scales = "free_x") +
  theme_minimal() +
  labs(
    title = "Manhattan Plot of DMRs in Alzheimer's Disease",
    x = "Genomic Position",
    y = "-log10(p-value)",
    color = "Significance (q < 0.05 & |meth.diff| > 15)"
  ) +
  scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
  theme(legend.position = "top")

# Creating volcano plot for visualization
dmr_data$significance <- ifelse(dmr_data$qvalue < 0.05 & abs(dmr_data$meth.diff) > 15,
                                "Significant", "Not Significant") #additional code to set meth.diff to 15%
volcano_plot <- ggplot(dmr_data, aes(x = meth.diff, y = -log10(pvalue), color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  theme_minimal() +
  labs(
    title = "Volcano Plot of DMRs in Alzheimer's Disease",
    x = "Methylation Difference (%)",
    y = "-log10(p-value)",
    color = "Significance (q < 0.05 & |meth.diff| > 15)"
  ) +
  scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
  geom_vline(xintercept = c(-15, 15), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  theme(legend.position = "top")

# Visualization of the plots in R
print(volcano_plot)
print(manhattan_plot)

# Filter significant DMRs for annotation
significant_dmrs <- dmr_data %>%
  filter(significant) %>%
  select(chr, start, end, strand, pvalue, qvalue, meth.diff)

# Create GRanges object for all significant DMRs
dmr_granges <- GRanges(
  seqnames = significant_dmrs$chr,
  ranges = IRanges(start = significant_dmrs$start, end = significant_dmrs$end),
  strand = significant_dmrs$strand,
  pvalue = significant_dmrs$pvalue,
  qvalue = significant_dmrs$qvalue,
  meth.diff = significant_dmrs$meth.diff
)

# Load TxDb object for gene annotation
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Perform gene annotation on all significant DMRs
annotation <- annotatePeak(dmr_granges, TxDb = txdb, annoDb = "org.Hs.eg.db")
annotation_df <- as.data.frame(annotation)

# Combine with original DMR data
annotated_dmrs <- data.frame(
  chr = significant_dmrs$chr,
  start = significant_dmrs$start,
  end = significant_dmrs$end,
  strand = significant_dmrs$strand,
  pvalue = significant_dmrs$pvalue,
  qvalue = significant_dmrs$qvalue,
  meth.diff = significant_dmrs$meth.diff,
  annotation = annotation_df$annotation,
  gene_id = annotation_df$geneId,
  gene_symbol = annotation_df$SYMBOL,
  distance_to_tss = annotation_df$distanceToTSS
)
fwrite(annotated_dmrs, "all_annotated_dmrs.csv")

# Filter DMRs for specific AD-related genes
ad_genes <- c("APP", "PSEN1", "PSEN2", "APOE", "MAPT", "TREM2")
filtered_dmrs <- annotated_dmrs[annotated_dmrs$gene_symbol %in% ad_genes | is.na(annotated_dmrs$gene_symbol), ]
fwrite(filtered_dmrs, "filtered_ad_dmrs.csv")
print("Filtered DMRs for AD-related genes:")
print(filtered_dmrs)  # Print the filtered data

# Create GRanges object for all significant DMRs
dmr_granges <- GRanges(
  seqnames = significant_dmrs$chr,
  ranges = IRanges(start = significant_dmrs$start, end = significant_dmrs$end),
  strand = significant_dmrs$strand,
  pvalue = significant_dmrs$pvalue,
  qvalue = significant_dmrs$qvalue,
  meth.diff = significant_dmrs$meth.diff
)

# Load TxDb object for gene annotation
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Perform gene annotation on all significant DMRs
annotation <- annotatePeak(dmr_granges, TxDb = txdb, annoDb = "org.Hs.eg.db")
annotation_df <- as.data.frame(annotation)

# Combine with original DMR data
annotated_dmrs <- data.frame(
  chr = significant_dmrs$chr,
  start = significant_dmrs$start,
  end = significant_dmrs$end,
  strand = significant_dmrs$strand,
  pvalue = significant_dmrs$pvalue,
  qvalue = significant_dmrs$qvalue,
  meth.diff = significant_dmrs$meth.diff,
  annotation = annotation_df$annotation,
  gene_id = annotation_df$geneId,
  gene_symbol = annotation_df$SYMBOL,
  distance_to_tss = annotation_df$distanceToTSS
)
fwrite(annotated_dmrs, "all_annotated_dmrs.csv")

# Filter DMRs for specific AD-related genes
ad_genes <- c("APP", "PSEN1", "PSEN2", "APOE", "MAPT", "TREM2")
filtered_dmrs <- annotated_dmrs[annotated_dmrs$gene_symbol %in% ad_genes | is.na(annotated_dmrs$gene_symbol), ]
fwrite(filtered_dmrs, "filtered_ad_dmrs.csv")

# Print the Filtered DMRs for AD-related genes
print("Filtered DMRs for AD-related genes:")

# Print the filtered data
print(filtered_dmrs)  

# Prepare gene list for GSEA
gene_list <- filtered_dmrs$meth.diff
names(gene_list) <- filtered_dmrs$gene_symbol

# Remove duplicates in gene names
gene_list <- gene_list[!duplicated(names(gene_list))]  # Remove duplicate gene symbols
gene_list <- gene_list[!is.na(names(gene_list))]  # Remove NA names
gene_list <- sort(gene_list, decreasing = TRUE)
entrez_ids <- mapIds(org.Hs.eg.db, keys = names(gene_list), keytype = "SYMBOL", column = "ENTREZID")
gene_list_entrez <- gene_list[!is.na(entrez_ids)]
names(gene_list_entrez) <- entrez_ids[!is.na(entrez_ids)]

# Perform KEGG pathway enrichment
> if (length(gene_list_entrez) >= 1) {
     background_genes <- keys(org.Hs.eg.db, keytype = "ENTREZID")
     kegg_result <- enrichKEGG(gene = names(gene_list_entrez), universe = background_genes, organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500)
     ad_kegg_result <- kegg_result[kegg_result$ID == "hsa05010" | grepl("hsa050", kegg_result$ID), ]
     if (nrow(ad_kegg_result) == 0) {
       cat("No AD-specific pathways found. Showing all significant results.\n")
       ad_kegg_result <- kegg_result
       }
     print(ad_kegg_result)
     save_path <- "C:/Users/YourUsername/Documents/AD_Analysis/Results/"
     file_name <- "kegg_all_ad_results.csv"
     if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
     full_path <- file.path(save_path, file_name)
     write.csv(as.data.frame(ad_kegg_result), full_path, row.names = FALSE)
     cat("\nKEGG results saved to:", full_path, "\n")
     if (nrow(ad_kegg_result) > 0) {
         dotplot(ad_kegg_result, showCategory = 10, title = "AD-Related KEGG Pathways (All Significant DMRs)")
     } else {
         cat("No significant KEGG results to plot.\n")
     }
 } else {
     cat("No valid genes for KEGG analysis.\n")
 }
# Create KEGG pathway plot using pathview for hsa05010
if (length(gene_list_entrez) >= 1) {
    # Prepare data for pathview (gene ID to meth.diff mapping)
    pv_data <- gene_list_entrez
    names(pv_data) <- names(gene_list_entrez)  # Ensure Entrez IDs as names
    
    # Specify the KEGG pathway (e.g., hsa05010 for Alzheimer's disease)
    pathway_id <- "hsa05010"
    
    # Generate pathway plot, saving to Downloads
    pv.out <- pathview(
        gene.data = pv_data,
        pathway.id = pathway_id,
        species = "hsa",
        out.suffix = "AD_DMRs",
        kegg.native = TRUE,
        limit = list(gene = max(abs(pv_data), na.rm = TRUE)),
        kegg.dir = "/path"
    )
    
    cat("\nKEGG pathway plot saved to Downloads as:", 
        file.path("/path", paste0(pathway_id, ".AD_DMRs.png")), "\n")
} else {
    cat("No valid genes for KEGG pathway plot.\n")
}

# Create KEGG pathway plot using pathview for hsa05022
if (length(gene_list_entrez) >= 1) {
    # Prepare data for pathview (gene ID to meth.diff mapping)
    pv_data <- gene_list_entrez
    names(pv_data) <- names(gene_list_entrez)  # Ensure Entrez IDs as names
    
    # Specify the KEGG pathway (hsa05022 for multiple neurodegeneration diseases)
    pathway_id <- "hsa05022"
    
    # Generate pathway plot, saving to Downloads
    pv.out <- pathview(
        gene.data = pv_data,
        pathway.id = pathway_id,
        species = "hsa",
        out.suffix = "AD_DMRs",
        kegg.native = TRUE,
        limit = list(gene = max(abs(pv_data), na.rm = TRUE)),
        kegg.dir = "/Users/admin/Downloads"
    )
    
    cat("\nKEGG pathway plot saved to Downloads as:", 
        file.path("/Users/admin/Downloads", paste0(pathway_id, ".AD_DMRs.png")), "\n")
} else {
    cat("No valid genes for KEGG pathway plot.\n")
}


