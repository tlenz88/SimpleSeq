#!/usr/bin/env Rscript

## Created: October 5, 2023
## Updated: March 8, 2024
## Author(s): Todd Lenz, tlenz001@ucr.edu

## Performs differential expression analysis via DESeq2.

args <- commandArgs(trailingOnly = TRUE)

ShowHelp <- function() {
  cat("differential_expression.R --help\n")
  cat("Usage: differential_expression.R -r READ_COUNTS -m METADATA -c CONTROL -q QVALUE\n")
  cat("---------------------------------------------------------------------\n")
  cat("  Input arguments:\n")
  cat("    -r|--read_counts READ_COUNTS   : Read count file.\n")
  cat("    -m|--metadata METADATA         : Metadata file.\n")
  cat("    -c|--control CONTROL           : String indicating control group.\n")
  cat("    -q|--qvalue QVALUE             : Float indicating q-value cutoff.\n")
  cat("    -h|--help HELP                 : Show help message.\n")
  cat("---------------------------------------------------------------------\n")
}

if ("--help" %in% args || "-h" %in% commandArgs() || length(commandArgs(trailingOnly = TRUE)) == 0) {
  ShowHelp()
  quit()
}


###########################################
## Load packages and install if missing. ##
###########################################
load_packages <- function(package_list) {
  for (pkg in package_list) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      eval(bquote(install.packages(.(pkg), dependencies = TRUE)))
    }
    suppressMessages(eval(bquote(library(.(pkg), quietly = TRUE))))
  }
}
load_packages(c("DESeq2", "ggplot2", "pheatmap", 
  "RColorBrewer", "apeglm", "dplyr", "EnhancedVolcano"))


############################
## Check input arguments. ##
############################
for (i in c(1, 3, 5, 7)) {
  if (!(args[i] %in% c("--read_counts", "-r", "--metadata", "-m", "--control", "-c", "--qval", "-q"))) {
    cat(paste("Error: ", args[i], "is an invalid argument.\n"))
    quit()
  } else if (args[i] == "--read_counts" || args[i] == "-r") {
    if (file.exists(args[i + 1])) {
      read_counts <- args[i + 1]
    } else {
      cat("Error: Invalid input for counts argument.\n")
    }
  } else if (args[i] == "--metadata" || args[i] == "-m") {
    if (file.exists(args[i + 1])) {
      metadata <- args[i + 1]
    } else {
      cat("Error: Invalid input for metadata argument.\n")
    }
  } else if (args[i] == "--control" || args[i] == "-c") {
    control <- args[i + 1]
  } else if (args[i] == "--qval" || args[i] == "-q") {
    qval <- as.numeric(args[i + 1])
  }
}


#################################################
## Load HTSeq read counts and sample metadata. ##
#################################################
read_counts <- normalizePath(read_counts)
countData <- read.csv(read_counts, header = TRUE, sep = "\t")
metaData <- read.csv(metadata, header = TRUE, sep = "\t")
if (!(control %in% metaData$group)) {
  cat(paste("Error: ", control, " not in 'Group' column of metadata\n"))
}


###############################################
## Perform differential expression analysis. ##
###############################################
findDiff <- function(dds, treatment, control, qval) {
  res <- results(dds, contrast = c("group", treatment, control))
  resSig <- res[which(res$padj < qval), ]
  deg <- resSig[order(resSig$log2FoldChange, decreasing = TRUE), ]
  deg <- cbind(rownames(deg), deg)
  colnames(deg)[colnames(deg) == "rownames(deg)"] <- "geneID"
  outfile <- file.path(dirname(read_counts), paste0(treatment, "_DEgenes.txt"))
  write.table(deg, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)
}

dds <- suppressWarnings(DESeqDataSetFromMatrix(countData, metaData, ~group, tidy = TRUE))
dds$group <- relevel(dds$group, ref = control)
dds <- DESeq(dds)

for (treatment in unique(metaData$group[!grepl(control, metaData$group)])) {
  findDiff(dds, treatment, control, qval)
}


###################################
## Plot output from DE analysis. ##
###################################
pdf(file.path(dirname(read_counts), "RNAseq_figures.pdf"))
color_map <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
  "#0072B2", "#D55E00", "#CC79A7", "#000000", "#999999"
)

# PCA plot
rld <- rlog(dds, blind = TRUE)
pca_rld <- plotPCA(rld, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pca_rld, "percentVar"))
ggplot(pca_rld, aes(PC1, PC2, color = group)) +
  labs(fill = "group") +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme(text = element_text(family = "serif", size = 8)) +
  scale_color_manual(values = color_map[1:length(unique(metaData$group))]) +
  coord_fixed()

# Sample correlation heatmap
normalized_counts <- counts(dds, normalized = TRUE)
spearman_corr <- cor(normalized_counts, method = "spearman")
colors <- colorRampPalette(brewer.pal(9, "Reds"))(255)
pheatmap(
  spearman_corr,
  display_numbers = TRUE,
  col = colors,
  number_color = "#000000",
  fontsize_number = 8,
  treeheight_row = 0,
  fontsize = 8
)

# MA plot and volcano plot
for (treatment in unique(metaData$group[!grepl(control, metaData$group)])) {
  resLFCapeglm <- lfcShrink(
    dds,
    coef = paste0("group_", treatment, "_vs_", control),
    type = "apeglm"
  )
  resLFCapeglm <- resLFCapeglm[is.na(resLFCapeglm$log2FoldChange) == "FALSE", ]
  resLFCapeglm <- resLFCapeglm[is.na(resLFCapeglm$padj) == "FALSE", ]
  min_log2FC <- min(resLFCapeglm$log2FoldChange, na.rm = TRUE)
  max_log2FC <- max(resLFCapeglm$log2FoldChange, na.rm = TRUE)
  min_padj <- -log10(max(resLFCapeglm$padj, na.rm = TRUE))
  max_padj <- -log10(min(resLFCapeglm$padj, na.rm = TRUE))
  plotMA(
    resLFCapeglm,
    alpha = qval,
    main = paste(treatment, " vs ", control),
    ylim = c(min_log2FC, max_log2FC),
    xaxt = "n",
    cex = 0.5,
    cex.lab = 1.2,
    colLine = "#000000",
    colNonSig = "#999999",
    colSig = "#000000"
  )
  MA_xlim <- c(1, 10, 100, 1000, 10000, 100000, 1000000)
  axis(1, at = MA_xlim, labels = MA_xlim)
  keyvals <- ifelse(between(resLFCapeglm$log2FoldChange, -1, 1), "#999999",
    ifelse(resLFCapeglm$padj > qval, "#999999", "#000000"))
  names(keyvals)[keyvals == "#000000"] <- "pass log2FC and p-value"
  names(keyvals)[keyvals == "#999999"] <- "fail log2FC or p-value"
  plot(EnhancedVolcano(
    resLFCapeglm,
    lab = "",
    title = paste(treatment, "vs", control),
    caption = paste(nrow(resLFCapeglm), "variables, padj <=", qval, "Log2FC >= 2"),
    colCustom = keyvals,
    pointSize = 1,
    subtitle = NULL,
    captionLabSize = 8,
    legendLabSize = 8,
    legendIconSize = 1.5,
    x = "log2FoldChange",
    y = "padj",
    xlim = c(min_log2FC, max_log2FC),
    ylim = c(min_padj, max_padj),
    axisLabSize = 8,
    pCutoff = qval,
    FCcutoff = 1.0,
    colAlpha = 1,
    border = "partial",
    borderWidth = 1.5,
    borderColour = "#000000"
  ))
}

# Heatmap
rld <- rlog(dds, blind = FALSE)
annoDF <- data.frame(colData(rld)[, "group"])
colnames(annoDF) <- "group"
rownames(annoDF) <- colnames(rld)
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 300)
sigMat <- assay(rld)[topVarGenes, ]
sigMat <- sigMat[, order(colnames(sigMat))]
ann_colors = list(group = color_map[1:length(unique(metaData$group))])
names(ann_colors$group) <- unique(metaData$group)
pheatmap(
  sigMat,
  main = "Top 300 variable genes",
  scale = "row",
  show_rownames = FALSE,
  show_colnames = TRUE,
  treeheight_col = 2,
  treeheight_row = 0,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annoDF,
  annotation_names_col = FALSE,
  annotation_colors = ann_colors,
  fontsize = 8,
  cellwidth = 18,
  cellheight = 1,
  border_color = NA
)

dev.off()

#var_genes <- read.csv("F:/organism_genome/Pfalciparum3D7/var.gff", sep = "\t", header = FALSE)
#var_genes <- subset(var_genes, V3 == "protein_coding_gene")
#var_genes <- var_genes["V9"]
#var_gene_list <- str_extract(var_genes$V9, "(?<=ID=)[^;]+")
#gene_colors <- ifelse(rownames(resLFCapeglm) %in% var_gene_list, "#FF0000", "#000000")
#plotMA(
#    resLFCapeglm,
#    alpha = qval,
#    main = paste(treatment, " vs ", control),
#    ylim = c(min_log2FC, max_log2FC),
#    xaxt = "n",
#    cex = 0.5,
#    cex.lab = 1.2,
#    colLine = "#000000",
#    colNonSig = "#999999",
#    col = gene_colors
#)
