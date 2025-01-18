#!/usr/bin/env Rscript

## Created: October 5, 2023
## Updated: August 7, 2024
## Author(s): Todd Lenz, tlenz001@ucr.edu

## Performs differential expression analysis via DESeq2.

args <- commandArgs(trailingOnly = TRUE)

ShowHelp <- function() {
  cat("differential_expression.R --help\n")
  cat("Usage: differential_expression.R -r READ_COUNTS -m METADATA -c CONTROL -g GENE_LIST -q QVALUE [OPTIONS]\n")
  cat("---------------------------------------------------------------------\n")
  cat("  Input arguments:\n")
  cat("    -r|--read_counts READ_COUNTS : Read count file.\n")
  cat("    -m|--metadata METADATA       : Metadata file.\n")
  cat("    -c|--control CONTROL         : String indicating control group.\n")
  cat("    -q|--qvalue QVALUE           : Float indicating q-value cutoff.\n")
  cat("    -g|--gene_list GENE_LIST     : List of genes for plot annotation.\n")
  cat("    -h|--help HELP               : Show help message.\n")
  cat("  Optional arguments:\n")
  cat("    --group_column GROUP_COL     : Column name for group in metadata (default: 'group').\n")
  cat("    --top_var_genes TOP_N        : Number of top variable genes for heatmap (default: 300).\n")
  cat("---------------------------------------------------------------------\n")
}

if ("--help" %in% args || "-h" %in% commandArgs() || length(commandArgs(trailingOnly = TRUE)) == 0) {
  ShowHelp()
  quit()
}


###########################################
## Load packages and install if missing. ##
###########################################
required_packages <- c("DESeq2", "ggplot2", "pheatmap", "RColorBrewer",
                       "apeglm", "dplyr", "EnhancedVolcano", "readr", "rlang")

missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
  cat("Error: The following required packages are missing:\n")
  cat(paste(missing_packages, collapse = ", "), "\n")
  cat("Please install them before running this script.\n")
  quit(status = 1)
}

# Load packages
suppressPackageStartupMessages({
  invisible(lapply(required_packages, library, character.only = TRUE))
})


############################
## Check input arguments. ##
############################
parse_args <- function(args) {
  parsed <- list()
  i <- 1
  while (i <= length(args)) {
    if (args[i] %in% c("--read_counts", "-r")) {
      parsed$read_counts <- args[i + 1]
      i <- i + 2
    } else if (args[i] %in% c("--metadata", "-m")) {
      parsed$metadata <- args[i + 1]
      i <- i + 2
    } else if (args[i] %in% c("--control", "-c")) {
      parsed$control <- args[i + 1]
      i <- i + 2
    } else if (args[i] %in% c("--qval", "-q")) {
      parsed$qval <- as.numeric(args[i + 1])
      i <- i + 2
    } else if (args[i] %in% c("--gene_list", "-g")) {
      parsed$gene_list <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--group_column") {
      parsed$group_column <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--top_var_genes") {
      parsed$top_var_genes <- as.integer(args[i + 1])
      i <- i + 2
    } else {
      cat("Error: Unknown argument", args[i], "\n")
      quit(status = 1)
    }
  }
  return(parsed)
}

parsed_args <- parse_args(args)

if (is.null(parsed_args$group_column)) parsed_args$group_column <- "group"
if (is.null(parsed_args$top_var_genes)) parsed_args$top_var_genes <- 300

#################################
## Function to load gene list. ##
#################################
read_gene_list <- function(file_path) {
  if (is.null(file_path)) return(NULL)
  tryCatch({
    raw_data <- readLines(file_path)
    combined_data <- paste(raw_data, collapse = " ")
    genes <- unlist(strsplit(combined_data, "[ \t,;\n]+"))
    genes <- genes[genes != ""]
    return(genes)
  }, error = function(e) {
    cat("Error reading gene list file:", file_path, "\n")
    cat("Error message:", e$message, "\n")
    return(NULL)
  })
}
annotateGenes <- read_gene_list(parsed_args$gene_list)

#################################################
## Load HTSeq read counts and sample metadata. ##
#################################################
safe_read_csv <- function(file_path, ...) {
  tryCatch(
    {
      read.csv(file_path, ...)
    },
    error = function(e) {
      cat("Error reading file:", file_path, "\n")
      cat("Error message:", e$message, "\n")
      quit(status = 1)
    }
  )
}

validate_input <- function(countData, metaData, control, group_column) {
  if (!group_column %in% colnames(metaData)) {
    cat("Error: '", group_column, "' column not found in metadata\n")
    quit(status = 1)
  }
  if (!control %in% metaData[[group_column]]) {
    cat("Error: '", control, "' not found in '", group_column, "' column of metadata\n")
    quit(status = 1)
  }
  count_samples <- colnames(countData)[-1]
  if (!all(count_samples %in% metaData$sample)) {
    cat("Error: Not all samples in count data match samples in metadata\n")
    cat("Samples in count data:", paste(count_samples, collapse = ", "), "\n")
    cat("Samples in metadata:", paste(metaData$sample, collapse = ", "), "\n")
    quit(status = 1)
  }
}

countData <- safe_read_csv(parsed_args$read_counts, header = TRUE, sep = "\t")
metaData <- safe_read_csv(parsed_args$metadata, header = TRUE, sep = "\t")
control <- parsed_args$control
group_column <- parsed_args$group_column
colnames(metaData) <- tolower(colnames(metaData))

validate_input(countData, metaData, control, group_column)

###############################################
## Perform differential expression analysis. ##
###############################################
findDiff <- function(dds, treatment, control, qval) {
  res <- results(dds, contrast = c("group", treatment, control))
  resSig <- res[which(res$padj < qval), ]
  deg <- resSig[order(resSig$log2FoldChange, decreasing = TRUE), ]
  deg <- cbind(rownames(deg), deg)
  colnames(deg)[colnames(deg) == "rownames(deg)"] <- "geneID"
  outfile <- file.path(dirname(parsed_args$read_counts), paste0(treatment, "_vs_", control, "_DEgenes.txt"))
  write.table(deg, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)
}

dds <- suppressWarnings(DESeqDataSetFromMatrix(countData, metaData, as.formula(paste0("~", group_column)), tidy = TRUE))
dds[[group_column]] <- relevel(dds[[group_column]], ref = control)
dds <- DESeq(dds)
normalized_counts <- as.data.frame(counts(dds, normalized = TRUE))
norm_outfile <- file.path(dirname(parsed_args$read_counts), paste0(control, "_control_normalized_read_counts.txt"))
write.table(normalized_counts, file = norm_outfile, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

for (treatment in unique(metaData[[group_column]][!grepl(parsed_args$control, metaData[[group_column]])])) {
  findDiff(dds, treatment, parsed_args$control, parsed_args$qval)
}


###################################
## Plot output from DE analysis. ##
###################################
pdf(file.path(dirname(parsed_args$read_counts), paste0(control, "_control_RNAseq_figures.pdf")))
color_map <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
  "#0072B2", "#D55E00", "#CC79A7", "#000000", "#999999"
)

# PCA plot
rld <- rlog(dds, blind = TRUE)
pca_rld <- plotPCA(rld, intgroup = group_column, returnData = TRUE)
percentVar <- round(100 * attr(pca_rld, "percentVar"))
ggplot(pca_rld, aes(x = PC1, y = PC2, color = !!sym(group_column))) +
  labs(fill = group_column) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme(text = element_text(family = "serif", size = 12)) +
  scale_color_manual(values = color_map[1:length(unique(metaData[[group_column]]))]) +
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
  fontsize_number = 16,
  treeheight_row = 0,
  fontsize = 16
)

# MA plot and volcano plot
for (treatment in unique(metaData[[group_column]][!grepl(parsed_args$control, metaData[[group_column]])])) {
  resLFCapeglm <- lfcShrink(
    dds,
    coef = paste0("group_", treatment, "_vs_", parsed_args$control),
    type = "apeglm"
  )
  resLFCapeglm <- resLFCapeglm[!is.na(resLFCapeglm$log2FoldChange) & !is.na(resLFCapeglm$padj), ]
  min_log2FC <- min(resLFCapeglm$log2FoldChange, na.rm = TRUE)
  max_log2FC <- max(resLFCapeglm$log2FoldChange, na.rm = TRUE)
  min_padj <- -log10(max(resLFCapeglm$padj, na.rm = TRUE))
  max_padj <- -log10(min(resLFCapeglm$padj, na.rm = TRUE))
  plotMA(
    resLFCapeglm,
    alpha = parsed_args$qval,
    main = paste(treatment, " vs ", parsed_args$control),
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
    ifelse(resLFCapeglm$padj > parsed_args$qval, "#999999", "#000000"))

  names(keyvals)[keyvals == "#000000"] <- "pass log2FC and p-value"
  names(keyvals)[keyvals == "#999999"] <- "fail log2FC or p-value"

  plot(EnhancedVolcano(
    resLFCapeglm,
    lab = "",
    title = paste(treatment, "vs", parsed_args$control),
    caption = paste(nrow(resLFCapeglm), "variables, padj <=", parsed_args$qval, "Log2FC >= 2"),
    colCustom = keyvals,
    pointSize = 1,
    subtitle = NULL,
    captionLabSize = 12,
    legendLabSize = 12,
    legendIconSize = 1.5,
    x = "log2FoldChange",
    y = "padj",
    xlim = c(min_log2FC, max_log2FC),
    ylim = c(min_padj, max_padj),
    axisLabSize = 12,
    pCutoff = parsed_args$qval,
    FCcutoff = 1.0,
    colAlpha = 1,
    border = "partial",
    borderWidth = 1.5,
    borderColour = "#000000",
  ))
}

# Heatmap
rld <- rlog(dds, blind = FALSE)
annoDF <- data.frame(colData(rld)[, group_column])
colnames(annoDF) <- "group"
rownames(annoDF) <- colnames(rld)
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), parsed_args$top_var_genes)
sigMat <- assay(rld)[topVarGenes, ]
sigMat <- sigMat[, order(colnames(sigMat))]
ann_colors = list(group = color_map[1:length(unique(metaData[[group_column]]))])
names(ann_colors$group) <- unique(metaData[[group_column]])
pheatmap(
  sigMat,
  main = paste("Top", parsed_args$top_var_genes, "variable genes"),
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
  fontsize = 12,
  cellwidth = 18,
  cellheight = 1,
  border_color = NA
)

dev.off()
