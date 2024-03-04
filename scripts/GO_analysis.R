#!/usr/bin/env Rscript

## Created: October 5, 2023
## Updated: February 24, 2024
## Author(s): Todd Lenz, tlenz001@ucr.edu

## Performs GO enrichment analysis on output from DESeq2.

suppressPackageStartupMessages({
  library(topGO)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
})

args <- commandArgs(trailingOnly = TRUE)

ShowHelp <- function() {
  cat("GO_analysis.R --help\n")
  cat("Usage: GO_analysis.R -d DEGenes -m MAP\n")
  cat("---------------------------------------------------------------------\n")
  cat("  Input arguments:\n")
  cat("    -d|--DEGenes DEGenes              : Differential expression data.\n")
  cat("    -m|--map MAP                      : GO term mapping file.\n")
  cat("    -h|--help HELP                    : Show help message.\n")
  cat("---------------------------------------------------------------------\n")
}

if ("--help" %in% args || "-h" %in% commandArgs() || length(commandArgs(trailingOnly = TRUE)) == 0) {
  ShowHelp()
  quit()
}


############################
## Check input arguments. ##
############################
for (i in c(1, 3)) {
  if (!(args[i] %in% c("--DEGenes", "-d", "--map", "-m"))) {
    cat(paste("Error: ", args[i], "is an invalid argument.\n"))
    quit()
  } else if (args[i] == "--DEGenes" || args[i] == "-d") {
    if (file.exists(args[i + 1])) {
      setwd(dirname(args[i + 1]))
      DEGenes <- args[i + 1]
      outdir <- dirname(DEGenes)
    } else {
      cat("Error: Invalid input for DEGenes argument.\n")
    }
  } else if (args[i] == "--map" || args[i] == "-m") {
    if (file.exists(args[i + 1])) {
      GO_map <- args[i + 1]
    } else {
      cat("Error: Invalid input for MAP argument.\n")
    }
  }
}


############################################################
## Function to perform gene ontology enrichment analysis. ##
############################################################
GOenrichment <- function(DEGenes, gene_names, GO_map) {
  gene_list <- factor(as.integer(gene_names %in% DEGenes[, 1]))
  names(gene_list) <- gene_names
  GO_data <- new("topGOdata",
    ontology = "BP",
    allGenes = gene_list,
    annot = annFUN.gene2GO,
    gene2GO = GO_map
  )
  fisher_test <- runTest(GO_data, statistic = "fisher")
  fisher_test_results <- GenTable(GO_data,
    weighted_Fisher = fisher_test,
    orderBy = "weighted_Fisher",
    ranksOf = "weighted_Fisher",
    topNodes = 50
  )
  fisher_terms = c(fisher_test_results$GO.ID)
  fisher_genes <- genesInTerm(GO_data, fisher_terms)
  fisher_pval <- as.vector(as.numeric(fisher_test_results$weighted_Fisher))
  fisher_test_results$logpval <- -log10(fisher_pval)
  return(fisher_test_results)
}


################################################
## Perform gene ontology enrichment analysis. ##
################################################
GO_map <- readMappings(GO_map)
gene_names <- names(GO_map)
DEGenes <- read.csv(DEGenes, head = TRUE, sep = "\t")
DEGenes_pos <- DEGenes[which(DEGenes$log2FoldChange > 0), ]
DEGenes_neg <- DEGenes[which(DEGenes$log2FoldChange < 0), ]
GO_results_pos <- GOenrichment(DEGenes_pos, gene_names, GO_map)
GO_results_neg <- GOenrichment(DEGenes_neg, gene_names, GO_map)
write.table(GO_results_pos,
  file = file.path(outdir, "upregulated_GO_enrichment.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(GO_results_neg,
  file = file.path(outdir, "downregulated_GO_enrichment.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


########################################################
## Plot results of gene ontology enrichment analysis. ##
########################################################
pdf(file.path(outdir, "GO_enrichment_barplot.pdf"))

max_val <- max(GO_results_pos[1:10, ]$logpval, GO_results_neg[1:10, ]$logpval)

ggplot(GO_results_pos[1:10, ]) +
  geom_bar(aes(x = reorder(Term, logpval), y = logpval),
    fill = "#F3766E",
    position = "dodge",
    stat = "identity",
    width = 0.8) +
  coord_flip() +
  ggtitle("GO enrichment of upregulated genes") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max_val)) +
  scale_fill_manual(values = "#F3766E") +
  xlab(element_blank()) +
  ylab("-log10 p-value") +
  theme(text = element_text(family = "serif"),
    axis.text.x = element_text(color = "#000000", size = 8),
    axis.text.y = element_text(color = "#000000", size = 8),
    axis.line.x = element_line(color = "#000000", linewidth = 0.1),
    axis.line.y = element_line(color = "#000000", linewidth = 0.1),
    panel.background = element_blank(),
    aspect.ratio = 1/2,
    panel.grid.major.x = element_line(linewidth = 0.1, 
      linetype = "dashed", 
      colour = "#666666"
    )
  )

ggplot(GO_results_neg[1:10, ]) +
  geom_bar(aes(x = reorder(Term, logpval), y = logpval),
    fill = "#1CBDC2",
    position = "dodge",
    stat = "identity",
    width = 0.8) +
  coord_flip() +
  ggtitle("GO enrichment of downregulated genes") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max_val)) +
  scale_fill_manual(values = "#1CBDC2") +
  xlab(element_blank()) +
  ylab("-log10 p-value") +
  theme(text = element_text(family = "serif"),
    axis.text.x = element_text(color = "#000000", size = 8),
    axis.text.y = element_text(color = "#000000", size = 8),
    axis.line.x = element_line(color = "#000000", linewidth = 0.1),
    axis.line.y = element_line(color = "#000000", linewidth = 0.1),
    panel.background = element_blank(),
    aspect.ratio = 1/2,
    panel.grid.major.x = element_line(linewidth = 0.1, 
      linetype = "dashed",
      colour = "#666666"
    )
  )

#figure <- ggarrange(p1, p2, ncol = 1, nrow = 2)

dev.off()