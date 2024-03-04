#!/usr/bin/env Rscript

## Created: June 19, 2022
## Updated: February 24, 2024
## Author(s): Todd Lenz, tlenz001@ucr.edu

## Performs differential peak calling via DiffBind.

suppressPackageStartupMessages({
  library(BiocParallel)
  library(DiffBind)
  library(tidyverse)
  library(RColorBrewer)
})

args <- commandArgs(trailingOnly = TRUE)

ShowHelp <- function() {
  cat("differential_peak_calling.R --help\n")
  cat("Usage: differential_peak_calling.R -m METADATA -c CONTROL\n")
  cat("\n")
  cat("-------------------------------------------------------------\n")
  cat("  Input arguments:\n")
  cat("    -m|--metadata METADATA : ChIP-seq metadata file name.\n")
  cat("    -c|--control CONTROL   : String indicating control group.\n")
  cat("    -h|--help HELP         : Show help message.\n")
  cat("-------------------------------------------------------------\n")
}

if ("--help" %in% args || "-h" %in% commandArgs() || length(commandArgs(trailingOnly = TRUE)) == 0) {
  ShowHelp()
  quit()
}

############################
## Check input arguments. ##
############################
for (i in c(1, 3, 5)) {
  if (!(args[i] %in% c("--metadata", "-m", "--control", "-c"))) {
    cat(paste("Error: ", args[i], " is an invalid argument.\n"))
    quit()
  } else if (args[i] == "--metadata" || args[i] == "-m") {
    if (file.exists(args[i + 1])) {
      metadata <- args[i + 1]
      outdir <- dirname(metadata)
    } else {
      cat("Error: Invalid input for metadata argument.\n")
    }
  } else if (args[i] == "--control" || args[i] == "-c") {
    control <- args[i + 1]
  }
}


############################################
## Load peaksets called by MACS callpeak. ##
############################################
samples <- read.csv(metadata, sep = "\t")
if (!(control %in% samples$Condition)) {
  cat(paste("Error: ", control, " not in 'Condition' column of metadata\n"))
}
chip_df <- dba(sampleSheet = samples)


#################################################
## Perform differential peak calling analysis. ##
#################################################
chip_consensus <- dba.peakset(chip_df,
  consensus = c(DBA_CONDITION),
  minOverlap = 1
)
chip_df <- dba.count(chip_df, bUseSummarizeOverlaps = TRUE, bParallel = FALSE)
chip_df <- dba.contrast(chip_df,
  reorderMeta = list(Condition = control),
  minMembers = 2
)
chip_df <- dba.analyze(chip_df, method = DBA_DESEQ2, bParallel = FALSE)
chip_df.report <- dba.report(chip_df,
  method = DBA_DESEQ2,
  th = 0.05,
  bUsePval = FALSE,
  fold = 1,
  bNormalized = TRUE,
  bFlip = FALSE,
  precision = 0,
  bCounts = TRUE
)
diff_peaks <- as.data.frame(chip_df.report)
write.table(diff_peaks,
  file = file.path(outdir, "diff_peaks.txt"),
  sep = "\t",
  quote = F,
  row.names = F
)

##########################################################
## Plot output from differential peak calling analysis. ##
##########################################################
pdf(file.path(outdir, "ChIPseq_figures.pdf"))

# Sample correlation heatmap
dba.plotHeatmap(chip_df,
  correlations = TRUE,
  report = chip_df.report,
  colScheme = "Reds"
)

#
if (nrow(diff_peaks) > 1) {
  hmap <- colorRampPalette(c("blue", "white", "red"))(n = 15)
  dba.plotHeatmap(chip_df,
    contrast = 1,
    correlations = FALSE,
    scale = "row",
    report = chip_df.report,
    colScheme = hmap
  )
  dba.plotPCA(chip_df,
    contrast = 1,
    th = 0.05,
    report = chip_df.report,
    vColors = c("#f8766d", "#00bfc4")
  )
}

# Venn diagram
dba.plotVenn(chip_consensus, chip_consensus$masks$Consensus)

# MA plot
dba.plotMA(chip_df,
  th = 0.05,
  fold = 1,
  bNormalized = TRUE,
  dotSize = 1,
  bSmooth = FALSE,
  bXY = FALSE
)

# Boxplot
dba.plotBox(chip_df, notch = FALSE, vColors = c("#f8766d", "#00bfc4"))

# Volcano plot
dba.plotVolcano(chip_df, th = 0.05, fold = 1)

dev.off()

# Custom volcano plot parameters
# How to customize volcano plot elements:
# debug(DiffBind:::pv.DBAplotVolcano)
# When you get the line line that says "plot(p)", enter the following,
# then save as PDF size 7" x 7":
#p + geom_point(aes(col = Legend), size = 1) + 
#  xlim(-7.6, 7.6) + 
#  labs(title = "NF54CSAh vs NF54", caption = "FDR <= 0.05, FC >= 2") +
#  ylab(expression(-log[10]~italic(FDR))) + 
#  xlab(expression(log[2]~fold~change)) + 
#  scale_color_manual(values=c("grey60", "black"), 
#    labels = c("fail log2FC or FDR", "pass log2FC and FDR")) + 
#  geom_hline(yintercept = -log10(0.05),
#    linetype = "longdash", 
#    color = "black",
#    size = 0.4) + 
#  geom_vline(xintercept = log2(2),
#    linetype = "longdash", 
#    color = "black",
#    size = 0.4) +
#  geom_vline(xintercept = -log2(2),
#    linetype = "longdash", 
#    color = "black",
#    size = 0.4) +
#  theme(plot.title = element_text(size = 18,
#      face = "bold",
#      family = "serif",
#      margin = margin(7, 0, 7, 0, "pt")),
#    axis.text.y = element_text(size = 10,
#      family = "serif",
#      margin = margin(0, 5, 0, 0, "pt")),
#    axis.text.x = element_text(size = 10,
#      family = "serif",
#      margin = margin(5, 0, 0, 0, "pt")),
#    axis.title.y = element_text(size = 10,
#      family = "serif",
#      margin = margin(0, 7, 0, 7, "pt")),
#    axis.title.x = element_text(size = 10,
#      family = "serif",
#      margin = margin(7, 0, 0, 0, "pt")),
#    axis.line = element_line(color = "black", size = 1),
#    axis.ticks = element_line(color = "black", linewidth = 1),
#    axis.ticks.length = unit(4, "pt"),
#    panel.grid = element_line(color = "grey92", size = 1),
#    panel.background = element_blank(),
#    legend.position = "top",
#    legend.title = element_blank(),
#    legend.text = element_text(color = "black",
#      size = 10,
#      family = "serif"),
#    legend.background = element_blank(),
#    legend.key = element_blank(),
#    legend.margin = margin(14, 0, 18, 0, "pt"),
#    plot.caption = element_text(size = 10,
#      family = "serif",
#      margin = margin(14, 10, 7, 0, "pt")),
#  plot.margin = margin(7, 14, 7, 7, "pt"))
