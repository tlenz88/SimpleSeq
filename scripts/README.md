# ChIPPy Scripts
There are several scripts used by ChIPPy to generate input files or intermediate files required by other tools within the pipeline.

## create_sizes.py

Finds the length of all chromosomes in a FASTA file and outputs a tab-delimited text file. The output will be stored in the same directory as the input FASTA file.

**Usage**: ```python3 create_sizes.py *.fasta```

## differential_peak_calling.R

Performs differential peak calling using a metadata file describing the experimental setup. See ```example_metadata.txt``` in the ```examples``` folder. The metadata file should contain the path to the BAM files with respect to the metadata file--e.g. if chip_metadata.txt is in the OUTPUT_DIR the path to the BAM files should be ```output_files/sample1_H3K9/sample1_H3K9.bam``` without OUTPUT_DIR in the path before output_files.

**Usage**: ```Rscript differential_peak_calling.R -o OUTPUT -m chip_metadata.txt -c CONTROL```

- **Input arguments**:

    --output [-o]: Path to parent directory containing output folder, logs folder and chip metadata file.

    --metadata [-m]: Path to the tab-delimited text file containing ChIP-seq experiment metadata. Columns of metadata file should have at least seven columns: SampleID, Factor, Condition, Replicate, bamReads, ControlID and bamControl. The metadata file will be automatically modified to include two additional columns--Peaks and PeakCaller--after peak calling. Use the ```example_metadata.txt``` file in the ```examples``` directory as a template. See the [DiffBind vignette](https://bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf) and the following section of this README for more information on setting up the metadata file.

    --control[-c]: String indicating the control samples in the experiment. This string should be under the 'Condition' column in the metadata file.

## peaks2genes.py

Finds peaks identified by ```macs2 callpeak``` whose summit is contained within genes and 5' promoter regions in a ```gff``` file. Accepts either ```narrowPeak``` or ```broadPeak``` input.

**Usage**: ```python3 peaks2genes.py *.gff *.narrowPeak```

## plot_chromosome_coverage.py

Plots genomewide coverage for ChIP-seq data. Track lengths are normalized with respect to the longest chromosome--i.e. the longest chromosome will fill the width of the figure and all other chromosomes are proportionally plotted against it. All samples will be grouped by chromosome and the desired genes will be plotted at the bottom.

**Usage**: ```plot_chromosome_coverage.py -b sample1.bed sample2.bed sample3.bed```

- **Required input arguments**:

    --bed [-b]: BED file(s) containing genome-wide per-base coverage. To generage this file from a BAM file: ```samtools depth -a -o input.bed input.bam```

- **Optional input arguments**:

    --gff [-g]: Standard GFF file for genome/organism of interest.

    --samples [-s]: List of sample names to annotate the y-axes of barplots. List should be the same length as the list of BED files and plots will be annotated in the same order as BED files. If no sample names are given, the BED file names are used to annotate plots.

    --output [-o]: Directory and name for the output PDF file. If no output is given, the output file is 'chromosome_coverage.pdf' and saved in the directory of the first BED input.

    --resolution [-r]: Integer indicating the binning resolution (default = 1000).

    --centromeres [-c]: Tab-delimited file containing centromere coordinates for each chromosome. File should have three columns: chromosome names, start coordinates and end coordinates. The values in the 'chromosome names' column should match the chromosome names in the BED files.

    --gene_list [-l]: List of genes to plot, separated by spaces, commas, newlines, or tabs. Genes will be extracted from the provided GFF file. If no list is provided, all genes will be annotated on plots.

    --normalization [-n]: If argument is provided, data is counts-per-million (CPM) normalized prior to plotting. This allows for direct comparison of samples regardless of sequencing depth.

    --ymax [-y]: If argument is provided, all chromosome plots will use this maximum y-value for the entire binned dataset. If not provided, each chromosome will use the maximum value for that chromosome.

## plot_gene_coverage.py

**Usage**: ```plot_gene_coverage.py -b sample1.bed sample2.bed sample3.bed -g genome.gff```

- **Required input arguments**:

    --bed [-b]: BED file(s) containing genome-wide per-base coverage. To generage this file from a BAM file: ```samtools depth -a -o input.bed input.bam```

    --gff [-g]: Standard GFF file for genome/organism of interest.

- **Optional input arguments**:

    --samples [-s]: List of sample names to annotate the y-axes of barplots. List should be the same length as the list of BED files and plots will be annotated in the same order as BED files. If no sample names are given, the BED file names are used to annotate plots.

    --output [-o]: Directory and name for the output PDF file. If no output is given, the output file is 'chromosome_coverage.pdf' and saved in the directory of the first BED input.

    --resolution [-r]: Integer indicating the binning resolution (default = 10).

    --gene_list [-l]: List of genes to plot, separated by spaces, commas, newlines, or tabs. Genes will be extracted from the provided GFF file. If no list is provided, all genes will be annotated on plots.

    --normalization [-n]: If argument is provided, data is counts-per-million (CPM) normalized prior to plotting. This allows for direct comparison of samples regardless of sequencing depth.

    --distance [-d]: Integer indicating length of 5' and 3' regions to plot. By default, only the coding region for the plotted genes is used.

    --ymax [-y]: If argument is provided, all chromosome plots will use the maximum y-value for the entire binned dataset. If not provided, each chromosome will use the maximum value for that chromosome.
