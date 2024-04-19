# SimpleSeq
This is a complete pipeline for RNA-seq and ChIP-seq data analysis. Pipeline can be started from any step, but specific files are required for the desired step.

It is important to note that SimpleSeq and any supplementary tools are merely templates to provide fast and simple analysis and may not work for all samples or experimental designs. If new to these analyses, or any of the tools utilized within the pipeline, please learn about each of these tools and cite their respective publications where appropriate.

## How do I organize my data?

Put all input FASTQ files in a single directory with one folder per sample. The sequences can be either single or paired-end and can be gzipped, though this is not a requirement. Paired-end files need to have '_R1' and '_R2' within the file names to specify forward and reverse reads. 

```
    + PATH_TO_INPUT
        + sample1
            ++ sample1_R1.fastq.gz
            ++ sample1_R2.fastq.gz
        + sample2
            ++ sample2_R1.fastq.gz
            ++ sample2_R2.fastq.gz
        + sample 3
            ++ sample3.fastq.gz
        + sample 4
            ++ sample4.fastq.gz
```

## What tools do I need?

Although SimpleSeq automatically downloads most of the necessary software/tools, depending on a user's system environment some packages may need to be installed manually. Below is a list of all necessary tools to run the complete pipeline on both RNA-seq and ChIP-seq data:

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [CutAdapt](https://cutadapt.readthedocs.io/en/stable/)
- [Hisat2](https://daehwankimlab.github.io/hisat2/)
- [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [Picard](https://broadinstitute.github.io/picard/)
- [Samtools](https://www.htslib.org/)
- [HTSeq](https://htseq.readthedocs.io/en/latest/index.html)

### Conda environment installation

To make installation of SimpleSeq dependencies easier, users can create a conda environment using the provided YAML file. Any dependencies for other utilities within this package will also be installed.

1. Download and install [miniconda](https://docs.conda.io/en/latest/miniconda.html) or [anaconda](https://www.anaconda.com/download).

2. Installing the package will automatically attempt to create the SimpleSeq conda environment. If conda was downloaded after installing SimpleSeq, the environment can be created using: ```conda env create -f environment.yml```

3. Activate environment: ```conda activate SimpleSeq```

## What else do I need?

Any genome can be used if properly formatted FASTA and GFF files are provided. Any additional files, such as alignment indexes, will be created automatically if not already present.

If performing ChIP-seq analysis, users may want to download [IGV](https://igv.org/) to visualize the output WIG file.

Additional scripts used to automate plotting are available to users in the ```scripts``` folder if a different pipeline or tools were used to perform the initial analysis.

## How do I run SimpleSeq?

Brief description of input arguments via help message:

```
    usage : SimpleSeq.sh -i INPUT -g GENOME [-o OUTPUT] [-s STEP] [-q QUALITY] [-t THREADS] [-c CHIP] [-r REMOVE] [-h]
    ---------------------------------------------------------------------
     Required inputs:
      -i|--input  INPUT       : Input data folder.
      -g|--genome GENOME      : Directory containing genome files.

     Optional inputs:
      -o|--output OUTPUT      : Output folder.
      -s|--step STEP          : Choose starting step.
             quality_check    : Initial quality check.
                  trimming    : Adapter trimming.
                 alignment    : Read alignment.
             deduplication    : Remove PCR duplicates (for ChIP-seq).
                 filtering    : Filtering low-quality reads.
                   sorting    : Sorting reads by coordinate.
                   mapping    : Mapping reads to each gene/base pair.
      -q|--quality QUALITY    : Phred quality score for filtering.
      -t|--threads THREADS    : Processor threads.
      -c|--chip CHIP          : Samples are from ChIP-seq experiment.
      -r|--remove REMOVE      : Remove intermediate files.
      -h|--help HELP          : Show help message.
    ---------------------------------------------------------------------
```

The following is a more detailed description of input arguments:

- **Required arguments**:

    --input [-i]: Directory containing a folder for each sample. Samples can be single or paired-end, but paired-end sequences need to contain '_R1' and '_R2' within the file names to determine forward and reverse reads. If starting from a step in the pipeline after genome alignment the input files are not required to contain '_R1' or '_R2'. RNA-seq and ChIP-seq data cannot be analyzed simultaneously due to the different tools used by each analysis pipeline.

    --genome [-g]: Directory containing genome files for the organism of interest. A FASTA file should be present in the directory so that any additional required genome files can be automatically created if missing. Alignment indexes should have the same basename as the FASTA file. There should only be a single FASTA file and GFF file in the ```genome``` directory.

- **Optional arguments**:

    --output [-o]: Output directory. If no output is given, the output files will be saved to the input folders for each sample.

    --step [-s]: Desired step at which to start the pipeline. If no step is entered, the entire pipeline will be run. 

    --quality [-q]: Integer indicating Phred quality score for trimming low-quality bases at the ends of reads during read pairing and for removing low-quality reads during filtering (default = 30).

    --threads [-t]: Integer indicating number of processor threads to use for tools that allow multithreading. If no value is given, the number of available threads will be determined automatically and will use half of available threads on the user's system.

    --chip [-c]: Indicates that the input data is from a ChIP-seq experiment. The ChIP-seq processing pipeline differs in two ways: 1. Bowtie2 is used for read alignment rather than Hisat2, and 2. Samtools depth is used to find read depth at each base pair rather than using HTSeq to find gene coverage.

    --remove [-r]: Removes intermediate files when no longer needed by pipeline to save hard drive space.

## What are the individual steps within the pipeline?

- **Quality check**:

    FastQC is used to check sequence quality prior to running the rest of the pipeline. Statistics such as per base sequence quality, per base sequence content, sequence duplication levels and adapter content are calculated from input FASTQ files. See the [project website](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for more details.

- **Trimming**:

    Cutadapt is used to perform initial trimming. Adapter sequences in the provided ```adapters.txt``` file and bases below the quality threshold set by the ```quality [-q]``` argument (default = 30) are removed from the ends of reads. The ```adapters.txt``` file contains common universal adapters for Illumina platforms, but can be modified to fit a specific protocol. Once the reads are trimmed, any reads that are shorter than 25 bp are removed. Single-end sequences are then output as trimmed FASTQ files. See the [project website](https://cutadapt.readthedocs.io/en/stable/) for more details.

- **Alignment**:

    The trimmed FASTQ files are then aligned to the parent genome using either Hisat2 or Bowtie2. Which read alignment tool is used is dependent on whether the data is from an RNA-seq or ChIP-seq experiment. Alignment indexes with extension '.ht2' or '.bt2', depending on which read aligner is used, are generated from the FASTA file found in the directory provided by the ```genome [-g]``` argument if not already available. To ensure that all possible alignments are identified, the --very-sensitive argument is used to specify the high-sensitivity mode. A single SAM file is output for each aligned sample. See the [Hisat2 project website](https://daehwankimlab.github.io/hisat2/) or [Bowtie2 project website](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) for more details.

- **Deduplication**:

    If the ```chip [-c]``` flag is utilized, PCR duplicates are flagged in the SAM files using ```picard markDuplicates```. Deduplication metrics are output in a separate text file. See the [project website](https://broadinstitute.github.io/picard/) for more details.

- **Filtering**:

    Low-quality reads with Phred score less than the value specifid by the ```quality [-q]``` argument are removed from the SAM files using ```samtools view```. The filtered BAM files only include properly paired aligned reads using the flags ```-f 0x02 -F 0x04```. See the [project website](https://www.htslib.org/doc/samtools.html) for more details.

- **Sorting**:

    Once quality filtered, the BAM files are sorted by genome coordinate using ```samtools sort``` to make mapping faster. See the [project website](https://www.htslib.org/doc/samtools.html) for more details.

- **Mapping**:

    If analyzing RNA-seq data, a GFF file in the ```genome [-g]``` directory is utilized to map reads in the sorted BAM file to find total read coverage for each gene using ```htseq-count```. The output BED file is a tab-delimited file with two columns: gene name and read count. See the [project website](https://htseq.readthedocs.io/en/latest/index.html) for more details.
    
    If the ```chip [-c]``` flag is utilized, the sorted BAM files are mapped to the genome to find the depth of coverage at each nucleotide using ```samtools depth```. The output BED file is a tab-delimited file with three columns: chromosome name, position and read count. See the [project website](https://www.htslib.org/doc/samtools.html) for more details.

# EasyDGE: RNA-seq differential gene expression analysis

EasyDGE is a supplementary tool that performs differential gene expression analysis. It is meant to be used after running SimpleSeq on a set of RNA-seq samples, but can be used on any dataset if the proper inputs are provided.

## What tools do I need?

Although EasyDGE attempts to automatically download most of the necessary software/tools, depending on a user's system environment some packages may need to be installed manually. Below is a list of all necessary tools to perform differential expression analysis and plotting of figures:

- [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
- [ggplot2](https://ggplot2.tidyverse.org/)
- [pheatmap](https://r-charts.com/correlation/pheatmap/)
- [RColorBrewer](https://r-graph-gallery.com/38-rcolorbrewers-palettes.html)
- [apeglm](https://bioconductor.org/packages/release/bioc/html/apeglm.html)
- [dplyr](https://dplyr.tidyverse.org/)
- [EnhancedVolcano](https://bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html)

The SimpleSeq conda environment described above also contains all of the necessary R packages and is therefore recommended while running EasyDGE.

## How do I run EasyDGE?

Brief description of input arguments via help message:

```
    EasyDGE.sh --help
    usage : EasyDGE.sh -i INPUT [-o OUTPUT] [-q QUALITY] [-h]

    ---------------------------------------------------------------------
     Required inputs:
      -i|--input INPUT          : Input read counts file.
      -m|--metadata METADATA    : Metadata file.
      -c|--control CONTROL      : Control group.

     Optional inputs:
      -o|--output OUTPUT         : Output folder.
      -q|--qvalue QVALUE         : q-value cutoff.
      -h|--help HELP             : Show help message.
    ---------------------------------------------------------------------
```

The following is a more detailed description of input arguments:

- **Required arguments**:

    --input [-i]: A tab-delimited file of gene read counts for each sample. If SimpleSeq was used to process samples a ```read_counts.txt``` file can be found in the ```output [-o]``` directory. If not, user's can manually generate a ```read_counts.txt``` file containing per-gene read coverage by running ```htseq-count``` on BAM files for each sample and then merge the outputs into a single file using the ```merge_htseq_counts.sh``` script contained in the ```scripts``` directory. The ```read_counts.txt``` file should have the following format:

    ```
    Gene_ID    Sample1    Sample2    Sample3 ... Sample4
    gene1
    gene2
    gene3
    ...
    geneN
    ```

    --metadata [-m]: A metadata file describing the experimental protocol. There should be at least two columns in the metadata file (see below), the first being a list of unique sample names corresponding to columns in the ```read_counts.txt``` file, and the second providing group names for the samples (e.g. KO and CTRL). The group names are used to differentiate a group of control samples from all other samples. Any number of conditional groups (i.e. non-control groups) can be specified.

    ```
    Sample     Group
    Sample1    KO
    Sample2    KO
    Sample3    CTRL
    Sample4    CTRL
    ```

    --control [-c]: String representing the control group in the ```metadata [-m] file```.

- **Optional arguments**:

    --qvalue [-q]: q-value cutoff for filtering results of [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) differential expression analysis. 

# EasyPeaks: ChIP-seq peak calling analysis

EasyPeaks is a supplementary tool that performs peak calling and differential peak calling. It is meant to be used after running SimpleSeq on a set of ChIP-seq samples, but can be used on any dataset if the proper inputs are provided.

## What tools do I need?

Although EasyPeaks attempts to automatically download most of the necessary software/tools, depending on a user's system environment some packages may need to be installed manually. Below is a list of all necessary tools to perform peak calling, differential peak calling and plotting of figures:

- [MACS2](https://github.com/macs3-project/MACS)
- [ggplot2](https://ggplot2.tidyverse.org/)
- [pheatmap](https://r-charts.com/correlation/pheatmap/)
- [RColorBrewer](https://r-graph-gallery.com/38-rcolorbrewers-palettes.html)
- [apeglm](https://bioconductor.org/packages/release/bioc/html/apeglm.html)
- [dplyr](https://dplyr.tidyverse.org/)
- [EnhancedVolcano](https://bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html)

The SimpleSeq conda environment described above also contains all of the necessary R packages and is therefore recommended while running EasyDGE.

## How do I run EasyPeaks?

Brief description of input arguments via help message:

```
    ChIPPeaks.sh --help
    usage : ChIPPeaks.sh -m METADATA -c CONTROL -g GENOME [-t THREADS] [-h]

    ---------------------------------------------------------------------
     Required inputs:
      -m|--metadata METADATA    : ChIP-seq sample metadata file.
      -g|--genome GENOME        : Directory containing genome files.

     Optional inputs:
      -c|--control CONTROL      : Control condition (and column).
      -q|--qval QVAL            : q-value cutoff for peak calling.
      -t|--threads THREADS      : Processor threads.
      -h|--help HELP            : Show help message.
    ---------------------------------------------------------------------
```

The following is a more detailed description of input arguments:

- **Input arguments**:

    --metadata [-m]: A tab-delimited text file containing ChIP-seq experiment metadata. Metadata file should have at least seven columns: SampleID, Factor, Condition, Replicate, bamReads, ControlID and bamControl. The metadata file will be automatically modified to include two additional columns--Peaks and PeakCaller--after peak calling. Use the ```example_metadata.txt``` file in the ```examples``` directory as a template. See the [DiffBind vignette](https://bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf) and the next section of this README for more information on setting up the metadata file.

    --control[-c]: String indicating the control samples in the experiment. This string should be under the 'Condition' column in the metadata file.

    --genome [-g]: Directory containing genome files for the organism of interest. A FASTA file should be present in the directory so that any additional required genome files can be automatically created if missing. All genome files should have the same basename.

## How should I set up the metadata?

The metadata file should contain a minimum of 7 columns used to describe the format of a user's ChIP-seq experiment, including sample names, relationships between sample groups, replicate numbers, and file paths. Although it isn't required for the metadata to be formatted exactly as shown in the ```example_chip_metadata.txt``` file, it's recommended that any additional columns be added at the end (to the right) of the existing columns. By default, the 'Condition' column is used to group samples for differential peak calling. File paths should be relative to the metadata file path.

The following are brief descriptions of the template columns:

    --SampleID: Short string unique to each target sample.

    --Factor: Antibody used for immunoprecipitation of each target sample. Used to determine peak calling parameters.

    --Condition: String used to differentiate sample groups (e.g. knockout (KO), control (CTRL), etc.). Can be the same as the 'Factor' column if differentiating binding sites of two (or more) antibodies.

    --Replicate: Integer indicating replicate number.

    --bamReads: Relative path to coordinate sorted BAM files for target samples.

    --ControlID: Short string unique to each control sample.

    --bamControl: Relative path to coordinate sorted BAM files for control samples.
