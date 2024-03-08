#!/usr/bin/env bash

## Created: September 22, 2022
## Updated: March 8, 2024
## Author(s): Todd Lenz, tlenz001@ucr.edu

## A complete pipeline for RNA-seq and ChIP-seq analysis.

function help {
    echo "SimpleSeq.sh --help"
    echo "usage : SimpleSeq.sh -i INPUT -g GENOME [-o OUTPUT] [-s STEP] [-q QUALITY] [-t THREADS] [-c] [-r] [-h]"
    echo
    echo "---------------------------------------------------------------------"
    echo " Required inputs:"
    echo "  -i|--input  INPUT       : Input data folder."
    echo "  -g|--genome GENOME      : Directory containing genome files."
    echo
    echo " Optional inputs:"
    echo "  -o|--output OUTPUT      : Output folder."
    echo "  -s|--step STEP          : Choose starting step."
    echo "         quality_check    : Initial quality check."
    echo "              trimming    : Adapter trimming."
    echo "             alignment    : Read alignment."
    echo "         deduplication    : Remove PCR duplicates (for ChIP-seq)."
    echo "             filtering    : Filtering low-quality reads."
    echo "               sorting    : Sorting reads by coordinate."
    echo "               mapping    : Mapping reads to each gene/base pair."
    echo "  -q|--quality QUALITY    : Phred quality score for filtering."
    echo "  -t|--threads THREADS    : Processor threads."
    echo "  -c|--chip CHIP          : Samples are from ChIP-seq experiment."
    echo "  -r|--remove REMOVE      : Remove intermediate files."
    echo "  -h|--help HELP          : Show help message."
    echo "---------------------------------------------------------------------"
    exit 0;
}


#############################
## Define input arguments. ##
#############################
for arg in "$@"; do
    case "$arg" in
        "--input") set -- "$@" "-i" ;;
        "--genome") set -- "$@" "-g" ;;
        "--output") set -- "$@" "-o" ;;
        "--step") set -- "$@" "-s" ;;
        "--quality") set -- "$@" "-q" ;;
        "--threads") set -- "$@" "-t" ;;
        "--chip") set -- "$@" "-c" ;;
        "--remove") set -- "$@" "-r" ;;
        "--help") set -- "$@" "-h" ;;
        *) set -- "$@" "$arg"
    esac
done

pipedir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
SCRIPTS="$pipedir"/scripts

while getopts ":i:g:o:s:q:t:c:r:h:" opt; do
    case $opt in
        i) INPUT="$OPTARG";;
        g) GENOME="$OPTARG";;
        o) OUTPUT="$OPTARG";;
        s) STEP="$OPTARG";;
        q) QUALITY="$OPTARG";;
        t) THREADS="$OPTARG";;
        c) CHIP=true;;
        r) REMOVE=true;;
        h) help;;
        *) echo "Error: '$OPTARG' is an invalid argument."
    esac
done


############################
## Check input arguments. ##
############################
if [[ -z "$INPUT" ]]; then
    echo "Error: Input argument required."
    exit 1
elif [[ ! -e "$INPUT" ]]; then
    echo "Error: Input folder not found."
    exit 1
elif [[ -z "$GENOME" ]]; then
    echo "Error: Genome argument required."
    exit 1
elif [[ ! -e "$GENOME" ]]; then
    echo "Error: Genome folder not found."
    exit 1
elif [[ -n "$THREADS" && "$THREADS" -gt $(nproc) ]]; then
    echo "Error: Invalid value for threads argument."
    exit 1
elif [[ -n "$THREADS" && ! "$THREADS" =~ ^[0-9]+$ ]]; then
    echo "Error: Invalid value for threads argument."
    exit 1
fi


####################
## Define output. ##
####################
if [[ -z "$OUTPUT" ]]; then
    echo "No output folder selected. Files will be saved to input directory."
    OUTPUT="$(dirname "$(readlink -f "$INPUT")")"/output_files
fi
if [[ -e "$OUTPUT" && "$OUTPUT" != "$INPUT" && -n "$STEP" && "$STEP" != "quality_check" ]]; then
    echo "$OUTPUT folder alreads exists. Do you want to overwrite it? (y/n)"
    read -r ans
    if [[ "$ans" = "y" ]]; then
        rm -rf "$OUTPUT"
        mkdir "$OUTPUT"
        for i in "$INPUT"/*; do
            mkdir "$OUTPUT"/"$(basename "$i")"
        done
    fi
elif [[ ! -e "$OUTPUT" ]]; then
    mkdir "$OUTPUT"
    for i in "$INPUT"/*; do
        mkdir "$OUTPUT"/"$(basename "$i")"
    done
fi


############################
## Define logs directory. ##
############################
if [[ ! -e "$(dirname "$OUTPUT")"/logs ]]; then
    log_dir="$(dirname "$OUTPUT")"/logs
    mkdir "$log_dir"
else
    log_dir="$(dirname "$OUTPUT")"/logs
fi


###############################################
## Define Phred score for quality filtering. ##
###############################################
if [[ -z $QUALITY ]]; then
    QUALITY=30
fi


#################################################
## Define number of processing threads to use. ##
#################################################
if [[ -z $THREADS ]]; then
    if command -v "$nproc" &> /dev/null; then
        THREADS=$(($(nproc) / 2))
    else
        THREADS=$($(sysctl -n hw.ncpu) / 2)
    fi
fi


##########################
## Check starting step. ##
##########################
STEP="${STEP/^ //}"
ALL_STEPS=("quality_check" "trimming" "alignment" "deduplication" "filtering" "sorting" "mapping")

if [[ -n $STEP ]]; then
    step="$(echo "$STEP" | tr '[:upper:]' '[:lower:]')"
    echo "Starting analysis pipeline at the $step step."
    case "$step" in
        "quality_check") FILTERED_STEPS=("${ALL_STEPS[@]}");;
        "trimming") FILTERED_STEPS=("${ALL_STEPS[@]:1}");;
        "alignment") FILTERED_STEPS=("${ALL_STEPS[@]:2}");;
        "deduplication") FILTERED_STEPS=("${ALL_STEPS[@]:3}");;
        "filtering") FILTERED_STEPS=("${ALL_STEPS[@]:4}");;
        "sorting") FILTERED_STEPS=("${ALL_STEPS[@]:5}");;
        "mapping") FILTERED_STEPS=("${ALL_STEPS[@]:6}");;
        *) echo "Invalid step: $step"; exit 1;;
    esac
else
    FILTERED_STEPS=("${ALL_STEPS[@]}")
fi

if [[ -z "$CHIP" ]]; then
    FILTERED_STEPS=("${FILTERED_STEPS[@]/deduplication}")
fi


##################################
## Functions/steps in pipeline. ##
##################################
function quality_check() {
    echo "Performing FastQC quality check."
    echo "---------------------------------------------------------------------"
    for i in "$INPUT"/*; do
        basename "$i"
        "$SCRIPTS"/QC_check.sh -i "$i" -o "$OUTPUT" -t "$THREADS" >> "$log_dir"/"$fs".log 2>&1
    done
    echo
}


function trimming() {
    if ! command -v cutadapt &> /dev/null; then
        "$SCRIPTS"/package_installer.sh cutadapt
    fi
    echo "Running Cutadapt to trim reads."
    echo "---------------------------------------------------------------------"
    for i in "$INPUT"/*; do
        basename "$i"
        "$SCRIPTS"/read_trimming.sh -i "$i" -g "$GENOME" -o "$OUTPUT" -q "$QUALITY" -t "$THREADS" >> "$log_dir"/"$fs".log 2>&1
    done
    echo
}


function alignment() {
    if [[ -n "$CHIP" ]]; then
        alignment_tool="bowtie2"
        alignment_flag="-c"
    else
        alignment_tool="hisat2"
        alignment_flag=""
    fi
    if ! command -v "$alignment_tool" &> /dev/null; then
        "$SCRIPTS"/package_installer.sh "$alignment_tool"
    fi
    if [[ "$STEP" != "$fs" ]]; then
        INPUT="$OUTPUT"
    fi
    echo "Running ${alignment_tool} to align reads to the genome."
    echo "---------------------------------------------------------------------"
    for i in "$INPUT"/*; do
        basename "$i"
        "$SCRIPTS"/read_alignment.sh -i "$i" -g "$GENOME" -o "$OUTPUT" -t "$THREADS" "$alignment_flag" >> "$log_dir"/"$fs".log 2>&1
        if [[ -n "$REMOVE" ]]; then
            rm "$(find -L "$i" -mindepth 1 -name "*_trimmed.fastq.gz")"
        fi
    done
    echo
}


function deduplication() {
    if [[ "$STEP" != "$fs" ]]; then
        INPUT="$OUTPUT"
    fi
    echo "Running Picartools to remove PCR duplicates."
    echo "---------------------------------------------------------------------"
    for i in "$INPUT"/*; do
        basename "$i"
        "$SCRIPTS"/pcr_deduplication.sh -i "$i" -o "$OUTPUT" -t "$THREADS" >> "$log_dir"/"$fs".log 2>&1
        if [[ -n "$REMOVE" ]]; then
            find -L "$i" -mindepth 1 -name "*.sam" | while IFS= read -r s; do
                if grep -q 'XT:A:U' "$s"; then
                    :
                else
                    rm "$s"
                fi
            done
        fi
    done
    echo
}


function filtering() {
    if [[ -n "$CHIP" ]]; then
        alignment_flag="-c"
    else
        alignment_flag=""
    fi
    if ! command -v samtools &> /dev/null; then
        "$SCRIPTS"/package_installer.sh samtools
    fi
    if [[ "$STEP" != "$fs" ]]; then
        INPUT="$OUTPUT"
    fi
    echo "Running Samtools to filter low quality alignments."
    echo "---------------------------------------------------------------------"
    for i in "$INPUT"/*; do
        basename "$i"
        samtools_step="filtering"
        "$SCRIPTS"/samtools_wrapper.sh -i "$i" -s "$samtools_step" -o "$OUTPUT" -q "$QUALITY" -t "$THREADS" >> "$log_dir"/"$fs".log 2>&1
        if [[ -n "$REMOVE" ]]; then
            rm "$(find -L "$i" -mindepth 1 -name "*.sam")"
        fi
    done
    echo
}


function sorting() {
    if [[ "$STEP" = "$fs" ]]; then
        INPUT="$OUTPUT"
        if ! command -v samtools &> /dev/null; then
            "$SCRIPTS"/package_installer.sh samtools
        fi
    fi
    echo "Running Samtools to sort and index reads."
    echo "---------------------------------------------------------------------"
    for i in "$INPUT"/*; do
        basename "$i"
        samtools_step="sorting"
        "$SCRIPTS"/samtools_wrapper.sh -i "$i" -s "$samtools_step" -o "$OUTPUT" -t "$THREADS" >> "$log_dir"/"$fs".log 2>&1
        if [[ -n "$REMOVE" ]]; then
            if samtools view -H "$i" | grep -q '^@HD.*SO:coordinate'; then
                :
            else
                rm "$(find -L "$i" -mindepth 1 -name "$i")"
            fi
        fi
    done
    echo
}


function mapping() {
    if [[ "$STEP" = "$fs" ]]; then
        INPUT="$OUTPUT"
    fi
    if [[ -z "$CHIP" ]]; then
        if ! command -v htseq-count &> /dev/null; then
            "$SCRIPTS"/package_installer.sh htseq
        fi
        echo "Running HTSeq to find read coverage for each gene."
        echo "---------------------------------------------------------------------"
        for i in "$INPUT"/*; do
            basename "$i"
            "$SCRIPTS"/read_mapping.sh -i "$i" -g "$GENOME" -o "$OUTPUT" -t "$THREADS" >> "$log_dir"/"$fs".log 2>&1
        done
        count_files=("$OUTPUT"/*/*_counts.txt)
        read_counts="${count_files[0]}"
        cut -f1 "$read_counts" > read_counts.txt
        for cf in "${count_files[@]:1}"; do
            paste read_counts.txt <(cut -f2- "$cf") > merged.tmp && mv merged.tmp read_counts.txt
        done
        header=$(printf "Gene_ID\t%s\n" "$(basename -a "${count_files[@]%%.*}")" | tr '\n' '\t' | sed 's/\t$/\n/')
        { echo "$header"; cat read_counts.txt; } > merged.tmp && mv merged.tmp read_counts.txt
    else
        if [[ "$STEP" = "$fs" ]]; then
            if ! command -v samtools &> /dev/null; then
                "$SCRIPTS"/package_installer.sh samtools
            fi
        fi
        echo "Running Samtools to find read coverage at each base pair."
        echo "---------------------------------------------------------------------"
        samtools_step="mapping"
        for i in "$INPUT"/*; do
            basename "$i"
            "$SCRIPTS"/samtools_wrapper.sh -i "$i" -s "$samtools_step" -o "$OUTPUT" -t "$THREADS" >> "$log_dir"/"$fs".log 2>&1
            bed_file=$(find -L "$i" -mindepth 1 -name "*.bed")
            "$SCRIPTS"/bed2wig.py -i "$bed_file"
        done
    fi
}


###################
## Run pipeline. ##
###################
for fs in "${FILTERED_STEPS[@]}"; do
    "$fs"
done
