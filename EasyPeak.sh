#!/usr/bin/env bash

## Created: August 25, 2023
## Updated: March 3, 2024
## Author(s): Todd Lenz, tlenz001@ucr.edu

## EasyPeaks: Performs peak calling and differential peak calling.

function help {
    echo "EasyPeaks.sh --help"
    echo "usage : EasyPeaks.sh -m METADATA -c CONTROL -g GENOME [-q QVALUE] [-t THREADS] [-h]"
    echo
    echo "---------------------------------------------------------------------"
    echo " Required inputs:"
    echo "  -m|--metadata METADATA    : ChIP-seq sample metadata file."
    echo "  -c|--control CONTROL      : Control group."
    echo "  -g|--genome GENOME        : Directory containing genome files."
    echo
    echo " Optional inputs:"
    echo "  -s|--step STEP            : Choose starting step."
    echo "            peak_calling    : Initial quality check."
    echo "       diff_peak_calling    : Adapter trimming."
    echo "           TSS_profiling    : Read alignment."
    echo "  -q|--qvalue QVALUE        : q-value cutoff for peak calling."
    echo "  -t|--threads THREADS      : Processor threads."
    echo "  -h|--help HELP            : Show help message."
    echo "---------------------------------------------------------------------"
    exit 0;
}


#############################
## Define input arguments. ##
#############################
for arg in "$@"; do
    case "$arg" in
        "--metadata") set -- "$@" "-m" ;;
        "--control") set -- "$@" "-c" ;;
        "--genome") set -- "$@" "-g" ;;
        "--qvalue") set -- "$@" "-q" ;;
        "--threads") set -- "$@" "-t" ;;
        "--help") set -- "$@" "-h" ;;
        *) set -- "$@" "$arg"
    esac
done

pipedir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
SCRIPTS="$pipedir"/scripts

while getopts ":m:c:g:q:t:h:" opt; do
    case $opt in
        m) METADATA="$OPTARG";;
        c) CONTROL="$OPTARG";;
        g) GENOME="$OPTARG";;
        q) QVALUE="$OPTARG";;
        t) THREADS="$OPTARG";;
        h) help;;
        *) echo "Error: '$OPTARG' is an invalid argument."
    esac
done


############################
## Check input arguments. ##
############################
if [[ -z "$METADATA" ]]; then
    echo "Error: Metadata file required."
    exit 1
elif [[ ! -e "$METADATA" ]]; then
    echo "Error: Metadata file not found."
    exit 1
elif [[ -z "$GENOME" ]]; then
    echo "Error: Genome argument required."
    exit 1
elif [[ ! -e "$GENOME" ]]; then
    echo "Error: Genome folder not found."
    exit 1
elif [[ -n "$QVALUE" && ! "$QVALUE" =~ ^-?[0-9]+(\.[0-9]+)?$ ]]; then
    echo "Error: Invalid value for qualue argument."
    exit 1
elif [[ -n "$QVALUE" && "$QVALUE" < 0 || "$QVALUE" > 1 ]]; then
    echo "Error: Invalid value for qualue argument."
    exit 1
elif [[ -n "$THREADS" && "$THREADS" -gt $(nproc) ]]; then
    echo "Error: Invalid value for threads argument."
    exit 1
elif [[ -n "$THREADS" && ! "$THREADS" =~ ^[0-9]+$ ]]; then
    echo "Error: Invalid value for threads argument."
    exit 1
fi


########################################
## Define figure and log directories. ##
########################################
data_dir="$(dirname "$(readlink -f "$METADATA")")"

if [[ ! -e "$data_dir"/figures ]]; then
    fig_dir="$data_dir"/figures
    mkdir "$fig_dir"
else
    fig_dir="$data_dir"/figures
fi

if [[ ! -e "$data_dir"/logs ]]; then
    log_dir="$data_dir"/logs
    mkdir "$log_dir"
else
    log_dir="$data_dir"/logs
fi


########################################################
## Define q-value cutoff for differential expression. ##
########################################################
if [[ -z $QVALUE ]]; then
    QVALUE=0.05
fi


#################################################
## Define number of processing threads to use. ##
#################################################
if [[ -z $THREADS ]]; then
    THREADS=$(($(nproc) / 2))
fi


#######################################
## Check for necessary genome files. ##
#######################################
function get_gtf() {
    gtf=$(find -L "$1" -mindepth 1 -name "*.gff")
    if [[ -z "$gtf" ]]; then
        echo "Error: No GTF file found."
        exit 1
    fi
}

if [[ -z "$GENOME" ]]; then
    echo "Error: Genome argument required."
    exit 1
elif [[ ! -e "$GENOME" ]]; then
    echo "Error: Genome folder not found."
    exit 1
else
    genome_dir="$(readlink -f "$GENOME")"
    get_gtf "$genome_dir"
fi


##########################
## Check starting step. ##
##########################
STEP="${STEP/^ //}"
ALL_STEPS=("peak_calling" "differential_peak_calling" "TSS_profiling")

if [[ -n $STEP ]]; then
    step="$(echo "$STEP" | tr '[:upper:]' '[:lower:]')"
    echo "Starting analysis pipeline at the $step step."
    case "$step" in
        "peak_calling") FILTERED_STEPS=("${ALL_STEPS[@]}");;
        "differential_peak_calling") FILTERED_STEPS=("${ALL_STEPS[@]:1}");;
        "TSS_profiling") FILTERED_STEPS=("${ALL_STEPS[@]:2}");;
        *) echo "Invalid step: $step"; exit 1;;
    esac
else
    FILTERED_STEPS=("${ALL_STEPS[@]}")
fi


##################################
## Functions/steps in pipeline. ##
##################################
function peak_calling() {
    pkgs=("macs2" "samtools")
    for pkg in "${pkgs[@]}"; do
        if ! command -v "macs2" &> /dev/null; then
            "$SCRIPTS"/package_installer.sh "macs2"
        fi
    done
    echo "Performing peak calling."
    echo "---------------------------------------------------------------------"
    "$SCRIPTS"/peak_calling.sh -m "$METADATA" -g "$GENOME" -q "$QVALUE" >> "$log_dir"/"$step".log 2>&1
}


function diff_peak_calling() {
    if [[ -z "$CONTROL" ]]; then
        echo "Error: Control argument required for differential peak calling."
    elif [[ $(awk '{ print $3 }' < "$METADATA" | grep -c "$CONTROL") -eq 0 ]]; then
        echo "Error: Control argument not found in metadata."
    else
        echo "Finding differential peaks."
        echo "-----------------------------------------------------------------"
        Rscript "$SCRIPTS"/diff_peaks.R -m "$METADATA" -c "$CONTROL" >> "$log_dir"/"$step".log 2>&1
    fi
}


function TSS_profiling() {
    echo "Performing TSS coverage profiling."
    echo "---------------------------------------------------------------------"
    unique_values=$(awk 'NR > 1 { print $3 }' "$METADATA" | sort -u)
    for value in $unique_values; do
        bamReads="$(awk -v prefix="$data_dir/" '$3 == $value { print prefix $5 }' "$METADATA")"
        bamControl="$(awk -v prefix="$data_dir/" '$3 == $value { print prefix $7 }' "$METADATA")"
        outdir="$OUTPUT/$(basename "$(dirname "$(readlink -f "$tfq_r1")")")/$(basename "${tfq_r1%*}")"
        {
            samtools merge -@ 18 "$output_dir"/"${value%%.*}"_target.bam "$bamReads"
            samtools index -@ 18 "$output_dir"/"${value%%.*}"_target.bam
            samtools merge -@ 18 "$output_dir"/"${value%%.*}"_control.bam "$bamControl"
            samtools index -@ 18 "$output_dir"/"${value%%.*}"_control.bam
            bamCompare -b1 "$output_dir"/"${value%%.*}"_target.bam -b2 "$output_dir"/"${value%%.*}"_control.bam --outFileName "$output_dir"/"${value%%.*}".bw --scaleFactorsMethod None --normalizeUsing RPKM --operation subtract -bs 1 --ignoreDuplicates -p "$THREADS"
        } >> "$log_dir"/"$step".log 2>&1
    done
    computeMatrix reference-point -S "$data_dir"/*.bw -R "$gtf" -o "$data_dir"/TSS_coverage.matrix -b 500 -a 1500 -bs 1 -p "$THREADS" >> "$log_dir"/TSS_coverage.log 2>&1
    plotProfile -m "$data_dir"/TSS_coverage.matrix -o "$fig_dir"/TSS_coverage_profile.pdf --perGroup >> "$log_dir"/TSS_coverage.log 2>&1
    echo
}


function plot_chromosomes() {
    echo
}


function plot_genes() {
    echo
}


###################
## Run pipeline. ##
###################
STEPS=("peak_calling") # "diff_peak_calling", "TSS_profiling", "plot_chromosomes", "plot_genes")
for step in "${STEPS[@]}"; do
    "$step"
done
