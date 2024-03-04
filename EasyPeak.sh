#!/usr/bin/env bash

## Created: August 25, 2023
## Updated: March 3, 2024
## Author(s): Todd Lenz, tlenz001@ucr.edu

## ChIPPeaks: Performs peak calling and differential peak calling using an input chip-seq metadata file describing samples.

function help {
    echo "ChIPPeaks.sh --help"
    echo "usage : ChIPPeaks.sh -m METADATA -g GENOME [-c CONTROL] [-q QUALITY] [-t THREADS] [-h]"
    echo
    echo "---------------------------------------------------------------------"
    echo " Required inputs:"
    echo "  -m|--metadata METADATA    : ChIP-seq sample metadata file."
    echo "  -g|--genome GENOME        : Directory containing genome files."
    echo
    echo " Optional inputs:"
    echo "  -c|--control CONTROL      : Control condition."
    echo "  -q|--qval QVAL            : q-value cutoff for peak calling."
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
        "--quality") set -- "$@" "-q" ;;
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
        q) QUALITY="$OPTARG";;
        t) THREADS="$OPTARG";;
        h) help;;
        *) echo "Error: '$OPTARG' is an invalid argument."
    esac
done


############################
## Check input arguments. ##
############################
if [[ -z "$METADATA" ]]; then
    echo "Error: Script must be run with a metadata file."
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
elif [[ -n "$QUALITY" && ! "$QUALITY" =~ ^[0-9]+$ ]]; then
    echo "Error: Invalid value for quality argument."
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


######################################################
## Check Python installation and required packages. ##
######################################################
if ! command -v "python3" &> /dev/null; then
    echo "Installing newest version of python3."
    if command -v apt &> /dev/null; then
        sudo apt update
        sudo apt upgrade
        sudo apt install -y python3
    elif command -v yum &> /dev/null; then
        sudo yum install -y python3
    elif command -v dnf &> /dev/null; then
        sudo dnf update
        sudo dnf install -y python3
    elif command -v zypper &> /dev/null; then
        sudo zypper refresh
        sudo zypper update
        sudo zypper --non-interactive install python3
    elif command -v pacman &> /dev/null; then
        sudo pacman -Syu
        sudo pacman install --noconfirm python3
    elif command -v brew &> /dev/null; then
        brew update
        brew upgrade
        brew install -y python3
    else
        echo "Can't determine package manager. Install python3 manually."
        exit 1
    fi
fi

if ! command -v pip3 &> /dev/null; then
    python -m ensurepip --upgrade
fi

if ! command -v macs2 &> /dev/null; then
    echo "Installing macs2."
    pip3 install macs2
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
function get_chrom_sizes() {
    fa=$(find -L "$1" -mindepth 1 -maxdepth 1 \( -name "*.fasta" -o -name "*.fa" \))
    chrom_sizes=$(find -L "$1" -mindepth 1 -maxdepth 1 \( -name "*.sizes" \))
    bname="${fa%%.*}"
    if [[ -e "$bname.chrom.sizes"  ]]; then
        gsize="$(awk -F '\t' '{ sum += $2 } END { print sum }' "$chrom_sizes")"
    elif [[ $fa = "" && ! -e "$bname.chrom.sizes"  ]]; then
        echo "Error: No FASTA or '.chrom.sizes' file found."
        exit 1
    elif [[ $fa != "" && ! -e "$bname.chrom.sizes" ]]; then
        echo "Finding chromosome lengths."
        python3 "$SCRIPTS"/create_sizes.py "$fa"
        gsize="$(awk -F '\t' '{ sum += $2 } END { print sum }' "$chrom_sizes")"
    fi
}

function get_gtf_file() {
    fa=$(find -L "$1" -mindepth 1 -maxdepth 1 \( -name "*.fasta" -o -name "*.fa" \))
    bname="${fa%%.*}"
    gtf="$bname".gtf
    if [[ $fa = "" && ! -e "$bname.gtf"  ]]; then
        echo "Error: No FASTA or GTF file found."
        exit 1
    elif [[ $fa != "" && ! -e "$bname.gtf" ]]; then
        echo "No GTF file found."
        exit 1
    fi
}

if [[ -z $GENOME ]]; then
    get_chrom_sizes "$pipedir"/genomes
    get_gtf_file "$pipedir"/genomes
else
    genome_dir="$(readlink -f "$GENOME")"
    get_chrom_sizes "$genome_dir"
    get_gtf_file "$genome_dir"
fi


#######################
## Run peak calling. ##
#######################
function check_bam_type() {
    if $(samtools view -f 0x1 "$bamControl" | head -n 1 | wc -l) -eq 1; then
        echo "BAMPE"
    else
        echo "BAM"
    fi
}

function check_alg_type() {
    Factor=$(echo "$1" | awk -F '\t' '{print $2}')
    if [[ $Factor =~ ^H2[AB]K || $Factor =~ ^H3K ]]; then
        echo "--broad"
    else
        echo ""
    fi
}

function peak_calling() {
    echo "Performing peak calling."
    while IFS= read -r line; do
        SampleID="$(echo "$line" | awk -F '\t' '{print $1}')"
        echo "$SampleID"
        bamReads="$(echo "$line" | awk -F '\t' '{print $5}')"
        bamReads="$data_dir/${bamReads#"$data_dir"}"
        bamControl="$(echo "$line" | awk -F '\t' '{print $7}')"
        bamControl="$data_dir/${bamControl#"$data_dir"}"
        bam_format="$(check_bam_type "$bamReads")"
        outdir="$(dirname "$bamReads")"
        alg="$(check_alg_type "$line")"
        macs2 callpeak -t "$bamReads" -c "$bamControl" -f "$bam_format" -g "$gsize" -n "$SampleID" -q "$QVAL" --outdir "$outdir" -B "$alg" >> "$log_dir"/peak_calling.log 2>&1
    done < <(tail -n +2 "$METADATA")
    while IFS= read -r line; do
        if [[ "$(echo "$line" | awk -F'\t' '{print NF}')" -eq 7 ]]; then
            if [[ "$(echo "$line" | awk -F '\t' '{print $1}')" = "SampleID" ]]; then
                echo -e "$line\tPeaks\tPeakCaller"
            else
                bamReads="$(echo "$line" | awk -F '\t' '{print $5}')"
                echo -e "$line\t${bamReads%_*}_peaks.xls\tmacs"
            fi
        else
            echo "$line"
        fi
    done < "$METADATA" > "${METADATA%.*}_updated.txt"
    mv "${METADATA%.*}_updated.txt" "$METADATA"
    echo
}


####################################
## Run differential peak calling. ##
####################################
function diff_peak_calling() {
    if [[ -z "$CONTROL" ]]; then
        echo "Warning: No control argument entered. Differential peak calling step can not be performed."
    elif [[ $(awk '{ print $3 }' < "$METADATA" | grep -c "$CONTROL") -eq 0 ]]; then
        echo "Warning: Control argument not found in metadata. Differential peak calling step can not be performed."
    else
        echo "Finding differential peaks."
        Rscript "$SCRIPTS"/differential_peak_calling.R -o "$data_dir" -m "$(basename "$METADATA")" -c "$CONTROL" >> "$log_dir"/diff_peak_calling.log 2>&1
    fi
}


#####################################################
## Run transcription start site coverage analysis. ##
#####################################################
function TSS_analysis() {
    echo "Performing TSS analysis."
    echo "----------------------------------------------------------------"
    unique_values=$(awk -F '\t' -v column="3" 'NR > 1 {print $column}' "$METADATA" | sort -u)
    for value in $unique_values; do
        bamReads="$(awk -F '\t' -v column="3" -v value="$value" -v prefix="$data_dir/" '$column == value {print prefix $5}' "$METADATA")"
        bamControl="$(awk -F '\t' -v column="3" -v value="$value" -v prefix="$data_dir/" '$column == value {print prefix $7}' "$METADATA")"
        output_dir="$data_dir"/"$value"
        mkdir "$output_dir"
        {
            samtools merge -@ 18 "$output_dir"/"${value%%.*}"_target.bam "$bamReads"
            samtools index -@ 18 "$output_dir"/"${value%%.*}"_target.bam
            samtools merge -@ 18 "$output_dir"/"${value%%.*}"_control.bam "$bamControl"
            samtools index -@ 18 "$output_dir"/"${value%%.*}"_control.bam
            bamCompare -b1 "$output_dir"/"${value%%.*}"_target.bam -b2 "$output_dir"/"${value%%.*}"_control.bam --outFileName "$output_dir"/"${value%%.*}".bw --scaleFactorsMethod None --normalizeUsing RPKM --operation subtract -bs 1 --ignoreDuplicates -p "$THREADS"
        } >> "$log_dir"/TSS_coverage.log 2>&1
    done
    computeMatrix reference-point -S "$data_dir"/*.bw -R "$gtf" -o "$data_dir"/TSS_coverage.matrix -b 500 -a 1500 -bs 1 -p "$THREADS" >> "$log_dir"/TSS_coverage.log 2>&1
    plotProfile -m "$data_dir"/TSS_coverage.matrix -o "$fig_dir"/TSS_coverage_profile.pdf --perGroup >> "$log_dir"/TSS_coverage.log 2>&1
    echo
}


###################
## Run pipeline. ##
###################
STEPS=("peak_calling" "diff_peak_calling")
for s in "${STEPS[@]}"; do
    "$s"
done
