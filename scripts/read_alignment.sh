#!/usr/bin/env bash

## Created: February 22, 2024
## Updated: April 19, 2024
## Author(s): Todd Lenz, tlenz001@ucr.edu

## Wrapper for running Hisat2/Bowtie2 on trimmed fastq files.

# Parse command-line options
while getopts ":i:g:o:t:c" opt; do
    case $opt in
        i) INPUT="$OPTARG";;
        g) GENOME="$OPTARG";;
        o) OUTPUT="$OPTARG";;
        t) THREADS="$OPTARG";;
        c) CHIP=true;;
        *) echo "Error: $OPTARG is an invalid argument." >&2
           exit 1
    esac
done

if [[ $# -eq 0 ]]; then
    echo "Check input arguments and retry."
    exit 1
fi

# Set alignment tool and extension based on CHIP option
if [[ -n "$CHIP" ]]; then
    alignment_tool="bowtie2"
    alignment_ext="bt2"
else
    alignment_tool="hisat2"
    alignment_ext="ht2"
fi

# Build genome indexes for alignment tool if needed
if [[ -n $GENOME ]]; then
    genome_fasta="$(find -L "$GENOME" -mindepth 1 \( -name "*.fasta" -o -name "*.fa" \))"
    bname="${genome_fasta%%.*}"
    if [[ ! -e "$bname".1."${alignment_ext}" ]]; then
        if [[ -z "$genome_fasta" ]]; then
            echo "Error: No FASTA file found in $GENOME directory to build" >&2
            echo "$alignment_tool indexes." >&2
        else
            echo "Creating ${alignment_tool} indexes."
            ${alignment_tool}-build --threads "$THREADS" "$genome_fasta" "$bname"
        fi
    else
        echo "${alignment_tool} indexes detected."
    fi
fi

# Find FASTQ files
read_files=$(find -L "$INPUT" -mindepth 1 \( -name "*_trimmed.fastq*" -o -name "*_trimmed.fq*" \))
num_files=$(echo "$read_files" | wc -l)

# Check the number of FASTQ files and run alignment tool
if [[ "$num_files" -eq 0 ]]; then
    echo "Error: Trimmed FASTQ files not found in $INPUT directory." >&2
    exit 1
elif [[ "$num_files" -eq 1 ]]; then
    fq_r1=$(echo "$read_files")
    out_prefix="$OUTPUT/$(basename "$(dirname "$fq_r1")")/$(basename "${fq_r1%%.*}")"
    "${alignment_tool}" -x "$bname" -1 "$fq_r1" -S "${out_prefix}_aligned.sam" -p "$THREADS" --very-sensitive
else
    fq_r1=$(echo "$read_files" | grep '_R1')
    fq_r2="${fq_r1/_R1/_R2}"
    out_prefix="$OUTPUT/$(basename "$(dirname "$fq_r1")")/$(basename "${fq_r1%_*}")"
    sam_prefix="${out_prefix/_R1/}"
    "${alignment_tool}" -x "$bname" -1 "$fq_r1" -2 "$fq_r2" -S "${sam_prefix}_aligned.sam" -p "$THREADS" --very-sensitive
fi
