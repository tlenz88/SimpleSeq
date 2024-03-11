#!/usr/bin/env bash

## Created: February 22, 2024
## Updated: March 1, 2024
## Author(s): Todd Lenz, tlenz001@ucr.edu

## Wrapper for running Cutadapt to remove low quality bases, adapter sequences
## and short reads from raw fastq files.

# Parse command-line options
while getopts ":i:g:o:q:t:" opt; do
    case $opt in
        i) INPUT="$OPTARG";;
        g) GENOME="$OPTARG";;
        o) OUTPUT="$OPTARG";;
        q) QUALITY="$OPTARG";;
        t) THREADS="$OPTARG";;
        *) echo "Error: $OPTARG is an invalid argument." >&2
           exit 1
    esac
done

# Find file containing adapter sequences
adapters=$(find -L "$GENOME" -mindepth 1 -name "*adapters*")
if [[ -z "$adapters" ]]; then
    echo "Error: No adapter files found in $GENOME directory." >&2
    exit 1
fi

# Check the number of FASTQ files and run cutadapt
read_files=$(find -L "$INPUT" -mindepth 1 \( -name "*.fastq*" -o -name "*.fq*" \))
num_files=$(echo "$read_files" | wc -l)
if [[ "$num_files" -eq 0 ]]; then
    echo "Error: No FASTQ files found in $INPUT directory." >&2
    exit 1
elif [[ "$num_files" -eq 1 ]]; then
    fq_r1=$(echo "$read_files")
    out_r1="$OUTPUT/$(basename "$(dirname "$fq_r1")")/$(basename "${fq_r1%%.*}")"
    cutadapt -j "$THREADS" -a file:"$adapters" -o "${out_r1}_trimmed.fastq.gz" -q "$QUALITY" -m 25 "$fq_r1"
else
    fq_r1=$(echo "$read_files" | grep '_R1')
    fq_r2="${fq_r1/_R1/_R2}"
    out_r1="$OUTPUT/$(basename "$(dirname "$fq_r1")")/$(basename "${fq_r1%%.*}")"
    out_r2="${out_r1/_R1/_R2}"
    cutadapt -j "$THREADS" -a file:"$adapters" -A file:"$adapters" -o "${out_r1}_trimmed.fastq.gz" -p "${out_r2}_trimmed.fastq.gz" -q "$QUALITY" -m 25 "$fq_r1" "$fq_r2"
fi
