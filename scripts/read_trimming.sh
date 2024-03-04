#!/usr/bin/env bash

## Created: February 22, 2024
## Updated: March 1, 2024
## Author(s): Todd Lenz, tlenz001@ucr.edu

## Wrapper for running Cutadapt to remove low quality bases, adapter sequences
## and short reads from raw fastq files.

while getopts ":i:g:o:q:t:" opt; do
    case $opt in
        i) INPUT="$OPTARG";;
        g) GENOME="$OPTARG";;
        o) OUTPUT="$OPTARG";;
        q) QUALITY="$OPTARG";;
        t) THREADS="$OPTARG";;
        *) echo "Error: $OPTARG is an invalid argument."
    esac
done

adapters=$(find -L "$GENOME" -mindepth 1 -name "*adapters*")
fq=$(find -L "$INPUT" -mindepth 1 -name "*.fastq*" -o -name "*.fq*" | wc -l)

if [[ "$fq" -eq 1 ]]; then
    fq_r1=$(find -L "$INPUT" -mindepth 1 -name "*.fastq*" -o -name "*.fq*")
    out_r1="$OUTPUT"/"$(basename "$(dirname "$fq_r1")")"/"$(basename "${fq_r1%%.*}")"
    cutadapt -j "$THREADS" -a file:"$adapters" -o "${out_r1}_trimmed.fastq.gz" -q "$QUALITY" -m 25 "$fq_r1"
elif [[ "$fq" -eq 2 ]]; then
    fq_r1=$(find -L "$INPUT" -mindepth 1 \( -name "*.fastq*" -o -name "*.fq*" \) -and -name "*_R1*")
    fq_r2="${fq_r1/_R1/_R2}"
    out_r1="$OUTPUT"/"$(basename "$(dirname "$fq_r1")")"/"$(basename "${fq_r1%%.*}")"
    out_r2="${out_r1/_R1/_R2}"
    cutadapt -j "$THREADS" -a file:"$adapters" -A file:"$adapters" -o "${out_r1}_trimmed.fastq.gz" -p "${out_r2}_trimmed.fastq.gz" -q "$QUALITY" -m 25 "$fq_r1" "$fq_r2"
elif [[ "$fq" -eq 0 ]]; then
    echo "Error: No FASTQ files found in $INPUT directory."
    exit 1
fi
