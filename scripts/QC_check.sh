#!/usr/bin/env bash

## Created: February 22, 2024
## Updated: March 1, 2024
## Author(s): Todd Lenz, tlenz001@ucr.edu

## Wrapper for running FastQC on raw fastq files.

while getopts ":i:o:t:" opt; do
    case $opt in
        i) INPUT="$OPTARG";;
        o) OUTPUT="$OPTARG";;
        t) THREADS="$OPTARG";;
        *) echo "Error: '$OPTARG' is an invalid argument."
    esac
done

fq=$(find -L "$INPUT" -mindepth 1 -name "*.fastq*" -o -name "*.fq*" | wc -l)

if [[ "$fq" -eq 1 ]]; then
    fq_r1=$(find -L "$INPUT" -mindepth 1 -name "*.fastq*" -o -name "*.fq*")
    out_r1="$OUTPUT"/"$(basename "$(dirname "$fq_r1")")"/"$(basename "${fq_r1%%.*}")"
    fastqc -t "$THREADS" "$fq_r1"
    if [[ "$OUTPUT" != "$INPUT" ]]; then
        mv "${fq_r1%%.*}"_fastqc.html "$out_r1"_fastqc.html
    fi
    rm "${fq_r1%%.*}"_fastqc.zip
elif [[ "$fq" -eq 2 ]]; then
    fq_r1=$(find -L "$INPUT" -mindepth 1 \( -name "*.fastq*" -o -name "*.fq*" \) -and -name "*_R1*")
    fq_r2="${fq_r1/_R1/_R2}"
    out_r1="$OUTPUT"/"$(basename "$(dirname "$fq_r1")")"/"$(basename "${fq_r1%%.*}")"
    out_r2="$OUTPUT"/"$(basename "$(dirname "$fq_r2")")"/"$(basename "${fq_r2%%.*}")"
    fastqc -t "$THREADS" "$fq_r1" "$fq_r2"
    if [[ "$OUTPUT" != "$INPUT" ]]; then
        mv "${fq_r1%%.*}"_fastqc.html "$out_r1"_fastqc.html
        mv "${fq_r2%%.*}"_fastqc.html "$out_r2"_fastqc.html
    fi
    rm "${fq_r1%%.*}"_fastqc.zip
    rm "${fq_r2%%.*}"_fastqc.zip
elif [[ "$fq" -eq 0 ]]; then
    echo "Error: No FASTQ files found in $INPUT directory."
    exit 1
fi
