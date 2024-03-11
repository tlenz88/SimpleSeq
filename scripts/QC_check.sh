#!/usr/bin/env bash

## Created: February 22, 2024
## Updated: March 10, 2024
## Author(s): Todd Lenz, tlenz001@ucr.edu

## Wrapper for running FastQC on raw fastq files

# Parse command-line options
while getopts ":i:o:t:" opt; do
    case $opt in
        i) INPUT="$OPTARG";;
        o) OUTPUT="$OPTARG";;
        t) THREADS="$OPTARG";;
        *) echo "Error: '$OPTARG' is an invalid argument." >&2
           exit 1
    esac
done

# Find FASTQ files
read_files=$(find -L "$INPUT" -mindepth 1 \( -name "*.fastq*" -o -name "*.fq*" \))
fq_count=$(echo "$read_files" | wc -l)

# Check the number of FASTQ files and run FastQC
case $fq_count in
    0) echo "Error: No FASTQ files found in $INPUT directory." >&2
       exit 1
       ;;
    1) # Single-end reads
       out_dir="$OUTPUT/$(basename "$(dirname "$read_files")")"
       out_file="$(basename "${read_files%%.*}")"
       fastqc -t "$THREADS" "$read_files"
       ;;
    2) # Paired-end reads
       out_dir="$OUTPUT/$(basename "$(dirname "$(echo "$read_files" | head -n 1)")")"
       out_file="$(basename "$(echo "$read_files" | head -n 1)")"
       fastqc -t "$THREADS" $read_files
       ;;
esac

# Move and remove FASTQC output files
if [ "$fq_count" -eq 1 ]; then
    if [ "$OUTPUT" != "$INPUT" ]; then
        mv "$INPUT"/"${out_file}"_fastqc.html "$out_dir/$out_file"_fastqc.html
    fi
    rm "${out_prefix}"_fastqc.zip
elif [ "$fq_count" -eq 2 ]; then
    if [ "$OUTPUT" != "$INPUT" ]; then
        out_file2="${out_file/_R1/_R2}"
        mv "$INPUT"/"${out_file%%.*}"_fastqc.html "$out_dir"/"${out_file%%.*}"_fastqc.html
        mv "$INPUT"/"${out_file2%%.*}"_fastqc.html "$out_dir"/"${out_file2%%.*}"_fastqc.html
    fi
    rm "$INPUT"/"${out_file%%.*}"_fastqc.zip
    rm "$INPUT"/"${out_file2%%.*}"_fastqc.zip
fi
