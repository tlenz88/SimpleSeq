#!/usr/bin/env bash

## Created: February 29, 2024
## Updated: April 19, 2024
## Author(s): Todd Lenz, tlenz001@ucr.edu

## Wrapper for running Samtools to perform quality filtering,
## coordinate sorting, indexing and per-base mapping.

# Parse command-line options
while getopts ":i:s:o:q:t:c" opt; do
    case $opt in
        i) INPUT="$OPTARG";;
        s) STEP="$OPTARG";;
        o) OUTPUT="$OPTARG";;
        q) QUALITY="$OPTARG";;
        t) THREADS="$OPTARG";;
        c) CHIP=true;;
        *) echo "Error: '$OPTARG' is an invalid argument." >&2
           exit 1
    esac
done

# Find SAM/BAM files based on the step
if [[ "$STEP" == "filtering" ]]; then
    input_files=$(find -L "$INPUT" -mindepth 1 -name "*.sam")
elif [[ "$STEP" == "sorting" || "$STEP" == "mapping" ]]; then
    input_files=$(find -L "$INPUT" -mindepth 1 -name "*.bam")
else
    echo "Error: Invalid step specified." >&2
    exit 1
fi

# Process SAM/BAM files based on the step
for input_file in "${input_files[@]}"; do
    bam_prefix="$OUTPUT/$(basename "$(dirname "$input_file")")/$(basename "${input_file%%.*}")"
    # Filters SAM file
    if [[ "$STEP" == "filtering" ]]; then
        if [[ -z "$CHIP" ]]; then
            samtools view -@ "$THREADS" -q "$QUALITY" -f 0X02 -F 0X04 -b "$input_file" -o "${bam_prefix}.bam"
        else
            if samtools view -H "$input_file" | grep -qE '^@PG.*ID:MarkDuplicates'; then
                samtools view -@ "$THREADS" -q "$QUALITY" -f 0X02 -F 0X04 -b "$input_file" -o "${bam_prefix}.bam"
            else
                echo "Error: SAM file must be deduplicated before filtering ChIP-seq data." >&2
                exit 1
            fi
        fi
    # Sorts BAM file by coordinate
    elif [[ "$STEP" == "sorting" ]]; then
        if samtools view -H "$input_file" | grep -E '^@PG.*-f\s*0X02\s*-F\s*0X04'; then
            if [[ -n "$CHIP" ]]; then
                if ! samtools view -H "$input_file" | grep -qE '^@PG.*ID:MarkDuplicates'; then
                    echo "Error: BAM files must be deduplicated and quality filtered before sorting ChIP-seq data." >&2
                    exit 1
                fi
            fi
            samtools sort -@ "$THREADS" "$input_file" -o "${bam_prefix}_sorted.bam"
            samtools stats -@ "$THREADS" -d -x "${bam_prefix}_sorted.bam" > "${bam_prefix}_stats.txt"
            samtools index -@ "$THREADS" -b "${bam_prefix}_sorted.bam"
        else
            echo "Error: BAM files must be quality filtered before sorting." >&2
            exit 1
        fi
    # Finds per-base coverage of BAM file
    elif [[ "$STEP" == "mapping" ]]; then
        if ! samtools view -H "$input_file" | grep -qE '^@HD.*SO:coordinate' || 
           ! samtools view -H "$input_file" | grep -qE '^@PG.*-f\s*0X02\s*-F\s*0X04' ||
           ! samtools view -H "$input_file" | grep -qE '^@PG.*ID:MarkDuplicates'; then
            echo "Error: BAM files must be deduplicated, quality filtered, and coordinate sorted prior to mapping." >&2
        else
            bed_file="$OUTPUT/$(basename "$(dirname "$input_file")")/$(basename "${input_file%%.*}")"
            samtools depth -a -o "${bed_file}.bed" "$input_file"
        fi
    else
        echo "Error: '$STEP' is an invalid step." >&2
        exit 1
    fi
done