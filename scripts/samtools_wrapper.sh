#!/usr/bin/env bash

## Created: February 29, 2024
## Updated: May 21, 2024
## Author(s): Todd Lenz, tlenz001@ucr.edu

## Wrapper for running Samtools to perform quality filtering,
## coordinate sorting, indexing and per-base mapping.

# Parse command-line options
while getopts ":i:s:o:q:t:c" opt; do
    case $opt in
        i) INPUT="${OPTARG%/}";;
        s) STEP="$OPTARG";;
        o) OUTPUT="${OPTARG%/}";;
        q) QUALITY="$OPTARG";;
        t) THREADS="$OPTARG";;
        c) CHIP=true;;
        *) echo "Error: '$OPTARG' is an invalid argument." >&2
           exit 1;;
    esac
done

if [[ $# -eq 0 ]]; then
    echo "Check input arguments and retry." >&2
    exit 1
fi

# Check required arguments
if [[ -z "$INPUT" || -z "$STEP" || -z "$OUTPUT" ]]; then
    echo "Error: Missing required arguments." >&2
    exit 1
fi

# Find SAM/BAM files based on the step
if [[ "$STEP" == "filtering" ]]; then
    file_ext=".sam"
elif [[ "$STEP" == "sorting" || "$STEP" == "deduplication" || "$STEP" == "mapping" ]]; then
    file_ext=".bam"
else
    echo "Error: Invalid step specified." >&2
    exit 1
fi

input_files=()
while IFS= read -r -d '' file; do
    input_files+=("$file")
done < <(find -L "$INPUT" -mindepth 1 -name "*$file_ext" -print0)

# Process SAM/BAM files based on the step
for input_file in "${input_files[@]}"; do
    bam_prefix="$OUTPUT/$(basename "$(dirname "$input_file")")/$(basename "${input_file%%.*}")"
    if [[ "$(samtools view -c -f 1 "$input_file")" -gt 0 ]]; then
        samflags="-f 0X02 -F 0X04"
    else
        samflags="-F 0X04"
    fi

    # Quality filters SAM file and converts to BAM
    if [[ "$STEP" == "filtering" ]]; then
        echo "Quality filtering $bam_prefix using samtools." >&2
        echo "samtools view -@ "$THREADS" -q "$QUALITY" "$samflags" -b "$input_file" -o "${bam_prefix}.bam""
        samtools view -@ "$THREADS" -q "$QUALITY" "$samflags" -b "$input_file" -o "${bam_prefix}.bam"

    # Sorts BAM file by coordinate
    elif [[ "$STEP" == "sorting" ]]; then
        echo "Sorting $bam_prefix using samtools." >&2
        if ! samtools view -H "$input_file" | grep -qE "^@PG.*${samflags}"; then
            echo "Error: BAM file must be quality filtered before sorting." >&2
        else
            if [ "$CHIP" != true ]; then
                echo "samtools sort -@ "$THREADS" -o "${bam_prefix}_sorted.bam" "$input_file""
                samtools sort -@ "$THREADS" -o "${bam_prefix}_sorted.bam" "$input_file"
                echo "samtools index -@ "$THREADS" -b "${bam_prefix}_sorted.bam""
                samtools index -@ "$THREADS" -b "${bam_prefix}_sorted.bam"
                echo "samtools stats -@ "$THREADS" -d -x "${bam_prefix}_sorted.bam" > "${bam_prefix}_sorted_stats.txt""
                samtools stats -@ "$THREADS" -d -x "${bam_prefix}_sorted.bam" > "${bam_prefix}_sorted_stats.txt"
            else
                echo "samtools sort -n -@ "$THREADS" -o "${bam_prefix}_sorted.bam" "$input_file""
                samtools sort -n -@ "$THREADS" -o "${bam_prefix}_sorted.bam" "$input_file"
            fi
        fi

    elif [[ "$STEP" == "deduplication" ]]; then
        echo "Deduplicating $bam_prefix using samtools." >&2
        if ! samtools view -H "$input_file" | grep -qE "^@PG.*${samflags}" ||
           ! samtools view -H "$input_file" | grep -qE '^@PG.*samtools sort'; then
            echo "Error: BAM file must be quality filtered and coordinate sorted prior to deduplication." >&2
        else
            echo "samtools fixmate -@ "$THREADS" -m "$input_file" "${bam_prefix}_fixmate.bam""
            samtools fixmate -@ "$THREADS" -m "$input_file" "${bam_prefix}_fixmate.bam"
            echo "samtools sort -@ "$THREADS" -o "${bam_prefix}_fixmate_sorted.bam" "${bam_prefix}_fixmate.bam""
            samtools sort -@ "$THREADS" -o "${bam_prefix}_fixmate_sorted.bam" "${bam_prefix}_fixmate.bam"
            echo "samtools markdup -@ "$THREADS" -r "${bam_prefix}_fixmate_sorted.bam" "${bam_prefix}_fixmate_sorted_dedup.bam""
            samtools markdup -@ "$THREADS" -r "${bam_prefix}_fixmate_sorted.bam" "${bam_prefix}_fixmate_sorted_dedup.bam"
            echo "samtools index -@ "$THREADS" -b "${bam_prefix}_fixmate_sorted_dedup.bam""
            samtools index -@ "$THREADS" -b "${bam_prefix}_fixmate_sorted_dedup.bam"
            echo "samtools stats -@ "$THREADS" -d -x "${bam_prefix}_fixmate_sorted_dedup.bam" > "${bam_prefix}_fixmate_sorted_dedup_stats.txt""
            samtools stats -@ "$THREADS" -d -x "${bam_prefix}_fixmate_sorted_dedup.bam" > "${bam_prefix}_fixmate_sorted_dedup_stats.txt"
        fi

    # Finds per-base coverage of BAM file
    elif [[ "$STEP" == "mapping" ]]; then
        echo "Mapping genome-wide coverage of $bam_prefix using samtools." >&2
        if ! samtools view -H "$input_file" | grep -qE "^@PG.*${samflags}" ||
           ! samtools view -H "$input_file" | grep -qE '^@PG.*samtools sort'; then
            echo "Error: BAM file must be quality filtered and coordinate sorted prior to mapping." >&2
        else
            if [ "$CHIP" != true ] && 
               ! samtools view -H "$input_file" | grep -qE '^@PG.*ID:MarkDuplicates'; then
                echo "Error: BAM file from ChIP-seq must be deduplicated prior to mapping." >&2
            fi
            bed_file="$OUTPUT/$(basename "$(dirname "$input_file")")/$(basename "${input_file%%.*}")"
            samtools depth -a -o "${bed_file}.bed" "$input_file"
        fi
    else
        echo "Error: '$STEP' is an invalid step." >&2
        exit 1
    fi
done
