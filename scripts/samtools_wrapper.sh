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
: << EOF
if [[ $# -eq 0 ]]; then
    echo "Check input arguments and retry."
    exit 1
fi

# Check required arguments
if [[ -z "$INPUT" || -z "$STEP" || -z "$OUTPUT" ]]; then
    echo "Error: Missing required arguments." >&2
    exit 1
fi
EOF
# Find SAM/BAM files based on the step
if [[ "$STEP" == "filtering" ]]; then
    file_ext=".sam"
elif [[ "$STEP" == "sorting" || "$STEP" == "mapping" ]]; then
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
    # Filters SAM file
    if [[ "$STEP" == "filtering" ]]; then
        echo "Quality filtering $input_file using samtools" >&2
        if [[ -z "$CHIP" ]]; then
            samtools view -@ "$THREADS" -q "$QUALITY" -f 0x02 -F 0x04 -b "$input_file" -o "${bam_prefix}.bam"
        else
            if samtools view -H "$input_file" | grep -qE '^@PG.*ID:MarkDuplicates'; then
                samtools view -@ "$THREADS" -q "$QUALITY" -f 0x02 -F 0x04 -b "$input_file" -o "${bam_prefix}.bam"
            else
                echo "Error: SAM file must be deduplicated before filtering ChIP-seq data." >&2
                continue
            fi
        fi
    
    # Sorts BAM file by coordinate
    elif [[ "$STEP" == "sorting" ]]; then
        echo "Sorting $input_file by coordinate using samtools" >&2
        if samtools view -H "$input_file" | grep -qE '^@PG.*-f\s*0x02\s*-F\s*0x04'; then
            if [[ -z "$CHIP" ]]; then
                samtools sort -@ "$THREADS" "$input_file" -o "${bam_prefix}_sorted.bam"
                samtools stats -@ "$THREADS" -d -x "${bam_prefix}_sorted.bam" > "${bam_prefix}_stats.txt"
                samtools index -@ "$THREADS" -b "${bam_prefix}_sorted.bam"
            else
                if ! samtools view -H "$input_file" | grep -qE '^@PG.*ID:MarkDuplicates'; then
                    echo "Error: BAM files must be deduplicated and quality filtered before sorting ChIP-seq data." >&2
                    continue
                fi
            fi
        else
            echo "Error: BAM files must be quality filtered before sorting." >&2
            exit 1
        fi
    
    # Finds per-base coverage of BAM file
    elif [[ "$STEP" == "mapping" ]]; then
        if ! samtools view -H "$input_file" | grep -qE '^@HD.*SO:coordinate' || 
           ! samtools view -H "$input_file" | grep -qE '^@PG.*-f\s*0x02\s*-F\s*0x04' ||
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
