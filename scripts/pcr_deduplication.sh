#!/usr/bin/env bash

## Created: February 22, 2024
## Updated: March 10, 2024
## Author(s): Todd Lenz, tlenz001@ucr.edu

## Wrapper for running picardtools on aligned SAM files to mark PCR duplicates.

# Parse command-line options
while getopts ":i:o:" opt; do
    case $opt in
        i) INPUT="$OPTARG";;
        o) OUTPUT="$OPTARG";;
        *) echo "Error: '$OPTARG' is an invalid argument." >&2
           exit 1
    esac
done

# Find script directory
SCRIPTS="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)/scripts"

# Determine the deduplication command
if ! command -v picard &> /dev/null; then
    ptjar=$(find -L "$SCRIPTS" -mindepth 1 -name "picard.jar")
    if [ -f "$ptjar" ]; then
        usemem=$(( $(free -m | awk '/Mem/{print $7}') / 2 ))
        dedup_command="java -Xmx${usemem}M -jar $ptjar"
    else
        echo "Error: Picard package and/or .jar file not found."
        echo "Download the Picard .jar file and place it in the scripts" >&2
        echo "folder or use the provided conda environment." >&2
        exit 1
    fi
else
    dedup_command="picard"
fi

# Find SAM file
sam_file=$(find -L "$INPUT" -mindepth 1 -name "*.sam")
num_files="${#sam_file[@]}"

# Check the number of SAM files and run picard
if [[ "$num_files" -eq 0 ]]; then
    echo "Error: A SAM file is required for PCR deduplication." >&2
    exit 1
elif [[ "$num_files" -eq 1 ]]; then
    out_prefix="$OUTPUT/$(basename "$(dirname "$sam_file")")/$(basename "${sam_file%%.*}")"
    echo "$dedup_command" MarkDuplicates -I "$sam_file" -O "${out_prefix}_dedup.sam" -M "${out_prefix}_dedup_metrics.txt" -ASO queryname
    "$dedup_command" MarkDuplicates -I "$sam_file" -O "${out_prefix}_dedup.sam" -M "${out_prefix}_dedup_metrics.txt" -ASO queryname
else
    echo "Error: Multiple SAM files found in the $INPUT directory." >&2
    echo "Unable to determine proper SAM file for PCR deduplication." >&2
    exit 1
fi
