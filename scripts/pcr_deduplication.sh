#!/usr/bin/env bash

## Created: February 22, 2024
## Updated: March 1, 2024
## Author(s): Todd Lenz, tlenz001@ucr.edu

## Wrapper for running picardtools on aligned SAM files to mark PCR duplicates.

while getopts ":i:o:" opt; do
    case $opt in
        i) INPUT="$OPTARG";;
        o) OUTPUT="$OPTARG";;
        *) echo "Error: '$OPTARG' is an invalid argument."
    esac
done

pipedir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
SCRIPTS="$pipedir"/scripts

if ! command -v picard &> /dev/null; then
    ptjar=$(find -L "$SCRIPTS" -mindepth 1 -name "picard.jar")
    if [ -f "$ptjar" ]; then
        usemem=$(($(free -m | grep Mem | awk '{print $7}') / 2))
        dedup_command="java -Xmx'$usemem'M -jar $ptjar"
    else
        echo "Picard package and/or .jar file not found."
        echo "Download picard .jar file and place in the scripts folder"
        echo "or use the provided conda environment."
    fi
else
    dedup_command="picard"
fi

sam=$(find -L "$INPUT" -mindepth 1 -name "*.sam" | wc -l)

if [[ "$sam" -eq 1 ]]; then
    asam=$(find -L "$INPUT" -mindepth 1 -name "*.sam")
    for samfile in $asam; do
        dsam="$OUTPUT/$(basename "$(dirname "$(readlink -f "$samfile")")")/$(basename "${samfile%%.*}")"
        "$dedup_command" MarkDuplicates -I "$samfile" -O "${dsam}_dedup.sam" -M "${dsam}_dedup_metrics.txt" -ASO queryname
    done
elif [[ "$sam" -eq 0 ]]; then
    echo "Error: A SAM file is required for read deduplication."
    exit 1
fi
