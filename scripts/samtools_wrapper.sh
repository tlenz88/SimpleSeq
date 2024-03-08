#!/usr/bin/env bash

## Created: February 29, 2024
## Updated: March 1, 2024
## Author(s): Todd Lenz, tlenz001@ucr.edu

## Wrapper for running Samtools to perform quality filtering,
## coordinate sorting, indexing and per-base mapping.

while getopts ":i:s:o:q:t:" opt; do
    case $opt in
        i) INPUT="$OPTARG";;
        s) STEP="$OPTARG";;
        o) OUTPUT="$OPTARG";;
        q) QUALITY="$OPTARG";;
        t) THREADS="$OPTARG";;
        *) echo "Error: '$OPTARG' is an invalid argument."
    esac
done

if [[ "$STEP" == "filtering" ]]; then
    sam_count=$(find -L "$INPUT" -mindepth 1 -name "*.sam" | wc -l)
    if [[ "$sam_count" -ge 1 ]]; then
        asam=$(find -L "$INPUT" -mindepth 1 -name "*.sam")
        for samfile in $asam; do
            abam="$OUTPUT/$(basename "$(dirname "$(readlink -f "$samfile")")")/$(basename "${samfile%%.*}")"
            samtools view -@ "$THREADS" -q "$QUALITY" -f 0X02 -F 0X04 -b "$samfile" -o "${abam}.bam"
        done
    elif [[ "$sam_count" -eq 0 ]]; then
        echo "Error: Genome aligned SAM file(s) required for read filtering."
        exit 1
    fi
elif [[ "$STEP" == "sorting" ]]; then
    bam_count=$(find -L "$INPUT" -mindepth 1 -name "*.bam" | wc -l)
    if [[ "$bam_count" -ge 1 ]]; then
        abam=$(find -L "$INPUT" -mindepth 1 -name "*.bam")
        for bfile in $abam; do
            if samtools view -H "$bfile" | grep -E '^@PG.*-f\s*0X02\s*-F\s*0X04'; then
                sbam="$OUTPUT/$(basename "$(dirname "$(readlink -f "$bfile")")")/$(basename "${bfile%%.*}")"
                samtools sort -@ "$THREADS" "$bfile" -o "${sbam}_sorted.bam"
                samtools stats -@ "$THREADS" -d -x "${sbam}_sorted.bam" > "${sbam}_stats.txt"
                samtools index -@ "$THREADS" -b "${sbam}_sorted.bam"
            else
                ((bam_count--))
            fi
        done
        if [[ "$bam_count" -eq 0 ]]; then
            echo "Error: Quality filtered BAM file(s) not found."
            echo "Start pipeline from the filtering step or run:"
            echo "samtools view -q $QUALITY -f 0X02 -F 0X04"
            exit 1
        fi
    elif [[ "$bam_count" -eq 0 ]]; then
        echo "Error: BAM file(s) required for sorting."
        exit 1
    fi
elif [[ "$STEP" == "mapping" ]]; then
    bam_count=$(find -L "$INPUT" -mindepth 1 -name "*.bam" | wc -l)
    if [[ "$bam_count" -ge 1 ]]; then
        sbam=$(find -L "$INPUT" -mindepth 1 -name "*.bam")
        for bfile in $abam; do
            if samtools view -H "$bfile" | grep -E '^@PG.*-f\s*0X02\s*-F\s*0X04' && samtools view -H "$bfile" | grep -E '^@HD.*SO:coordinate'; then
                sbed="$OUTPUT/$(basename "$(dirname "$(readlink -f "$sbam")")")/$(basename "${bfile%%.*}")"
                samtools depth -a -o "${sbed}.bed" "$bfile"
            else
                ((bam_count--))
            fi
        done
        if [[ "$bam_count" -eq 0 ]]; then
            echo "Error: Quality filtered & coordinate sorted BAM file(s) not"
            echo "found. Start pipeline from the filtering step or run samtools"
            echo "view using the args -q $QUALITY -f 0X02 -F 0X04, and then run"
            echo "samtools sort to sort by coordinate (default)."
            exit 1
        fi
    elif [[ "$bam_count" -eq 0 ]]; then
        echo "Error: BAM file(s) required for sorting."
        exit 1
    fi
fi
