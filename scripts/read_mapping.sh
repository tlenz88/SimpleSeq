#!/usr/bin/env bash

## Created: February 22, 2024
## Updated: May 25, 2024
## Author(s): Todd Lenz, tlenz001@ucr.edu

## Wrapper for running HTSeq to map BAM reads to genes.

while getopts ":i:g:o:t:" opt; do
    case $opt in
        i) INPUT="$OPTARG";;
        g) GENOME="$OPTARG";;
        o) OUTPUT="$OPTARG";;
        t) THREADS="$OPTARG";;
        *) echo "Error: '$OPTARG' is an invalid argument."
    esac
done

GENOME="$(readlink -f "$GENOME")"
gff="$(find -L "$GENOME" -mindepth 1 -name "*.gff")"

if [[ -z "$gff" ]]; then
    echo "Error: GFF file not found in $GENOME directory."
    exit 1
fi

if [[ -n $(awk -v search="protein_coding_gene" '$3 == search {print $3}' "$gff") ]]; then
    gene_string="protein_coding_gene"
else
    gene_string="gene"
fi

bam_count=$(find -L "$INPUT" -mindepth 1 -name "*.bam" | wc -l)
if [[ "$bam_count" -ge 1 ]]; then
    sbam=$(find -L "$INPUT" -mindepth 1 -name "*.bam")
    for bfile in $sbam; do
        if samtools view -H "$bfile" | grep -E '^@PG.*-f\s*0x02\s*-F\s*0x04' && samtools view -H "$bfile" | grep -E '^@HD.*SO:coordinate'; then
            stxt="$OUTPUT/$(basename "$(dirname "$(readlink -f "$bfile")")")/$(basename "${bfile%%.*}")"
            htseq-count -s reverse -t "$gene_string" -i ID "$bfile" "$gff" -n "$THREADS" >> "${stxt}_counts.txt"
        else
            ((bam_count--))
        fi
    done
    if [[ "$bam_count" -eq 0 ]]; then
        echo "Error: Quality filtered & coordinate sorted BAM file(s) not."
        echo "found. Start pipeline from the filtering step or run samtools"
        echo "view using the args -q $QUALITY -f 0x02 -F 0x04, and then run"
        echo "samtools sort to sort by coordinate (default)."
        exit 1
    fi
elif [[ "$bam_count" -eq 0 ]]; then
    echo "Error: Quality filtered and coordinate sorted BAM file(s) required"
    echo "for gene mapping."
    exit 1
fi
