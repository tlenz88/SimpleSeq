#!/usr/bin/env bash

## Created: March 8, 2024
## Updated: March 8, 2024
## Author(s): Todd Lenz, tlenz001@ucr.edu

## Wrapper for running deeptools.

print_help() {
    echo "Usage: TSS_coverage.sh [options]"
    echo ""
    echo "Options:"
    echo "  -m <metadata_file>    Specify the metadata file."
    echo "  -t <threads>          Specify the number of threads to use (default=1)."
    echo "  -g <gtf_file>         Specify the GTF file." 
    echo "  -b <binsize>          Specify the bin size for bamCompare (default=1000)."
    echo "  -h, --help            Show this help message and exit."
}

THREADS=1
BINSIZE=1000

while getopts ":m:t:g:b:h:" opt; do
    case $opt in
        m) METADATA="$OPTARG";; ## ChIP metadata file containing file names/paths
        t) THREADS="$OPTARG";; ## Number of processor threads
        g) GTF="$OPTARG";; ## GTF file containing genes of interest
        b) BINSIZE="$OPTARG";; ## Binsize for plotting
        h|*) 
            print_help
            exit 0
            ;;
        *) echo "Error: '$OPTARG' is an invalid argument."
    esac
done
: '
unique_values=$(awk 'NR > 1 { print $3 }' "$METADATA" | sort -u)
echo "$unique_values" | while IFS= read -r value; do
    bamReads=$(awk -v val="$value" '$3 == val { print $5 }' "$METADATA" | tr ' ' '\n' | sort -u | tr '\n' ' ')
    read -r -a bamReadsArray <<< "$bamReads"
    bamControl=$(awk -v val="$value" '$3 == val { print $7 }' "$METADATA" | tr ' ' '\n' | sort -u | tr '\n' ' ')
    read -r -a bamControlArray <<< "$bamControl"
    outdir=$(dirname "$(echo "$bamReads" | awk '{print $1}')")

    target_bam="${value%%.*}_target.bam"
    control_bam="${value%%.*}_control.bam"

    if [ ! -f "$target_bam" ]; then
        if [ "${#bamReadsArray[@]}" -gt 0 ]; then
            samtools merge -@ "$THREADS" "$target_bam" "${bamReadsArray[@]}"
            samtools index -@ "$THREADS" "$target_bam"
        fi
    else
        echo "$target_bam already exists, skipping merging for target."
    fi

    if [ ! -f "$control_bam" ]; then
        if [ "${#bamControlArray[@]}" -gt 0 ]; then
            samtools merge -@ "$THREADS" "$control_bam" "${bamControlArray[@]}"
            samtools index -@ "$THREADS" "$control_bam"
        fi
    else
        echo "$control_bam already exists, skipping merging for control."
    fi

    bigwig_file="${value%%.*}.bw"
    if [ -f "$target_bam" ] && [ -f "$control_bam" ]; then
        if [ ! -f "$bigwig_file" ]; then
            bamCompare -b1 "$target_bam" \
            -b2 "$control_bam" \
            --outFileName "$bigwig_file" \
            --scaleFactorsMethod None \
            --normalizeUsing RPKM \
            --operation subtract \
            -bs "$BINSIZE" \
            --ignoreDuplicates \
            -p "$THREADS"
        else
            echo "$bigwig_file already exists."
        fi
    else
        echo "Skipping bamCompare because one or both BAM files do not exist."
    fi
done
: '
computeMatrix reference-point -S *.bw -R "$GTF" -o AP2_TSS_coverage.matrix -b 1000 -a 1000 -bs "$BINSIZE" -p "$THREADS"
plotHeatmap -m AP2_TSS_coverage.matrix -o AP2_TSS_coverage_heatmap.pdf
plotProfile -m AP2_TSS_coverage.matrix -o AP2_TSS_coverage_profile.pdf --perGroup
