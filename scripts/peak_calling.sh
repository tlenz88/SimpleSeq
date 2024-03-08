#!/usr/bin/env bash

## Created: March, 2024
## Updated: March 4, 2024
## Author(s): Todd Lenz, tlenz001@ucr.edu

## Wrapper for running MACS callpeak.

while getopts ":m:g:q:" opt; do
    case $opt in
        m) METADATA="$OPTARG";;
        g) GENOME="$OPTARG";;
        q) QVALUE="$OPTARG";;
        *) echo "Error: '$OPTARG' is an invalid argument."
    esac
done


############################
## Calculate genome size. ##
############################
function get_chrom_sizes() {
    fa=$(find -L "$1" -mindepth 1 -name "*.fasta" -o -name "*.fa")
    chrom_sizes=$(find -L "$1" -mindepth 1 -name "*.sizes")
    if [[ -n "$chrom_sizes" ]]; then
        gsize="$(awk '{ sum += $2 } END { print sum }' "$chrom_sizes")"
    elif [[ -z "$fa" && -z "$chrom_sizes" ]]; then
        echo "Error: No FASTA file found. Can't generate chromosome sizes file."
        exit 1
    elif [[ -n "$fa" && -z "$chrom_sizes" ]]; then
        echo "Using FASTA file to find chromosome lengths."
        python3 "$SCRIPTS"/create_sizes.py "$fa"
        gsize="$(awk '{ sum += $2 } END { print sum }' "$chrom_sizes")"
    fi
}

if [[ -z "$GENOME" ]]; then
    echo "Error: Genome argument required."
    exit 1
elif [[ ! -e "$GENOME" ]]; then
    echo "Error: Genome folder not found."
    exit 1
else
    genome_dir="$(readlink -f "$GENOME")"
    get_chrom_sizes "$genome_dir"
fi


###########################
## Perform peak calling. ##
###########################
function check_bam_type() {
    if [[ $(samtools view -f 0x1 "$1" | head -n 1 | wc -l) -eq 1 ]]; then
        echo "BAMPE"
    else
        echo "BAM"
    fi
}


function check_alg_type() {
    Factor=$(echo "$1" | awk '{ print $2 }')
    if [[ $Factor =~ ^H2[AB]K || $Factor =~ ^H3K ]]; then
        echo "--broad"
    else
        echo ""
    fi
}


function check_BAM() {
    bamfile="$data_dir/$1"
    if [[ $(samtools view --rf 0x400  "$bamfile" | head -n 1 | wc -l) -eq 0 ]]; then
        echo "The BAM file $(basename $bamfile) has not been PCR deduplicated."
        echo "It is recommended that users run picard MarkDuplicates prior to"
        echo "performing peak calling."
    elif [[ $(samtools view -H "$bamfile" | grep -E '^@PG.*-f\s*0X02\s*-F\s*0X04' | head -n 1 | wc -l) -eq 0 ]]; then
        echo "The BAM file $(basename $bamfile) has not been quality filtered."
        echo "It is recommended that users run samtools view using the args"
        echo "-q 30 -f 0X02 -F 0X04 prior to performing peak calling."
    elif [[ $(samtools view -H "$bamfile" | grep -E '^@HD.*SO:coordinate' | head -n 1 | wc -l) -eq 0 ]]; then
        echo "Error: The BAM file $(basename $bamfile) has not been coordinate"
        echo "sorted, which is required to run macs3 callpeak. Run samtools"
        echo "sort to sort by coordinate (default)."
        exit 1
    fi
}


data_dir="$(dirname "$(readlink -f "$METADATA")")"
while IFS= read -r line; do
    SampleID="$(echo "$line" | awk '{ print $1 }')"
    bamReads="$(echo "$line" | awk '{ print $5 }')"
    bamReads="$(check_BAM "$bamReads")"
    bamControl="$(echo "$line" | awk '{ print $7 }')"
    bamControl="$(check_BAM "$bamControl")"
    bam_format="$(check_bam_type "$bamReads")"
    alg="$(check_alg_type "$line")"
    outdir="$(dirname "$(readlink -f "$bamReads")")"
    macs2 callpeak -t "$bamReads" -c "$bamControl" -f "$bam_format" -g "$gsize" -n "$SampleID" -q "$QVALUE" --outdir "$outdir" -B "$alg"
done < <(tail -n +2 "$METADATA")


while IFS= read -r line; do
    if [[ "$(echo "$line" | awk '{ print NF }')" -eq 7 ]]; then
        if [[ "$(echo "$line" | awk '{ print $1 }')" = "SampleID" ]]; then
            echo -e "$line\tPeaks\tPeakCaller"
        else
            bamReads="$(echo "$line" | awk '{ print $5 }')"
            echo -e "$line\t${bamReads%_*}_peaks.xls\tmacs"
        fi
    else
        echo "$line"
    fi
done < "$METADATA" > "${METADATA%.*}_updated.txt"
mv "${METADATA%.*}_updated.txt" "$METADATA"
