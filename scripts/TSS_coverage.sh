#!/usr/bin/env bash

## Created: March 8, 2024
## Updated: March 8, 2024
## Author(s): Todd Lenz, tlenz001@ucr.edu

## Wrapper for running deeptools.

while getopts ":m:t:q:" opt; do
    case $opt in
        m) METADATA="$OPTARG";;
        t) THREADS="$OPTARG";;
        *) echo "Error: '$OPTARG' is an invalid argument."
    esac
done


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


unique_values=$(awk 'NR > 1 { print $3 }' "$METADATA" | sort -u)
for value in $unique_values; do
    bamReads="$(awk -v prefix="$data_dir/" '$3 == $value { print prefix $5 }' "$METADATA")"
    bamControl="$(awk -v prefix="$data_dir/" '$3 == $value { print prefix $7 }' "$METADATA")"
    outdir="$OUTPUT/$(basename "$(dirname "$(readlink -f "$tfq_r1")")")/$(basename "${tfq_r1%*}")"
    {
        samtools merge -@ 18 "$output_dir"/"${value%%.*}"_target.bam "$bamReads"
        samtools index -@ 18 "$output_dir"/"${value%%.*}"_target.bam
        samtools merge -@ 18 "$output_dir"/"${value%%.*}"_control.bam "$bamControl"
        samtools index -@ 18 "$output_dir"/"${value%%.*}"_control.bam
        bamCompare -b1 "$output_dir"/"${value%%.*}"_target.bam \
		-b2 "$output_dir"/"${value%%.*}"_control.bam \
		--outFileName "$output_dir"/"${value%%.*}".bw \
		--scaleFactorsMethod None \
		--normalizeUsing RPKM \
		--operation subtract \
		-bs 1 \
		--ignoreDuplicates \
		-p "$THREADS"
    } 
done
computeMatrix reference-point -S "$data_dir"/*.bw -R "$gtf" -o "$data_dir"/TSS_coverage.matrix -b 500 -a 1500 -bs 1 -p "$THREADS" >> "$log_dir"/TSS_coverage.log 2>&1
plotProfile -m "$data_dir"/TSS_coverage.matrix -o "$fig_dir"/TSS_coverage_profile.pdf --perGroup >> "$log_dir"/TSS_coverage.log 2>&1
