#!/usr/bin/env bash

## Created: February 22, 2024
## Updated: March 1, 2024
## Author(s): Todd Lenz, tlenz001@ucr.edu

## Wrapper for running Hisat2/Bowtie2 on trimmed fastq files.

while getopts ":i:g:o:t:c:" opt; do
    case $opt in
        i) INPUT="$OPTARG";;
        g) GENOME="$OPTARG";;
        o) OUTPUT="$OPTARG";;
        t) THREADS="$OPTARG";;
        c) CHIP=true;;
        *) echo "Error: '$OPTARG' is an invalid argument."
    esac
done

if [[ -n "$CHIP" ]]; then
    alignment_tool="bowtie2"
    alignment_ext="bt2"
else
    alignment_tool="hisat2"
    alignment_ext="ht2"
fi

if [[ -n $GENOME ]]; then
    GENOME="$(readlink -f "$GENOME")"
    fa="$(find -L "$(readlink -f "$GENOME")" -mindepth 1 -name "*.fasta" -o -name "*.fa")"
    bname="${fa%%.*}"
    if [[ $fa = "" && ! -e "$bname.1.${alignment_ext}" ]]; then
        echo "Error: No FASTA found to build ${alignment_tool} indexes."
    elif [[ ! -e "$bname.1.${alignment_ext}" ]]; then
        echo "Creating ${alignment_tool} indexes."
        ${alignment_tool}-build --threads "$THREADS" "$fa" "$bname"
    else
        echo "${alignment_tool} indexes detected."
    fi
fi

fq=$(find -L "$INPUT" -mindepth 1 -name "*_trimmed.fastq*" -o -name "*_trimmed.fq*" | wc -l)

if [[ "$fq" -eq 1 ]]; then
    tfq_r1=$(find -L "$INPUT" -mindepth 1 -name "*_trimmed.fastq*" -o -name "*_trimmed.fq*")
    out="$OUTPUT/$(basename "$(dirname "$(readlink -f "$tfq_r1")")")/$(basename "${tfq_r1%*}")"
    "${alignment_tool}" -x "$bname" -1 "$tfq_r1" -S "${out}_aligned.sam" -p "$THREADS" --very-sensitive
elif [[ $fq -eq 2 ]]; then
    pfq_r1="$(find -L "$INPUT" -mindepth 1 \( -name "*_trimmed.fastq*" -o -name "*_trimmed.fq*" \) -and -name "*_R1*")"
    pfq_r2="${pfq_r1/_R1/_R2}"
    out="$OUTPUT/$(basename "$(dirname "$(readlink -f "$pfq_r1")")")/$(basename "${pfq_r1%_*}")"
    asam="${out/_R1/}"
    "${alignment_tool}" -x "$bname" -1 "$pfq_r1" -2 "$pfq_r2" -S "${asam}_aligned.sam" -p "$THREADS" --very-sensitive
elif [[ "$fq" -eq 0 ]]; then
    echo "Error: Trimmed FASTQ file(s) required for alignment."
    echo "Start pipeline from the trimming step or run Cutadapt"
    echo "on the raw FASTQ file(s)."
    exit 1
fi
