#!/usr/bin/env bash

## Created: March 3, 2024
## Updated: March 3, 2024
## Author(s): Todd Lenz, tlenz001@ucr.edu

## Merges the outputs of htseq-count.

count_files=("$@")
read_counts="${count_files[0]}"
cut -f1 "$read_counts" > read_counts.txt
for cf in "${count_files[@]:1}"; do
    paste read_counts.txt <(cut -f2- "$cf") > merged.tmp && mv merged.tmp read_counts.txt
done
header=$(printf "Gene_ID\t%s\n" "$(basename -a "${count_files[@]%%.*}")" | tr '\n' '\t' | sed 's/\t$/\n/')
{ echo "$header"; cat read_counts.txt; } > merged.tmp && mv merged.tmp read_counts.txt
