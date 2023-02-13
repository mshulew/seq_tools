#!/bin/bash
set -e

# renames files across 4 directories

for a in {1..13}; do
    sample=$a
    if [ "$a" == 1 ]; then
        library=7
        index=ATTACTCG+AGGCTATA;
    elif [ "$a" == 2 ]; then
        library=8
        index=TCCGGAGA+AGGCTATA;
    elif [ "$a" == 3 ]; then
        library=9
        index=CGCTCATT+AGGCTATA;
    elif [ "$a" == 4 ]; then
        library=10
        index=GAGATTCC+AGGCTATA;
    elif [ "$a" == 5 ]; then
        library=11
        index=ATTCAGAA+GCCTCTAT;
    elif [ "$a" == 6 ]; then
        library=12
        index=GAATTCGT+GCCTCTAT;
    elif [ "$a" == 7 ]; then
        library=13
        index=CTGAAGCT+GCCTCTAT;
    elif [ "$a" == 8 ]; then
        library=14
        index=TAATGCGC+AGGATAGG;
    elif [ "$a" == 9 ]; then
        library=15
        index=TCTCGCGC+AGGATAGG;
    elif [ "$a" == 10 ]; then
        library=16
        index=CGGCTATG+AGGATAGG;
    elif [ "$a" == 11 ]; then
        library=17
        index=TCCGCGAA+AGGATAGG;
    elif [ "$a" == 12 ]; then
        library=18
        index=TAATGCGC+TCAGAGCC;
    elif [ "$a" == 13 ]; then
        library=19
        index=CGGCTATG+TCAGAGCC;
    fi
    echo processing "$dir"/"$index"_S1_L001_R1_001.fastq
    for b in {1..4}; do
        dir=processed_files_"$b";
        if [ -f "$dir"/"$index"_S1_L001_R1_001.fastq ]; then
            mv "$dir"/"$index"_S1_L001_R1_001.fastq "$dir"/LibV"$library"_S"$sample"_L00"$b"_R1_001.fastq;
        fi
        if [ -f "$dir"/"$index"_S1_L001_R2_001.fastq ]; then
         mv "$dir"/"$index"_S1_L001_R2_001.fastq "$dir"/LibV"$library"_S"$sample"_L00"$b"_R2_001.fastq;
        fi
    done
done
