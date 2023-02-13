#!/bin/bash
set -e

# truncates 4 nextseq fastq file and concatenates to a single file
# options: --I [dir with input files] --name a --S s --R r --size c --O [dir for output file]
# a = file name, s = S #, r = 1 or 2 for R1 or R2, c = # of reads in final file in millions

while [ "$1" != "" ]; do
        case "$1" in
                --I | --i ) inputdir=$2; shift;;
                --O | --o ) outputdir=$2; shift;;
                --name ) name=$2; shift;;
                --S | --s ) s=$2; shift;;
                --R | --r ) r=$2; shift;;
                --size ) size=$2; shift;;
        esac
        shift
done

lines=$((size*1000000))
cat <(head -"$lines" "$inputdir"/"$name"_S"$s"_L001_R"$r"_001.fastq) <(head -"$lines" "$inputdir"/"$name"_S"$s"_L002_R"$r"_001.fastq) <(head -"$lines" "$inputdir"/"$name"_S"$s"_L003_R"$r"_001.fastq) <(head -"$lines" "$inputdir"/"$name"_S"$s"_L004_R"$r"_001.fastq) > "$outputdir"/"$name"_R"$r"_"$size"M.fastq;
