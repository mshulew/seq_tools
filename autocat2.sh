#!/bin/bash
set -e

# concatenates individual NextSeq files into a single file

while [ "$1" != "" ]; do
        case "$1" in
                --i ) inputdir=$2; shift;;
                --o ) outputdir=$2; shift;;
        esac
        shift
done

if [ ${outputdir: -1} != "/" ]; then
    outputdir=${outputdir}/
fi

if ! [ -d "$outputdir" ]; then
    echo "output directory ${outputdir} does not exist; making directory"
    mkdir $outputdir
fi
        

for r1file1 in $inputdir/*; do
    if [[ $r1file1 == *"L001_R1"* ]]; then
        r1file1=${r1file1//\/\//\/}
        r1file2=${r1file1/L001/L002}
        r1file3=${r1file1/L001/L003}
        r1file4=${r1file1/L001/L004}
        r2file1=${r1file1/R1/R2}
        r2file2=${r2file1/L001/L002}
        r2file3=${r2file1/L001/L003}
        r2file4=${r2file1/L001/L004}
        basename="${r1file1##*/}"
        basename=${basename/_S*_L001_R1_001/}
        basename=${basename/.fastq.gz/}
        echo "cat $r1file1 $r1file2 $r1file3 $r1file4 > ${outputdir}${basename}_R1.fastq.gz"
        cat $r1file1 $r1file2 $r1file3 $r1file4 > ${outputdir}${basename}_R1.fastq.gz
        echo "cat $r2file1 $r2file2 $r2file3 $r2file4 > ${outputdir}${basename}_R2.fastq.gz"
        cat $r2file1 $r2file2 $r2file3 $r2file4 > ${outputdir}${basename}_R2.fastq.gz
    fi
done


