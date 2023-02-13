#!/bin/bash
set -e

# removes rRNA reads from fastq

while [ "$1" != "" ]; do
    case "$1" in
        --fastq ) r1=$2; if [[ "--bam" != *"$3"* ]] || [[ "--outDir" != *"$3"* ]] || [[ "--genomes_bases" != *"$3"* ]] ; then r2=$3; fi; shift;;
        --star ) stardir=$2; shift;;
        --genomes_base ) refs=$2; shift;;
        --outDir ) outdir=$2; shift;;
    esac
    shift
done

echo ' '
echo 'Removing ribosomal reads from FASTQ files(s)'
echo ' '

if [ ${#outdir} == 0 ]; then
    echo 'no outDir selected'
    exit
fi

if [ ${outdir: -1} != "/" ]; then
    outputdir=${outputdir}/
fi

if ! [ -d "$outdir" ]; then
    echo "output directory ${outdir} does not exist; making directory"
    mkdir $outdir
fi
        
echo ' '
echo '***********************'
echo "FASTQ files: " ${r1} ${r2}
echo "STAR directory: " ${stardir}
echo "Genome Reference: " ${refs}
echo "Output directory: " ${outdir}
echo '***********************'
echo ' '
echo 'Preparing ribosomal interval bed file...'
tail -1776 ${refs}/hg38/anno_ncbi/hg38_ribosomal_intervals.txt | awk 'BEGIN {FS="\t";OFS="\t"} {print $1,$2,$3,"rRNA",".",$4}' > ${outdir}/ribosomal.bed
#grep -e 'gbkey "rRNA";' ${refs}/hg38/anno_ncbi/hg38_longRNA_annotation.bed > ${outdir}/ribosomal.bed
echo 'Removing ribosomal reads from STAR BAM file...   (can take a few minutes)'
bedtools intersect -a $stardir/star/Aligned.sortedByCoord.out.bam -b ${outdir}/ribosomal.bed -v -s > ${outdir}/depleted.bam
echo 'Generate list of reads remaining after ribosomsal removal...'
samtools view ${outdir}/depleted.bam | cut -f 1 | sed -e 's/^.*_//' | sort | uniq > ${outdir}/reads_remaining
echo 'Generating R1 without ribosomal reads...'
subr1="${r1##*/}"
r1prefix="${subr1%.fastq*}"
zcat $r1 | sed -e 's/\r\$//g' | paste - - - - > ${outdir}/r1
echo "input reads: " $(wc -l ${outdir}/r1 | cut -d' ' -f 1)
grep -w -F -f ${outdir}/reads_remaining ${outdir}/r1 | tr '\t' '\n' > ${outdir}/${r1prefix}_depleted.fastq 
echo "output reads: " $(( $(wc -l ${outdir}/${r1prefix}_depleted.fastq | cut -d' ' -f 1)/4 ))
rm ${outdir}/r1
if [ ${#r2} > 0 ]; then
    echo 'Generating R2 without ribosomal reads...'
    subr2="${r2##*/}"
    r2prefix="${subr2%.fastq*}"
    zcat $r2 | sed -e 's/\r\$//g' | paste - - - - > ${outdir}/r2
    echo "input reads: " $(wc -l ${outdir}/r2 | cut -d' ' -f 1)
    grep -w -F -f ${outdir}/reads_remaining ${outdir}/r2 | tr '\t' '\n' > ${outdir}/${r2prefix}_depleted.fastq 
    echo "output reads: " $(( $(wc -l ${outdir}/${r2prefix}_depleted.fastq | cut -d' ' -f 1)/4 ))
    rm ${outdir}/r2
fi
echo ' '
echo 'Compressing FASTQ files...'
gzip ${outdir}/*.fastq







