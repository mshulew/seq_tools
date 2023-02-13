#!/bin/bash
set -e

# batch process directory of NEB R2 & R3 files > single R2

inDir=$1
outDir=$2

mkdir $outDir -p

for a in {1..4};do
    for file in $inDir/*L00${a}_R2_*;do
        echo "processing ${file} ${file/R2/R3} > ${outDir}/${file##*/}..."
        python3.6 ~/seq_tools/sequoia-express-neb/get-umi.py --read ${file/R2/R3} --umi ${file} --output ${outDir}/${file##*/}
    done
done
