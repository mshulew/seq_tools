#!/bin/bash
set -e

# create reverse complement of fastq file

outputfile=revcomp.fastq;
while [ "$1" != "" ]; do
        case "$1" in
                --input ) inputfile=$2; shift;;
                --output ) outputfile=$2; shift;;
        esac
        shift
done
if [ -f "$outputfile" ];then
     rm "$outputfile";
fi;
touch "$outputfile";
line_count=0;
while read line; do
        line_count=$(( $line_count + 1 ));
        if [ $line_count == 1 ]; then
                echo "$line" >> "$outputfile";
        elif [ $line_count == 2 ]; then
                for((i=${#line}-1;i>=0;i--));do
                        complement=$( echo "${line:$i:1}" | tr "[ACGTacgt]" "[TGCAtgca]" )        
                        rev="$rev$complement";
                done
                echo "$rev" >> "$outputfile";
                unset rev;
        elif [ $line_count == 3 ]; then
                echo "$line" >> "$outputfile";
        elif [ $line_count == 4 ]; then
                for((i=${#line}-1;i>=0;i--));do
                        rev="$rev${line:$i:1}";
                done
                echo "$rev" >> "$outputfile";
                unset rev;
                line_count=0;
         fi
done <"$inputfile"
