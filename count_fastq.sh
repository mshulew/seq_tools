#!/bin/bash
set -e

# counts FASTQ entries in a file

while [ "$1" != "" ]; do
        case "$1" in
                --i ) inputfile=$2; shift;;
        esac
        shift
done

fastqlines="$(wc -l "$inputfile" | cut -d ' ' -f1)";
echo $(( fastqlines / 4 ))


