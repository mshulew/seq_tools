#!/bin/bash
set -e

# moves all files in directories to a single directory

while [ "$1" != "" ]; do
        case "$1" in
                --d ) destination=$2; shift;;
        esac
        shift
done

for folder in *; do
    for file in $folder/*; do
        mv "$file" "$destination"
    done
done

