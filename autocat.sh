#!/bin/bash
set -e

# concatenates individual NextSeq files into a single file

while [ "$1" != "" ]; do
        case "$1" in
                --i ) inputdir=$2; shift;;
        esac
        shift
done

# concatenate R1
for file1 in $inputdir/*; do
    if [[ $file1 == *"L001_R1"* ]]; then
        for file2 in $inputdir/*; do
            temp=${file2/L002/L001}
            if [[ $file1 == $temp  ]] && [[ $file1 != $file2 ]]; then
                for file3 in $inputdir/*; do
                    temp=${file3/L003/L001}
                    if [[ $file1 == $temp  ]] && [[ $file1 != $file3 ]]; then
                        for file4 in $inputdir/*; do
                            temp=${file4/L004/L001}
                            if [[ $file1 == $temp  ]] && [[ $file1 != $file4 ]]; then
                                temp=${file1/_S*_L001/}
                                output=${temp/_001/}
                                echo "cat $file1 $file2 $file3 $file4 > $output"
                                cat $file1 $file2 $file3 $file4 > $output
                                break
                            fi
                        done
                        break
                    fi
                done
                break
            fi
        done
    fi
done

# concatenate R2
for file1 in $inputdir/*; do
    if [[ $file1 == *"L001_R2"* ]]; then
        for file2 in $inputdir/*; do
            temp=${file2/L002/L001}
            if [[ $file1 == $temp  ]] && [[ $file1 != $file2 ]]; then
                for file3 in $inputdir/*; do
                    temp=${file3/L003/L001}
                    if [[ $file1 == $temp  ]] && [[ $file1 != $file3 ]]; then
                        for file4 in $inputdir/*; do
                            temp=${file4/L004/L001}
                            if [[ $file1 == $temp  ]] && [[ $file1 != $file4 ]]; then
                                temp=${file1/_S*_L001/}
                                output=${temp/_001/}
                                echo "cat $file1 $file2 $file3 $file4 > $output"
                                cat $file1 $file2 $file3 $file4 > $output
                                break
                            fi
                        done
                        break
                    fi
                done
                break
            fi
        done
    fi
done
