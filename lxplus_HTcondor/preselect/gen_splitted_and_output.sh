#!/usr/bin/env bash

set -eu -o pipefail

originalListSpace=ntuple_filelist
splittedListSpace=ntuple_filelist_splitted.tmp
outputSpace=/eos/user/y/yuehshun/ntuple_filelist_output.tmp
outputLogSpace=ntuple_filelist_outputlog.tmp
NO_WRITE_OUTPUT_DIRS=0
if [ ! -e "$(dirname "$outputSpace")" ]; then
    echo "$(dirname "$outputSpace") not found. Not writing output directories"
    NO_WRITE_OUTPUT_DIRS=1
fi

declare -i nFilesPerJob=50

while IFS= read -d $'\0' file; do
    echo "working on filelist $file" >&2
    subdir="${file:$(( ${#originalListSpace} + 1)):-4}"
    splittedListDir="$splittedListSpace/$subdir"
    if [[ -e  "$splittedListDir" ]]; then
        echo "removing:$splittedListDir" >&2
        rm -rf "$splittedListDir"
    fi
    mkdir -vp "$splittedListDir"
    i=1;
    j=0;
    while read line; do
        if [ -n "$line" ]; then
            fileToCreate="$splittedListDir/splitted_$(basename "$splittedListDir")_$i.txt"
            if ! [ -e "$fileToCreate" ]; then
                echo "creating $fileToCreate" >&2
            fi
            echo "$line" >> "$fileToCreate"
            j=$(( j + 1 ))
        fi
        if (( j >= nFilesPerJob )); then
            j=0
            i=$(( i + 1 ))
        fi
    done < "$file"
    if ! (( NO_WRITE_OUTPUT_DIRS )); then
      outputDir="$outputSpace/$subdir"
      [[ -e "$outputDir" ]] || mkdir -vp "$outputDir";
    fi
    outputLogDir="$outputLogSpace/$subdir"
    [[ -e "$outputLogDir" ]] || mkdir -p "$outputLogDir"
    for suf in out err; do
      [[ -e "$outputLogDir/$suf" ]] || mkdir -v "$outputLogDir/$suf"
    done
done < <(find "$originalListSpace" -name "*.txt" -print0)
