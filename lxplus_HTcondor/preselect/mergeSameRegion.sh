#!/usr/bin/env bash

source ./common_header.sh

fillToArray() {
    local __arrayName="$1"
    local __file=""
    local -a __subArray=()
    shift
    while read -r -d $'\0' __file; do
        __subArray+=( "$__file" )
    done < <(find "$@" -print0)
    declare -a "${__arrayName}=( \"\${${__arrayName}[@]}\" \"\${__subArray[@]}\" )"
}

# originalListSpace="/afs/cern.ch/work/y/yuehshun/private/Projects/ShamrockLee/ctau-proper/lxplus_HTcondor/preselect/ntuple_filelist.tmp"
originalListSpace=./ntuple_filelist
outputSpace=/eos/user/y/yuehshun/ntuple_filelist_output.tmp
outputMergedSpace=/eos/user/y/yuehshun/ntuple_filelist_outputmerged.tmp

subDirRoot="";

clusterID="";

dryRun=0

while [[ "${#queuedArgsArray[@]}" -gt 0 ]]; do
    case "${queuedArgsArray[0]}" in
        --help)
            cat <<END_OF_HELP
END_OF_HELP
            exit 0
            ;;
        --cluster)
            clusterID="${queuedArgsArray[1]}"
            shiftArgs 2
            ;;
        --dry-run)
            dryRun=1
            shiftArgs
            ;;
        --maxdepth)
            maxDepth="${queuedArgsArray[1]}"
            shiftArgs 2
            ;;
        --mindepth)
            minDepth="${queuedArgsArray[1]}"
            shiftArgs 2
            ;;
        --original-list-space)
            originalListSpace="${queuedArgsArray[1]}"
            shiftArgs 2
            ;;
        --output-merged-space)
            outputMergedSpace="${queuedArgsArray[1]}"
            shiftArgs 2
            ;;
        --output-space)
            outputSpace="${queuedArgsArray[1]}"
            shiftArgs 2
            ;;
        --*)
            echo "Unexpected flag ${queuedArgsArray[0]}"
            exit 1
            ;;
    esac
done

declare -a searchCommand=()
searchCommand+=( find "$originalListSpace/$subDirRoot" )
if [[ -n "${minDepth:-}" ]]; then
    searchCommand+=( "-mindepth" "$minDepth" )
fi
if [[ -n "${maxDepth:-}" ]]; then
    searchCommand+=( "-maxdepth" "$maxDepth" )
fi
searchCommand+=( "-name" "*.txt" )
searchCommand+=( "-print0" )

while read -r -d $'\0' pathOriginalList; do
    # Change the prefix and remove the .txt file extension
    subDirPlusPrefSlash="${pathOriginalList:${#originalListSpace}:-4}"
    declare -a filesToMerge=()
    if [[ ! -d "$outputSpace$subDirPlusPrefSlash" && ! ( -L "$outputSpace$subDirPlusPrefSlash" && -d "$(realpath "$outputSpace$subDirPlusPrefSlash")" ) ]]; then
        echo "$outputSpace$subDirPlusPrefSlash is not a directory / a symlink to a directory" >&2
        exit 1
    fi
    fillToArray filesToMerge "$outputSpace$subDirPlusPrefSlash" -mindepth 1 -maxdepth 1 -name "_$clusterID.root"
    wrapDryRun mkdir -p 
    wrapDryRun hadd -f "${outputMergedSpace}${subDirPlusPrefSlash}_merged.root" "${filesToMerge[@]}"
done < <("${searchCommand[@]}")
