#!/usr/bin/env bash

source ./common_header.sh

# # Works only for bash version >4.4
# fillToArray() {
#     local __arrayName="$1"
#     shift
#     readarray -d $'\0' "$__arrayName" < <(find "$@" -print0)
# }

originalListSpace="/afs/cern.ch/work/y/yuehshun/private/Projects/ShamrockLee/ctau-proper/lxplus_HTcondor/preselect/ntuple_filelist"
# originalListSpace=./ntuple_filelist
outputSpace=/eos/user/y/yuehshun/ntuple_filelist_output.tmp
outputMergedSpace=/eos/user/y/yuehshun/ntuple_filelist_outputmerged.tmp

subDirRoot="";

clusterID="";

dryRun=0
force=0

while [[ "${#queuedArgsArray[@]}" -gt 0 ]]; do
    case "${queuedArgsArray[0]}" in
        -h|--help)
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
        -f|--force)
            force=1
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
        --subdir)
            subDirRoot="${queuedArgsArray[1]}"
            shiftArgs 2
            ;;
        --?*)
            echo "Unexpected flag ${queuedArgsArray[0]}"
            exit 1
            ;;
        --)
            echo "Unexpected flag ${queuedArgsArray[0]}"
            exit 1
            ;;
        -??*)
            manageShorthands
            ;;

    esac
done

if [[ -z "$clusterID" ]]; then
    echo "Expected --cluster clusterID" >&2
    exit 1
fi

if [[ -z "$subDirRoot" ]]; then
  echo "subDirRoot is usually not set as empty." >&2
  echo "Set it with --subdir subDirRoot" >&2
fi

declare -a forceArgArray=()
if (( force )); then
    forceArgArray=( "--force" )
fi

declare -a searchCommand=()
# If subDirRoot is empty, it will find "$originalListSpace"
searchCommand+=( find "$originalListSpace${subDirRoot:+/$subDirRoot}" )
if [[ -n "${minDepth:-}" ]]; then
    searchCommand+=( "-mindepth" "$minDepth" )
fi
if [[ -n "${maxDepth:-}" ]]; then
    searchCommand+=( "-maxdepth" "$maxDepth" )
fi
searchCommand+=( "-name" "*.txt" )
searchCommand+=( "-print0" )

while IFS= read -r -d $'\0' pathOriginalList; do
    # Change the prefix and remove the .txt file extension
    subDirPlusPrefSlash="${pathOriginalList:${#originalListSpace}:-4}"
    declare -a filesToMerge=()
    if [[ ! -d "$outputSpace$subDirPlusPrefSlash" && ! ( -L "$outputSpace$subDirPlusPrefSlash" && -d "$(realpath "$outputSpace$subDirPlusPrefSlash")" ) ]]; then
        echo "$outputSpace$subDirPlusPrefSlash is not a directory / a symlink to a directory" >&2
        exit 1
    fi
    filesToMerge=()
    while IFS= read -r -d $'\0' file; do
        filesToMerge+=( "$file" )
    done < <(find "$outputSpace$subDirPlusPrefSlash" -mindepth 1 -maxdepth 1 -name "*_$clusterID.root" -print0)
    wrapDryRun mkdir -p "$(dirname "${outputMergedSpace}${subDirPlusPrefSlash}_merged.root")"
    wrapDryRun hadd "${forceArgArray[@]}" "${outputMergedSpace}${subDirPlusPrefSlash}_merged_$clusterID.root" "${filesToMerge[@]}"
done < <("${searchCommand[@]}")
