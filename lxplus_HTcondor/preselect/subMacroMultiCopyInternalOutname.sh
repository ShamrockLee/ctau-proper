#!/usr/bin/env bash

SECONDS=0

sendSpace="$1"
singularityImage="$2"
macroName="$3"
datagroupName="$4"
clusterId="$5"
inputListFile="$6"

splittedSpace=ntuple_filelist_splitted.tmp
outputSpace=/eos/user/y/yuehshun/ntuple_filelist_output.tmp
outputFile="$(dirname "$outputSpace${inputListFile:${#splittedSpace}}")/$(basename -s .txt "$inputListFile")_$clusterId.root"

# sandboxPath=$(mktemp -u "$(basename -s .sif "$singularityImage")_XXXXXX")

# singularity build --sandbox "$sandboxPath" "$singularityImage"

export X509_USER_PROXY=/afs/cern.ch/user/y/yuehshun/private/grid.proxy
echo "X509_USER_PROXY=$X509_USER_PROXY"

infileDir="infiles_${datagroupName}"
mkdir "$infileDir"

# singularity exec --bind /pool/condor --bind /afs "$sandboxPath" parallel -j 4 --retry 5 gfal-copy -v "{}" "$infileDir/" :::: "$sendSpace/$inputListFile"
# shellcheck disable=SC2046
xrdcp -v --parallel 4 --retry 5 $(cat "$sendSpace/$inputListFile") "$infileDir/"
ls -s "$infileDir"

# singularity exec --bind /pool/condor "$sandboxPath" "$macroName" -j 8 "$outputFile" "$infileDir/*.root"
singularity exec --bind /pool/condor --bind /eos "$singularityImage" "$macroName" -vj 8 "$outputFile" "$infileDir/*.root"

# chmod -R u+w "$sandboxPath"
# rm -rf "$sandboxPath"

duration="$SECONDS"
echo -e "RUN TIME: $((duration / 60)) minutes and $((duration % 60)) seconds"
echo "Done"
