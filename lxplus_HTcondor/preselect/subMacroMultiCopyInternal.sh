#!/usr/bin/env bash

SECONDS=0

singularityImage="$1"
macroName="$2"
datagroupName="$3"
outputFile="$4"
inputListFile="$5"

sandboxPath=$(mktemp -d "$(basename -s .sif "$singularityImage")_XXXXXX")

singularity build --sandbox "$sandboxPath" "$singularityImage"

infileDir="infiles_${datagroupName}"
mkdir "$infileDir"

singularity exec "$sandboxPath" parallel -j 4 --retry 5 gfal-cp :::: "inputListFile" ::: "$infileDir/"

singularity exec -v "$sandboxPath" "$macroName" -j 8 "$outputFile" "$infileDir/*.root"

chmod -R u+w "$sandboxPath"
rm -rf "$sandboxPath"

duration="$SECONDS"
echo -e "RUN TIME: $(($duration / 60)) minutes and $(($duration % 60)) seconds"
echo "Done"
