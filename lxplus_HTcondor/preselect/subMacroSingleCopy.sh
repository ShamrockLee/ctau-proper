#!/usr/bin/env bash

SECONDS=0

singularityImage="$1"
macroName="$2"
datagroupName="$3"
outputFile="$4"
inputFile="$5"

infileDir="infiles_${datagroupName}"
mkdir "$infileDir"
if (echo "$inputFile" | grep -q -e "^root://"); then
  gfal-cp "$inputFile" "infileDir/"
  inputFileNew="${infileDir}/$(basename "$inputFile")"
  inputFile="$inputFileNew"
  unset inputFileNew
fi

singularity exec "$singularityImage" -v "$macroName" -j 8 "$outputFile" "$inputFile"

duration="$SECONDS"
echo -e "RUN TIME: $(($duration / 60)) minutes and $(($duration % 60)) seconds"
echo "Done"
