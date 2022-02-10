#!/usr/bin/env bash

SECONDS=0

singularityImage="$1"
macroName="$2"
datagroupName="$3"
inputFile="$4"
outputFile="$5"

if (echo "$inputFile" | grep -q -e "^root://"); then
  gfal-cp "$inputFile" .
  inputFileNew="$(dirname "$datagroupName_$inputFile")"
  mv "$inputFile" "$inputFileNew"
  inputFile="$inputFileNew"
  unset inputFileNew
fi

singularity run "$singularityImage" "$macroName" "$inputFile" "$outputFile"

duration=$SECONDS
echo -e "RUN TIME: $(($duration / 60)) minutes and $(($duration % 60)) seconds"
echo "Done"

