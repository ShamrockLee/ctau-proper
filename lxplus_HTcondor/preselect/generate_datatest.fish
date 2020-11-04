#!/usr/bin/env fish

set datagroupName $argv[1]
if test (count $argv) -ge 2
    set redirector $argv[2]
else
    set redirector ""
end

set redirectorGlobal "root://cms-xrd-global.cern.ch"
if test -z "$redirector"
    set redirector $redirectorGlobal
end
if test ""(string sub -s (string length $redirector) -l 1 $redirector) = "/"
    set redirector (string sub -s 1 -l (math (string length $redirector) - 1) $redirector)
end

set pathOutput "dataTest_"$datagroupName".txt"
set pathOutputNumbers "fileNumbers_"$datagroupName".txt"
echo pathOutputNumbers: $pathOutputNumbers
set dirDatasetLists "dataset_list"

echo -n "" > $pathOutput
echo -n "" > $pathOutputNumbers

for dataset in (cat $dirDatasetLists"/"$datagroupName"_list")
    set filelistSub (dasgoclient --query="file dataset="$dataset)
    for filePath in $filelistSub
        echo $filePath >> $pathOutput
    end
    echo writing (count $filelistSub) files
    count $filelistSub >> $pathOutputNumbers
end

sed -i -r "s/(^.*\$)/"(echo $redirector | sed "s/\\//\\\\\\//g")\\/\\1"/g" $pathOutput
