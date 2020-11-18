#!/usr/bin/env fish

function rid_end_slash
    set strIn $argv[1]
    if test (count (string match "*/" $strIn)) -gt 0
        string sub -s 1 -l (math (string length $strIn) - 1) $strIn
    else
        echo $strIn
    end
end

set datagroupName $argv[1]
if test (count $argv) -ge 2
    set redirector $argv[2]
else
    set redirector ""
end
if test (count $argv) -ge 3
    set dirWork $argv[3]
else
    set dirWork "."
end

set redirectorGlobal "root://cms-xrd-global.cern.ch"
if test -z "$redirector"
    set redirector $redirectorGlobal
end
# if test ""(string sub -s (string length $redirector) -l 1 $redirector) = "/"
#     set redirector (string sub -s 1 -l (math (string length $redirector) - 1) $redirector)
# end
set redirector (rid_end_slash (string trim $redirector))


set dirWork (rid_end_slash (string trim $dirWork))
set pathOutput $dirWork"/dataTest_"$datagroupName".txt"
set pathOutputNumbers $dirWork"/fileNumbers_"$datagroupName".txt"
echo pathOutputNumbers: $pathOutputNumbers
set subdirDatasetLists "dataset_list"
set subdirDatasetLists (rid_end_slash (string trim $subdirDatasetLists))
set dirDatasetLists $dirWork"/"$subdirDatasetLists

echo -n "" > $pathOutput
echo -n "" > $pathOutputNumbers

for dataset in (cat $dirDatasetLists"/"$datagroupName"_list")
    if test (count (echo $dataset | grep "_ext[0-9]")) -gt 0
        echo Opt out extended: $dataset
        continue
    end
    set filelistSub (dasgoclient --query="file dataset="$dataset)
    for filePath in $filelistSub
        echo $filePath >> $pathOutput
    end
    echo writing (count $filelistSub) files
    count $filelistSub >> $pathOutputNumbers
end

sed -i -r "s/(^.*\$)/"(echo $redirector | sed "s/\\//\\\\\\//g")\\/\\1"/g" $pathOutput
