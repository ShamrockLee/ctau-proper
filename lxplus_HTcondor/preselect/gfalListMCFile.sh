#!/bin/bash

scriptname=`basename $0`
EXPECTED_ARGS=4

userid="syu"

if [ $# -eq $(( EXPECTED_ARGS - 1 )) ]
then
    echo "user ID is set to "$userid
else if [ $# -ne $EXPECTED_ARGS ]
then
    echo "Usage: ./$scriptname remoteDirectoryString string grid userid"
    echo "Example: ./$scriptname QCD NCUGlobal ncu/nchc syu"
    exit 1
else 
    userid=$4
fi
fi

string=$2
grid=$3


if [ "$grid" == "ncu" ]; then
    dpmprefix="root://grid71.phy.ncu.edu.tw//dpm/phy.ncu.edu.tw/home/cms/store/user/"$userid

else if [ "$grid" == "nchc" ]; then
    dpmprefix="root://se01.grid.nchc.org.tw//dpm/grid.nchc.org.tw/home/cms/store/user/"$userid
else
    echo "Invalid grid! You need to modify the script by yourself!"
    exit 1
fi
fi


echo $prefix
echo $dpmprefix


tempfile=dir1.txt ## this file determines the directory name
gfal-ls $dpmprefix | grep -a $1 > $tempfile

iteration=0
lastfile=`cat $tempfile | wc -l`
echo "There are "$lastfile" datasets"

## Now check the subdirectories, if more than one is created, all directories will be included
while [ $iteration -lt $lastfile ]; 
do
    iteration=$(( iteration + 1 ))
    input=(`head -n $iteration $tempfile  | tail -1`)
    echo $input
    output=${input}.txt
    rm -rf $output
    echo "Output file "$output" is cleared"
    dir=$dpmprefix/$input
    file2=dir2\_${iteration}.txt
    gfal-ls $dir > $file2
    ## Now check the subdirectories, if more than one is created, all directories will be included
    iter2=0
    nfiles2=`cat $file2 | wc -l`
    while [ $iter2 -lt $nfiles2 ]; 
    do
	iter2=$(( iter2 + 1 ))
	input2=(`head -n $iter2 $file2  | tail -1`)
	echo $input2
	dir2=$dir/$input2
	file3=dir3\_${iteration}\_${iter2}.txt
	gfal-ls $dir2 > $file3

    ## Now check the subdirectories, if more than one is created, all directories will be included
	iter3=0
	nfiles3=`cat $file3 | wc -l`
	while [ $iter3 -lt $nfiles3 ]; 
	do
	    iter3=$(( iter3 + 1 ))
	    input3=(`head -n $iter3 $file3  | tail -1`)
	    echo $input3
	    dir3=$dir2/$input3
	    file4=dir4\_${iteration}\_${iter2}\_${iter3}.txt
	    gfal-ls $dir3 > $file4

       ## Now List the files
	    iter4=0
	    nfiles4=`cat $file4 | wc -l`
	    while [ $iter4 -lt $nfiles4 ]; 
	    do
		iter4=$(( iter4 + 1 ))
		input4=(`head -n $iter4 $file4  | tail -1`)
		echo $input4
		dir4=$dir3/$input4
		gfal-ls $dir4 ##  This line is added so that the response of gfal-ls will be fast when writing output to a text file
		gfal-ls $dir4 | grep -a $string | awk -v my_var=$dir4 '{print my_var"/"$1}' >> $output
	    done
	
	done

    done

done

rm -rf dir*txt

