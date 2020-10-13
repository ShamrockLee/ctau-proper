#!/bin/bash

#mode
compile=1

if [ $# -lt 2 ]
then
    echo "2 parameters are required"
    echo "$0 macroName inputFile"
    exit 1
fi

# set environment variables
SECONDS=0
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc530
export X509_USER_PROXY=/afs/cern.ch/user/s/shuxiao/private/grid.proxy
echo "X509_USER_PROXY=$X509_USER_PROXY"
# setting cmsenv

# RUNNING
if [[ $1 =~ ".so" ]]
then
    echo "generate script to run share library"
cat > run.C << EOF
void run() {
     gROOT->ProcessLine(".L $1");
     gROOT->ProcessLine("${1%_C*}(\"$2\",\"$3\")");
}
EOF
    echo ""
    cat run.C
    echo ""
    echo "running script..."
    root -b -q run.C

elif [[ $1 =~ ".C" ]] || [[ $1 =~ ".c" ]]  
then
    echo "running macro..."
    if [ $compile -eq 1 ]
    then
        root -b -q ${1}+\(\"$2\",\"$3\"\) 
    else
        root -b -q ${1}\(\"$2\",\"$3\"\) 
    fi
else
    echo "File name extension is wrong"
    echo "Check your macro name!"
fi

duration=$SECONDS
echo -e "RUN TIME: $(($duration / 60)) minutes and $(($duration % 60)) seconds"
echo "Done"


