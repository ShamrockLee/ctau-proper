#/bin/bash

scriptname=`basename $0`
EXPECTED_ARGS=2
if [ $# -ne $EXPECTED_ARGS ]
then
    echo "Usage: $scriptname outputTextFileName rootFileDirPath"
    echo "Example: ./$scriptname JetHT_Run2016B.txt root://se01.grid.nchc.org.tw//dpm/grid.nchc.org.tw/home/cms/store/user/syu/JetHT/crab_JetHT-Run2016B/170205_120657/0000"
    exit 1
fi

fileName=$1
rootFileDirPath=$2
nfile=$( gfal-ls ${rootFileDirPath} | grep -h ".root" | wc -l)
echo "There are ${nfile} root files. "
gfal-ls $rootFileDirPath | grep -h ".root" > $fileName
if [[ ${rootFileDirPath:(-1)} != "/" ]]; then
    sed -i "s#^#/#" $fileName
fi
sed -i "s#^#${rootFileDirPath}#" $fileName
echo "Finish! root file pathes are recorded in ${fileName}"
