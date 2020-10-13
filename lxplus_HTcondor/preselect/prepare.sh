#!/bin/bash
for f in *.txt
do
    mkdir -p ${f%.txt}/err
    mkdir -p ${f%.txt}/out
done
echo "preparing directories..."

export X509_USER_PROXY=$HOME/private/grid.proxy
sed -i "s#export X509.*#export X509_USER_PROXY=$X509_USER_PROXY#g" subMacro.sh
echo "execute 'voms-proxy-init --voms cms'"
voms-proxy-init --voms cms

# comment block
<<multi-line-comment-1
exefile=skimTree.C
if [ ! -f $exefile ]
then
    echo "cannot find the target macro"
    echo "please compile your macro"
    echo 'root -b -q -e ".L yourMacro.C++"'
elif ! find . -name ${exefile%.*}_*.so -printf 1 -quit|grep -q 1
then 
    echo "don't find share library"
    echo "compiling..."
    root -b -q -e ".L $exefile++"
fi
multi-line-comment-1

echo "Done"
