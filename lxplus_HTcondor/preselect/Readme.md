# How to submit lxplus batch jobs to skim data ntuples located at NCU Tier 3 or Taiwan Tier 2
## Compiling of skimTree

```
setenv SCRAM_ARCH slc7_amd64_gcc530  ## tch
export SCRAM_ARCH=slc7_amd64_gcc530  ## bash
cmsrel CMSSW_8_0_27
cd CMSSW_8_0_27/src
cmsenv
git init scripts
cd scripts
git remote add origin https://github.com/syuvivida/xtozh_common
git config core.sparsecheckout true
echo "lxplus_HTcondor/skimming/*" >> .git/info/sparse-checkout
git pull --depth=1 origin 80X_analysis
cd lxplus_HTcondor/skimming/
chmod 755 *sh
```
Set the environment variables and some required directories and files
```
source prepare.sh
```

## Prepare inputfile that contains a list of NCU ntuples at NCU Tier 3 or Taiwan Tier 2

You can find the full list of input files here: https://github.com/syuvivida/xtozh_common/tree/80X_analysis/2016data


If the files are at NCU Tier 3:
```
gfal-ls root://grid71.phy.ncu.edu.tw:1094//dpm/phy.ncu.edu.tw/home/cms/store/user/syu/SingleMuon/
```
Note, you need to include all the files that are in different subfolders!!
An example input file of ntuples from /JetHT/Run2016B-23Sep2016-v3/MINIAOD is "JetMET_Run2016B"

If the files are at Taiwan Tier 2:
```
gfal-ls root://se01.grid.nchc.org.tw//dpm/grid.nchc.org.tw/home/cms/store/user/syu/SingleMuon
```

To prepare the data root files list, you can use [gfalListDataFile.sh](gfalListDataFile.sh) to make it. 
By giving the keyword of the directory name, it will generates a file which records root file pathes. BUT NOTE, if you have directories from the old jobs under the same dataset name, PLEASE REMEMBER to remove them via "gfal-rm -r xxx"!!
```
./gfalListDataFile.sh JetHT NCUGlobal ncu syu
```

To prepare the MC root files list, you can use [gfalListMCFile.sh](gfalListMCFile.sh) to make it. 
By giving the keyword of the directory name, it will generates a file which records root file pathes. BUT NOTE, if you have directories from the old jobs under the same dataset name, PLEASE REMEMBER to remove them via "gfal-rm -r xxx"!!
```
./gfalListMCFile.sh QCD NCUGlobal ncu syu
```


## Submit the jobs at lxplus

### Submit jobs with HTcondor
Before submit the jobs, please check you have run prepare.sh
```
source prepare.sh
```

you will have a grid.proxy at $HOME/private. Then, you can submit jobs if you want to process all input files of a certain dataset at the same time. There are two method to submit the jobs. In the first method, you SHOULD edit the configuration in main.sub and then submit the jobs. In the second method, condor system reads the input file from your submit option so that the file edition is not required. 


```
# I. Change variable 'listFile' in main.sub and run the command
# listFile  = yourDataList
condor_submit main.sub

# II. Run condor_submit with '- append' option
# If you want to run multiple files, use this way
condor_submit main.sub -append "listFile = yourDataList"
```
## Useful commands of HTcondor
### Check status of job
Using the following command to check the status of jobs
```
condor_q
condor_q <ClusterID>
```
### Problem handling
If the status of job is "hold" which means something is wrong, using following command to study the reason.
```
condor_q -analyze <ClusterID.jobID>
```
After fixing it, using following command to go on

```
condor_release <ClusterID>   ## release jobs from hold state
```
Or you can remove jobs
```
condor_rm <ClusterID>        ## used to remove specific job
condor_rm -all               ## used to remove all jobs
```
### Common troubleshooting
```
-- Failed to fetch ads from: <... : bigbirdxx.cern.ch
SECMAN:2007:Failed to end classad message.
```
If you see this error, wait a while and submit your jobs again. This is usually due to heavy load, either on a specific schedd host or the central manager. In exceptional cases this might be caused by a centralized outage causing delays to the system.
#### solution
To avoid the problem of heavy load of machine, please switch remote machine
```
# show the status of machine, check which machine is not busy
condor_status -schedd
# e.g. using bigbird15.cern.ch
export _CONDOR_SCHEDD_HOST=bigbird15.cern.ch
export _CONDOR_CREDD_HOST=bigbird15.cern.ch 
condor_submit main.sub
```
## How to run with other macro?
Please see [guide.md](Guide/guide.md)

the detail of this framework is writen in this file.

## Useful Links
[CERN Batch Service User Guide](http://batchdocs.web.cern.ch/batchdocs/tutorial/introduction.html)

[HTCondor Manuals](http://research.cs.wisc.edu/htcondor/manual/)
