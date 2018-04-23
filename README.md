# NtupleMaker

### Installing
1. `curl https://raw.githubusercontent.com/cmstas/NtupleMaker/master/install.sh > install.sh`.
2. Specify the CMSSW_release and CMS3Tag (branch name or tag name) you want to use at the top of `install.sh`.
3. `source install.sh` will check out the CMSSW release and NtupleMaker repository, and build everything.

### Running

1. If you're not there already, `cd $CMSSW_BASE/src/CMS3/NtupleMaker/test/`.
2. The main pset is, unsuprisingly, `main_pset.py`, and it supports running on data (PromptReco or a Re-Reco), MC (FullSim or FastSim) with various options.
3. Check out the last block of code in the pset to see what parameters you need to modify. Typically, the parameters to modify are:
```
process.GlobalTag.globaltag = "80X_dataRun2_2016SeptRepro_v7"
process.out.fileName = cms.untracked.string('ntuple.root') # output
process.source.fileNames = cms.untracked.vstring('file:DataDoubleEG2016C.root') # input
process.eventMaker.CMS3tag = cms.string('V08-00-18') # doesn't affect ntupling, only for bookkeeping later on
process.eventMaker.datasetName = cms.string('/DoubleEG/Run2016C-03Feb2017-v1/MINIAOD') # doesn't affect ntupling, only for bookkeeping later on
process.maxEvents.input = cms.untracked.int32(3000) # max number of events; note that crab overrides this to -1
```
4. Finally, 
  * `cmsRun main_pset.py data=True prompt=True` to run on _prompt_ data
  * `cmsRun main_pset.py data=False` to run on _FullSim_ MC
  * `cmsRun main_pset.py fastsim=True` to run on _FastSim_ MC


### Some quickstart parameters
In `install.sh`, use `CMS3Tag=CMS4_V00-00-07` and `CMSSW_release=CMSSW_9_4_0` to run on the RunII 2017 re-reco sample for `/DoubleEG/Run2017F-17Nov2017-v1/MINIAOD`

And paste the following at the end of `main_pset.py`.

```
process.GlobalTag.globaltag = "94X_dataRun2_ReReco_EOY17_v2"
process.out.fileName = cms.untracked.string('ntuple.root') # output
process.source.fileNames = cms.untracked.vstring(['/store/data/Run2017F/DoubleEG/MINIAOD/17Nov2017-v1/60000/EAED912B-F7DE-E711-8E9B-0242AC1C0500.root']) # input
process.eventMaker.CMS3tag = cms.string('CMS4_V00-00-07') # doesn't affect ntupling, only for bookkeeping later on
process.eventMaker.datasetName = cms.string('/DoubleEG/Run2017F-17Nov2017-v1/MINIAOD') # doesn't affect ntupling, only for bookkeeping later on
process.maxEvents.input = cms.untracked.int32(3000) # max number of events; note that crab overrides this to -1
```

Run it with `cmsRun main_pset.py data=True`.

### ProjectMetis details
To make a tarfile for use with Metis, the current string looks something like this
```bash
mtarfile lib_CMS4_V09-04-13_946p1.tar.gz -e $CMSSW_BASE/external/$SCRAM_ARCH/lib/libmxnet_predict.so $CMSSW_BASE/config/toolbox/$SCRAM_ARCH/tools/selected/mxnet_predict.xml
```
due to the presence of the heavy object tagger, which unfortunately requires us to copy .so files :( So messy.
