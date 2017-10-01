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
In `install.sh`, use `CMS3Tag=CMS4_V00-00-06` and `CMSSW_release=CMSSW_9_2_8` to run on the RunIISummer17 sample for `/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer17MiniAOD-92X_upgrade2017_realistic_v10-v1/MINIAODSIM`

And paste the following at the end of `main_pset.py`.

```
process.GlobalTag.globaltag = "92X_upgrade2017_realistic_v10"
process.out.fileName = cms.untracked.string('ntuple.root') # output
process.source.fileNames = cms.untracked.vstring(['/store/mc/RunIISummer17MiniAOD/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/92X_upgrade2017_realistic_v10-v1/90000/DAC03321-BD91-E711-BACC-24BE05C636E1.root']) # input
process.eventMaker.CMS3tag = cms.string('CMS4_V00-00-06') # doesn't affect ntupling, only for bookkeeping later on
process.eventMaker.datasetName = cms.string('/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer17MiniAOD-92X_upgrade2017_realistic_v10-v1/MINIAODSIM') # doesn't affect ntupling, only for bookkeeping later on
process.maxEvents.input = cms.untracked.int32(3000) # max number of events; note that crab overrides this to -1
```

Run it with `cmsRun main_pset.py data=False`.
