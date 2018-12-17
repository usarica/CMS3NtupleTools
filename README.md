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
  * `cmsRun main_pset.py data=False year=2018` to run on _FullSim_ MC; note the year is needed for genMaker to pick the right weights
  * `cmsRun main_pset.py fastsim=True year=2018` to run on _FastSim_ MC; note the year is needed for genMaker to pick the right weights
  * There are more options explained inside `main_pset.py`

### Testing
With `./run_tests.sh`, one can test the following campaigns
* 2016
   * 94X Re-reco Data (`/*/*17Jul2018*/MINIAOD`)
   * 94X (MiniAODv3) Re-reco MC (`/*/*RunIISummer16MiniAODv3*/MINIAODSIM`)
* 2017
   * Re-reco Data (`/*/*31Mar2018*/MINIAOD`)
   * MiniAODv2 MC (`/*/*RunIIFall17MiniAODv2*/MINIAODSIM`)
* 2018
   * Prompt-reco Data (`/*/Run2018D-PromptReco-v2/MINIAOD`)
   * Re-reco Data (`/*/*17Sep2018*/MINIAOD`)
   * MiniAODv1 MC (`/*/*RunIIAutumn18MiniAOD*/MINIAODSIM`)
When I wrote this sentence, they all worked.

### Some quickstart parameters
In this `install.sh`, point to this tag or branch and `CMSSW_release=CMSSW_10_2_5` to run on the RunII 2018 data re-reco sample for `/DoubleMuon/Run2018A-17Sep2018-v2/MINIAOD`

And paste the following at the end of `main_pset.py`.

```
process.GlobalTag.globaltag = "102X_dataRun2_Sep2018Rereco_v1"
process.out.fileName = cms.untracked.string('ntuple.root') # output
process.source.fileNames = cms.untracked.vstring(['/store/data/Run2018A/DoubleMuon/MINIAOD/17Sep2018-v2/00000/7B954B49-BE06-B64C-89DC-F568513B41A3.root']) # input
process.eventMaker.CMS3tag = cms.string('CMS4_V10-02-01') # doesn't affect ntupling, only for bookkeeping later on
process.eventMaker.datasetName = cms.string('/DoubleMuon/Run2018A-17Sep2018-v2/MINIAOD') # doesn't affect ntupling, only for bookkeeping later on
process.maxEvents.input = cms.untracked.int32(3000) # max number of events; note that crab overrides this to -1
```

Run it with `cmsRun main_pset.py data=True prompt=True`. Yes, you did see
`prompt=True` when this is a re-reco data file. `edmDumpEventContent` shows
that all the collections have `RECO` process instead of `PAT` (which is what we
should expect for a re-reco campaign), so we tell the pset this is PromptReco
to use `RECO` instead. 

In an ideal world, we would be free from hardcoded process names
(i.e., don't specify that parameter in all the InputTag objects), and then the last
(and hopefully only) process name will be used. If you choose to sacrifice a few hours
of time, I'll throw in some cookies as a reward.

### ProjectMetis details
To make a tarfile for use with Metis, the current command (for 10X+) to execute after setup is
```bash
mtarfile lib_CMS4_V10-02-03_1025.tar.xz --xz --xz_level 9 -x "ZZMatrixElement/MELA/data/Pdfdata" "*ZZMatrixElement/MELA/data/*.root"
# note, extract with `tar xf` to detect the compression algorithm automatically
```
...ignoring some files due to their filesize.

### Misc notes
* `cmsRun` options (`cmsRun -h`) must be put before the pset or else the VarParsing library intercepts it and crashes
* Locally, 4 threads (`cmsRun -n 4 main_pset.py ...`) gives best results
