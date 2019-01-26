# NtupleMaker

### Installing
1. Copy the install.sh locally.
2. Specify the CMSSW_release and CMS3Tag (branch name or tag name) you want to use at the top of `install.sh`.
3. `source install.sh` will check out the CMSSW release and NtupleMaker repository, and build everything.

### Running

1. If you're not there already, `cd $CMSSW_BASE/src/CMS3/NtupleMaker/test/`.
2. The main pset is, unsuprisingly, `main_pset.py`
3. Input files, event counts, global tag, output file name, and others are specified as command line arguments
4. Finally, 
  * `cmsRun main_pset.py data=True prompt=True` to run on _prompt_ data
  * `cmsRun main_pset.py data=False year=2018` to run on _FullSim_ MC; note the year is needed for genMaker to pick the right weights
  * `cmsRun main_pset.py fastsim=True year=2018` to run on _FastSim_ MC; note the year is needed for genMaker to pick the right weights
  * There are more options explained inside `main_pset.py`, and examples for different campaigns in the testing scripts.

### Testing
With `./run_tests.sh` (*actually* that's outdated, and you should
use `python py_run_tests.py`), one can test the following campaigns:
* 2016
   * 80X MiniAODv2 Re-reco Data (`/*/*03Feb2017*/MINIAOD`)
   * 80X MiniAODv2 MC (`/*/*RunIISummer16MiniAODv2*/MINIAODSIM`)
   * 80X FastSim MC (`/*/*Spring16Fast*/MINIAODSIM`)
   * 94X MiniAODv3 Re-reco Data (`/*/*17Jul2018*/MINIAOD`)
   * 94X MiniAODv3 MC (`/*/*RunIISummer16MiniAODv3*/MINIAODSIM`)
   * 94X MiniAODv3 FastSim MC (`/*/*Summer16v3Fast*/MINIAODSIM`)
* 2017
   * Re-reco Data (`/*/*31Mar2018*/MINIAOD`)
   * Re-reco Data EraF (`/*/*09May2018*/MINIAOD`)
   * MiniAODv2 MC (`/*/*RunIIFall17MiniAODv2*/MINIAODSIM`)
   * MiniAODv2 FastSim MC (`/*/*Fall17Fast*/MINIAODSIM`)
* 2018
   * Prompt-reco Data (`/*/Run2018D-PromptReco-v2/MINIAOD`)
   * Re-reco Data (`/*/*17Sep2018*/MINIAOD`)
   * MiniAODv1 MC (`/*/*RunIIAutumn18MiniAOD*/MINIAODSIM`)

*When I wrote this sentence, all 13 campaigns worked* -- though, with the caveat(s) in the GitHub issues page.

Also, the test script just runs the ntuples and checks that they don't crash. But you should still check branch outputs.
Adding another function to diff two sets is a todo.

Finally, for your own mental safety, I highly recommend locally downloading files first, rather than using xrootd while testing, via [this script](test/profiling/copy_to_local_hadoop.sh).

### Profiling
Details in [here](test/profiling/README.md).

### ProjectMetis details
To make a tarfile for use with Metis, the current command (for 10X+) to execute after setup is
```bash
mtarfile lib_CMS4_V10-02-03_1025.tar.xz --xz --xz_level 9 -x "ZZMatrixElement/MELA/data/Pdfdata" "*ZZMatrixElement/MELA/data/*.root"
# note, extract with `tar xf` to detect the compression algorithm automatically
```
...ignoring some files due to their filesize.

### Misc notes
* `cmsRun` options (`cmsRun -h`) must be put before the pset or else the VarParsing library intercepts it and crashes
* Don't be alarmed if it takes a while to finish processing the first event (see notes in [here](test/profiling/warmup_cache.sh))
* Locally, 4 threads (`cmsRun -n 4 main_pset.py ...`) gives best results, but be wary that this can sometimes cause crashes in the subJetMaker (deep tagger stuff)
