# Miscellaneous NtupleMaker profiling tools

## Time profiling (CPU) - IgProf
### First time setup
* IgProf comes with a nice web GUI which requires some setup. In your `~/public_html/.htaccess` file, ensure you have the following
```
AddHandler cgi-script .cgi .py
Options +ExecCGI
```
and do
```bash
mkdir -p ~/public_html/cgi-bin/data/
cp igprof-navigator.py ~/public_html/cgi-bin/
chmod 755 -R ~/public_html/cgi-bin/
```
Note that permissions for cgi scripts and folders they reside in are _very_ annoying to get right. If you see issues with the web site, try making the permissions on the data files within `data/` 644.

### Running
* Modify the pset to run over the desired number of events. If you run over less than a few thousand, your profiling will be dominated by startup overhead.
* Modify the parameters in `run_igprof.sh` to reproduce your local running. Be sure the path to the pset is correct, and it can't hurt to read the script before running it...
* `./run_igprof.sh`
* Visit the link that gets printed out.

## Time profiling (user time)
The poor man's profiling is to run on many events and calculate the event rate from the `Begin processing` lines. Use the script `print_timing.py` to do this after you've redirected such output into a log file. The script has other features to be discovered.

## Memory profiling (Valgrind)
A good reference for Valgrind with CMSSW is [here](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideValgrindMemcheckParser#valgrind_MemcheckGraph_pl_graphi).

### Instructions
* Modify the pset to run over only a _handful_ of events, as this kind of profiling is quite intensive. Then, put this into the pset:
```python
process.ProfilerService = cms.Service (
        "ProfilerService",
        firstEvent = cms.untracked.int32(3),
        lastEvent = cms.untracked.int32(12),
        paths = cms.untracked.vstring('p1')
        )
```
* Execute the following, tweaking the `cmsRun` statement and arguments to match your local setup
```bash
valgrind --leak-check=yes  cmsRun main_pset.py data=False >& log.txt
valgrindMemcheckParser.pl --preset=prod,-prod1+ log.txt  > memory_profiling.html
cp memory_profiling.html ~/public_html/
```
* Finally, check the output in your browser.
