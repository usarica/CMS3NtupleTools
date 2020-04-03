# How to run a data file

Say the csv you would like to run is samples_Data_2018.csv. The first thing would be to create the list of samples out of it.
Make sure you have your Grid proxy fresh.

Do

```
createCMS3SampleListFromCSV.py --csv=samples_Data_2018.csv --outfile=samples_Data_2018.txt
```

This creates a txt file with each line corresponding to one ROOT file. If you want to have N ROOT files per job, do

```
createCMS3SampleListFromCSV.py --csv=samples_Data_2018.csv --outfile=samples_Data_2018.txt --ninputsperjob=[N]
```


Then, you need to do

```
submitCMS3NtupleProduction.sh infile=samples_Data_2018.txt
```

to create the condor.sub scripts in output/[YYMMDD].

Afterward, the workflow to submit jobs is the same as resubmitting them:

```
resubmitCMS3NtupleProduction.sh output/[YYMMDD]
```


In order to check the status of the jobs, do

```
checkCMS3NtupleProduction.sh output/[YYMMDD]
```

If the jobs fail consistently and they process multiple input files, you can split them into invidual input files by doing

```
splitCMS3NtupleProduction.sh output/[YYMMDD]
```


Sometimes, you may want to put the submission scripts at a different location. In this case, you should do

```
submitCMS3NtupleProduction.sh infile=samples_Data_2018.txt outdir=[your directory] date=[your subdirectory]
```

This would create [your directory]/[your subdirectory] to put the scripts there.

# Submission with ProjectMetis

1. Checkout ProjectMetis and make sure it is in the PATH (via the `setup.sh` script):
```
cmsenv
git clone https://github.com/aminnj/ProjectMetis
cd ProjectMetis
source setup.sh
cd ..
```
2. Make a tarball for the worker node:
```bash
mtarfile tarball_v1.tar.xz --xz --xz_level 3 -x "ZZMatrixElement/MELA/data/Pdfdata" "*ZZMatrixElement/MELA/data/*.root"
```
3. Submit jobs in a GNU screen with `python submit_jobs.py /home/users/usarica/work/public/for200313/*.csv samples_Data_*.csv --tarfile tarball_v0.tar.xz --tag OFFSHELL_v0` (after appropriate
edits of the arguments).
Note: I recommend using `DO_TEST=True` with a dummy/different `tag` to submit a handful of events
for the first sample in each csv file. If those jobs succeed, switch to the actual `tag` for production and turn off `DO_TEST`.
4. Visit the monitoring page to view progress and output location.
