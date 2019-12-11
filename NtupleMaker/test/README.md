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


Sometimes, you may want to put the submission scripts at a different location. In this case, you should do


```
submitCMS3NtupleProduction.sh infile=samples_Data_2018.txt outdir=[your directory] date=[your subdirectory]
```

This would create [your directory]/[your subdirectory] to put the scripts there.