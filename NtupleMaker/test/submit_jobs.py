import os
import sys
import csv
import glob
import re
from pprint import pprint

from metis.Sample import DBSSample, DirectorySample
from metis.CMSSWTask import CMSSWTask
from metis.StatsParser import StatsParser
from metis.Utils import send_email, interruptible_sleep
from metis.Optimizer import Optimizer

DO_TEST = True
tarfile = "tarball_v0.tar.xz"
tag = "OFFSHELL_v0"
csvs = glob.glob("/home/users/usarica/work/public/for200313/*.csv")

def get_tasks():


    samples = []
    for fname in csvs:
        with open(fname) as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                dataset = row["#Dataset"]
                if dataset.strip().startswith("#"): continue
                sample = DBSSample(dataset=dataset, xsec=row["xsec"], efact=row["BR"])
                sample.info["options"] = row["options"]
                samples.append(sample)
                if DO_TEST:
                    break # only do one sample per file... FIXME delete this after testing

    tasks = []
    for sample in samples:
        isdata = "Run201" in sample.info["dataset"]
        events_per_output = (200e3 if isdata else 200e3)
        pset_args = sample.info["options"]
        global_tag = re.search("globaltag=(\w+)",pset_args).groups()[0]
        extra = dict()
        if DO_TEST:
            extra["max_jobs"] = 2 # 2 condor jobs per sample... FIXME delete this after testing
            extra["max_nevents_per_job"] = 200 # 200 events per job... FIXME delete this after testing
        task = CMSSWTask(
                sample=sample,
                tarfile=tarfile,
                tag=tag,
                global_tag=global_tag,
                pset_args=pset_args,
                scram_arch="slc6_amd64_gcc700",
                cmssw_version="CMSSW_10_2_18",
                output_name="ntuple.root",
                events_per_output=events_per_output,
                pset="main_pset.py",
                is_tree_output=False,
                dont_check_tree=True,
                **extra
                )
        tasks.append(task)

    return tasks

if __name__ == "__main__":

    total_summary = {}
    optimizer = Optimizer()

    for i in range(10000):

        for task in get_tasks():
            task.process(optimizer=optimizer)
            total_summary[task.get_sample().get_datasetname()] = task.get_task_summary()

        StatsParser(data=total_summary, webdir="~/public_html/dump/metis_offshell/").do()
        interruptible_sleep(4*3600)
