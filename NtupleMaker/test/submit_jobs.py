import os
import sys
import csv
import glob
import re
from pprint import pprint
import argparse

from metis.Sample import DBSSample, DirectorySample
from metis.CMSSWTask import CMSSWTask
from metis.StatsParser import StatsParser
from metis.Utils import send_email, interruptible_sleep
from metis.Optimizer import Optimizer

DO_TEST = False

# To make tarfile use command:  mtarfile tarball_v1.tar.xz --xz --xz_level 3 -x "ZZMatrixElement/MELA/data/Pdfdata" "*ZZMatrixElement/MELA/data/*.root"

def get_tasks(csvs, tarfile, tag):

    if not os.path.exists(tarfile):
        raise RuntimeError("{} doesn't exist!".format(tarfile))

    scram_arch = os.getenv("SCRAM_ARCH")
    cmssw_version = os.getenv("CMSSW_VERSION")

    # Make list of sample objects, one per CSV line
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
        pset_args = "xsec={} BR={} ".format(sample.info["xsec"], sample.info["efact"]) + sample.info["options"]
        global_tag = re.search("globaltag=(\w+)",pset_args).groups()[0]
        extra = dict()
        if DO_TEST:
            extra["max_jobs"] = 2 # 2 condor jobs per sample... FIXME delete this after testing
            extra["max_nevents_per_job"] = 200 # 200 events per job... FIXME delete this after testing

        # build output directory
        hadoop_user = os.environ.get("GRIDUSER","").strip()  # Set by Metis. Can be different from $USER for some people.
        if not hadoop_user: hadoop_user = os.environ.get("USER")
        part1 = sample.get_datasetname().split("/")[1]
        part2 = sample.get_datasetname().split("/")[2]
        output_dir = "/hadoop/cms/store/user/{}/Offshell_2L2Nu/Ntuples/{}/{}/{}/".format(
                hadoop_user, tag, part1, part2
                )

        task = CMSSWTask(
                sample=sample,
                tarfile=tarfile,
                tag=tag,
                global_tag=global_tag,
                pset_args=pset_args,
                scram_arch=scram_arch,
                cmssw_version=cmssw_version,
                output_name="ntuple.root",
                output_dir=output_dir,
                events_per_output=events_per_output,
                pset="main_pset.py",
                is_tree_output=False,
                dont_check_tree=True,
                # no_load_from_backup=True,
                **extra
                )
        tasks.append(task)

    return tasks

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("csvs", help="csv files with samples", nargs="+")
    parser.add_argument("--tarfile", help="path to tarball", type=str, required=True)
    parser.add_argument("--tag", help="production tag", type=str, required=True)
    args = parser.parse_args()

    tasks = get_tasks(csvs=args.csvs, tarfile=args.tarfile, tag=args.tag)

    total_summary = {}
    optimizer = Optimizer()

    for i in range(10000):

        for task in tasks:
            task.process(optimizer=optimizer)
            total_summary[task.get_sample().get_datasetname()] = task.get_task_summary()

        StatsParser(data=total_summary, webdir="~/public_html/dump/metis_offshell/").do()
        interruptible_sleep(4*3600)
