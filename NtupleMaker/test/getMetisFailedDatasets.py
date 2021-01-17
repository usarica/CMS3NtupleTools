#!/bin/env python

import sys
import imp
import copy
import os
import filecmp
import shutil
import pickle
import math
import pprint
import subprocess
import csv
import json
from datetime import date
from optparse import OptionParser


class BatchManager:
   def __init__(self):
      # define options and arguments ====================================
      self.parser = OptionParser()

      self.parser.add_option("--metis_json", type="string", help="'web_summary.json' or equivalent")
      self.parser.add_option("--outfile", type="string", help="Output file to contain a rucio command")
      self.parser.add_option("--dry", action="store_true", default=False, help="Test run without creation of jobs")

      (self.opt,self.args) = self.parser.parse_args()

      optchecks=[
         "metis_json",
         "outfile"
         ]
      for theOpt in optchecks:
         if not hasattr(self.opt, theOpt) or getattr(self.opt, theOpt) is None:
            sys.exit("Need to set --{} option".format(theOpt))

      self.json_sample_nfiles_list = []
      fjson = open(self.opt.metis_json, "r")
      djs = json.load(fjson)
      for task in djs["tasks"]:
         strsf = str(task["general"]["dataset"])
         nfiles = int(task["general"]["njobs_total"])
         outdir = str(task["general"]["output_dir"])
         self.json_sample_nfiles_list.append([strsf,nfiles,outdir])
      fjson.close()

      self.run()


   def run(self):
      fout = open(self.opt.outfile, "w")
      samples_incomplete = []
      for sample_nfiles_pair in self.json_sample_nfiles_list:
         strsample = sample_nfiles_pair[0]
         ffoutcore = strsample
         ffoutcore = ffoutcore.replace('/MINIAODSIM','')
         ffoutcore = ffoutcore.replace('/MINIAOD','')
         ffoutcore = ffoutcore.lstrip('/')

         print("Checking {}...".format(strsample))

         spath = sample_nfiles_pair[2]
         isIncomplete = False
         if not os.path.isdir(spath):
            samples_incomplete.append(strsample)
            isIncomplete = True
            print(" - Sample directory {} does not exist.".format(spath))
         else:
            filelist = [f for f in os.listdir(spath) if (os.path.isfile(os.path.join(spath, f)) and '.root' in f)]
            nfiles = len(filelist)
            isIncomplete = (nfiles < sample_nfiles_pair[1])
            if isIncomplete:
               print(" - Number of files ({}) < number of total files to run ({})".format(nfiles, sample_nfiles_pair[1]))

         if isIncomplete:
            samples_incomplete.append(strsample)

      for strsample in samples_incomplete:
         fout.write("cms:"+strsample+" ")
      fout.write("\n")
      fout.close()

if __name__ == '__main__':
   batchManager = BatchManager()
