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
from datetime import date
from optparse import OptionParser
from IvyFramework.IvyDataTools.TranslateStringBetweenPythonAndShell import *


class BatchManager:
   def __init__(self):
      # define options and arguments ====================================
      self.parser = OptionParser()

      self.parser.add_option("--csv", type="string", help="CSV file to expand")
      self.parser.add_option("--tree_req", type="string", help="Comma-separated list of trees")
      self.parser.add_option("--data", action="store_true", default=False, help="List data samples")
      self.parser.add_option("--sim", action="store_true", default=False, help="List simulation samples")

      (self.opt,self.args) = self.parser.parse_args()

      optchecks=[
         "csv",
         "tree_req"
         ]
      for theOpt in optchecks:
         if not hasattr(self.opt, theOpt) or getattr(self.opt, theOpt) is None:
            sys.exit("Need to set the --{} option".format(theOpt))

      self.infile = self.opt.csv
      if not os.path.isfile(self.infile):
         sys.exit("File {} does not exist.".format(self.infile))

      if not self.opt.data and not self.opt.sim:
         sys.exit("Need to set either the --data or the --sim option")

      self.opt.tree_req = ''.join(self.opt.tree_req.split()) # Remove all whitespaces first...
      treelist = self.opt.tree_req.split(',') # ...then split by comma
      self.reqDilepton = False
      self.reqDileptonControl = False
      self.reqSingleLepton = False
      self.reqSinglePhoton = False
      for tree in treelist:
         if tree == "Dilepton":
            self.reqDilepton = True
         if tree == "Dilepton_Control":
            self.reqDileptonControl = True
         if tree == "SingleLepton":
            self.reqSingleLepton = True
         if tree == "SinglePhoton":
            self.reqSinglePhoton = True

      self.run()


   def run(self):
      with open(self.infile,"rb") as csvfile:
         csvreader = csv.DictReader(csvfile)
         for row in csvreader:
            strsample = row["Dataset"]
            if strsample.strip().startswith("#"): continue

            isMC = ('MINIAODSIM' in strsample)
            hasDilepton = (row["doDilepton"] == 'true')
            hasDileptonControl = (row["doDilepton_Control"] == 'true')
            hasSingleLepton = (row["doSingleLepton"] == 'true')
            hasSinglePhoton = (row["doSinglePhoton"] == 'true')

            if (not isMC and not self.opt.data) or (isMC and not self.opt.sim):
               continue

            if (not hasDilepton and self.reqDilepton):
               continue
            if (not hasDileptonControl and self.reqDileptonControl):
               continue
            if (not hasSingleLepton and self.reqSingleLepton):
               continue
            if (not hasSinglePhoton and self.reqSinglePhoton):
               continue

            print(strsample)


if __name__ == '__main__':
   batchManager = BatchManager()
