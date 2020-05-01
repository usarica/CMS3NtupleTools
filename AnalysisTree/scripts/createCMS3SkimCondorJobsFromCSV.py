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
from CMSDataTools.AnalysisTree.TranslateStringBetweenPythonAndShell import *
from CMSDataTools.AnalysisTree.eostools import listFiles


class BatchManager:
   def __init__(self):
      # define options and arguments ====================================
      self.parser = OptionParser()

      self.parser.add_option("--csv", type="string", help="CSV file to expand")
      self.parser.add_option("--date", type="string", help="Skim date")
      self.parser.add_option("--production_tag", type="string", help="Production tag")
      self.parser.add_option("--production_dir", type="string", help="Production directory")
      self.parser.add_option("--nfilesperjob_MC", type="int", default=35, help="Approximate number of input files per MC skim job")
      self.parser.add_option("--nfilesperjob_data", type="int", default=60, help="Approximate number of input files per data skim job")
      self.parser.add_option("--recreate", action="store_true", default=False, help="Recreate the job directories")

      (self.opt,self.args) = self.parser.parse_args()

      optchecks=[
         "csv",
         "date",
         "production_tag",
         "production_dir"
         ]
      for theOpt in optchecks:
         if not hasattr(self.opt, theOpt) or getattr(self.opt, theOpt) is None:
            sys.exit("Need to set --{} option".format(theOpt))

      self.infile = self.opt.csv
      self.extracmd = ""
      if self.opt.recreate:
         self.extracmd += " recreate"

      self.run()


   def run(self):
      cmdlinebase="\"<strSampleSet>\",\"<period>\",\"<prodVersion>\",\"<strdate>\",<ichunk>,<nchunks>,<doDilepton>,<doDilepton_Control>,<doSingleLepton>,<doGJets>"

      with open(self.infile,"rb") as csvfile:
         csvreader = csv.DictReader(csvfile)
         cmdlist = []
         for row in csvreader:
            strsample = row["Dataset"]
            if strsample.strip().startswith("#"): continue

            isMC = ('MINIAODSIM' in strsample)
            nfilesperjob = self.opt.nfilesperjob_data
            if isMC:
               nfilesperjob = self.opt.nfilesperjob_MC

            ffoutcore = strsample
            ffoutcore = ffoutcore.replace('/MINIAODSIM','')
            ffoutcore = ffoutcore.replace('/MINIAOD','')
            ffoutcore = ffoutcore.lstrip('/')
            cmdline = cmdlinebase
            cmdline = cmdline.replace("<strSampleSet>",strsample)
            cmdline = cmdline.replace("<period>",row["Period"])
            cmdline = cmdline.replace("<prodVersion>",self.opt.production_tag)
            cmdline = cmdline.replace("<strdate>",self.opt.date)
            cmdline = cmdline.replace("<doDilepton>",row["doDilepton"])
            cmdline = cmdline.replace("<doDilepton_Control>",row["doDilepton_Control"])
            cmdline = cmdline.replace("<doSingleLepton>",row["doSingleLepton"])
            cmdline = cmdline.replace("<doGJets>",row["doGJets"])

            spath = self.opt.production_dir+'/'+self.opt.production_tag+'/'+ffoutcore
            print "=========="
            print "Checking {}:{}".format(strsample,spath)
            filelist = [f for f in os.listdir(spath) if (os.path.isfile(os.path.join(spath, f)) and '.root' in f)]
            nfiles = len(filelist)
            nchunks = 0
            if nfilesperjob > 0:
               nchunks = nfiles/nfilesperjob
               if nfiles > nchunks*nfilesperjob:
                  nchunks = nchunks+1
               if nchunks == 1:
                  nchunks = 0
            for ichunk in range(0,max(1,nchunks)):
               cmdstr = cmdline
               cmdstr = cmdstr.replace("<ichunk>",str(ichunk))
               cmdstr = cmdstr.replace("<nchunks>",str(nchunks))
               print cmdstr
               cmdlist.append(cmdstr)
            print "=========="
         # Configure jobs
         for cmdarg in cmdlist:
            os.system(
               "submitCMS3AnalysisProduction.sh{} script=produceSkims.cc function=produceSkims arguments=\'{}\' date=skimProduction_{}".format(
                  self.extracmd, cmdarg, self.infile.replace('.csv','')
               )
            )



if __name__ == '__main__':
   batchManager = BatchManager()
