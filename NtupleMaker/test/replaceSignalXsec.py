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
import re
from datetime import date
from optparse import OptionParser

class BatchManager:
   def __init__(self):
      # define options and arguments ====================================
      self.parser = OptionParser()

      self.parser.add_option("--csv", type="string", help="CSV file to expand")
      self.parser.add_option("--xseccsv", type="string", help="CSV file to read xsec")

      (self.opt,self.args) = self.parser.parse_args()

      optchecks=[
         "csv"
         ]

      for theOpt in optchecks:
         if not hasattr(self.opt, theOpt) or getattr(self.opt, theOpt) is None:
            sys.exit("Need to set the --{} option".format(theOpt))

      self.infile = self.opt.csv
      if not os.path.isfile(self.infile):
         sys.exit("File {} does not exist.".format(self.infile))
      self.outfile = self.infile.replace(".csv","_new.csv")

      self.xseccsv = None
      if hasattr(self.opt, "xseccsv") and getattr(self.opt, "xseccsv") is not None:
         self.xseccsv=self.opt.xseccsv

      self.run()


   def run(self):
      with open(self.outfile,"w") as outfile:
         csvwriter = csv.writer(outfile)
         with open(self.infile,"rb") as csvfile:
            csvreader = csv.DictReader(csvfile)
            fields = csvreader.fieldnames
            print(','.join(fields))
            csvwriter.writerow(fields)
            for row in csvreader:
               outrow=[]
               dset = row["#Dataset"]
               xsecline = row["xsec"]
               if dset and (dset.startswith("#/") or dset.startswith("/")):
                  xsec = str(row["xsec"])
                  xsecline=xsec

                  if self.xseccsv is not None:
                     with open(self.xseccsv,'r') as xsecfile:
                        xsecreader = csv.DictReader(xsecfile)
                        for xsecrow in xsecreader:
                           xsdset = xsecrow["#Dataset"]
                           xsrep = str(xsecrow["xsec"])
                           if xsdset == dset:
                              xsecline=xsrep
                              break
                  else:
                     matchres = re.search("_M[0-9]*_",dset)
                     massval=None
                     if matchres is not None:
                        matchstr = matchres.group()
                        massval = matchstr
                        massval = massval.replace("_M","")
                        massval = massval.replace("_","")

                     refdir="/home/users/usarica/work/GenStudies/Offshell2020_Gridpacks/"
                     year=""
                     short_process=""
                     if "ZZTo" in dset:
                        refdir = refdir + "ZZ2L2Nu"
                     else:
                        refdir = refdir + "WW2L2Nu"
                     if "RunIISummer16" in dset:
                        year="2016"
                     if "RunIIFall17" in dset:
                        year="2017"
                     if "RunIIAutumn18" in dset:
                        year="2018"
                     if "GluGlu" in dset and "minlo" in dset.lower():
                        short_process="GGH_minloHJJ"
                     elif "GluGlu" in dset:
                        short_process="GGH"
                     elif "VBF" in dset:
                        short_process="VBFH"
                     elif "Wminus" in dset:
                        short_process="WminusH"
                     elif "Wplus" in dset:
                        short_process="WplusH"
                     else:
                        short_process="ZH"
                     refdir = refdir+'/'+year
                     refdir = refdir+'/'+short_process+'/'+massval
                     fname = refdir+"/pwg-stat.dat"

                     with open(fname,"rb") as reffile:
                        for line in reffile:
                           if "total (btilde+remnants) cross section in pb" in line:
                              line=line.replace("total (btilde+remnants) cross section in pb","")
                              largs=line.split()
                              xsecline=largs[0]
                              break

                  print "{}: {} -> {}".format(dset,row["xsec"],xsecline)

               for key in fields:
                  if key == "xsec":
                     outrow.append(xsecline)
                  else:
                     outrow.append(row[key])
               csvwriter.writerow(outrow)


if __name__ == '__main__':
   batchManager = BatchManager()
