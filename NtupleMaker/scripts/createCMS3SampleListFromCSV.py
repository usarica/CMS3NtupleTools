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
from IvyFramework.IvyDataTools.cmseostools import listFiles


class BatchManager:
   def __init__(self):
      # define options and arguments ====================================
      self.parser = OptionParser()

      self.parser.add_option("--csv", type="string", help="CSV file to expand")
      self.parser.add_option("--outfile", type="string", help="Output file to write")
      self.parser.add_option("--filterfile", type="string", help="List of files to select (useful for rerunning jobs manually)")
      self.parser.add_option("--method", type="string", default="dbs", help="Method to list the data files")
      self.parser.add_option("--options", type="string", default=None, help="Other options specific to each method")
      self.parser.add_option("--nfiles", type="int", default=-1, help="Limit on the number of files per process")
      self.parser.add_option("--ninputsperjob", type="int", default=1, help="Number of input files per job")

      (self.opt,self.args) = self.parser.parse_args()

      optchecks=[
         "csv",
         "outfile"
         ]
      for theOpt in optchecks:
         if not hasattr(self.opt, theOpt) or getattr(self.opt, theOpt) is None:
            sys.exit("Need to set --{} option".format(theOpt))

      if hasattr(self.opt, "filterfile") and getattr(self.opt, "filterfile") is not None:
         self.filterfile = self.opt.filterfile
      else:
         self.filterfile = None

      self.infile = self.opt.csv
      self.outfile = self.opt.outfile

      self.run()


   def run(self):
      filterlist = []
      if self.filterfile:
         with open(self.filterfile,"rb") as filterfile:
            for filterline in filterfile:
               filterline=filterline.lstrip()
               filterline=filterline.rstrip()
               filterline=filterline.replace(',','')
               filterline=filterline.replace('"','')
               if filterline.startswith('#'): continue
               filterlist.append(filterline)

      firstLine = True
      indices = []
      with open(self.outfile,"wb") as outfile:
         with open(self.infile,"rb") as csvfile:
            csvreader = csv.reader(csvfile)
            for row in csvreader:
               rowstr = ','.join(row)
               rowstr.lstrip()
               if firstLine:
                  for el in row:
                     el=el.lstrip()
                     el=el.rstrip()
                     el=el.replace('#','')
                     indices.append(el)
                  firstLine = False
               elif rowstr.startswith('#'):
                  continue
               elif (len(row)==0):
                  continue
               else:
                  slines = []
                  strsample=row[0]
                  strsample = strsample.lstrip()
                  strsample = strsample.rstrip()
                  row[-1] = row[-1].split('#')[0].rstrip()

                  ffoutcore = strsample
                  ffoutcore = ffoutcore.replace('/MINIAODSIM','')
                  ffoutcore = ffoutcore.replace('/MINIAOD','')
                  ffoutcore = ffoutcore.lstrip('/')
                  condorffout = ffoutcore
                  ffoutcore = ffoutcore.replace('/','_')

                  print("Checking {}".format(strsample))
                  filelist = listFiles(
                     sample = strsample,
                     path = self.opt.method,
                     rec = True,
                     other_options = self.opt.options
                     )
                  if len(filterlist)>0:
                     tmplist = []
                     for fl in filelist:
                        if fl in filterlist:
                           tmplist.append(fl)
                     filelist = tmplist

                  groupedfilelist = []
                  if self.opt.ninputsperjob > 0:
                     ifilecount=0
                     infilestr=""
                     for ff in filelist:
                        if infilestr:
                           infilestr = "{},{}".format(infilestr,ff)
                        else:
                           infilestr=ff
                        ifilecount = ifilecount+1
                        if ifilecount == self.opt.ninputsperjob:
                           groupedfilelist.append(infilestr)
                           infilestr=""
                           ifilecount=0
                  else:
                     groupedfilelist = filelist

                  index_ff=0
                  for ff in groupedfilelist:
                     stroutlist=[]
                     ffout = "{}_{}.root".format(ffoutcore,index_ff)
                     for ix in range(len(row)):
                        if ix == 0:
                           stroutlist.append('inputs={}'.format(ff))
                           stroutlist.append('output={}'.format(ffout))
                        elif ix == len(row)-1:
                           stroutlist.append(row[ix])
                        elif indices[ix] == "condoroutdir":
                           stroutlist.append('{}={}/{}'.format(indices[ix],row[ix],condorffout))
                        else:
                           if row[ix]:
                              stroutlist.append('{}={}'.format(indices[ix],row[ix]))
                     strout = " ".join(stroutlist)
                     print(strout)
                     outfile.write(strout+'\n')
                     index_ff = index_ff+1
                     if self.opt.nfiles>0 and index_ff == self.opt.nfiles:
                        break




if __name__ == '__main__':
   batchManager = BatchManager()
