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

      self.parser.add_option("--indir", type="string", help="Directory to search for files")
      self.parser.add_option("--outfile", type="string", help="Name of the output file")

      (self.opt,self.args) = self.parser.parse_args()

      self.indir = self.opt.indir

      self.run()


   def run(self):
      with open(self.opt.outfile,'w') as fout:
         fout.write("#Dataset,xsec,xsecerr\n")
         dirs=os.listdir(self.indir)
         dirs.sort()
         for d1 in dirs:
            dd1='/'.join([self.indir,d1])
            subdirs=os.listdir(dd1)
            for d2 in subdirs:
               dd2='/'.join([dd1,d2])
               sfiles=[ '/'.join([dd2,ff]) for ff in os.listdir(dd2) if ff.endswith('.txt') ]

               dset="/{}/{}/MINIAODSIM".format(d1,d2)
               sum_xsec=0.
               sum_dxsec=0.

               for fname in sfiles:
                  with open(fname,"r") as fin:
                     for ll in fin:
                        llist=ll.split()
                        for ipos in range(0,len(llist)):
                           if llist[ipos]=='+-':
                              xsec = float(llist[ipos-1])
                              dxsec = float(llist[ipos+1])
                              sum_xsec += xsec/(dxsec**2)
                              sum_dxsec += 1./(dxsec**2)

               sum_xsec /= sum_dxsec
               sum_dxsec = 1./math.sqrt(sum_dxsec)
               print("{}: {} +- {}".format(dset,sum_xsec,sum_dxsec))
               fout.write("{},{},{}\n".format(dset,sum_xsec,sum_dxsec))


if __name__ == '__main__':
   batchManager = BatchManager()
