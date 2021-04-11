#!/bin/env python

import os
import sys
import socket
import glob
import re
import subprocess
from pprint import pprint
import argparse
import multiprocessing as mp


def run_single(njdir, refcsubs):
   if not os.path.isdir(njdir):
      os.makedirs(njdir)
   print("Creating new job file {}".format(njdir + "/condor.sub"))
   with open(njdir + "/condor.sub", 'w') as fout:
      for refcsub in refcsubs:
         refjobdir = refcsub.replace("/condor.sub","")
         fout.write("# Job from {}\n".format(os.path.abspath(refjobdir)))
         with open(refcsub, 'r') as fin:
            fout.writelines(fin.readlines())
         fout.write("#####\n\n")


def run(args):
   oldjobdir = args.oldjobdir
   newjobdir = args.newjobdir
   nsubjobs = args.nsubjobs
   reset_dirs = args.reset_dirs
   nthreads = min(args.nthreads, mp.cpu_count())

   if os.path.abspath(newjobdir) == os.getcwd():
      raise RuntimeError("You may not use the new job directory to be the current directory. Please specify something else.")

   if reset_dirs:
      os.system("rm -rf {}".format(newjobdir))
   if not os.path.isdir(newjobdir):
      os.makedirs(newjobdir)

   csubfiles = glob.glob(oldjobdir+"/**/condor.sub")
   noldjobs = len(csubfiles)
   nnewjobs = noldjobs / nsubjobs
   if noldjobs % nsubjobs != 0:
      nnewjobs = nnewjobs + 1

   print("Number of new jobs: ".format(nnewjobs))

   sub_data = []
   job_count = 0
   for ijob in range(nnewjobs):
      njdir = newjobdir + "/job_{}".format(ijob)
      jbegin = job_count
      jend = min(jbegin+nsubjobs, noldjobs)
      job_count = jend
      refcsubs = [ csubfiles[jjob] for jjob in range(jbegin, jend) ]
      sub_data.append([njdir, refcsubs])

   pool = mp.Pool(nthreads)
   #starmap does not exist until Python3...
   #pool.starmap_async(run_single, [(njdir, refcsubs) for njdir, refcsubs in sub_data])
   [ pool.apply_async(run_single, args=(njdir, refcsubs)) for njdir, refcsubs in sub_data ]
   pool.close()
   pool.join()



if __name__ == "__main__":
   parser = argparse.ArgumentParser()
   parser.add_argument("--oldjobdir", help="Unparallelized job directory", type=str, required=True)
   parser.add_argument("--newjobdir", help="Directory for the collected jobs", type=str, required=True)
   parser.add_argument("--nsubjobs", help="Number of subjobs grouped together", type=int, required=True)
   parser.add_argument("--nthreads", help="Number of threads", type=int, required=False, default=4)
   parser.add_argument("--reset_dirs", help="Reset the main collected job directory", action='store_true', required=False)
   args = parser.parse_args()
   run(args)
