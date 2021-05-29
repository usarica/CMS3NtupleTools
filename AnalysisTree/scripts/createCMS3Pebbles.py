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
from datetime import date
import argparse
import subprocess
import commands
import multiprocessing as mp


def makeDirectory(subDirName):
   if not os.path.exists(subDirName):
      cmd = 'mkdir -p ' + subDirName
      status, output = commands.getstatusoutput(cmd)
      if status != 0:
         print('makeDirectory: Error in creating',subDirName,'.')
         sys.exit()


def run(args):
   currentdir = os.getcwd()
   CMSSWBASE = os.getenv("CMSSW_BASE")
   SCRAMARCH = os.getenv("SCRAM_ARCH")
   batchscript = CMSSWBASE + '/bin/{}/runCMS3Pebbles.condor.sh'.format(SCRAMARCH)

   if not os.path.isfile(batchscript):
      sys.exit("Script {} does not exist. Exiting...".format(batchscript))

   jobmaindir=currentdir+"/output/Pebbles_{}_{}".format(args.npebbles,args.date)
   makeDirectory(jobmaindir+"/Logs")
   tarfile="cms3analysistree.tar"
   if not os.path.isfile(jobmaindir+'/'+tarfile):
      os.system("cd {}; createCMS3AnalysisTreeTarball.sh; cd -;".format(jobmaindir))

   currentCMSSWBASESRC = CMSSWBASE+"/src/" # Need the trailing '/'
   currendir_noCMSSWsrc = currentdir.replace(currentCMSSWBASESRC,'')

   singularityver = "cms:rhel6"
   if "slc7" in SCRAMARCH:
      singularityver = "cms:rhel7"

   scriptargs = {
      "batchScript" : batchscript,

      "CMSSWVERSION" : os.getenv("CMSSW_VERSION"),
      "SCRAMARCH" : SCRAMARCH,
      "SUBMITDIR" : currendir_noCMSSWsrc,
      "TARFILE" : tarfile,
      "home" : os.path.expanduser("~"),
      "uid" : os.getuid(),
      "outDir" : jobmaindir,
      "outLog" : "log_job",
      "errLog" : "err_job",
      "QUEUE" : "vanilla",
      "SINGULARITYVERSION" : singularityver,
      "REQMEM" : "2048M",
      "JOBFLAVOR" : "tomorrow"
   }

   scriptcontents = """
universe={QUEUE}
+DESIRED_Sites="T2_US_UCSD"
executable              = {batchScript}
arguments               = {CMSSWVERSION} {SCRAMARCH} {SUBMITDIR} {TARFILE}
Initialdir              = {outDir}
output                  = Logs/{outLog}.$(ClusterId).$(ProcId).txt
error                   = Logs/{errLog}.$(ClusterId).$(ProcId).err
log                     = Logs/$(ClusterId).$(ProcId).log
request_memory          = {REQMEM}
+JobFlavour             = "{JOBFLAVOR}"
x509userproxy           = {home}/x509up_u{uid}
#https://www-auth.cs.wisc.edu/lists/htcondor-users/2010-September/msg00009.shtml
periodic_remove         = JobStatus == 5
transfer_executable=True
transfer_input_files    = {TARFILE}
transfer_output_files = ""
+Owner = undefined
+project_Name = "cmssurfandturf"
notification=Never
should_transfer_files = YES
when_to_transfer_output = ON_EXIT_OR_EVICT
Requirements = (HAS_SINGULARITY=?=True) || (regexp("(uaf-[0-9]{{1,2}}|uafino)\.", TARGET.Machine) && !(TARGET.SlotID>(TotalSlots<14 ? 3:7) && regexp("uaf-[0-9]", TARGET.Machine)))
+SingularityImage = "/cvmfs/singularity.opensciencegrid.org/cmssw/{SINGULARITYVERSION}"

"""
   scriptcontents = scriptcontents.format(**scriptargs)

   ijob=0
   make_new=True
   condorScriptFile=None
   for it in range(0,args.npebbles):
      if make_new:
         condorScriptDir = "{}/pebbles_{}".format(jobmaindir, ijob)
         makeDirectory(condorScriptDir)
         condorScriptName = condorScriptDir + "/condor.sub"
         print(" - Opening {}...".format(condorScriptName))
         condorScriptFile = open(condorScriptName,'w')
         condorScriptFile.write(scriptcontents)
         make_new = False
      condorScriptFile.write("queue\n")
      if ((it+1)%100)==0 or it==args.npebbles-1:
         ijob = ijob+1
         make_new = True
         condorScriptFile.close()


if __name__ == "__main__":
   parser = argparse.ArgumentParser()
   parser.add_argument("--date", type=str, help="Tag for the pebbles", required=True)
   parser.add_argument("--npebbles", help="Number of pebbles", type=int, required=True)
   args = parser.parse_args()
   run(args)
