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
from datetime import date
from optparse import OptionParser
from IvyFramework.IvyDataTools.TranslateStringBetweenPythonAndShell import *


class BatchManager:
   def __init__(self):
      # define options and arguments ====================================
      self.parser = OptionParser()

      self.parser.add_option("--batchqueue", type="string", help="Batch queue")
      self.parser.add_option("--batchscript", type="string", help="Name of the HTCondor script")
      self.parser.add_option("--tarfile", type="string", help="Name of the tar file to upload")
      self.parser.add_option("--outdir", type="string", help="Name of the output directory")

      self.parser.add_option("--condorsite", type="string", help="Name of the HTCondor site")
      self.parser.add_option("--condoroutdir", type="string", help="Name of the HTCondor output directory")

      self.parser.add_option("--outlog", type="string", help="Name of the output log file")
      self.parser.add_option("--errlog", type="string", help="Name of the output error file")

      self.parser.add_option("--script", type="string", help="Name of the script")
      self.parser.add_option("--fcn", type="string", help="Name of the function")
      self.parser.add_option("--fcnargs", type="string", help="Arguments to the function")

      self.parser.add_option("--required_memory", type="string", default="2048M", help="Required RAM for the job")
      self.parser.add_option("--job_flavor", type="string", default="tomorrow", help="Time limit for job (tomorrow = 1 day, workday = 8 hours, see https://batchdocs.web.cern.ch/local/submit.html#job-flavours for more)")

      self.parser.add_option("--dry", dest="dryRun", action="store_true", default=False, help="Do not submit jobs, just set up the files")

      (self.opt,self.args) = self.parser.parse_args()
      optchecks=[
         "batchqueue",
         "batchscript",
         "tarfile",
         "outdir",
         "condorsite",
         "condoroutdir",
         "outlog",
         "errlog",
         "script",
         "fcn",
         "fcnargs"
      ]
      for theOpt in optchecks:
         if not hasattr(self.opt, theOpt) or getattr(self.opt, theOpt) is None:
            sys.exit("Need to set --{} option".format(theOpt))

      if self.opt.outdir.startswith("./"):
         self.opt.outdir = self.opt.outdir.replace(".",os.getcwd(),1)

      if not os.path.isfile(self.opt.batchscript):
         print("Batch script does not exist in current directory, will search for CMSSW_BASE/bin")
         if os.path.isfile(os.getenv("CMSSW_BASE")+"/bin/"+os.getenv("SCRAM_ARCH")+"/"+self.opt.batchscript):
            self.opt.batchscript = os.getenv("CMSSW_BASE")+"/bin/"+os.getenv("SCRAM_ARCH")+"/"+self.opt.batchscript
            print("\t- Found the batch script")
         else:
            sys.exit("Batch script {} does not exist. Exiting...".format(self.opt.batchscript))

      if not os.path.isfile(self.opt.script):
         sys.exit("Script {} does not exist. Exiting...".format(self.opt.script))

      for theOpt in optchecks:
         print("Option {}={}".format(theOpt,getattr(self.opt, theOpt)))

      self.submitJobs()


   def produceCondorScript(self):
      currentdir = os.getcwd()
      currentCMSSWBASESRC = os.getenv("CMSSW_BASE")+"/src/" # Need the trailing '/'
      currendir_noCMSSWsrc = currentdir.replace(currentCMSSWBASESRC,'')
      if self.opt.fcnargs is not None:
         self.opt.fcnargs = translateFromPythonToShell(self.opt.fcnargs)

      scramver = os.getenv("SCRAM_ARCH")
      singularityver = "cms:rhel6"
      if "slc7" in scramver:
         singularityver = "cms:rhel7"

      scriptargs = {
         "home" : os.path.expanduser("~"),
         "uid" : os.getuid(),
         "batchScript" : self.opt.batchscript,
         "CONDORSITE" : self.opt.condorsite,
         "CONDOROUTDIR" : self.opt.condoroutdir,
         "outDir" : self.opt.outdir,
         "outLog" : self.opt.outlog,
         "errLog" : self.opt.errlog,
         "QUEUE" : self.opt.batchqueue,
         "SINGULARITYVERSION" : singularityver,
         "CMSSWVERSION" : os.getenv("CMSSW_VERSION"),
         "SCRAMARCH" : scramver,
         "SUBMITDIR" : currendir_noCMSSWsrc,
         "TARFILE" : self.opt.tarfile,
         "SCRIPT" : self.opt.script,
         "FCN" : self.opt.fcn,
         "FCNARGS" : self.opt.fcnargs,
         "REQMEM" : self.opt.required_memory,
         "JOBFLAVOR" : self.opt.job_flavor
      }

      scriptcontents = """
universe={QUEUE}
+DESIRED_Sites="T2_US_UCSD"
executable              = {batchScript}
arguments               = {CMSSWVERSION} {SCRAMARCH} {SUBMITDIR} {TARFILE} {SCRIPT} {FCN} {FCNARGS} {CONDORSITE} {CONDOROUTDIR}
Initialdir              = {outDir}
output                  = {outLog}.$(ClusterId).$(ProcId).txt
error                   = {errLog}.$(ClusterId).$(ProcId).err
log                     = $(ClusterId).$(ProcId).log
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
Requirements = (HAS_SINGULARITY=?=True) && !(( regexp("(mh-epyc7662-1)\..*",TARGET.Machine) || regexp("(mh-epyc7662-5)\..*",TARGET.Machine) || regexp("(mh-epyc7662-6)\..*",TARGET.Machine) || regexp("(mh-epyc7662-9)\..*",TARGET.Machine) || regexp("(mh-epyc7662-10)\..*",TARGET.Machine) || regexp("(sdsc-84)\..*",TARGET.Machine) || regexp("(sdsc-3)\..*",TARGET.Machine) || regexp("(cabinet-0-0-29)\..*",TARGET.Machine) || regexp("(cabinet-0-0-23)\..*",TARGET.Machine) || regexp("(cabinet-0-0-21)\..*",TARGET.Machine) || regexp("(cabinet-11-11-3)\..*",TARGET.Machine) )=?=True)
+SingularityImage = "/cvmfs/singularity.opensciencegrid.org/cmssw/{SINGULARITYVERSION}"


queue

"""
      scriptcontents = scriptcontents.format(**scriptargs)

      self.condorScriptName = "condor.sub"
      condorScriptFile = open(self.opt.outdir+"/"+self.condorScriptName,'w')
      condorScriptFile.write(scriptcontents)
      condorScriptFile.close()


   def submitJobs(self):
      self.produceCondorScript()

      jobcmd = "cd {}; condor_submit {}; cd -".format(self.opt.outdir, self.condorScriptName)
      if self.opt.dryRun:
         print("Job command: '{}'".format(jobcmd))
      else:
         ret = os.system( jobcmd )



if __name__ == '__main__':
   batchManager = BatchManager()
