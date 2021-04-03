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
import commands
from optparse import OptionParser


class BatchManager:
   def __init__(self):
      # define options and arguments ====================================
      self.parser = OptionParser()

      self.parser.add_option("--date", type="string", help="Tag for the fits")
      self.parser.add_option("--tag", type="string", help="Tag for the workspaces")
      self.parser.add_option("--years", type="string", default="2016,2017,2018", help="Years to run")
      self.parser.add_option("--dry", action="store_true", default=False, help="Test run without creation of jobs")
      self.parser.add_option("--direct_submit", action="store_true", default=False, help="Submut jobs directly")

      (self.opt,self.args) = self.parser.parse_args()

      optchecks=[
         "date",
         "tag"
         ]
      for theOpt in optchecks:
         if not hasattr(self.opt, theOpt) or getattr(self.opt, theOpt) is None:
            sys.exit("Need to set --{} option".format(theOpt))

      self.years = self.opt.years.split(',')

      self.run()


   def makeDirectory(self,subDirName):
      if not os.path.exists(subDirName):
         cmd = 'mkdir -p ' + subDirName
         status, output = commands.getstatusoutput(cmd)
         if status != 0:
            print('makeDirectory: Error in creating',subDirName,'.')
            sys.exit()


   def getVOMSProxy(self):
      gridproxy = None
      if os.getenv("X509_USER_PROXY") is None or not os.getenv("X509_USER_PROXY"):
         gridproxycheckfiles = [
            "{home}/x509up_u{uid}".format(home=os.path.expanduser("~"), uid=os.getuid()),
            "/tmp/x509up_u{uid}".format(uid=os.getuid())
            ]
         for gridproxycheckfile in gridproxycheckfiles:
            if os.path.exists(gridproxycheckfile):
               gridproxy = gridproxycheckfile
               break
      else:
         gridproxy = os.getenv("X509_USER_PROXY")
      if gridproxy is None or not os.path.exists(gridproxy):
         sys.exit("Cannot find a valid grid proxy")
      return gridproxy


   def run(self):
      CMSSWBASE = os.getenv("CMSSW_BASE")
      SCRAMARCH = os.getenv("SCRAM_ARCH")
      batchscript = CMSSWBASE + '/bin/{}/runCMS3TnPFit.condor.sh'.format(SCRAMARCH)
      condorsite="t2.ucsd.edu"

      gridproxy = self.getVOMSProxy()
      grid_user = subprocess.check_output("voms-proxy-info -identity -file={} | cut -d '/' -f6 | cut -d '=' -f2".format(gridproxy), shell=True)
      if not grid_user:
         grid_user = os.environ.get("USER")
      grid_user = grid_user.strip()

      curdir=os.getcwd()
      jobmaindir=curdir+"/output/{}_LeptonTnP_CombineFits".format(self.opt.date)
      self.makeDirectory(jobmaindir)
      os.system("cd {}; createCMS3TnPTarball.sh; cd -;".format(jobmaindir))
      for year in self.years:
         indir=curdir+"/output/LeptonEfficiencies/WSandDCs/{}/{}".format(self.opt.tag, year)
         if not os.path.isdir(indir):
            raise RuntimeError("{} does not exist.".format(indir))
         wsnames=[ xx for xx in os.listdir(indir) if '.tar' in xx ]
         for wsname in wsnames:
            fname=indir+'/'+wsname
            jobdir=jobmaindir+"/{}/{}".format(year,wsname.replace('.tar',''))
            self.makeDirectory(jobdir)
            os.system("ln -sf {} {}/".format(batchscript, jobdir))
            os.system("ln -sf {}/cms3tnpfit.tar {}/".format(jobmaindir, jobdir))
            if os.path.exists(jobdir+"/wstarfile.tar"):
               os.unlink(jobdir+"/wstarfile.tar")
            os.symlink(fname,jobdir+"/wstarfile.tar")
            condoroutdir="/hadoop/cms/store/user/{}/Offshell_2L2Nu/Worker/output/LeptonEfficiencies/DataFits/{}/{}/{}".format(grid_user,self.opt.date,year,wsname)

            jobargs = {
               "batchqueue" : "vanilla",
               "batchscript" : batchscript,
               "tarfile" : "cms3tnpfit.tar",
               "wstarfile" : "wstarfile.tar",
               "outdir" : jobdir,
               "condorsite" : condorsite,
               "condoroutdir" : condoroutdir,
               "outlog" : "Logs/log_job",
               "errlog" : "Logs/err_job",
               "required_memory" : "2048M",
               "job_flavor" : "tomorrow"
            }
            runCmd = "configureCMS3TnPCondorJob.py"
            if not self.opt.direct_submit:
               runCmd = runCmd + " --dry"
            for key,val in jobargs.iteritems():
               runCmd = runCmd + " --{}={}".format(key,val)
            if self.opt.dry:
               print("Job command: '{}'".format(runCmd))
            else:
               os.system(runCmd)



if __name__ == '__main__':
   batchManager = BatchManager()
