#!/bin/env python

import os
import pprint
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


def getVOMSProxy():
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


def run_single(args,grid_user,jobmaindir,batchscript,condorsite,year,indir,wsname):
   fname=indir+'/'+wsname
   jobdir=jobmaindir+"/{}/{}".format(year,wsname.replace('.tar',''))
   if (not args.remake and os.path.exists(jobdir+".tar")) or (args.skip_existing and os.path.exists(jobdir)):
      return

   condoroutdir="/ceph/cms/store/user/{}/Offshell_2L2Nu/Worker/output/LeptonEfficiencies/DataFits/{}/{}/{}".format(grid_user,args.date,year,wsname.replace('.tar',''))
   if args.robust_fit:
      # This means the directory should exist...
      if not os.path.exists(condoroutdir):
         print("Directory {} does not exist.".format(condoroutdir))
         return
      elif not os.path.exists(condoroutdir+"/combined_withSnapshot.root"):
         print("Directory {} does not contain combined_withSnapshot.root.".format(condoroutdir))
         return
      elif os.path.exists(condoroutdir+"/combined_withSnapshot_withRobustFit.root"):
         print("Directory {} already contains combined_withSnapshot_withRobustFit.root.".format(condoroutdir))
         return

   print("Creating {}".format(jobdir))
   makeDirectory(jobdir+"/Logs")
   ret = os.system("ln -sf {} {}/".format(batchscript, jobdir))
   if ret != 0:
      print("Could not link the batch script in {}".format(jobdir))

   ret = os.system("ln -sf {}/cms3tnpfit.tar {}/".format(jobmaindir, jobdir))
   if ret != 0:
      print("Could not link the main tar file in {}".format(jobdir))

   if not args.robust_fit:
      ret = os.system("ln -sf {} {}/wstarfile.tar".format(fname, jobdir))
      if ret != 0:
         print("Could not link the workspace tar file in {}".format(jobdir))

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
   if not args.direct_submit:
      runCmd = runCmd + " --dry"
   for key,val in jobargs.iteritems():
      runCmd = runCmd + " --{}={}".format(key,val)
   if args.robust_fit:
      runCmd = runCmd + " --robust_fit"
   if args.dry:
      print("Job command: '{}'".format(runCmd))
   else:
      os.system(runCmd + " > /dev/null")


def run(args):
   CMSSWBASE = os.getenv("CMSSW_BASE")
   SCRAMARCH = os.getenv("SCRAM_ARCH")
   batchscript = CMSSWBASE + '/bin/{}/runCMS3TnPFit.condor.sh'.format(SCRAMARCH)
   condorsite="t2.ucsd.edu"
   nthreads = min(args.nthreads, mp.cpu_count())

   gridproxy = getVOMSProxy()
   grid_user = subprocess.check_output("voms-proxy-info -identity -file={} | cut -d '/' -f6 | cut -d '=' -f2".format(gridproxy), shell=True)
   if not grid_user:
      grid_user = os.environ.get("USER")
   grid_user = grid_user.strip()

   curdir=os.getcwd()
   jobmaindir=curdir+"/output/{}_LeptonTnP_CombineFits".format(args.date)
   if args.robust_fit:
      jobmaindir = jobmaindir + "_RobustFit"
   makeDirectory(jobmaindir)
   os.system("cd {}; createCMS3TnPTarball.sh; cd -;".format(jobmaindir))

   sub_data = []
   for year in args.years.split(','):
      indir=curdir+"/output/LeptonEfficiencies/WSandDCs/{}/{}".format(args.tag, year)
      if not os.path.isdir(indir):
         raise RuntimeError("{} does not exist.".format(indir))
      wsnames=[ xx for xx in os.listdir(indir) if '.tar' in xx ]
      for wsname in wsnames:
         sub_data.append([year,indir,wsname])

   pool = mp.Pool(nthreads)
   #starmap does not exist until Python3...
   #pool.starmap_async(run_single, [(njdir, refcsubs) for njdir, refcsubs in sub_data])
   [ pool.apply_async(run_single, args=(args, grid_user, jobmaindir, batchscript, condorsite, year, indir, wsname)) for year, indir, wsname in sub_data ]
   pool.close()
   pool.join()


if __name__ == "__main__":
   parser = argparse.ArgumentParser()
   parser.add_argument("--date", type=str, help="Tag for the fits", required=True)
   parser.add_argument("--tag", type=str, help="Tag for the workspaces", required=True)
   parser.add_argument("--nthreads", help="Number of threads", type=int, required=False, default=4)
   parser.add_argument("--years", type=str, default="2016,2017,2018", help="Years to run", required=False)
   parser.add_argument("--remake", action="store_true", help="Remake all jobs", required=False)
   parser.add_argument("--skip_existing", action="store_true", help="Skip existing folders in addition to tar files to save time", required=False)
   parser.add_argument("--robust_fit", action="store_true", help="Recover robust fit", required=False)
   parser.add_argument("--dry", action="store_true", help="Test run without creation of jobs", required=False)
   parser.add_argument("--direct_submit", action="store_true", help="Submut jobs directly", required=False)
   args = parser.parse_args()
   run(args)
