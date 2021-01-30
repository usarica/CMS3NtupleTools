#!/bin/env python

import os
import sys
import csv
import glob
import re
from pprint import pprint
import argparse


def getNFilesRequired(sname):
   if "GluGlu" in sname or "VBF" in sname or ("ZZTo2L2Nu" in sname and ("Wminus" in sname or "Wplus" in sname)) or "ZH_LNuQQ" in sname or "HZJ_HToWW2L2NuInclusiveZ" in sname:
      return 500
   elif "WminusH_HToZZTo2L2Q" in sname or "WplusH_HToZZTo2L2Q" in sname:
      return 50
   elif "ZH_HTo2L2Q" in sname:
      return 560
   elif "ZH_HTo4Q" in sname:
      if "RunIISummer16" in sname:
         return 2000
      else:
         return 1000
   elif "ZH_HTo2Nu2X" in sname:
      return 2500
   else:
      raise RuntimeError("Cannot identify the number of required files for {}".format(sname))


def getCSVLine(sname):
   line_tpl="{SNAME},/hadoop/cms/store/user/<USER>/Offshell_2L2Nu/Production/<DATE>,-1,-1,{VVMODE},{VVDECAYMODE},{MEFILE},globaltag={GLOBALTAG} nevents=-1 year={YEAR} {METRECIPE}doTrigObjMatching=True triggerListFromFile=OffshellTriggerFilterList_{YEAR}.lst {KFACTOROPT}data=False minNleptons=1 minNphotons=1 includeLJetsSelection=True keepGenParticles=reducedfinalstatesandhardprocesses keepGenJets=True LHEInputTag=source disableDuplicateCheck=True"

   kfactoropt=""
   if "GluGlu" in sname:
      kfactoropt="applyKFactorQCDNLOtoNNLOggVVSig=True " # Space is important

   year=None
   globaltag=None
   metrecipe=""
   vvmode=None
   vvdecaymode=None
   mefile=None

   if "RunIISummer16" in sname:
      year=2016
      globaltag="94X_mcRun2_asymptotic_v3"
   elif "RunIIFall17" in sname:
      year=2017
      globaltag="94X_mc2017_realistic_v17"
      metrecipe="metrecipe=True " # Space is important
   elif "RunIIAutumn18" in sname:
      year=2018
      globaltag="102X_upgrade2018_realistic_v21"
   else:
      raise RuntimeError("Year cannot be identified for {}".format(sname))

   if any(fsname in sname for fsname in [ "ZZTo2L2Nu", "ZZTo2L2Q", "ZZTo2Nu2X", "ZZTo4Q", "HTo2L2Q", "HTo2Nu2X", "HTo4Q" ]):
      vvmode="ZZ"
      if "2L2Nu" in sname:
         vvdecaymode=3
      else:
         vvdecaymode=-1
   else:
      vvmode="WW"
      if "2L2Nu" in sname:
         vvdecaymode=0
      else:
         vvdecaymode=-1

   if "GluGlu" in sname:
      mefile="GG_SIG_{}_0PM_H-HMASS-_POWHEG.me".format(vvmode)
   elif "VBF" in sname:
      mefile="VBF_SIG_{}_0PM_H-HMASS-_POWHEG.me".format(vvmode)
   elif "Wminus" in sname or "Wplus" in sname:
      mefile="WH_SIG_{}_0PM_H-HMASS-_POWHEG.me".format(vvmode)
   elif "ZH" in sname or "Wplus" in sname:
      mefile="ZH_SIG_{}_0PM_H-HMASS-_POWHEG.me".format(vvmode)
   else:
      raise RuntimeError("Cannot identify ME file for {}".format(sname))

   lineargs = {
      "SNAME" : sname,
      "GLOBALTAG" : globaltag,
      "MEFILE" : mefile,
      "YEAR" : year,
      "METRECIPE" : metrecipe,
      "VVMODE" : vvmode,
      "VVDECAYMODE" : vvdecaymode,
      "KFACTOROPT" : kfactoropt,
   }

   sline = line_tpl.format(**lineargs)
   if '{' in sname or '}' in sname:
      raise RuntimeError("There are leftover options.")

   return [ year, sline ]


def run(args):
   localSearchDir=args.localSearchDir
   sample_line_list = []
   valid_years = []
   for search_dir in localSearchDir:
      if not search_dir.endswith('/'):
         search_dir = search_dir + '/'
      subdirs=glob.glob(search_dir+"*/*/MINIAODSIM")
      for subdir in subdirs:
         sname = subdir.replace(search_dir,'/')
         nfiles = len(glob.glob(subdir+"/*.root"))
         nreq = getNFilesRequired(sname)
         if nfiles>=nreq:
            print("{}: {} / {}".format(sname, nfiles, nreq))
            csv_line = getCSVLine(sname)
            sample_line_list.append(csv_line)
            if csv_line[0] not in valid_years:
               valid_years.append(csv_line[0])
   sample_line_list.sort()
   valid_years.sort()
   for year in valid_years:
      csvname = args.csv.replace(".csv","_{}.csv".format(year))
      with open(csvname, 'w') as fout:
         fout.write('#Dataset,condoroutdir,xsec,BR,VVMode,VVDecayMode,lheMEfragment,options\n')
         for ss in sample_line_list:
            if ss[0]==year:
               fout.write(ss[1]+'\n')
      os.system("python replaceSignalXsec.py --csv={}".format(csvname))


if __name__ == "__main__":
   parser = argparse.ArgumentParser()
   parser.add_argument("--csv", help="Output csv file", type=str, required=True)
   parser.add_argument("localSearchDir", help="Search directory for local files. Can specify multiple.", nargs="+")
   args = parser.parse_args()

   run(args)