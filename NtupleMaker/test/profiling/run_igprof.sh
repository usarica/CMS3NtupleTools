#!/usr/bin/env bash

tag=myTag_Oct30
pset=../main_pset.py
psetflags="data=False"

igprof -d -pp -z -o igprof_${tag}.pp.gz cmsRun ${pset} ${psetflags} >& igtest_${tag}.pp.log
# make sql file to put in web area
igprof-analyse --sqlite -d -v -g igprof_${tag}.pp.gz | sqlite3 igreport_${tag}_perf.sql3 >& /dev/null
# copy to uaf
mkdir -p ~/public_html/cgi-bin/data/
cp igreport_${tag}_perf.sql3 ~/public_html/cgi-bin/data/
echo "Navigate to http://uaf-8.t2.ucsd.edu/~${USER}/cgi-bin/igprof-navigator.py/igreport_${tag}_perf/"

