#!/usr/bin/env python

import sys
import os
from tqdm import tqdm
import time
import datetime

def linfit(xs,ys):
    import math
    n = len(xs)
    sumx, sumy = sum(xs), sum(ys)
    sumxx, sumyy = sum([x*x for x in xs]), sum([y*y for y in ys])
    sumxy = sum([xs[i]*ys[i] for i in range(len(xs))])
    avgx, avgy = 1.0*sumx/n, 1.0*sumy/n
    ssxx, ssyy = sumxx-n*avgx*avgx, sumyy-n*avgy*avgy
    ssxy = sumxy-n*avgx*avgy
    m = 1.0*ssxy/ssxx
    b = avgy-m*avgx
    try:
        s = math.sqrt((ssyy-m*ssxy)/(n-2))
        errorm = s/math.sqrt(ssxx)
        errorb = s/math.sqrt(1.0/n+avgx*avgx/ssxx)
        if(n == 2): errorm, errorb = 0.0, 0.0
        return m,b, errorm, errorb
    except:
        print "ERROR:",m,b,ssxx,ssyy,ssxy,avgx,n
        return m,b,-1,-1

if __name__ == "__main__":

    if len(sys.argv) < 2: 
        print "give me a filename of cmsRun output with"
        print "    process.Timing = cms.Service(\"Timing\")"
        print "in the pset!"
        sys.exit()
    fname = sys.argv[-1]

    # fname = "timing.txt"

    # fname = "timing_1k.txt"

    d_times = {}
    d_events = set([])
    eventtimes = []
    processingpairs = []
    with open(fname, "r") as fhin:
        for line in tqdm(fhin):
            if line.startswith("TimeEvent> "): eventtimes.append(float(line.split()[-1]))
            if line.startswith("Begin processing the"):
                record = float("".join([b for b in line.split("record")[0].split("the")[-1] if b in "1234567890"]))
                dtobj = datetime.datetime.strptime( line.split()[-2], "%H:%M:%S.%f" ).replace(year=2016)
                ts = time.mktime(dtobj.timetuple())+(dtobj.microsecond/1.e6)
                processingpairs.append([record,ts])
            if not line.startswith("TimeModule> "): continue
            parts = line.split()
            module = parts[-3]
            t = parts[-1]
            event = parts[1]
            if module not in d_times: d_times[module] = []
            d_events.add(event)
            d_times[module].append(float(t))

    # print d_times.keys()

    # drop first and last two points
    processingpairs = processingpairs[2:-2]
    # print processingpairs
    mint = min([pp[-1] for pp in processingpairs])
    processingpairs = map(lambda x: [x[0],x[1]-mint], processingpairs)
    m, b, merr, berr = linfit(*zip(*processingpairs))

    # nevents = len(d_events)
    nevents = len(eventtimes)

    # drop first and last 10 events
    for mname in d_times.keys():
        d_times[mname] = d_times[mname][10:-10]
    eventtimes = eventtimes[10:-10]
    nevents -= 20

    # print d_times["hltMaker"]


    tot_avg_time = sum(eventtimes)/nevents
    tot_avg_time_weird = sum([sum(d_times[mod]) for mod in d_times.keys()])/nevents

    avgs = []
    for maker in d_times.keys():
        avg = sum(d_times[maker])/tot_avg_time_weird/len(d_times[maker])
        avgs.append( [avg,maker] )

    print "Total events: {}".format(int(processingpairs[-1][0]))
    print "Average event rate: {0:.1f}Hz (linear fit to Begin Processing lines)".format(1.0/m)
    if eventtimes:
        print "Average event rate: {0:.1f}Hz (reported by TimeEvent)".format(tot_avg_time**-1.0)
        print "Average event rate: {0:.1f}Hz (manual sum of module times)".format(tot_avg_time_weird**-1.0)
        avgs = sorted(avgs,reverse=True)
        print "{0:50s} {1:10s}".format("Module", "frac real time")
        print "-"*70
        for tavg,module in avgs:
            print "{0:50s} {1:10.2f}".format(module, tavg)

