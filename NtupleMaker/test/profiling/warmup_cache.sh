#!/usr/bin/env bash

# Problem: Magnetic field data is stored in ~11k files (with avg size of 4KB
# each) on cvmfs, so sometimes cmsRun will hang for 5 minutes-10 minutes on the
# first event, when setting up Bfield stuff It gets better after this happens
# once and things are in the cache, but it's excruciatingly annoying when
# that's not the case. Hopefully I can remember to warmup the cache in the
# background after checking out NtupleMaker so I don't have to wait when
# testing.

# Based on https://github.com/cms-sw/cmssw/blob/CMSSW_10_2_X/MagneticField/Engine/python/volumeBasedMagneticField_160812_cfi.py
# or getting PID of cmsRun process and 
# strace -s 1024 -p <pid> 2>&1 | grep grid
# to see all the files/paths that need to be ls'd and opened. It's pretty dumb...

# Usage: ./warmup_cache.sh &
# If all the files are successfully put into the cache, then the Bfield stuff
# in the first event should only take 20-25 seconds.

# 94X 2016, 2017 and 2018
tag=grid_160812_3_8t
cat ${CMSSW_RELEASE_BASE}/external/${SCRAM_ARCH}/data/MagneticField/Interpolation/data/${tag}/*/*.bin >& /dev/null

# 80X 2016
tag=grid_120812_3_8t_v7_large
cat ${CMSSW_RELEASE_BASE}/external/${SCRAM_ARCH}/data/MagneticField/Interpolation/data/${tag}/*/*.bin >& /dev/null
