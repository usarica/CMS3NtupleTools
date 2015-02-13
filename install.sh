#!/bin/bash

#USER INPUTS
GlobalTag=PHYS14_25_V2::All
CMS3Tag=CMS3_V07-02-05

curl https://raw.githubusercontent.com/cmstas/NtupleTools/master/AutoTuple/setup.sh > setup2.sh

sed -i '2,6d' setup2.sh
sed -i "3s/.*/gtag=$GlobalTag/" setup2.sh
sed -i "4s/.*/tag=$CMS3Tag/" setup2.sh

head -n -20 setup2.sh > setup.sh
rm setup2.sh
sed -i 's,git clone git@github.com:cmstas/NtupleMaker.git CMS3/NtupleMaker,mv ../../* CMSSW_7_2_0/src/,' setup.sh
