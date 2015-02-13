#!/bin/bash

#USER INPUTS
GlobalTag=PHYS14_25_V2::All
CMS3Tag=CMS3_V07-02-05

curl https://raw.githubusercontent.com/cmstas/NtupleTools/master/AutoTuple/setup.sh > setup.sh

sed -i '2,6d' setup.sh
sed -i "3s/.*/gtag=$GlobalTag/" setup.sh
sed -i "4s/.*/tag=$CMS3Tag/" setup.sh
sed -i '$d' setup.sh
sed -i '$d' setup.sh
