#!/usr/bin/env bash

# cmd=python
cmd=cmsRun

$cmd main_pset.py \
data=False \
setup=2017 \
metrecipe=True \
globaltag=94X_mc2017_realistic_v14 \
inputs=/store/mc/RunIIFall17MiniAODv2/ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/50000/00436CBC-6B70-E811-850C-00259075D70C.root

# $cmd main_pset.py \
# data=True prompt=True \
# globaltag=102X_dataRun2_Sep2018Rereco_v1 \
# inputs=/store/data/Run2018A/DoubleMuon/MINIAOD/17Sep2018-v2/00000/7B954B49-BE06-B64C-89DC-F568513B41A3.root
