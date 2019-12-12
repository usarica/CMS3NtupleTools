#!/usr/bin/env bash

cmd=cmsRun
pset=main_pset.py
nevents=100
outputdir=outputs/

mkdir -p $outputdir

# Handy PPD/NanoAOD info/links
# https://indico.cern.ch/event/778476/contributions/3239228/attachments/1770992/2877831/18-12-13_News_PPD.pdf
# https://indico.cern.ch/event/778476/contributions/3239274/attachments/1771132/2878044/xpog_ppd_13dec18.pdf

# NOTE 0. Use the python version of this instead...

# NOTE 1. If you want to use the "central" files, remove "/user/namin/localcache".
# I just copied them to my hadoop because xrootd somehow became the suckiest
# mechanism in the world, and I don't have 30 minutes to wait for a file to open

# NOTE 2. It's better to be explicit about all the arguments like the example below.
# I just put in some simple checks based on the input names for convenience.

# cmsRun main_pset.py \
#     data=False \
#     year=2016 \
#     is80x=True \
#     globaltag=80X_mcRun2_asymptotic_2016_TrancheIV_v8 \
#     inputs=/store/user/namin/localcache/mc/RunIISummer16MiniAODv2/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/1C8BFCBE-26B6-E611-80E5-A0000420FE80.root \
#     nevents=-1

# 2016 MC 80X MiniAODv2 -- /TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM
$cmd $pset dumpProcess=True globaltag=94X_mcRun2_asymptotic_v3 nevents=$nevents \
    inputs=/store/user/namin/localcache/mc/RunIISummer16MiniAODv2/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/1C8BFCBE-26B6-E611-80E5-A0000420FE80.root \
    output=$outputdir/ntuple_2016_mc_80xminiaodv2.root >& $outputdir/log_2016_mc_80xminiaodv2.txt &

# 2016 MC Fastsim 80X MiniAODv2 -- /TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM
#$cmd $pset dumpProcess=True globaltag=94X_mcRun2_asymptotic_v3 nevents=$nevents fastsim=True \
#    inputs=/store/user/namin/localcache/mc/RunIISpring16MiniAODv2/SMS-T1tttt_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16Fast_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/5A5E2A72-4F2D-E611-B2F5-02163E017638.root \
#    output=$outputdir/ntuple_2016_mc_80xfastsim.root >& $outputdir/log_2016_mc_80xfastsim.txt &

# 2016 Re-reco Data 94X MiniAODv3 -- /DoubleMuon/Run2016C-17Jul2018-v1/MINIAOD
$cmd $pset dumpProcess=True globaltag=94X_dataRun2_v10 nevents=$nevents \
    inputs=/store/user/namin/localcache/data/Run2016C/DoubleMuon/MINIAOD/17Jul2018-v1/50000/D229CC30-1E8B-E811-844A-A0369FD0B228.root \
    output=$outputdir/ntuple_2016_data_94xminiaodv3.root >& $outputdir/log_2016_data_94xminiaodv3.txt &

# 2016 MC 94X MiniAODv3 -- /DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2/MINIAODSIM
$cmd $pset dumpProcess=True globaltag=94X_mcRun2_asymptotic_v3 nevents=$nevents \
    inputs=/store/user/namin/localcache/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2/100000/6053E770-34E9-E811-98D2-246E96D10C28.root \
    output=$outputdir/ntuple_2016_mc_94xminiaodv3.root >& $outputdir/log_2016_mc_94xminiaodv3.txt &

# 2017 MC (with MET recipe) -- /ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM
$cmd $pset dumpProcess=True globaltag=94X_mc2017_realistic_v17 nevents=$nevents metrecipe=True \
    inputs=/store/user/namin/localcache/mc/RunIIFall17MiniAODv2/ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/50000/00436CBC-6B70-E811-850C-00259075D70C.root \
    output=$outputdir/ntuple_2017_mc.root >& $outputdir/log_2017_mc.txt &

# 2017 Re-reco Data (with MET recipe v2) -- /DoubleMuon/Run2017C-31Mar2018-v1/MINIAOD
$cmd $pset dumpProcess=True globaltag=94X_dataRun2_v11 nevents=$nevents metrecipe=True \
    inputs=/store/user/namin/localcache/data/Run2017C/DoubleMuon/MINIAOD/31Mar2018-v1/80000/04FCFB0D-FF39-E811-94C7-AC162DA6D2F8.root \
    output=$outputdir/ntuple_2017_data.root >& $outputdir/log_2017_data.txt &

# 2017 Re-reco Data (Re-Re-Reco of 2017F, so no MET recipe needed) -- /DoubleMuon/Run2017C-31Mar2018-v1/MINIAOD
$cmd $pset dumpProcess=True globaltag=94X_dataRun2_v11 nevents=$nevents \
    inputs=/store/user/namin/localcache/data/Run2017F/DoubleEG/MINIAOD/09May2018-v1/10000/444E03EB-B75F-E811-AFBA-F01FAFD8F16A.root \
    output=$outputdir/ntuple_2017_dataf.root >& $outputdir/log_2017_dataf.txt &

# 2018 Re-reco Data -- /DoubleMuon/Run2018A-17Sep2018-v2/MINIAOD
$cmd $pset dumpProcess=True globaltag=102X_dataRun2_v11 nevents=$nevents \
    inputs=/store/user/namin/localcache/data/Run2018A/DoubleMuon/MINIAOD/17Sep2018-v2/00000/7B954B49-BE06-B64C-89DC-F568513B41A3.root \
    output=$outputdir/ntuple_2018_data_rereco.root >& $outputdir/log_2018_data_rereco.txt &

# 2018 Prompt Data -- /EGamma/Run2018D-PromptReco-v2/MINIAOD
$cmd $pset dumpProcess=True globaltag=102X_dataRun2_Prompt_v14 nevents=$nevents \
    inputs=/store/user/namin/localcache/data/Run2018D/EGamma/MINIAOD/PromptReco-v2/000/322/204/00000/F09A218A-71B3-E811-9A04-02163E013E33.root \
    output=$outputdir/ntuple_2018_data_prompt.root >& $outputdir/log_2018_data_prompt.txt &

# 2018 MC -- /DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM
$cmd $pset dumpProcess=True globaltag=102X_upgrade2018_realistic_v19 nevents=$nevents \
    inputs=/store/user/namin/localcache/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/100000/042C8EE9-9431-5443-88C8-77F1D910B3A5.root \
    output=$outputdir/ntuple_2018_mc.root >& $outputdir/log_2018_mc.txt &

echo "After they are done, check the appropriate logs and root files manually if you"
echo "want echo to be sure that your change worked. For simple checks, just do the"
echo "following when they finish to check just the event counts:"
echo '   for output in $(ls '${outputdir}'/*.root); do edmFileUtil $output; done'
echo "Monitor with:"
echo '   tail -f '${outputdir}'/*.txt'
echo "Note that 80X probably failed. Still requires some debugging of jet collection updates"
