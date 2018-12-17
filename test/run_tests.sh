#!/usr/bin/env bash

cmd=cmsRun
# cmd=python
pset=main_pset.py
nevents=500
outputdir=outputs/

mkdir -p $outputdir

# 2016 Re-reco Data MiniAODv3
# /DoubleMuon/Run2016C-17Jul2018-v1/MINIAOD
$cmd $pset \
    data=True \
    prompt=False \
    name=DQM \
    globaltag=94X_dataRun2_v10 \
    inputs=/store/data/Run2016C/DoubleMuon/MINIAOD/17Jul2018-v1/50000/D229CC30-1E8B-E811-844A-A0369FD0B228.root \
    nevents=$nevents \
    output=$outputdir/ntuple_2016_data_miniaodv3.root >& $outputdir/log_2016_data_miniaodv3.txt &


# 2016 MC MiniAODv3
# /DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2/MINIAODSIM
$cmd $pset \
    data=False \
    setup=2016 \
    globaltag=94X_mcRun2_asymptotic_v3 \
    inputs=/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2/100000/6053E770-34E9-E811-98D2-246E96D10C28.root \
    nevents=$nevents \
    output=$outputdir/ntuple_2016_mc_miniaodv3.root >& $outputdir/log_2016_mc_miniaodv3.txt &

# 2017 MC (with MET recipe)
# /ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM
$cmd $pset \
    data=False \
    setup=2017 \
    metrecipe=True \
    globaltag=94X_mc2017_realistic_v14 \
    inputs=/store/mc/RunIIFall17MiniAODv2/ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/50000/00436CBC-6B70-E811-850C-00259075D70C.root \
    nevents=$nevents \
    output=$outputdir/ntuple_2017_mc.root >& $outputdir/log_2017_mc.txt &

# 2017 Re-reco Data (with MET recipe)
# /DoubleMuon/Run2017C-31Mar2018-v1/MINIAOD
$cmd $pset \
    data=True \
    prompt=False \
    setup=2017 \
    metrecipe=True \
    globaltag=94X_dataRun2_ReReco_EOY17_v6 \
    inputs=/store/data/Run2017C/DoubleMuon/MINIAOD/31Mar2018-v1/80000/04FCFB0D-FF39-E811-94C7-AC162DA6D2F8.root \
    nevents=$nevents \
    output=$outputdir/ntuple_2017_data.root >& $outputdir/log_2017_data.txt &

# 2018 Re-reco Data
# /DoubleMuon/Run2018A-17Sep2018-v2/MINIAOD
$cmd $pset \
    data=True \
    prompt=True \
    globaltag=102X_dataRun2_Sep2018Rereco_v1 \
    inputs=/store/data/Run2018A/DoubleMuon/MINIAOD/17Sep2018-v2/00000/7B954B49-BE06-B64C-89DC-F568513B41A3.root \
    nevents=$nevents \
    output=$outputdir/ntuple_2018_data_rereco.root >& $outputdir/log_2018_data_rereco.txt &

# 2018 Prompt Data
# /EGamma/Run2018D-PromptReco-v2/MINIAOD
$cmd $pset \
    data=True \
    prompt=True \
    globaltag=102X_dataRun2_Prompt_v1 \
    inputs=/store/data/Run2018D/EGamma/MINIAOD/PromptReco-v2/000/320/500/00000/F6727085-F895-E811-A2FE-FA163ED11D1B.root \
    nevents=$nevents \
    output=$outputdir/ntuple_2018_data_prompt.root >& $outputdir/log_2018_data_prompt.txt &

# 2018 MC
# /DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM
$cmd $pset \
    data=False \
    setup=2018 \
    globaltag=102X_upgrade2018_realistic_v15 \
    inputs=/store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/100000/042C8EE9-9431-5443-88C8-77F1D910B3A5.root \
    nevents=$nevents \
    output=$outputdir/ntuple_2018_mc.root >& $outputdir/log_2018_mc.txt &


    # goldenjson=Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt

echo "After they are done, check the appropriate logs and root files manually if you"
echo "want echo to be sure that your change worked. For simple checks, just do the"
echo "following when they finish to check just the event counts:"
echo '   for output in $(ls '${outputdir}'/*.root); do edmFileUtil $output; done'
