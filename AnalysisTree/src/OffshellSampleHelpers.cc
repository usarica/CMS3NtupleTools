// For the messier stuff
#include <cassert>
#include "OffshellSampleHelpers.h"
#include "HelperFunctions.h"
#include "HostHelpersCore.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;
using namespace HelperFunctions;


void SampleHelpers::constructSamplesList(TString const& sname, SystematicsHelpers::SystematicVariationTypes syst, std::vector<TString>& samples){
  assert(runConfigure);

  if (theDataYear == 2018){
    // Data
    if (sname == "SingleMuon"){
      constructSamplesList("SingleMuon_2018A", syst, samples);
      constructSamplesList("SingleMuon_2018B", syst, samples);
      constructSamplesList("SingleMuon_2018C", syst, samples);
      constructSamplesList("SingleMuon_2018D", syst, samples);
    }
    if (sname == "SingleMuon_2018A") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SingleMuon/Run2018A-17Sep2018-v2/MINIAOD" });
    if (sname == "SingleMuon_2018B") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SingleMuon/Run2018B-17Sep2018-v1/MINIAOD" });
    if (sname == "SingleMuon_2018C") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SingleMuon/Run2018C-17Sep2018-v1/MINIAOD" });
    if (sname == "SingleMuon_2018D") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SingleMuon/Run2018D-PromptReco-v2/MINIAOD" });
    if (sname == "DoubleMuon"){
      constructSamplesList("DoubleMuon_2018A", syst, samples);
      constructSamplesList("DoubleMuon_2018B", syst, samples);
      constructSamplesList("DoubleMuon_2018C", syst, samples);
      constructSamplesList("DoubleMuon_2018D", syst, samples);
    }
    if (sname == "DoubleMuon_2018A") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DoubleMuon/Run2018A-17Sep2018-v2/MINIAOD" });
    if (sname == "DoubleMuon_2018B") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DoubleMuon/Run2018B-17Sep2018-v1/MINIAOD" });
    if (sname == "DoubleMuon_2018C") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DoubleMuon/Run2018C-17Sep2018-v1/MINIAOD" });
    if (sname == "DoubleMuon_2018D") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DoubleMuon/Run2018D-PromptReco-v2/MINIAOD" });
    if (sname == "MuonEG"){
      constructSamplesList("MuonEG_2018A", syst, samples);
      constructSamplesList("MuonEG_2018B", syst, samples);
      constructSamplesList("MuonEG_2018C", syst, samples);
      constructSamplesList("MuonEG_2018D", syst, samples);
    }
    if (sname == "MuonEG_2018A") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MuonEG/Run2018A-17Sep2018-v1/MINIAOD" });
    if (sname == "MuonEG_2018B") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MuonEG/Run2018B-17Sep2018-v1/MINIAOD" });
    if (sname == "MuonEG_2018C") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MuonEG/Run2018C-17Sep2018-v1/MINIAOD" });
    if (sname == "MuonEG_2018D") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MuonEG/Run2018D-PromptReco-v2/MINIAOD" });
    if (sname == "EGamma"){
      constructSamplesList("EGamma_2018A", syst, samples);
      constructSamplesList("EGamma_2018B", syst, samples);
      constructSamplesList("EGamma_2018C", syst, samples);
      constructSamplesList("EGamma_2018D", syst, samples);
    }
    if (sname == "EGamma_2018A") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/EGamma/Run2018A-17Sep2018-v2/MINIAOD" });
    if (sname == "EGamma_2018B") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/EGamma/Run2018B-17Sep2018-v1/MINIAOD" });
    if (sname == "EGamma_2018C") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/EGamma/Run2018C-17Sep2018-v1/MINIAOD" });
    if (sname == "EGamma_2018D") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/EGamma/Run2018D-22Jan2019-v2/MINIAOD" });
    // Group by runs
    if (sname == "Run2018A"){
      constructSamplesList("DoubleMuon_2018A", syst, samples);
      constructSamplesList("MuonEG_2018A", syst, samples);
      constructSamplesList("EGamma_2018A", syst, samples);
      constructSamplesList("SingleMuon_2018A", syst, samples);
    }
    if (sname == "Run2018B"){
      constructSamplesList("DoubleMuon_2018B", syst, samples);
      constructSamplesList("MuonEG_2018B", syst, samples);
      constructSamplesList("EGamma_2018B", syst, samples);
      constructSamplesList("SingleMuon_2018B", syst, samples);
    }
    if (sname == "Run2018C"){
      constructSamplesList("DoubleMuon_2018C", syst, samples);
      constructSamplesList("MuonEG_2018C", syst, samples);
      constructSamplesList("EGamma_2018C", syst, samples);
      constructSamplesList("SingleMuon_2018C", syst, samples);
    }
    if (sname == "Run2018D"){
      constructSamplesList("DoubleMuon_2018D", syst, samples);
      constructSamplesList("MuonEG_2018D", syst, samples);
      constructSamplesList("EGamma_2018D", syst, samples);
      constructSamplesList("SingleMuon_2018D", syst, samples);
    }

    // Simulation for the main signals
    // GGH POWHEG
    if (sname == "GGH_M200_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToZZTo2L2Nu_M200_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "GGH_M300_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToZZTo2L2Nu_M300_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "GGH_M400_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToZZTo2L2Nu_M400_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "GGH_M500_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToZZTo2L2Nu_M500_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "GGH_M600_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToZZTo2L2Nu_M600_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "GGH_M700_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToZZTo2L2Nu_M700_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "GGH_M800_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToZZTo2L2Nu_M800_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "GGH_M900_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToZZTo2L2Nu_M900_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "GGH_M1000_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToZZTo2L2Nu_M1000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "GGH_M1500_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToZZTo2L2Nu_M1500_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "GGH_M2000_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToZZTo2L2Nu_M2000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "GGH_M2500_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToZZTo2L2Nu_M2500_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "GGH_M3000_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToZZTo2L2Nu_M3000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    // VBF POWHEG
    if (sname == "VBF_M200_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/VBF_HToZZTo2L2Nu_M200_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "VBF_M300_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/VBF_HToZZTo2L2Nu_M300_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "VBF_M400_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/VBF_HToZZTo2L2Nu_M400_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "VBF_M500_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/VBF_HToZZTo2L2Nu_M500_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "VBF_M600_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/VBF_HToZZTo2L2Nu_M600_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "VBF_M700_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/VBF_HToZZTo2L2Nu_M700_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "VBF_M800_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/VBF_HToZZTo2L2Nu_M800_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "VBF_M900_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/VBF_HToZZTo2L2Nu_M900_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "VBF_M1000_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/VBF_HToZZTo2L2Nu_M1000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "VBF_M1500_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/VBF_HToZZTo2L2Nu_M1500_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "VBF_M2000_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/VBF_HToZZTo2L2Nu_M2000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "VBF_M2500_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/VBF_HToZZTo2L2Nu_M2500_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "VBF_M3000_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/VBF_HToZZTo2L2Nu_M3000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });

    // Simulation for the main backgrounds
    if (sname == "DY_2l_M_10to50") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "DY_2l_M_50") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "DY_2l"){
      constructSamplesList("DY_2l_M10-50", syst, samples);
      constructSamplesList("DY_2l_M50", syst, samples);
    }
    if (sname == "TT_2l2nu") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "qqZZ_2l2nu") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZZTo2L2Nu_TuneCP5_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v2/MINIAODSIM" });
    if (sname == "qqWW_2l2nu") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "qqZZ_4l") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZZTo4L_TuneCP5_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v2/MINIAODSIM" });
    if (sname == "qqZZ_2q2nu") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZZTo2Q2Nu_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "qqZZ_2l2q") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "qqWZ_2l2q") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "qqVV"){
      constructSamplesList("qqZZ_2l2nu", syst, samples);
      constructSamplesList("qqWW_2l2nu", syst, samples);
      constructSamplesList("qqZZ_4l", syst, samples);
      constructSamplesList("qqZZ_2q2nu", syst, samples);
      constructSamplesList("qqZZ_2l2q", syst, samples);
      constructSamplesList("qqWZ_2l2q", syst, samples);
    }

    // Simulation for the photon CR processes
    if (sname == "GJets_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{
        "/GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "/GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-4cores5k_102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "/GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "/GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "/GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM"
      }
    );
    if (sname == "TGJets") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM" });
    if (sname == "TTGJets") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "WG_lnu") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "WZG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WZG_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "WJets_lnu") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{
        "/WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
      }
    );
    if (sname == "ZJets_nunu") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{
        "/ZJetsToNuNu_HT-100To200_13TeV-madgraph/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "/ZJetsToNuNu_HT-200To400_13TeV-madgraph/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "/ZJetsToNuNu_HT-400To600_13TeV-madgraph/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
        "/ZJetsToNuNu_HT-600To800_13TeV-madgraph/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "/ZJetsToNuNu_HT-800To1200_13TeV-madgraph/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "/ZJetsToNuNu_HT-1200To2500_13TeV-madgraph/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "/ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
      }
    );
    if (sname == "QCD_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{
        "/QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "/QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "/QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "/QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "/QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "/QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "/QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "/QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
      }
    );
    if (sname == "TTJets") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v1/MINIAODSIM" });
    if (sname == "GJets_topology"){
      constructSamplesList("WJets_lnu", syst, samples);
      constructSamplesList("ZJets_nunu", syst, samples);
      constructSamplesList("GJets_HT", syst, samples);
      constructSamplesList("QCD_HT", syst, samples);
      constructSamplesList("TTJets", syst, samples);
      constructSamplesList("TGJets", syst, samples);
      constructSamplesList("TTGJets", syst, samples);
      constructSamplesList("WG_lnu", syst, samples);
      constructSamplesList("WZG", syst, samples);
    }
  }
}
void SampleHelpers::getSamplesList(std::vector<TString> const& s, std::vector<TString>& vs, SystematicsHelpers::SystematicVariationTypes syst, std::vector<size_t>* ns){
  for (auto const& ss:s){
    std::vector<TString> dumappend; constructSamplesList(ss, syst, dumappend);
    HelperFunctions::appendVector<TString>(vs, dumappend);
    if (ns) ns->push_back(dumappend.size());
  }
}
