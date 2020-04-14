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

  if (sname.Contains("MINIAOD")){
    TString stasample=sname;
    if (!sname.BeginsWith("/")) stasample = Form("/%s", stasample.Data());
    MELAerr << "SampleHelpers::constructSamplesList: Warning! Sample " << stasample << " is already a standalone sample. Use at your own risk!" << endl;
    samples.push_back(stasample);
    return;
  }

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
    if (sname == "MET"){
      constructSamplesList("MET_2018A", syst, samples);
      constructSamplesList("MET_2018B", syst, samples);
      constructSamplesList("MET_2018C", syst, samples);
      constructSamplesList("MET_2018D", syst, samples);
    }
    if (sname == "MET_2018A") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MET/Run2018A-17Sep2018-v1/MINIAOD" });
    if (sname == "MET_2018B") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MET/Run2018B-17Sep2018-v1/MINIAOD" });
    if (sname == "MET_2018C") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MET/Run2018C-17Sep2018-v1/MINIAOD" });
    if (sname == "MET_2018D") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MET/Run2018D-PromptReco-v2/MINIAOD" });
    if (sname == "JetHT"){
      constructSamplesList("JetHT_2018A", syst, samples);
      constructSamplesList("JetHT_2018B", syst, samples);
      constructSamplesList("JetHT_2018C", syst, samples);
      constructSamplesList("JetHT_2018D", syst, samples);
    }
    if (sname == "JetHT_2018A") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/JetHT/Run2018A-17Sep2018-v1/MINIAOD" });
    if (sname == "JetHT_2018B") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/JetHT/Run2018B-17Sep2018-v1/MINIAOD" });
    if (sname == "JetHT_2018C") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/JetHT/Run2018C-17Sep2018-v1/MINIAOD" });
    if (sname == "JetHT_2018D") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/JetHT/Run2018D-PromptReco-v2/MINIAOD" });
    // Group by runs
    if (sname == "Run2018A"){
      constructSamplesList("EGamma_2018A", syst, samples);
      constructSamplesList("MuonEG_2018A", syst, samples);
      constructSamplesList("DoubleMuon_2018A", syst, samples);
      constructSamplesList("SingleMuon_2018A", syst, samples);
      constructSamplesList("MET_2018A", syst, samples);
      constructSamplesList("JetHT_2018A", syst, samples);
    }
    if (sname == "Run2018B"){
      constructSamplesList("EGamma_2018B", syst, samples);
      constructSamplesList("MuonEG_2018B", syst, samples);
      constructSamplesList("DoubleMuon_2018B", syst, samples);
      constructSamplesList("SingleMuon_2018B", syst, samples);
      constructSamplesList("MET_2018B", syst, samples);
      constructSamplesList("JetHT_2018B", syst, samples);
    }
    if (sname == "Run2018C"){
      constructSamplesList("EGamma_2018C", syst, samples);
      constructSamplesList("MuonEG_2018C", syst, samples);
      constructSamplesList("DoubleMuon_2018C", syst, samples);
      constructSamplesList("SingleMuon_2018C", syst, samples);
      constructSamplesList("MET_2018C", syst, samples);
      constructSamplesList("JetHT_2018C", syst, samples);
    }
    if (sname == "Run2018D"){
      constructSamplesList("EGamma_2018D", syst, samples);
      constructSamplesList("MuonEG_2018D", syst, samples);
      constructSamplesList("DoubleMuon_2018D", syst, samples);
      constructSamplesList("SingleMuon_2018D", syst, samples);
      constructSamplesList("MET_2018D", syst, samples);
      constructSamplesList("JetHT_2018D", syst, samples);
    }

    // Simulation for the main signals
    // GGH ZZ 2l2nu POWHEG
    if (sname == "GGH_ZZ2L2Nu_M200_POWHEG" || sname == "GGH_ZZ2L2Nu_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToZZTo2L2Nu_M200_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "GGH_ZZ2L2Nu_M300_POWHEG" || sname == "GGH_ZZ2L2Nu_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToZZTo2L2Nu_M300_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "GGH_ZZ2L2Nu_M400_POWHEG" || sname == "GGH_ZZ2L2Nu_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToZZTo2L2Nu_M400_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "GGH_ZZ2L2Nu_M500_POWHEG" || sname == "GGH_ZZ2L2Nu_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToZZTo2L2Nu_M500_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "GGH_ZZ2L2Nu_M600_POWHEG" || sname == "GGH_ZZ2L2Nu_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToZZTo2L2Nu_M600_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "GGH_ZZ2L2Nu_M700_POWHEG" || sname == "GGH_ZZ2L2Nu_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToZZTo2L2Nu_M700_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "GGH_ZZ2L2Nu_M800_POWHEG" || sname == "GGH_ZZ2L2Nu_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToZZTo2L2Nu_M800_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "GGH_ZZ2L2Nu_M900_POWHEG" || sname == "GGH_ZZ2L2Nu_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToZZTo2L2Nu_M900_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "GGH_ZZ2L2Nu_M1000_POWHEG" || sname == "GGH_ZZ2L2Nu_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToZZTo2L2Nu_M1000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "GGH_ZZ2L2Nu_M1500_POWHEG" || sname == "GGH_ZZ2L2Nu_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToZZTo2L2Nu_M1500_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "GGH_ZZ2L2Nu_M2000_POWHEG" || sname == "GGH_ZZ2L2Nu_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToZZTo2L2Nu_M2000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "GGH_ZZ2L2Nu_M2500_POWHEG" || sname == "GGH_ZZ2L2Nu_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToZZTo2L2Nu_M2500_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "GGH_ZZ2L2Nu_M3000_POWHEG" || sname == "GGH_ZZ2L2Nu_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToZZTo2L2Nu_M3000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    // VBF ZZ 2l2nu POWHEG
    if (sname == "VBF_ZZ2L2Nu_M200_POWHEG" || sname == "VBF_ZZ2L2Nu_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/VBF_HToZZTo2L2Nu_M200_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "VBF_ZZ2L2Nu_M300_POWHEG" || sname == "VBF_ZZ2L2Nu_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/VBF_HToZZTo2L2Nu_M300_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "VBF_ZZ2L2Nu_M400_POWHEG" || sname == "VBF_ZZ2L2Nu_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/VBF_HToZZTo2L2Nu_M400_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "VBF_ZZ2L2Nu_M500_POWHEG" || sname == "VBF_ZZ2L2Nu_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/VBF_HToZZTo2L2Nu_M500_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "VBF_ZZ2L2Nu_M600_POWHEG" || sname == "VBF_ZZ2L2Nu_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/VBF_HToZZTo2L2Nu_M600_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "VBF_ZZ2L2Nu_M700_POWHEG" || sname == "VBF_ZZ2L2Nu_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/VBF_HToZZTo2L2Nu_M700_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "VBF_ZZ2L2Nu_M800_POWHEG" || sname == "VBF_ZZ2L2Nu_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/VBF_HToZZTo2L2Nu_M800_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "VBF_ZZ2L2Nu_M900_POWHEG" || sname == "VBF_ZZ2L2Nu_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/VBF_HToZZTo2L2Nu_M900_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "VBF_ZZ2L2Nu_M1000_POWHEG" || sname == "VBF_ZZ2L2Nu_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/VBF_HToZZTo2L2Nu_M1000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "VBF_ZZ2L2Nu_M1500_POWHEG" || sname == "VBF_ZZ2L2Nu_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/VBF_HToZZTo2L2Nu_M1500_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "VBF_ZZ2L2Nu_M2000_POWHEG" || sname == "VBF_ZZ2L2Nu_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/VBF_HToZZTo2L2Nu_M2000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "VBF_ZZ2L2Nu_M2500_POWHEG" || sname == "VBF_ZZ2L2Nu_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/VBF_HToZZTo2L2Nu_M2500_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "VBF_ZZ2L2Nu_M3000_POWHEG" || sname == "VBF_ZZ2L2Nu_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/VBF_HToZZTo2L2Nu_M3000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });

    // GGH WW 2l2nu POWHEG
    if (sname == "GGH_WW2L2Nu_M125_POWHEG" || sname == "GGH_WW2L2Nu_POWHEG")  HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "GGH_WW2L2Nu_M160_POWHEG" || sname == "GGH_WW2L2Nu_POWHEG")  HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToWWTo2L2Nu_M160_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "GGH_WW2L2Nu_M165_POWHEG" || sname == "GGH_WW2L2Nu_POWHEG")  HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToWWTo2L2Nu_M165_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "GGH_WW2L2Nu_M170_POWHEG" || sname == "GGH_WW2L2Nu_POWHEG")  HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToWWTo2L2Nu_M170_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "GGH_WW2L2Nu_M175_POWHEG" || sname == "GGH_WW2L2Nu_POWHEG")  HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToWWTo2L2Nu_M175_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "GGH_WW2L2Nu_M180_POWHEG" || sname == "GGH_WW2L2Nu_POWHEG")  HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToWWTo2L2Nu_M180_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "GGH_WW2L2Nu_M190_POWHEG" || sname == "GGH_WW2L2Nu_POWHEG")  HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToWWTo2L2Nu_M190_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "GGH_WW2L2Nu_M200_POWHEG" || sname == "GGH_WW2L2Nu_POWHEG")  HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToWWTo2L2Nu_M200_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "GGH_WW2L2Nu_M210_POWHEG" || sname == "GGH_WW2L2Nu_POWHEG")  HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToWWTo2L2Nu_M210_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "GGH_WW2L2Nu_M230_POWHEG" || sname == "GGH_WW2L2Nu_POWHEG")  HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToWWTo2L2Nu_M230_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "GGH_WW2L2Nu_M250_POWHEG" || sname == "GGH_WW2L2Nu_POWHEG")  HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToWWTo2L2Nu_M250_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "GGH_WW2L2Nu_M270_POWHEG" || sname == "GGH_WW2L2Nu_POWHEG")  HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToWWTo2L2Nu_M270_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "GGH_WW2L2Nu_M300_POWHEG" || sname == "GGH_WW2L2Nu_POWHEG")  HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToWWTo2L2Nu_M300_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "GGH_WW2L2Nu_M450_POWHEG" || sname == "GGH_WW2L2Nu_POWHEG")  HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToWWTo2L2Nu_M450_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "GGH_WW2L2Nu_M500_POWHEG" || sname == "GGH_WW2L2Nu_POWHEG")  HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToWWTo2L2Nu_M500_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "GGH_WW2L2Nu_M550_POWHEG" || sname == "GGH_WW2L2Nu_POWHEG")  HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToWWTo2L2Nu_M550_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "GGH_WW2L2Nu_M600_POWHEG" || sname == "GGH_WW2L2Nu_POWHEG")  HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToWWTo2L2Nu_M600_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "GGH_WW2L2Nu_M650_POWHEG" || sname == "GGH_WW2L2Nu_POWHEG")  HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToWWTo2L2Nu_M650_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "GGH_WW2L2Nu_M700_POWHEG" || sname == "GGH_WW2L2Nu_POWHEG")  HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToWWTo2L2Nu_M700_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "GGH_WW2L2Nu_M800_POWHEG" || sname == "GGH_WW2L2Nu_POWHEG")  HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToWWTo2L2Nu_M800_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "GGH_WW2L2Nu_M900_POWHEG" || sname == "GGH_WW2L2Nu_POWHEG")  HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToWWTo2L2Nu_M900_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "GGH_WW2L2Nu_M1000_POWHEG" || sname == "GGH_WW2L2Nu_POWHEG")  HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToWWTo2L2Nu_M1000_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "GGH_WW2L2Nu_M1500_POWHEG" || sname == "GGH_WW2L2Nu_POWHEG")  HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToWWTo2L2Nu_M1500_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "GGH_WW2L2Nu_M2000_POWHEG" || sname == "GGH_WW2L2Nu_POWHEG")  HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToWWTo2L2Nu_M2000_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "GGH_WW2L2Nu_M2500_POWHEG" || sname == "GGH_WW2L2Nu_POWHEG")  HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToWWTo2L2Nu_M2500_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "GGH_WW2L2Nu_M3000_POWHEG" || sname == "GGH_WW2L2Nu_POWHEG")  HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GluGluHToWWTo2L2Nu_M3000_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });

    // VBF ZZ 2l2nu POWHEG
    if (sname == "VBF_WW2L2Nu_M125_POWHEG" || sname == "VBF_WW2L2Nu_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/VBFHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "VBF_WW2L2Nu_M1000_POWHEG" || sname == "VBF_WW2L2Nu_POWHEG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/VBFHToWWTo2L2Nu_M1000_13TeV_powheg2_JHUGenV710_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });

    // Simulation for the main backgrounds
    if (sname == "DY_2l_M_10to50" || sname == "DY_2l") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "DY_2l_M_10to50_ext" || sname == "DY_2l") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM" });
    if (sname == "DY_2l_M_50" || sname == "DY_2l") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "DY_2l_M_50_ext" || sname == "DY_2l") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v1/MINIAODSIM" });
    if (sname == "DY_2l_M_50_HT_70-100" || sname == "DY_2l_M_50_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DYJetsToLL_M-50_HT-70to100_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "DY_2l_M_50_HT_100-200" || sname == "DY_2l_M_50_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "DY_2l_M_50_HT_200-400" || sname == "DY_2l_M_50_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "DY_2l_M_50_HT_400-600" || sname == "DY_2l_M_50_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v3/MINIAODSIM" });
    if (sname == "DY_2l_M_50_HT_600-800" || sname == "DY_2l_M_50_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "DY_2l_M_50_HT_800-1200" || sname == "DY_2l_M_50_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "DY_2l_M_50_HT_1200-2500" || sname == "DY_2l_M_50_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "DY_2l_M_50_HT_2500-inf" || sname == "DY_2l_M_50_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    // Bkgs. with tops
    if (sname == "TT_2l2nu") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "ST_t-channel_top_5f") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ST_t-channel_top_5f_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "ST_t-channel_antitop_5f") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ST_t-channel_antitop_5f_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "ST_s-channel_top_leptonDecays") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ST_s-channel_top_leptonDecays_13TeV-PSweights_powheg-pythia/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "ST_s-channel_antitop_leptonDecays") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ST_s-channel_antitop_leptonDecays_13TeV-PSweights_powheg-pythia/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "ST_tW_top_5f_NoFullyHadronicDecays") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v3/MINIAODSIM" });
    if (sname == "ST_tW_antitop_5f_NoFullyHadronicDecays") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v3/MINIAODSIM" });
    if (sname == "TTZ_2l2nu") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM" });
    if (sname == "TTW_lnu") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM" });
    // qqVV
    if (sname == "qqZZ_2l2nu" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZZTo2L2Nu_TuneCP5_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v2/MINIAODSIM" });
    if (sname == "qqWW_lnu2q" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "qqWW_2l2nu" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "qqZZ_4l" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZZTo4L_TuneCP5_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v2/MINIAODSIM" });
    if (sname == "qqZZ_2q2nu" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZZTo2Q2Nu_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "qqZZ_2l2q" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "qqWZ_lnu2q" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "qqWZ_2l2q" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "qqWZ_3lnu_POWHEG" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WZTo3LNu_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM" });
    if (sname == "qqWZ_3lnu_MG" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM" }); // Has more stats
    if (sname == "qqWZ_13nu" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "qqWG_lnu" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    // VVV
    if (sname == "WWW_4f" || sname == "VVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM" });
    if (sname == "WWZ" || sname == "VVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WWZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM" });
    if (sname == "WZG" || sname == "VVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WZG_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "WZZ" || sname == "VVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM" });
    if (sname == "ZZZ" || sname == "VVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM" });

    // Simulation for the photon CR processes
    // G + jets
    if (sname == "GJets_HT_40-100" || sname == "GJets_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "GJets_HT_100-200" || sname == "GJets_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-4cores5k_102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "GJets_HT_200-400" || sname == "GJets_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "GJets_HT_400-600" || sname == "GJets_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "GJets_HT_600-inf" || sname == "GJets_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM" });
    // T/TT + G + jets
    if (sname == "TGJets") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM" });
    if (sname == "TTGJets") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    // W + jets
    if (sname == "WJets_lnu_HT_70-100" || sname == "WJets_lnu" || sname == "WJets_lnu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "WJets_lnu_HT_100-200" || sname == "WJets_lnu" || sname == "WJets_lnu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "WJets_lnu_HT_200-400" || sname == "WJets_lnu" || sname == "WJets_lnu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "WJets_lnu_HT_400-600" || sname == "WJets_lnu" || sname == "WJets_lnu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "WJets_lnu_HT_600-800" || sname == "WJets_lnu" || sname == "WJets_lnu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "WJets_lnu_HT_800-1200" || sname == "WJets_lnu" || sname == "WJets_lnu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "WJets_lnu_HT_1200-2500" || sname == "WJets_lnu" || sname == "WJets_lnu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "WJets_lnu_HT_2500-inf" || sname == "WJets_lnu" || sname == "WJets_lnu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    //if (sname == "ZJets_nunu") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{
    if (sname == "ZJets_nunu_HT_100-200" || sname == "ZJets_nunu" || sname == "ZJets_nunu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZJetsToNuNu_HT-100To200_13TeV-madgraph/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "ZJets_nunu_HT_200-400" || sname == "ZJets_nunu" || sname == "ZJets_nunu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZJetsToNuNu_HT-200To400_13TeV-madgraph/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "ZJets_nunu_HT_400-600" || sname == "ZJets_nunu" || sname == "ZJets_nunu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZJetsToNuNu_HT-400To600_13TeV-madgraph/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "ZJets_nunu_HT_600-800" || sname == "ZJets_nunu" || sname == "ZJets_nunu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZJetsToNuNu_HT-600To800_13TeV-madgraph/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "ZJets_nunu_HT_800-1200" || sname == "ZJets_nunu" || sname == "ZJets_nunu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZJetsToNuNu_HT-800To1200_13TeV-madgraph/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "ZJets_nunu_HT_1200-2500" || sname == "ZJets_nunu" || sname == "ZJets_nunu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZJetsToNuNu_HT-1200To2500_13TeV-madgraph/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "ZJets_nunu_HT_2500-inf" || sname == "ZJets_nunu" || sname == "ZJets_nunu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "ZGJets_nunu_pTG_40-130" || sname == "ZGJets_nunu" || sname == "ZGJets_nunu_pTG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZNuNuGJets_MonoPhoton_PtG-40to130_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "ZGJets_nunu_pTG_130-inf" || sname == "ZGJets_nunu" || sname == "ZGJets_nunu_pTG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZNuNuGJets_MonoPhoton_PtG-130_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "ZGJets_nunu_nlo_pTG_130-inf" || sname == "ZGJets_nunu_nlo" || sname == "ZGJets_nunu_nlo_pTG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZGTo2NuG_PtG-130_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    // if (sname == "ZGJets_nunu_nlo_inclusive") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZGTo2NuG_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "ZGJets_ll_pTG_40-130" || sname == "ZGJets_ll" || sname == "ZGJets_ll_pTG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZLLGJets_MonoPhoton_PtG-40to130_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "ZGJets_ll_pTG_130-inf" || sname == "ZGJets_ll" || sname == "ZGJets_ll_pTG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZLLGJets_MonoPhoton_PtG-130_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "ZGJets_ll_nlo_pTG_130-inf" || sname == "ZGJets_ll_nlo" || sname == "ZGJets_ll_nlo_pTG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "ZGTo2LG_PtG-130_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    //if (sname == "QCD_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{
    if (sname == "QCD_HT_100-200" || sname == "QCD_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "QCD_HT_200-300" || sname == "QCD_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "QCD_HT_300-500" || sname == "QCD_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "QCD_HT_500-700" || sname == "QCD_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "QCD_HT_700-1000" || sname == "QCD_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "QCD_HT_1000-1500" || sname == "QCD_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "QCD_HT_1500-2000" || sname == "QCD_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "QCD_HT_2000-inf" || sname == "QCD_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "TTJets") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v1/MINIAODSIM" });
    if (sname == "GJets_topology"){
      constructSamplesList("GJets_HT", syst, samples);
      constructSamplesList("QCD_HT", syst, samples);
      constructSamplesList("WJets_lnu", syst, samples);
      constructSamplesList("ZJets_nunu", syst, samples);
      constructSamplesList("ZGJets_nunu_pTG_40-130", syst, samples); // FIXME: LO sample at the moment
      constructSamplesList("ZGJets_nunu_nlo_pTG_130-inf", syst, samples);
      constructSamplesList("ZGJets_ll_pTG_40-130", syst, samples); // FIXME: LO sample at the moment
      constructSamplesList("ZGJets_ll_nlo_pTG_130-inf", syst, samples);
      constructSamplesList("TTJets", syst, samples);
      constructSamplesList("TGJets", syst, samples);
      constructSamplesList("TTGJets", syst, samples);
      constructSamplesList("qqWG_lnu", syst, samples);
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
