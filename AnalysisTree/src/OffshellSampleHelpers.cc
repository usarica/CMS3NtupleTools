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

  using namespace SystematicsHelpers;

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
    // Z + jets, Z->ll
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
    if (sname == "TZ_2l_4f") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/tZq_ll_4f_ckm_NLO_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM" });
    if (sname == "TTZ_2l2nu") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM" });
    if (sname == "TTW_lnu") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM" });
    // qqVV
    if (sname == "qqZZ_2l2nu_mZ_18-inf" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZZTo2L2Nu_mZMin-18_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "qqZZ_2l2nu" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZZTo2L2Nu_TuneCP5_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v2/MINIAODSIM" });
    if (sname == "qqZZ_2l2nu_ext" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZZTo2L2Nu_TuneCP5_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM" });
    if (sname == "qqWW_lnu2q" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "qqWW_2l2nu" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "qqZZ_4l" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZZTo4L_TuneCP5_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v2/MINIAODSIM" });
    if (sname == "qqZZ_2q2nu" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZZTo2Q2Nu_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "qqZZ_2l2q" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "qqWZ_lnu2q" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "qqWZ_2l2q" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "qqWZ_3lnu_POWHEG_mll_0p1-inf" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WZTo3LNu_mllmin01_NNPDF31_TuneCP5_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
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
    if (sname == "WJets_lnu" || sname == "WJets_lnu_Njets" || sname == "WJets_lnu_0j") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "WJets_lnu" || sname == "WJets_lnu_Njets" || sname == "WJets_lnu_1j") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "WJets_lnu" || sname == "WJets_lnu_Njets" || sname == "WJets_lnu_2j") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "WJets_lnu_HT_70-100" || sname == "WJets_lnu" || sname == "WJets_lnu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "WJets_lnu_HT_100-200" || sname == "WJets_lnu" || sname == "WJets_lnu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "WJets_lnu_HT_200-400" || sname == "WJets_lnu" || sname == "WJets_lnu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "WJets_lnu_HT_400-600" || sname == "WJets_lnu" || sname == "WJets_lnu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "WJets_lnu_HT_600-800" || sname == "WJets_lnu" || sname == "WJets_lnu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "WJets_lnu_HT_800-1200" || sname == "WJets_lnu" || sname == "WJets_lnu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "WJets_lnu_HT_1200-2500" || sname == "WJets_lnu" || sname == "WJets_lnu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "WJets_lnu_HT_2500-inf" || sname == "WJets_lnu" || sname == "WJets_lnu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    // Z + jets, Z->nunu
    if (sname == "ZJets_nunu_HT_100-200" || sname == "ZJets_nunu" || sname == "ZJets_nunu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZJetsToNuNu_HT-100To200_13TeV-madgraph/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "ZJets_nunu_HT_200-400" || sname == "ZJets_nunu" || sname == "ZJets_nunu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZJetsToNuNu_HT-200To400_13TeV-madgraph/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "ZJets_nunu_HT_400-600" || sname == "ZJets_nunu" || sname == "ZJets_nunu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZJetsToNuNu_HT-400To600_13TeV-madgraph/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "ZJets_nunu_HT_600-800" || sname == "ZJets_nunu" || sname == "ZJets_nunu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZJetsToNuNu_HT-600To800_13TeV-madgraph/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "ZJets_nunu_HT_800-1200" || sname == "ZJets_nunu" || sname == "ZJets_nunu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZJetsToNuNu_HT-800To1200_13TeV-madgraph/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "ZJets_nunu_HT_1200-2500" || sname == "ZJets_nunu" || sname == "ZJets_nunu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZJetsToNuNu_HT-1200To2500_13TeV-madgraph/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "ZJets_nunu_HT_2500-inf" || sname == "ZJets_nunu" || sname == "ZJets_nunu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    // ZG + jets, Z->nunu
    if (sname == "ZGJets_nunu_nlo_inclusive") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZGTo2NuG_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "ZGJets_nunu_pTG_40-130" || sname == "ZGJets_nunu" || sname == "ZGJets_nunu_pTG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZNuNuGJets_MonoPhoton_PtG-40to130_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "ZGJets_nunu_pTG_130-inf" || sname == "ZGJets_nunu" || sname == "ZGJets_nunu_pTG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZNuNuGJets_MonoPhoton_PtG-130_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "ZGJets_nunu_nlo_pTG_130-inf" || sname == "ZGJets_nunu_nlo" || sname == "ZGJets_nunu_nlo_pTG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZGTo2NuG_PtG-130_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    // ZG + jets, Z->ll
    if (sname == "ZGJets_ll_pTG_40-130" || sname == "ZGJets_ll" || sname == "ZGJets_ll_pTG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZLLGJets_MonoPhoton_PtG-40to130_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });

    if (sname == "ZGJets_ll_pTG_130-inf" || sname == "ZGJets_ll" || sname == "ZGJets_ll_pTG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZLLGJets_MonoPhoton_PtG-130_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "ZGJets_ll_nlo_pTG_130-inf" || sname == "ZGJets_ll_nlo" || sname == "ZGJets_ll_nlo_pTG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZGTo2LG_PtG-130_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    // QCD mutijets
    if (sname == "QCD_HT_50-100" || sname == "QCD_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/QCD_HT50to100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "QCD_HT_100-200" || sname == "QCD_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "QCD_HT_200-300" || sname == "QCD_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "QCD_HT_300-500" || sname == "QCD_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "QCD_HT_500-700" || sname == "QCD_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "QCD_HT_700-1000" || sname == "QCD_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "QCD_HT_1000-1500" || sname == "QCD_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "QCD_HT_1500-2000" || sname == "QCD_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "QCD_HT_2000-inf" || sname == "QCD_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    // tt + jets
    if (sname == "TTJets") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v1/MINIAODSIM" });
  }
  else if (theDataYear == 2017){
    // Data
    if (
      sname == "DoubleEG" || sname == "DoubleMuon" || sname == "MuonEG"
      ||
      sname == "SingleElectron" || sname == "SingleMuon" || sname == "SinglePhoton"
      ||
      sname == "JetHT" || sname == "HTMHT" || sname == "MET"
      ||
      sname == "EGamma"
      ){
      for (auto const& dp:SampleHelpers::getValidDataPeriods()) constructSamplesList(sname+"_"+dp, syst, samples);
    }
    for (auto const& dp:SampleHelpers::getValidDataPeriods()){
      if (sname == Form("EGamma_%s", dp.Data())){
        constructSamplesList(Form("DoubleEG_%s", dp.Data()), syst, samples);
        constructSamplesList(Form("SingleElectron_%s", dp.Data()), syst, samples);
        constructSamplesList(Form("SinglePhoton_%s", dp.Data()), syst, samples);
      }
    }
    // Group by runs
    if (sname == "Run2017B"){
      constructSamplesList("DoubleEG_2017B", syst, samples);
      constructSamplesList("DoubleMuon_2017B", syst, samples);
      constructSamplesList("MuonEG_2017B", syst, samples);
      constructSamplesList("SingleElectron_2017B", syst, samples);
      constructSamplesList("SingleMuon_2017B", syst, samples);
      constructSamplesList("SinglePhoton_2017B", syst, samples);
      constructSamplesList("JetHT_2017B", syst, samples);
      constructSamplesList("HTMHT_2017B", syst, samples);
      constructSamplesList("MET_2017B", syst, samples);
    }
    if (sname == "Run2017C"){
      constructSamplesList("DoubleEG_2017C", syst, samples);
      constructSamplesList("DoubleMuon_2017C", syst, samples);
      constructSamplesList("MuonEG_2017C", syst, samples);
      constructSamplesList("SingleElectron_2017C", syst, samples);
      constructSamplesList("SingleMuon_2017C", syst, samples);
      constructSamplesList("SinglePhoton_2017C", syst, samples);
      constructSamplesList("JetHT_2017C", syst, samples);
      constructSamplesList("HTMHT_2017C", syst, samples);
      constructSamplesList("MET_2017C", syst, samples);
    }
    if (sname == "Run2017D"){
      constructSamplesList("DoubleEG_2017D", syst, samples);
      constructSamplesList("DoubleMuon_2017D", syst, samples);
      constructSamplesList("MuonEG_2017D", syst, samples);
      constructSamplesList("SingleElectron_2017D", syst, samples);
      constructSamplesList("SingleMuon_2017D", syst, samples);
      constructSamplesList("SinglePhoton_2017D", syst, samples);
      constructSamplesList("JetHT_2017D", syst, samples);
      constructSamplesList("HTMHT_2017D", syst, samples);
      constructSamplesList("MET_2017D", syst, samples);
    }
    if (sname == "Run2017E"){
      constructSamplesList("DoubleEG_2017E", syst, samples);
      constructSamplesList("DoubleMuon_2017E", syst, samples);
      constructSamplesList("MuonEG_2017E", syst, samples);
      constructSamplesList("SingleElectron_2017E", syst, samples);
      constructSamplesList("SingleMuon_2017E", syst, samples);
      constructSamplesList("SinglePhoton_2017E", syst, samples);
      constructSamplesList("JetHT_2017E", syst, samples);
      constructSamplesList("HTMHT_2017E", syst, samples);
      constructSamplesList("MET_2017E", syst, samples);
    }
    if (sname == "Run2017F"){
      constructSamplesList("DoubleEG_2017F", syst, samples);
      constructSamplesList("DoubleMuon_2017F", syst, samples);
      constructSamplesList("MuonEG_2017F", syst, samples);
      constructSamplesList("SingleElectron_2017F", syst, samples);
      constructSamplesList("SingleMuon_2017F", syst, samples);
      constructSamplesList("SinglePhoton_2017F", syst, samples);
      constructSamplesList("JetHT_2017F", syst, samples);
      constructSamplesList("HTMHT_2017F", syst, samples);
      constructSamplesList("MET_2017F", syst, samples);
    }
    if (sname == "DoubleEG_2017B") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DoubleEG/Run2017B-31Mar2018-v1/MINIAOD" });
    if (sname == "DoubleEG_2017C") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DoubleEG/Run2017C-31Mar2018-v1/MINIAOD" });
    if (sname == "DoubleEG_2017D") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DoubleEG/Run2017D-31Mar2018-v1/MINIAOD" });
    if (sname == "DoubleEG_2017E") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DoubleEG/Run2017E-31Mar2018-v1/MINIAOD" });
    if (sname == "DoubleEG_2017F") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DoubleEG/Run2017F-31Mar2018-v1/MINIAOD" });
    if (sname == "DoubleMuon_2017B") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DoubleMuon/Run2017B-31Mar2018-v1/MINIAOD" });
    if (sname == "DoubleMuon_2017C") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DoubleMuon/Run2017C-31Mar2018-v1/MINIAOD" });
    if (sname == "DoubleMuon_2017D") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DoubleMuon/Run2017D-31Mar2018-v1/MINIAOD" });
    if (sname == "DoubleMuon_2017E") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DoubleMuon/Run2017E-31Mar2018-v1/MINIAOD" });
    if (sname == "DoubleMuon_2017F") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DoubleMuon/Run2017F-31Mar2018-v1/MINIAOD" });
    if (sname == "MuonEG_2017B") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MuonEG/Run2017B-31Mar2018-v1/MINIAOD" });
    if (sname == "MuonEG_2017C") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MuonEG/Run2017C-31Mar2018-v1/MINIAOD" });
    if (sname == "MuonEG_2017D") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MuonEG/Run2017D-31Mar2018-v1/MINIAOD" });
    if (sname == "MuonEG_2017E") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MuonEG/Run2017E-31Mar2018-v1/MINIAOD" });
    if (sname == "MuonEG_2017F") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MuonEG/Run2017F-31Mar2018-v1/MINIAOD" });
    if (sname == "SingleElectron_2017B") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SingleElectron/Run2017B-31Mar2018-v1/MINIAOD" });
    if (sname == "SingleElectron_2017C") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SingleElectron/Run2017C-31Mar2018-v1/MINIAOD" });
    if (sname == "SingleElectron_2017D") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SingleElectron/Run2017D-31Mar2018-v1/MINIAOD" });
    if (sname == "SingleElectron_2017E") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SingleElectron/Run2017E-31Mar2018-v1/MINIAOD" });
    if (sname == "SingleElectron_2017F") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SingleElectron/Run2017F-31Mar2018-v1/MINIAOD" });
    if (sname == "SingleMuon_2017B") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SingleMuon/Run2017B-31Mar2018-v1/MINIAOD" });
    if (sname == "SingleMuon_2017C") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SingleMuon/Run2017C-31Mar2018-v1/MINIAOD" });
    if (sname == "SingleMuon_2017D") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SingleMuon/Run2017D-31Mar2018-v1/MINIAOD" });
    if (sname == "SingleMuon_2017E") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SingleMuon/Run2017E-31Mar2018-v1/MINIAOD" });
    if (sname == "SingleMuon_2017F") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SingleMuon/Run2017F-31Mar2018-v1/MINIAOD" });
    if (sname == "SinglePhoton_2017B") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SinglePhoton/Run2017B-31Mar2018-v1/MINIAOD" });
    if (sname == "SinglePhoton_2017C") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SinglePhoton/Run2017C-31Mar2018-v1/MINIAOD" });
    if (sname == "SinglePhoton_2017D") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SinglePhoton/Run2017D-31Mar2018-v1/MINIAOD" });
    if (sname == "SinglePhoton_2017E") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SinglePhoton/Run2017E-31Mar2018-v1/MINIAOD" });
    if (sname == "SinglePhoton_2017F") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SinglePhoton/Run2017F-31Mar2018-v1/MINIAOD" });
    if (sname == "JetHT_2017B") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/JetHT/Run2017B-31Mar2018-v1/MINIAOD" });
    if (sname == "JetHT_2017C") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/JetHT/Run2017C-31Mar2018-v1/MINIAOD" });
    if (sname == "JetHT_2017D") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/JetHT/Run2017D-31Mar2018-v1/MINIAOD" });
    if (sname == "JetHT_2017E") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/JetHT/Run2017E-31Mar2018-v1/MINIAOD" });
    if (sname == "JetHT_2017F") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/JetHT/Run2017F-31Mar2018-v1/MINIAOD" });
    if (sname == "HTMHT_2017B") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/HTMHT/Run2017B-31Mar2018-v1/MINIAOD" });
    if (sname == "HTMHT_2017C") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/HTMHT/Run2017C-31Mar2018-v1/MINIAOD" });
    if (sname == "HTMHT_2017D") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/HTMHT/Run2017D-31Mar2018-v1/MINIAOD" });
    if (sname == "HTMHT_2017E") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/HTMHT/Run2017E-31Mar2018-v1/MINIAOD" });
    if (sname == "HTMHT_2017F") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/HTMHT/Run2017F-31Mar2018-v1/MINIAOD" });
    if (sname == "MET_2017B") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MET/Run2017B-31Mar2018-v1/MINIAOD" });
    if (sname == "MET_2017C") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MET/Run2017C-31Mar2018-v1/MINIAOD" });
    if (sname == "MET_2017D") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MET/Run2017D-31Mar2018-v1/MINIAOD" });
    if (sname == "MET_2017E") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MET/Run2017E-31Mar2018-v1/MINIAOD" });
    if (sname == "MET_2017F") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MET/Run2017F-31Mar2018-v1/MINIAOD" });

    if (sname == "DY_2l_M_10to50" || sname == "DY_2l") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v2/MINIAODSIM" });
    if (sname == "DY_2l_M_10to50_ext" || sname == "DY_2l") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2/MINIAODSIM" });
    if (sname == "DY_2l_M_50" || sname == "DY_2l") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "DY_2l_M_50_ext" || sname == "DY_2l") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM" });
    if (sname == "DY_2l_M_50_ext2" || sname == "DY_2l") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext2-v1/MINIAODSIM" });
    if (sname == "TT_2l2nu") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v2/MINIAODSIM" });
    if (sname == "TTJets") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "ST_t-channel_top_5f") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ST_t-channel_top_5f_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM" });
    if (sname == "ST_t-channel_antitop_5f") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ST_t-channel_antitop_5f_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM" });
    if (sname == "ST_s-channel_top_leptonDecays") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ST_s-channel_top_leptonDecays_13TeV-PSweights_powheg-pythia/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "ST_s-channel_antitop_leptonDecays") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ST_s-channel_antitop_leptonDecays_13TeV-PSweights_powheg-pythia/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "ST_s-channel_4f_leptonDecays") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ST_s-channel_4f_leptonDecays_TuneCP5_PSweights_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "ST_tW_top_5f_NoFullyHadronicDecays") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "ST_tW_antitop_5f_NoFullyHadronicDecays") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "TZ_2l_4f") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/tZq_ll_4f_ckm_NLO_TuneCP5_PSweights_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v2/MINIAODSIM" });
    if (sname == "TTZ_2l2nu") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM" });
    if (sname == "TTW_lnu") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM" });
    if (sname == "qqZZ_2l2nu_mZ_18-inf" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZZTo2L2Nu_mZMin-18_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "qqZZ_2l2nu" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZZTo2L2Nu_13TeV_powheg_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "qqZZ_4l" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZZTo4L_13TeV_powheg_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "qqZZ_4l_ext" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZZTo4L_13TeV_powheg_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext2-v1/MINIAODSIM" });
    if (sname == "qqZZ_2q2nu" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZZTo2Q2Nu_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "qqZZ_2l2q" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "qqWZ_lnu2q" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM" });
    if (sname == "qqWZ_2l2q" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "qqWZ_13nu" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8_v2/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "qqWZ_3lnu_POWHEG_mll_0p1-inf" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WZTo3LNu_mllmin01_NNPDF31_TuneCP5_13TeV_powheg_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v4/MINIAODSIM" });
    if (sname == "qqWZ_3lnu_POWHEG" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WZTo3LNu_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM" });
    if (sname == "qqWZ_3lnu_MG" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (
      syst == SystematicsHelpers::tPDFScaleDn || syst == SystematicsHelpers::tPDFScaleUp
      ||
      syst == SystematicsHelpers::tQCDScaleDn || syst == SystematicsHelpers::tQCDScaleUp
      ||
      syst == SystematicsHelpers::tAsMZDn || syst == SystematicsHelpers::tAsMZUp
      ||
      syst == SystematicsHelpers::tPDFReplicaDn || syst == SystematicsHelpers::tPDFReplicaUp
      ){
      if (sname == "qqWW_lnu2q" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WWToLNuQQ_NNPDF31_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM" });
      if (sname == "qqWW_2l2nu" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WWTo2L2Nu_NNPDF31_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM" });
    }
    else if (syst == SystematicsHelpers::tPythiaTuneDn){
      if (sname == "qqWW_lnu2q" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
      if (sname == "qqWW_lnu2q_ext" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WWToLNuQQ_NNPDF31_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM" });
      if (sname == "qqWW_2l2nu" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WWTo2L2Nu_NNPDF31_TuneCP5Down_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
      if (sname == "qqWW_2l2nu_ext" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WWTo2L2Nu_NNPDF31_TuneCP5Down_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
      if (sname == "qqWW_2l2nu_ext2" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WWTo2L2Nu_NNPDF31_TuneCP5Down_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2/MINIAODSIM" });
    }
    else if (syst == SystematicsHelpers::tPythiaTuneUp){
      if (sname == "qqWW_lnu2q" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
      if (sname == "qqWW_lnu2q_ext" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WWToLNuQQ_NNPDF31_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM" });
      if (sname == "qqWW_2l2nu" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WWTo2L2Nu_NNPDF31_TuneCP5Up_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM" });
      if (sname == "qqWW_2l2nu_ext" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WWTo2L2Nu_NNPDF31_TuneCP5Up_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM" });
      if (sname == "qqWW_2l2nu_ext2" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WWTo2L2Nu_NNPDF31_TuneCP5Up_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2/MINIAODSIM" });
    }
    else{
      if (sname == "qqWW_lnu2q" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
      if (sname == "qqWW_lnu2q_ext" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WWToLNuQQ_NNPDF31_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM" });
      if (sname == "qqWW_2l2nu" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
      if (sname == "qqWW_2l2nu_ext" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WWTo2L2Nu_NNPDF31_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM" });
    }
    if (sname == "qqWG_lnu" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "WWW_4f" || sname == "VVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM" });
    if (sname == "WWZ_4f" || sname == "VVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM" });
    if (sname == "WZG" || sname == "VVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WZG_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "WZZ" || sname == "VVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "ZZZ" || sname == "VVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM" });

    if (sname == "GJets_HT_40-100" || sname == "GJets_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14_ext1-v2/MINIAODSIM" });
    if (sname == "GJets_HT_100-200" || sname == "GJets_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "GJets_HT_200-400" || sname == "GJets_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v2/MINIAODSIM" });
    if (sname == "GJets_HT_400-600" || sname == "GJets_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "GJets_HT_600-inf" || sname == "GJets_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM" });
    if (sname == "TGJets") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/TGJets_leptonDecays_TuneCP5_PSweights_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "TTGJets") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "WJets_lnu" || sname == "WJets_lnu_Njets" || sname == "WJets_lnu_0j") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM" });
    if (sname == "WJets_lnu" || sname == "WJets_lnu_Njets" || sname == "WJets_lnu_1j") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "WJets_lnu" || sname == "WJets_lnu_Njets" || sname == "WJets_lnu_1j_ext") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM" });
    if (sname == "WJets_lnu" || sname == "WJets_lnu_Njets" || sname == "WJets_lnu_2j") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "WJets_lnu" || sname == "WJets_lnu_Njets" || sname == "WJets_lnu_2j_ext") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2/MINIAODSIM" });
    if (sname == "WJets_lnu_HT_100-200" || sname == "WJets_lnu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM" });
    if (sname == "WJets_lnu_HT_200-400" || sname == "WJets_lnu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "WJets_lnu_HT_400-600" || sname == "WJets_lnu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "WJets_lnu_HT_600-800" || sname == "WJets_lnu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "WJets_lnu_HT_800-1200" || sname == "WJets_lnu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "WJets_lnu_HT_1200-2500" || sname == "WJets_lnu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "WJets_lnu_HT_2500-inf" || sname == "WJets_lnu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/MINIAODSIM" });
    if (sname == "ZJets_nunu_HT_100-200" || sname == "ZJets_nunu" || sname == "ZJets_nunu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZJetsToNuNu_HT-100To200_13TeV-madgraph/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "ZJets_nunu_HT_200-400" || sname == "ZJets_nunu" || sname == "ZJets_nunu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZJetsToNuNu_HT-200To400_13TeV-madgraph/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "ZJets_nunu_HT_400-600" || sname == "ZJets_nunu" || sname == "ZJets_nunu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZJetsToNuNu_HT-400To600_13TeV-madgraph/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "ZJets_nunu_HT_600-800" || sname == "ZJets_nunu" || sname == "ZJets_nunu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZJetsToNuNu_HT-600To800_13TeV-madgraph/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "ZJets_nunu_HT_800-1200" || sname == "ZJets_nunu" || sname == "ZJets_nunu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZJetsToNuNu_HT-800To1200_13TeV-madgraph/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "ZJets_nunu_HT_1200-2500" || sname == "ZJets_nunu" || sname == "ZJets_nunu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZJetsToNuNu_HT-1200To2500_13TeV-madgraph/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "ZJets_nunu_HT_2500-inf" || sname == "ZJets_nunu" || sname == "ZJets_nunu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "ZGJets_nunu_pTG_40-130" || sname == "ZGJets_nunu" || sname == "ZGJets_nunu_pTG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZNuNuGJets_MonoPhoton_PtG-40to130_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "ZGJets_nunu_pTG_130-inf" || sname == "ZGJets_nunu" || sname == "ZGJets_nunu_pTG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZNuNuGJets_MonoPhoton_PtG-130_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "ZGJets_nunu_nlo_inclusive") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZGTo2NuG_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "ZGJets_nunu_nlo_pTG_130-inf" || sname == "ZGJets_nunu_nlo" || sname == "ZGJets_nunu_nlo_pTG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZGTo2NuG_PtG-130_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "ZGJets_ll_pTG_40-130" || sname == "ZGJets_ll" || sname == "ZGJets_ll_pTG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZLLGJets_MonoPhoton_PtG-40to130_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "ZGJets_ll_pTG_130-inf" || sname == "ZGJets_ll" || sname == "ZGJets_ll_pTG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZLLGJets_MonoPhoton_PtG-130_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "ZGJets_ll_nlo_pTG_130-inf" || sname == "ZGJets_ll_nlo" || sname == "ZGJets_ll_nlo_pTG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZGTo2LG_PtG-130_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "QCD_HT_50-100" || sname == "QCD_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/QCD_HT50to100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "QCD_HT_100-200" || sname == "QCD_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/QCD_HT100to200_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM" });
    if (sname == "QCD_HT_200-300" || sname == "QCD_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/QCD_HT200to300_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v2/MINIAODSIM" });
    if (sname == "QCD_HT_300-500" || sname == "QCD_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/QCD_HT300to500_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "QCD_HT_500-700" || sname == "QCD_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/QCD_HT500to700_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM" });
    if (sname == "QCD_HT_700-1000" || sname == "QCD_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/QCD_HT700to1000_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v2/MINIAODSIM" });
    if (sname == "QCD_HT_1000-1500" || sname == "QCD_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/QCD_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
    if (sname == "QCD_HT_1500-2000" || sname == "QCD_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/QCD_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM" });
    if (sname == "QCD_HT_2000-inf" || sname == "QCD_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/QCD_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM" });
  }
  else if (theDataYear == 2016){
    // Data
    if (
      sname == "DoubleEG" || sname == "DoubleMuon" || sname == "MuonEG"
      ||
      sname == "SingleElectron" || sname == "SingleMuon" || sname == "SinglePhoton"
      ||
      sname == "JetHT" || sname == "HTMHT" || sname == "MET"
      ||
      sname == "EGamma"
      ){
      for (auto const& dp:SampleHelpers::getValidDataPeriods()) constructSamplesList(sname+"_"+dp, syst, samples);
    }
    for (auto const& dp:SampleHelpers::getValidDataPeriods()){
      if (sname == Form("EGamma_%s", dp.Data())){
        constructSamplesList(Form("DoubleEG_%s", dp.Data()), syst, samples);
        constructSamplesList(Form("SingleElectron_%s", dp.Data()), syst, samples);
        constructSamplesList(Form("SinglePhoton_%s", dp.Data()), syst, samples);
      }
    }
    // Group by runs
    if (sname == "Run2016B"){
      constructSamplesList("DoubleEG_2016B", syst, samples);
      constructSamplesList("DoubleMuon_2016B", syst, samples);
      constructSamplesList("MuonEG_2016B", syst, samples);
      constructSamplesList("SingleElectron_2016B", syst, samples);
      constructSamplesList("SingleMuon_2016B", syst, samples);
      constructSamplesList("SinglePhoton_2016B", syst, samples);
      constructSamplesList("JetHT_2016B", syst, samples);
      constructSamplesList("HTMHT_2016B", syst, samples);
      constructSamplesList("MET_2016B", syst, samples);
    }
    if (sname == "Run2016C"){
      constructSamplesList("DoubleEG_2016C", syst, samples);
      constructSamplesList("DoubleMuon_2016C", syst, samples);
      constructSamplesList("MuonEG_2016C", syst, samples);
      constructSamplesList("SingleElectron_2016C", syst, samples);
      constructSamplesList("SingleMuon_2016C", syst, samples);
      constructSamplesList("SinglePhoton_2016C", syst, samples);
      constructSamplesList("JetHT_2016C", syst, samples);
      constructSamplesList("HTMHT_2016C", syst, samples);
      constructSamplesList("MET_2016C", syst, samples);
    }
    if (sname == "Run2016D"){
      constructSamplesList("DoubleEG_2016D", syst, samples);
      constructSamplesList("DoubleMuon_2016D", syst, samples);
      constructSamplesList("MuonEG_2016D", syst, samples);
      constructSamplesList("SingleElectron_2016D", syst, samples);
      constructSamplesList("SingleMuon_2016D", syst, samples);
      constructSamplesList("SinglePhoton_2016D", syst, samples);
      constructSamplesList("JetHT_2016D", syst, samples);
      constructSamplesList("HTMHT_2016D", syst, samples);
      constructSamplesList("MET_2016D", syst, samples);
    }
    if (sname == "Run2016E"){
      constructSamplesList("DoubleEG_2016E", syst, samples);
      constructSamplesList("DoubleMuon_2016E", syst, samples);
      constructSamplesList("MuonEG_2016E", syst, samples);
      constructSamplesList("SingleElectron_2016E", syst, samples);
      constructSamplesList("SingleMuon_2016E", syst, samples);
      constructSamplesList("SinglePhoton_2016E", syst, samples);
      constructSamplesList("JetHT_2016E", syst, samples);
      constructSamplesList("HTMHT_2016E", syst, samples);
      constructSamplesList("MET_2016E", syst, samples);
    }
    if (sname == "Run2016F"){
      constructSamplesList("DoubleEG_2016F", syst, samples);
      constructSamplesList("DoubleMuon_2016F", syst, samples);
      constructSamplesList("MuonEG_2016F", syst, samples);
      constructSamplesList("SingleElectron_2016F", syst, samples);
      constructSamplesList("SingleMuon_2016F", syst, samples);
      constructSamplesList("SinglePhoton_2016F", syst, samples);
      constructSamplesList("JetHT_2016F", syst, samples);
      constructSamplesList("HTMHT_2016F", syst, samples);
      constructSamplesList("MET_2016F", syst, samples);
    }
    if (sname == "Run2016G"){
      constructSamplesList("DoubleEG_2016G", syst, samples);
      constructSamplesList("DoubleMuon_2016G", syst, samples);
      constructSamplesList("MuonEG_2016G", syst, samples);
      constructSamplesList("SingleElectron_2016G", syst, samples);
      constructSamplesList("SingleMuon_2016G", syst, samples);
      constructSamplesList("SinglePhoton_2016G", syst, samples);
      constructSamplesList("JetHT_2016G", syst, samples);
      constructSamplesList("HTMHT_2016G", syst, samples);
      constructSamplesList("MET_2016G", syst, samples);
    }
    if (sname == "Run2016H"){
      constructSamplesList("DoubleEG_2016H", syst, samples);
      constructSamplesList("DoubleMuon_2016H", syst, samples);
      constructSamplesList("MuonEG_2016H", syst, samples);
      constructSamplesList("SingleElectron_2016H", syst, samples);
      constructSamplesList("SingleMuon_2016H", syst, samples);
      constructSamplesList("SinglePhoton_2016H", syst, samples);
      constructSamplesList("JetHT_2016H", syst, samples);
      constructSamplesList("HTMHT_2016H", syst, samples);
      constructSamplesList("MET_2016H", syst, samples);
    }
    if (sname == "DoubleEG_2016B") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DoubleEG/Run2016B-17Jul2018_ver2-v1/MINIAOD" });
    if (sname == "DoubleEG_2016C") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DoubleEG/Run2016C-17Jul2018-v1/MINIAOD" });
    if (sname == "DoubleEG_2016D") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DoubleEG/Run2016D-17Jul2018-v1/MINIAOD" });
    if (sname == "DoubleEG_2016E") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DoubleEG/Run2016E-17Jul2018-v1/MINIAOD" });
    if (sname == "DoubleEG_2016F") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DoubleEG/Run2016F-17Jul2018-v1/MINIAOD" });
    if (sname == "DoubleEG_2016G") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DoubleEG/Run2016G-17Jul2018-v1/MINIAOD" });
    if (sname == "DoubleEG_2016H") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DoubleEG/Run2016H-17Jul2018-v1/MINIAOD" });
    if (sname == "DoubleMuon_2016B") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DoubleMuon/Run2016B-17Jul2018_ver2-v1/MINIAOD" });
    if (sname == "DoubleMuon_2016C") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DoubleMuon/Run2016C-17Jul2018-v1/MINIAOD" });
    if (sname == "DoubleMuon_2016D") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DoubleMuon/Run2016D-17Jul2018-v1/MINIAOD" });
    if (sname == "DoubleMuon_2016E") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DoubleMuon/Run2016E-17Jul2018-v1/MINIAOD" });
    if (sname == "DoubleMuon_2016F") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DoubleMuon/Run2016F-17Jul2018-v1/MINIAOD" });
    if (sname == "DoubleMuon_2016G") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DoubleMuon/Run2016G-17Jul2018-v1/MINIAOD" });
    if (sname == "DoubleMuon_2016H") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DoubleMuon/Run2016H-17Jul2018-v1/MINIAOD" });
    if (sname == "MuonEG_2016B") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MuonEG/Run2016B-17Jul2018_ver2-v1/MINIAOD" });
    if (sname == "MuonEG_2016C") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MuonEG/Run2016C-17Jul2018-v1/MINIAOD" });
    if (sname == "MuonEG_2016D") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MuonEG/Run2016D-17Jul2018-v1/MINIAOD" });
    if (sname == "MuonEG_2016E") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MuonEG/Run2016E-17Jul2018-v2/MINIAOD" });
    if (sname == "MuonEG_2016F") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MuonEG/Run2016F-17Jul2018-v1/MINIAOD" });
    if (sname == "MuonEG_2016G") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MuonEG/Run2016G-17Jul2018-v1/MINIAOD" });
    if (sname == "MuonEG_2016H") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MuonEG/Run2016H-17Jul2018-v1/MINIAOD" });
    if (sname == "SingleElectron_2016B") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SingleElectron/Run2016B-17Jul2018_ver2-v1/MINIAOD" });
    if (sname == "SingleElectron_2016C") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SingleElectron/Run2016C-17Jul2018-v1/MINIAOD" });
    if (sname == "SingleElectron_2016D") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SingleElectron/Run2016D-17Jul2018-v1/MINIAOD" });
    if (sname == "SingleElectron_2016E") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SingleElectron/Run2016E-17Jul2018-v1/MINIAOD" });
    if (sname == "SingleElectron_2016F") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SingleElectron/Run2016F-17Jul2018-v1/MINIAOD" });
    if (sname == "SingleElectron_2016G") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SingleElectron/Run2016G-17Jul2018-v1/MINIAOD" });
    if (sname == "SingleElectron_2016H") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SingleElectron/Run2016H-17Jul2018-v1/MINIAOD" });
    if (sname == "SingleMuon_2016B") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SingleMuon/Run2016B-17Jul2018_ver2-v1/MINIAOD" });
    if (sname == "SingleMuon_2016C") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SingleMuon/Run2016C-17Jul2018-v1/MINIAOD" });
    if (sname == "SingleMuon_2016D") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SingleMuon/Run2016D-17Jul2018-v1/MINIAOD" });
    if (sname == "SingleMuon_2016E") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SingleMuon/Run2016E-17Jul2018-v1/MINIAOD" });
    if (sname == "SingleMuon_2016F") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SingleMuon/Run2016F-17Jul2018-v1/MINIAOD" });
    if (sname == "SingleMuon_2016G") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SingleMuon/Run2016G-17Jul2018-v1/MINIAOD" });
    if (sname == "SingleMuon_2016H") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SingleMuon/Run2016H-17Jul2018-v1/MINIAOD" });
    if (sname == "SinglePhoton_2016B") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SinglePhoton/Run2016B-17Jul2018_ver2-v1/MINIAOD" });
    if (sname == "SinglePhoton_2016C") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SinglePhoton/Run2016C-17Jul2018-v1/MINIAOD" });
    if (sname == "SinglePhoton_2016D") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SinglePhoton/Run2016D-17Jul2018-v1/MINIAOD" });
    if (sname == "SinglePhoton_2016E") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SinglePhoton/Run2016E-17Jul2018-v1/MINIAOD" });
    if (sname == "SinglePhoton_2016F") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SinglePhoton/Run2016F-17Jul2018-v1/MINIAOD" });
    if (sname == "SinglePhoton_2016G") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SinglePhoton/Run2016G-17Jul2018-v1/MINIAOD" });
    if (sname == "SinglePhoton_2016H") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SinglePhoton/Run2016H-17Jul2018-v1/MINIAOD" });
    if (sname == "JetHT_2016B") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/JetHT/Run2016B-17Jul2018_ver2-v2/MINIAOD" });
    if (sname == "JetHT_2016C") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/JetHT/Run2016C-17Jul2018-v1/MINIAOD" });
    if (sname == "JetHT_2016D") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/JetHT/Run2016D-17Jul2018-v1/MINIAOD" });
    if (sname == "JetHT_2016E") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/JetHT/Run2016E-17Jul2018-v1/MINIAOD" });
    if (sname == "JetHT_2016F") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/JetHT/Run2016F-17Jul2018-v1/MINIAOD" });
    if (sname == "JetHT_2016G") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/JetHT/Run2016G-17Jul2018-v1/MINIAOD" });
    if (sname == "JetHT_2016H") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/JetHT/Run2016H-17Jul2018-v1/MINIAOD" });
    if (sname == "HTMHT_2016B") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/HTMHT/Run2016B-17Jul2018_ver2-v1/MINIAOD" });
    if (sname == "HTMHT_2016C") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/HTMHT/Run2016C-17Jul2018-v1/MINIAOD" });
    if (sname == "HTMHT_2016D") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/HTMHT/Run2016D-17Jul2018-v1/MINIAOD" });
    if (sname == "HTMHT_2016E") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/HTMHT/Run2016E-17Jul2018-v1/MINIAOD" });
    if (sname == "HTMHT_2016F") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/HTMHT/Run2016F-17Jul2018-v1/MINIAOD" });
    if (sname == "HTMHT_2016G") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/HTMHT/Run2016G-17Jul2018-v1/MINIAOD" });
    if (sname == "HTMHT_2016H") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/HTMHT/Run2016H-17Jul2018-v1/MINIAOD" });
    if (sname == "MET_2016B") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MET/Run2016B-17Jul2018_ver2-v1/MINIAOD" });
    if (sname == "MET_2016C") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MET/Run2016C-17Jul2018-v1/MINIAOD" });
    if (sname == "MET_2016D") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MET/Run2016D-17Jul2018-v1/MINIAOD" });
    if (sname == "MET_2016E") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MET/Run2016E-17Jul2018-v1/MINIAOD" });
    if (sname == "MET_2016F") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MET/Run2016F-17Jul2018-v1/MINIAOD" });
    if (sname == "MET_2016G") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MET/Run2016G-17Jul2018-v1/MINIAOD" });
    if (sname == "MET_2016H") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MET/Run2016H-17Jul2018-v2/MINIAOD" });

    // Simulation for the main signals
    // GGH ZZ 2l2nu POWHEG

    // VBF ZZ 2l2nu POWHEG

    // GGH WW 2l2nu POWHEG

    // VBF ZZ 2l2nu POWHEG

    // Simulation for the main backgrounds
    // Z + jets, Z->ll
    if (sname == "DY_2l_M_10to50" || sname == "DY_2l") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM" });
    if (sname == "DY_2l_M_50" || sname == "DY_2l") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/MINIAODSIM" });
    // Bkgs. with tops
    if (sname == "TT_2l2nu") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM" });
    if (sname == "ST_t-channel_top_4f") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM" });
    if (sname == "ST_t-channel_antitop_4f") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM" });
    if (sname == "ST_s-channel_4f_leptonDecays") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM" });
    if (sname == "ST_tW_top_5f_NoFullyHadronicDecays") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM" });
    if (sname == "ST_tW_antitop_5f_NoFullyHadronicDecays") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM" });
    if (sname == "TZ_2l_4f") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/tZq_ll_4f_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/MINIAODSIM" });
    if (sname == "TTZ_2q") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM" });
    if (sname == "TTZ_2l2nu_M_1to10" || sname == "TTZ_2l2nu") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/TTZToLL_M-1to10_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-94X_mcRun2_asymptotic_v3-v1/MINIAODSIM" });
    if (sname == "TTZ_2l2nu_M_10" || sname == "TTZ_2l2nu") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/MINIAODSIM" });
    if (sname == "TTW_lnu") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/MINIAODSIM" });
    // qqVV
    if (sname == "qqZZ_2l2nu" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZZTo2L2Nu_13TeV_powheg_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM" });
    if (sname == "qqZZ_2l2nu_ext" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZZTo2L2Nu_13TeV_powheg_pythia8_ext1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM" });
    if (sname == "qqWW_lnu2q" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WWToLNuQQ_13TeV-powheg/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM" });
    if (sname == "qqWW_lnu2q_MG" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM" });
    if (sname == "qqWW_2l2nu" || sname == "qqVV"){
      if (syst == tPythiaTuneDn) HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WWTo2L2Nu_13TeV-powheg-CUETP8M1Down/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM" });
      else if (syst == tPythiaTuneUp) HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WWTo2L2Nu_13TeV-powheg-CUETP8M1Up/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM" });
      else HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WWTo2L2Nu_13TeV-powheg/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM" });
    }
    if (sname == "qqZZ_4l") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZZTo4L_13TeV_powheg_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM" });
    if (sname == "qqZZ_4l_ext" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZZTo4L_13TeV_powheg_pythia8_ext1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM" });
    if (sname == "qqZZ_2q2nu" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZZTo2Q2Nu_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM" });
    if (sname == "qqZZ_2l2q" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZZTo2L2Q_13TeV_powheg_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM" });
    if (sname == "qqWZ_lnu2q" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM" });
    if (sname == "qqWZ_2l2q" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM" });
    if (sname == "qqWZ_3lnu_POWHEG_mll_0p1-inf" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WZTo3LNu_mllmin01_13TeV-powheg-pythia8_ext1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM" });
    if (sname == "qqWZ_3lnu_POWHEG" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/MINIAODSIM" });
    if (sname == "qqWZ_3lnu_MG" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM" });
    if (sname == "qqWZ_13nu" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM" });
    if (sname == "qqWG_lnu" || sname == "qqVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/MINIAODSIM" });
    // VVV
    if (sname == "WWW_4f" || sname == "VVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM" });
    if (sname == "WWZ" || sname == "VVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM" });
    if (sname == "WZG" || sname == "VVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WZG_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM" });
    if (sname == "WZZ" || sname == "VVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM" });
    if (sname == "ZZZ" || sname == "VVV") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM" });

    // Simulation for the photon CR processes
    // G + jets
    if (sname == "GJets_HT_40-100" || sname == "GJets_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM" });
    if (sname == "GJets_HT_100-200" || sname == "GJets_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM" });
    if (sname == "GJets_HT_200-400" || sname == "GJets_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM" });
    if (sname == "GJets_HT_400-600" || sname == "GJets_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM" });
    if (sname == "GJets_HT_600-inf" || sname == "GJets_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM" });
    // T/TT + G + jets
    if (sname == "TGJets") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/TGJets_leptonDecays_13TeV_amcatnlo_madspin_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM" });
    if (sname == "TTGJets") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM" });
    // W + jets
    if (sname == "WJets_lnu" || sname == "WJets_lnu_inclusive") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM" });
    if (sname == "WJets_lnu" || sname == "WJets_lnu_inclusive_ext") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/MINIAODSIM" });
    if (sname == "WJets_lnu_HT_100-200" || sname == "WJets_lnu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM" });
    if (sname == "WJets_lnu_HT_200-400" || sname == "WJets_lnu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2/MINIAODSIM" });
    if (sname == "WJets_lnu_HT_400-600" || sname == "WJets_lnu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM" });
    if (sname == "WJets_lnu_HT_600-800" || sname == "WJets_lnu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM" });
    if (sname == "WJets_lnu_HT_800-1200" || sname == "WJets_lnu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM" });
    if (sname == "WJets_lnu_HT_1200-2500" || sname == "WJets_lnu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM" });
    if (sname == "WJets_lnu_HT_2500-inf" || sname == "WJets_lnu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM" });
    // Z + jets, Z->nunu
    if (sname == "ZJets_nunu_HT_100-200" || sname == "ZJets_nunu" || sname == "ZJets_nunu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZJetsToNuNu_HT-100To200_13TeV-madgraph/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM" });
    if (sname == "ZJets_nunu_HT_200-400" || sname == "ZJets_nunu" || sname == "ZJets_nunu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZJetsToNuNu_HT-200To400_13TeV-madgraph/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM" });
    if (sname == "ZJets_nunu_HT_400-600" || sname == "ZJets_nunu" || sname == "ZJets_nunu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZJetsToNuNu_HT-400To600_13TeV-madgraph/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM" });
    if (sname == "ZJets_nunu_HT_600-800" || sname == "ZJets_nunu" || sname == "ZJets_nunu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZJetsToNuNu_HT-600To800_13TeV-madgraph/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM" });
    if (sname == "ZJets_nunu_HT_800-1200" || sname == "ZJets_nunu" || sname == "ZJets_nunu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZJetsToNuNu_HT-800To1200_13TeV-madgraph/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM" });
    if (sname == "ZJets_nunu_HT_1200-2500" || sname == "ZJets_nunu" || sname == "ZJets_nunu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZJetsToNuNu_HT-1200To2500_13TeV-madgraph/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM" });
    if (sname == "ZJets_nunu_HT_2500-inf" || sname == "ZJets_nunu" || sname == "ZJets_nunu_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM" });
    // ZG + jets, Z->nunu
    if (sname == "ZGJets_nunu_pTG_40-130" || sname == "ZGJets_nunu" || sname == "ZGJets_nunu_pTG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZNuNuGJets_MonoPhoton_PtG-40to130_TuneCUETP8M1_13TeV-madgraph/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM" });
    if (sname == "ZGJets_nunu_nlo_inclusive" || sname == "ZGJets_nunu_nlo") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZGTo2NuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM" });
    if (sname == "ZGJets_nunu_nlo_pTG_130-inf" || sname == "ZGJets_nunu_nlo" || sname == "ZGJets_nunu_nlo_pTG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZGTo2NuG_PtG-130_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM" });
    // ZG + jets, Z->ll
    if (sname == "ZGJets_ll_nlo_inclusive" || sname == "ZGJets_ll_nlo") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/MINIAODSIM" });
    if (sname == "ZGJets_ll_nlo_pTG_130-inf" || sname == "ZGJets_ll_nlo" || sname == "ZGJets_ll_nlo_pTG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZGTo2LG_PtG-130_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM" });
    // QCD mutijets
    if (sname == "QCD_HT_50-100" || sname == "QCD_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/QCD_HT50to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM" });
    if (sname == "QCD_HT_100-200" || sname == "QCD_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM" });
    if (sname == "QCD_HT_200-300" || sname == "QCD_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM" });
    if (sname == "QCD_HT_300-500" || sname == "QCD_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM" });
    if (sname == "QCD_HT_500-700" || sname == "QCD_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM" });
    if (sname == "QCD_HT_700-1000" || sname == "QCD_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM" });
    if (sname == "QCD_HT_1000-1500" || sname == "QCD_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM" });
    if (sname == "QCD_HT_1500-2000" || sname == "QCD_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM" });
    if (sname == "QCD_HT_2000-inf" || sname == "QCD_HT") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM" });
    // tt + jets
    if (sname == "TTJets") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM" });
  }


}
void SampleHelpers::getSamplesList(std::vector<TString> const& s, std::vector<TString>& vs, SystematicsHelpers::SystematicVariationTypes syst, std::vector<size_t>* ns){
  for (auto const& ss:s){
    std::vector<TString> dumappend; constructSamplesList(ss, syst, dumappend);
    HelperFunctions::appendVector<TString>(vs, dumappend);
    if (ns) ns->push_back(dumappend.size());
  }
}
