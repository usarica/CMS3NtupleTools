// For the messier stuff
#include <cassert>
#include "OffshellSampleHelpers.h"
#include "HelperFunctions.h"
#include "HostHelpersCore.h"
#include "HiggsXSBRReader.h"
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

  if (sname == Form("Run%i", theDataYear)){
    for (auto const& pp:getValidDataPeriods()) constructSamplesList(Form("Run%s", pp.Data()), syst, samples);
  }

  if (theDataYear == 2018){
    // Data
    if (
      sname == "SingleMuon" || sname == "DoubleMuon" || sname == "MuonEG"
      ||
      sname == "EGamma"
      ||
      sname == "JetHT" || sname == "MET"
      ||
      sname == "JetMET" || sname == "SingleLepton" // Composite names
      ){
      for (auto const& dp:SampleHelpers::getValidDataPeriods()) constructSamplesList(sname+"_"+dp, syst, samples);
    }
    for (auto const& dp:SampleHelpers::getValidDataPeriods()){
      if (sname == Form("JetMET_%s", dp.Data())){
        constructSamplesList(Form("JetHT_%s", dp.Data()), syst, samples);
        constructSamplesList(Form("MET_%s", dp.Data()), syst, samples);
      }
      else if (sname == Form("SingleLepton_%s", dp.Data())){
        constructSamplesList(Form("SingleMuon_%s", dp.Data()), syst, samples);
        constructSamplesList(Form("EGamma_%s", dp.Data()), syst, samples);
      }
    }
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
    if (sname == "SingleMuon_2018A") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SingleMuon/Run2018A-17Sep2018-v2/MINIAOD" });
    if (sname == "SingleMuon_2018B") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SingleMuon/Run2018B-17Sep2018-v1/MINIAOD" });
    if (sname == "SingleMuon_2018C") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SingleMuon/Run2018C-17Sep2018-v1/MINIAOD" });
    if (sname == "SingleMuon_2018D") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/SingleMuon/Run2018D-22Jan2019-v2/MINIAOD" });
    if (sname == "DoubleMuon_2018A") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DoubleMuon/Run2018A-17Sep2018-v2/MINIAOD" });
    if (sname == "DoubleMuon_2018B") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DoubleMuon/Run2018B-17Sep2018-v1/MINIAOD" });
    if (sname == "DoubleMuon_2018C") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DoubleMuon/Run2018C-17Sep2018-v1/MINIAOD" });
    if (sname == "DoubleMuon_2018D") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/DoubleMuon/Run2018D-PromptReco-v2/MINIAOD" });
    if (sname == "MuonEG_2018A") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MuonEG/Run2018A-17Sep2018-v1/MINIAOD" });
    if (sname == "MuonEG_2018B") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MuonEG/Run2018B-17Sep2018-v1/MINIAOD" });
    if (sname == "MuonEG_2018C") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MuonEG/Run2018C-17Sep2018-v1/MINIAOD" });
    if (sname == "MuonEG_2018D") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MuonEG/Run2018D-PromptReco-v2/MINIAOD" });
    if (sname == "EGamma_2018A") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/EGamma/Run2018A-17Sep2018-v2/MINIAOD" });
    if (sname == "EGamma_2018B") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/EGamma/Run2018B-17Sep2018-v1/MINIAOD" });
    if (sname == "EGamma_2018C") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/EGamma/Run2018C-17Sep2018-v1/MINIAOD" });
    if (sname == "EGamma_2018D") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/EGamma/Run2018D-22Jan2019-v2/MINIAOD" });
    if (sname == "MET_2018A") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MET/Run2018A-17Sep2018-v1/MINIAOD" });
    if (sname == "MET_2018B") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MET/Run2018B-17Sep2018-v1/MINIAOD" });
    if (sname == "MET_2018C") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MET/Run2018C-17Sep2018-v1/MINIAOD" });
    if (sname == "MET_2018D") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/MET/Run2018D-PromptReco-v2/MINIAOD" });
    if (sname == "JetHT_2018A") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/JetHT/Run2018A-17Sep2018-v1/MINIAOD" });
    if (sname == "JetHT_2018B") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/JetHT/Run2018B-17Sep2018-v1/MINIAOD" });
    if (sname == "JetHT_2018C") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/JetHT/Run2018C-17Sep2018-v1/MINIAOD" });
    if (sname == "JetHT_2018D") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/JetHT/Run2018D-PromptReco-v2/MINIAOD" });

    // Simulation for the main signals
    // POWHEG signal samples
    if (sname.Contains("POWHEG")){
      std::vector<TString> samplelist;

      if (sname.Contains("GGH")){
        // gg->H
        if (sname.Contains("ZZTo2L2Nu")){
          // H->ZZ->2l2nu
          if (syst==tHardJetsDn || syst==tHardJetsUp){
            samplelist = std::vector<TString>{
              "/GluGluHToZZTo2L2Nu_M125_TuneCP5_13TeV_powheg2_minloHJJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M300_TuneCP5_13TeV_powheg2_minloHJJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/GluGluHToZZTo2L2Nu_M125_TuneCP5_13TeV_tunedown_powheg2_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M300_TuneCP5_13TeV_tunedown_powheg2_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/GluGluHToZZTo2L2Nu_M125_TuneCP5_13TeV_tuneup_powheg2_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M300_TuneCP5_13TeV_tuneup_powheg2_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/GluGluHToZZTo2L2Nu_M125_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M160_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M170_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M180_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M190_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M200_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M210_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M230_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M250_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M270_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M300_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M350_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M400_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M450_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M500_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M550_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M600_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M700_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M800_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M900_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M1000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M1500_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M2000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M2500_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M3000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM"
            };
          }
        } // End H->ZZ->2l2nu
        else if (sname.Contains("ZZTo4L")){
          // H->ZZ->4l
          if (syst==tHardJetsDn || syst==tHardJetsUp){
            samplelist = std::vector<TString>{
              "/GluGluHToZZTo4L_M125_13TeV_powheg2_minloHJJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToZZTo4L_M300_13TeV_powheg2_minloHJJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/GluGluHToZZTo4L_M160_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToZZTo4L_M170_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToZZTo4L_M180_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToZZTo4L_M190_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToZZTo4L_M200_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToZZTo4L_M210_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToZZTo4L_M230_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToZZTo4L_M250_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToZZTo4L_M270_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToZZTo4L_M300_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToZZTo4L_M350_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToZZTo4L_M400_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToZZTo4L_M450_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToZZTo4L_M500_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToZZTo4L_M550_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToZZTo4L_M600_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToZZTo4L_M700_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToZZTo4L_M800_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToZZTo4L_M900_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToZZTo4L_M1000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToZZTo4L_M1500_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToZZTo4L_M2000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToZZTo4L_M2500_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToZZTo4L_M3000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
            };
          }
        } // End H->ZZ->4l
        else if (sname.Contains("WWTo2L2Nu")){
          // H->WW->2l2nu
          if (syst==tHardJetsDn || syst==tHardJetsUp){
            samplelist = std::vector<TString>{
              "/GluGluHToWWTo2L2N_M125_NNPDF31_TuneCP5_PSweights_13TeV_powheg_minloHJJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2N_M300_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_minloHJJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/GluGluHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8-CP5Down/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/GluGluHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8-CP5Up/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/GluGluHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M160_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M170_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M180_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M190_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M200_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M210_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M230_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M250_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M270_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M300_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M350_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M450_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M500_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M550_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M600_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M700_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M800_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M900_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M1000_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M1500_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M2000_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M2500_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M3000_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
            };
          }
        } // End H->WW->2l2nu
      } // End gg->H
      else if (sname.Contains("VBF")){
        // VBF
        if (sname.Contains("ZZTo2L2Nu")){
          // H->ZZ->2l2nu
          if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/VBF_HToZZTo2L2Nu_M125_TuneCP5_13TeV_tunedown_powheg2_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M300_TuneCP5_13TeV_tunedown_powheg2_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/VBF_HToZZTo2L2Nu_M125_TuneCP5_13TeV_tuneup_powheg2_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M300_TuneCP5_13TeV_tuneup_powheg2_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/VBF_HToZZTo2L2Nu_M125_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M160_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M170_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M180_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M190_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M200_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M210_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M230_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M250_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M300_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M350_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M400_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M450_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M500_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M550_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M600_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M700_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M800_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M900_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M1000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M1500_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M2000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M2500_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M3000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
            };
          }
        } // End H->ZZ->2l2nu
        else if (sname.Contains("ZZTo4L")){
          // H->ZZ->4l
          if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/VBF_HToZZTo4L_M125_13TeV_tunedown_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/VBF_HToZZTo4L_M125_13TeV_tuneup_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/VBF_HToZZTo4L_M160_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/VBF_HToZZTo4L_M170_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/VBF_HToZZTo4L_M180_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/VBF_HToZZTo4L_M190_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/VBF_HToZZTo4L_M200_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/VBF_HToZZTo4L_M210_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/VBF_HToZZTo4L_M230_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/VBF_HToZZTo4L_M250_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBF_HToZZTo4L_M270_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBF_HToZZTo4L_M300_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBF_HToZZTo4L_M350_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBF_HToZZTo4L_M400_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBF_HToZZTo4L_M450_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBF_HToZZTo4L_M500_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBF_HToZZTo4L_M550_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBF_HToZZTo4L_M600_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBF_HToZZTo4L_M700_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBF_HToZZTo4L_M800_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBF_HToZZTo4L_M900_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBF_HToZZTo4L_M1000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBF_HToZZTo4L_M1500_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBF_HToZZTo4L_M2000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBF_HToZZTo4L_M2500_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBF_HToZZTo4L_M3000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
            };
          }
        } // End H->ZZ->4l
        else if (sname.Contains("WWTo2L2Nu")){
          // H->WW->2l2nu
          if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/VBFHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8_CP5Down/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBF_HToWWTo2L2Nu_M300_TuneCP5_13TeV_tunedown_powheg2_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/VBFHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8_CP5Up/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBF_HToWWTo2L2Nu_M300_TuneCP5_13TeV_tuneup_powheg2_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              // This sample is a special one to check for the dipole recoil adjustment in VBF
              //"/VBFHToWWTo2L2Nu_M125_TuneCP5_13TeV_powheg2_JHUGenV714_pythia8_withDipoleRecoil/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M160_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M170_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M180_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M190_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M200_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M210_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M230_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M270_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M300_13TeV_powheg2_JHUGenV710_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M350_13TeV_powheg2_JHUGenV710_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M450_13TeV_powheg2_JHUGenV710_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M500_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M550_13TeV_powheg2_JHUGenV710_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M600_13TeV_powheg2_JHUGenV710_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M700_13TeV_powheg2_JHUGenV710_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M800_13TeV_powheg2_JHUGenV710_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M900_13TeV_powheg2_JHUGenV710_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M1000_13TeV_powheg2_JHUGenV710_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M1500_13TeV_powheg2_JHUGenV710_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M2000_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M2500_13TeV_powheg2_JHUGenV710_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M3000_13TeV_powheg2_JHUGenV710_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
            };
          }
        } // End H->WW->2l2nu
      } // End VBF
      else if (sname.Contains("ZH")){
        // ZH
        if (sname.Contains("To2Nu2X_2LFilter")){
          // H->ZZ->2nu2x
          if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/ZH_HToZZTo2Nu2X_2LFilter_M125_TuneCP5_13TeV_tunedown_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_M300_TuneCP5_13TeV_tunedown_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/ZH_HToZZTo2Nu2X_2LFilter_M125_TuneCP5_13TeV_tuneup_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_M300_TuneCP5_13TeV_tuneup_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M125_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M160_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M170_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M180_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M190_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M200_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HTo2Nu2X_2LFilter_4LVetoFilter_M210_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo2Nu2X_2LFilter_4LVetoFilter_M230_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M250_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HTo2Nu2X_2LFilter_4LVetoFilter_M270_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M300_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M350_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M400_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M450_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M500_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M550_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M600_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M700_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M800_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M900_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M1000_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HTo2Nu2X_2LFilter_4LVetoFilter_M1500_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M2000_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M2500_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M3000_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
            };
          }
        } // End H->ZZ->2nu2x
        else if (sname.Contains("To2L2Q_2LFilter")){
          // H->ZZ->2l2q
          if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M125_TuneCP5_13TeV_tunedown_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M300_TuneCP5_13TeV_tunedown_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M125_TuneCP5_13TeV_tuneup_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M300_TuneCP5_13TeV_tuneup_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M125_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M160_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M180_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M190_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M200_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M210_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M230_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M250_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M270_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M300_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M350_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M400_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M450_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M500_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M550_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M600_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M700_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M800_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M900_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              // MH=1000???
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M1500_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M2000_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M2500_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M3000_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM"
            };
          }
        } // End H->ZZ->2l2q
        else if (sname.Contains("To4Q_2LFilter") && !sname.Contains("WWTo4Q_2LFilter")){
          // H->ZZ->4q
          if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              // 125?
              "/ZH_HTo4Q_2LFilter_M300_TuneCP5_13TeV_tunedown_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/ZH_HTo4Q_2LFilter_M125_TuneCP5_13TeV_tuneup_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M300_TuneCP5_13TeV_tuneup_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/ZH_HTo4Q_2LFilter_M125_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M160_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M170_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M180_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M190_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M210_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M230_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M250_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M300_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M350_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M450_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M500_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M550_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M600_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M700_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M900_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M1000_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M1500_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M2000_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M2500_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M3000_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM"
            };
          }
        } // End H->ZZ->4q
        else if (sname.Contains("ZZ_4LFilter")){
          // H->ZZ->2l2x, 4l filter
          if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/ZH_HToZZ_4LFilter_M125_13TeV_tunedown_powheg2-minlo-HZJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/ZH_HToZZ_4LFilter_M125_13TeV_tuneup_powheg2-minlo-HZJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZ_4LFilter_M160_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZ_4LFilter_M170_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZ_4LFilter_M180_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/ZH_HToZZ_4LFilter_M190_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/ZH_HToZZ_4LFilter_M200_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/ZH_HToZZ_4LFilter_M210_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZ_4LFilter_M230_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZ_4LFilter_M250_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZ_4LFilter_M270_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZ_4LFilter_M300_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZ_4LFilter_M350_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZ_4LFilter_M400_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZ_4LFilter_M450_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/ZH_HToZZ_4LFilter_M500_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZ_4LFilter_M550_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZ_4LFilter_M600_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZ_4LFilter_M700_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZ_4LFilter_M800_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZ_4LFilter_M900_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZ_4LFilter_M1000_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZ_4LFilter_M1500_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZ_4LFilter_M2000_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZ_4LFilter_M2500_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/ZH_HToZZ_4LFilter_M3000_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
            };
          }
        } // End H->ZZ->2l2x, 4l filter
        else if (sname.Contains("WWTo2L2Nu")){
          // H->WW->2l2nu, better to use the lnuqq sample for the 3l contribution MC
          if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/HZJ_HToWW2L2NuInclusiveZ_M125_NNPDF31_TuneCP5Down_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M300_CPS_NNPDF31_TuneCP5Down_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/HZJ_HToWW2L2NuInclusiveZ_M125_NNPDF31_TuneCP5Up_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M300_CPS_NNPDF31_TuneCP5Up_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/HZJ_HToWW2L2NuInclusiveZ_M125_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M160_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M170_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M180_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M190_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M200_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M210_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M230_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M250_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M270_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M300_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M350_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M400_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M450_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M500_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M550_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M600_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M700_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M800_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M900_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M1000_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M1500_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M2000_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M2500_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M3000_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
            };
          }
        } // End H->WW->2l2nu
        else if (sname.Contains("LNuQQ_2LFilter")){
          // H->WW->lnuqq
          if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/ZH_HToWWToLNuQQ_2LFilter_M125_TuneCP5_13TeV_tunedown_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/ZH_HToWWToLNuQQ_2LFilter_M125_TuneCP5_13TeV_tuneup_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M300_TuneCP5_13TeV_tuneup_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/ZH_HToWWToLNuQQ_2LFilter_M125_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M160_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M170_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M180_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M190_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M200_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M210_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M230_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M250_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M270_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M300_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M1000_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M2500_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M3000_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM"
            };
          }
        } // End H->WW->lnuqq
      } // End ZH
      else if (sname.Contains("WminusH")){
        // WminusH
        if (sname.Contains("ZZTo2L2Nu")){
          // H->ZZ->2l2nu
          if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/WminusH_HToZZTo2L2Nu_M125_TuneCP5_13TeV_tunedown_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/WminusH_HToZZTo2L2Nu_M125_TuneCP5_13TeV_tuneup_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M300_TuneCP5_13TeV_tuneup_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/WminusH_HToZZTo2L2Nu_M125_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M160_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M170_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M180_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M190_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M200_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M210_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M230_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M250_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M270_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M300_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M350_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M400_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M450_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M500_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M550_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M600_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M700_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M800_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M900_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M1000_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M1500_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M2000_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M2500_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M3000_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
            };
          }
        } // End H->ZZ->2l2nu
        else if (sname.Contains("ZZTo2L2Q")){
          // H->ZZ->2l2q
          if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/WminusH_HToZZTo2L2Q_M125_TuneCP5_13TeV_tunedown_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M300_TuneCP5_13TeV_tunedown_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/WminusH_HToZZTo2L2Q_M125_TuneCP5_13TeV_tuneup_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M300_TuneCP5_13TeV_tuneup_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/WminusH_HToZZTo2L2Q_M125_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M160_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M170_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M180_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M190_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M200_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M210_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M230_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M250_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M270_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M300_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M350_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M400_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M450_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M500_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M550_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M600_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M700_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M800_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M900_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M2500_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM"
            };
          }
        } // End H->ZZ->2l2q
        else if (sname.Contains("ZZTo4L")){
          // H->ZZ->4l
          if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/WminusH_HToZZTo4L_M125_13TeV_tunedown_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/WminusH_HToZZTo4L_M125_13TeV_tuneup_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/WminusH_HToZZTo4L_M160_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/WminusH_HToZZTo4L_M170_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/WminusH_HToZZTo4L_M180_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/WminusH_HToZZTo4L_M190_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/WminusH_HToZZTo4L_M200_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/WminusH_HToZZTo4L_M210_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/WminusH_HToZZTo4L_M230_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/WminusH_HToZZTo4L_M250_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/WminusH_HToZZTo4L_M270_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/WminusH_HToZZTo4L_M300_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/WminusH_HToZZTo4L_M350_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/WminusH_HToZZTo4L_M400_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/WminusH_HToZZTo4L_M450_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/WminusH_HToZZTo4L_M500_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/WminusH_HToZZTo4L_M550_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/WminusH_HToZZTo4L_M600_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/WminusH_HToZZTo4L_M700_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/WminusH_HToZZTo4L_M800_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/WminusH_HToZZTo4L_M900_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/WminusH_HToZZTo4L_M1000_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/WminusH_HToZZTo4L_M1500_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/WminusH_HToZZTo4L_M2000_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/WminusH_HToZZTo4L_M2500_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/WminusH_HToZZTo4L_M3000_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM"
            };
          }
        } // End H->ZZ->4l
        else if (sname.Contains("WW_2LOSFilter")){
          // H->WW->lnu2x with 2l OS filter
          if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/WminusHToWW_2LOSFilter_M125_NNPDF31_TuneCP5Down_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M300_CPS_NNPDF31_TuneCP5Down_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/WminusHToWW_2LOSFilter_M125_NNPDF31_TuneCP5Up_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M300_CPS_NNPDF31_TuneCP5Up_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/WminusHToWW_2LOSFilter_M125_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M160_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M170_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M180_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M190_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M200_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M210_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M230_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M250_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M270_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M300_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M350_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M400_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M450_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M500_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M550_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M600_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M700_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M800_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M900_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M1000_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M1500_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M2000_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M2500_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M3000_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
            };
          }
        } // End H->WW->lnu2x
      } // End WminusH
      else if (sname.Contains("WplusH")){
        // WplusH
        if (sname.Contains("ZZTo2L2Nu")){
          // H->ZZ->2l2nu
          if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/WplusH_HToZZTo2L2Nu_M125_TuneCP5_13TeV_tunedown_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M300_TuneCP5_13TeV_tunedown_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/WplusH_HToZZTo2L2Nu_M125_TuneCP5_13TeV_tuneup_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M300_TuneCP5_13TeV_tuneup_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/WplusH_HToZZTo2L2Nu_M125_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M160_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M170_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M180_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M190_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M200_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M210_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M230_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M250_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M270_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M300_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M350_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M400_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M450_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M500_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M550_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M600_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M700_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M800_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M900_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M1000_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M1500_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M2000_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M2500_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M3000_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
            };
          }
        } // End H->ZZ->2l2nu
        else if (sname.Contains("ZZTo2L2Q")){
          // H->ZZ->2l2q
          if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/WplusH_HToZZTo2L2Q_M125_TuneCP5_13TeV_tunedown_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M300_TuneCP5_13TeV_tunedown_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/WplusH_HToZZTo2L2Q_M125_TuneCP5_13TeV_tuneup_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M300_TuneCP5_13TeV_tuneup_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/WplusH_HToZZTo2L2Q_M125_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M160_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M170_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M180_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M190_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M200_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M210_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M230_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M250_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M270_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M300_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M350_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M400_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M450_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M500_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M550_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M600_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M700_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M800_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M900_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M1000_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M1500_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M2000_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_private/MINIAODSIM"
            };
          }
        } // End H->ZZ->2l2q
        else if (sname.Contains("ZZTo4L")){
          // H->ZZ->4l
          if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/WplusH_HToZZTo4L_M125_13TeV_tunedown_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/WplusH_HToZZTo4L_M125_13TeV_tuneup_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/WplusH_HToZZTo4L_M160_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
              "/WplusH_HToZZTo4L_M170_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo4L_M180_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo4L_M190_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo4L_M200_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo4L_M210_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo4L_M230_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo4L_M250_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo4L_M270_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo4L_M300_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo4L_M350_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo4L_M400_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo4L_M450_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo4L_M500_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo4L_M550_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo4L_M600_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo4L_M700_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo4L_M800_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo4L_M900_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo4L_M1000_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo4L_M1500_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo4L_M2000_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo4L_M2500_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusH_HToZZTo4L_M3000_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
            };
          }
        } // End H->ZZ->4l
        else if (sname.Contains("WW_2LOSFilter")){
          // H->WW->lnu2x with 2l OS filter
          if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/WplusHToWW_2LOSFilter_M125_NNPDF31_TuneCP5Down_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M300_CPS_NNPDF31_TuneCP5Down_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/WplusHToWW_2LOSFilter_M125_NNPDF31_TuneCP5Up_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M300_CPS_NNPDF31_TuneCP5Up_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/WplusHToWW_2LOSFilter_M125_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M160_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M170_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M180_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M190_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M200_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M210_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M230_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M250_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M270_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M300_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M350_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M400_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M450_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M500_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M550_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M600_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M700_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M800_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M900_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M1000_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M1500_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M2000_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M2500_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M3000_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
            };
          }
        } // End H->WW->lnu2x
      } // End WplusH

      float const mh_req = SampleHelpers::findPoleMass(sname);
      for (unsigned int isample=0; isample<samplelist.size(); isample++){
        auto const& ss = samplelist.at(isample);
        float const mh_sample = SampleHelpers::findPoleMass(ss);
        if (isample>0){
          float const mh_sample_prev = SampleHelpers::findPoleMass(samplelist.at(isample-1));
          if (mh_sample<=mh_sample_prev){
            MELAerr << "SampleHelpers::constructSamplesList: Mass = " << mh_sample << " <= previous mass = " << mh_sample_prev << " during the request for " << sname << endl;
            assert(0);
          }
        }
        if (mh_req<0.f || mh_sample==mh_req) HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ ss });
      }
    }


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
    if (sname == "ZGJets_nunu_nlo_inclusive" || sname == "ZGJets_nunu_nlo") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZGTo2NuG_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "ZGJets_nunu_nlo_pTG_130-inf" || sname == "ZGJets_nunu_nlo" || sname == "ZGJets_nunu_nlo_pTG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZGTo2NuG_PtG-130_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" });
    if (sname == "ZGJets_nunu_pTG_40-130" || sname == "ZGJets_nunu" || sname == "ZGJets_nunu_pTG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZNuNuGJets_MonoPhoton_PtG-40to130_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
    if (sname == "ZGJets_nunu_pTG_130-inf" || sname == "ZGJets_nunu" || sname == "ZGJets_nunu_pTG") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZNuNuGJets_MonoPhoton_PtG-130_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM" });
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
      sname == "EGamma" || sname == "JetMET" || sname == "SingleLepton" // Composite names
      ){
      for (auto const& dp:SampleHelpers::getValidDataPeriods()) constructSamplesList(sname+"_"+dp, syst, samples);
    }
    for (auto const& dp:SampleHelpers::getValidDataPeriods()){
      if (sname == Form("EGamma_%s", dp.Data())){
        constructSamplesList(Form("DoubleEG_%s", dp.Data()), syst, samples);
        constructSamplesList(Form("SingleElectron_%s", dp.Data()), syst, samples);
        constructSamplesList(Form("SinglePhoton_%s", dp.Data()), syst, samples);
      }
      else if (sname == Form("JetMET_%s", dp.Data())){
        constructSamplesList(Form("JetHT_%s", dp.Data()), syst, samples);
        constructSamplesList(Form("HTMHT_%s", dp.Data()), syst, samples);
        constructSamplesList(Form("MET_%s", dp.Data()), syst, samples);
      }
      else if (sname == Form("SingleLepton_%s", dp.Data())){
        constructSamplesList(Form("SingleMuon_%s", dp.Data()), syst, samples);
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

    // Simulation for the main signals
    // POWHEG signal samples
    if (sname.Contains("POWHEG")){
      std::vector<TString> samplelist;

      if (sname.Contains("GGH")){
        // gg->H
        if (sname.Contains("ZZTo2L2Nu")){
          // H->ZZ->2l2nu
          if (syst==tHardJetsDn || syst==tHardJetsUp){
            samplelist = std::vector<TString>{
              "/GluGluHToZZTo2L2Nu_M125_TuneCP5_13TeV_powheg2_minloHJJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M300_TuneCP5_13TeV_powheg2_minloHJJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/GluGluHToZZTo2L2Nu_M125_TuneCP5_13TeV_tunedown_powheg2_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M300_TuneCP5_13TeV_tunedown_powheg2_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/GluGluHToZZTo2L2Nu_M125_TuneCP5_13TeV_tuneup_powheg2_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M300_TuneCP5_13TeV_tuneup_powheg2_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/GluGluHToZZTo2L2Nu_M125_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M160_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M170_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M180_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M200_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M210_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M230_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M250_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M270_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M300_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M350_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M400_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M450_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M500_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M550_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M600_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M700_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M800_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M900_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M1000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M1500_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M2000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M2500_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M3000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM"
            };
          }
        } // End H->ZZ->2l2nu
        else if (sname.Contains("WWTo2L2Nu")){
          // H->WW->2l2nu
          if (syst==tHardJetsDn || syst==tHardJetsUp){
            samplelist = std::vector<TString>{
            };
          }
          else if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/GluGluHToWWTo2L2Nu_M125_NNPDF31_TuneCP5Down_13TeV_powheg_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/GluGluHToWWTo2L2Nu_M125_NNPDF31_TuneCP5Up_13TeV_powheg_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/GluGluHToWWTo2L2Nu_M125_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M160_13TeV_powheg2_JHUGenV714_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M170_13TeV_powheg2_JHUGenV714_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M180_13TeV_powheg2_JHUGenV714_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M190_13TeV_powheg2_JHUGenV714_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M200_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M210_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M230_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M250_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M270_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M300_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGenV714_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M350_13TeV_powheg2_JHUGenV714_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M400_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGenV727_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M450_13TeV_powheg2_JHUGenV714_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M500_13TeV_powheg2_JHUGenV714_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M550_13TeV_powheg2_JHUGenV714_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M600_13TeV_powheg2_JHUGenV714_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M700_13TeV_powheg2_JHUGenV714_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M800_13TeV_powheg2_JHUGenV714_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M900_13TeV_powheg2_JHUGenV714_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M1000_13TeV_powheg2_JHUGenV714_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M1500_13TeV_powheg2_JHUGenV714_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M2000_13TeV_powheg2_JHUGenV714_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M2500_13TeV_powheg2_JHUGenV714_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M3000_13TeV_powheg2_JHUGenV714_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
            };
          }
        } // End H->WW->2l2nu
      } // End gg->H
      else if (sname.Contains("VBF")){
        // VBF
        if (sname.Contains("ZZTo2L2Nu")){
          // H->ZZ->2l2nu
          if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/VBF_HToZZTo2L2Nu_M125_TuneCP5_13TeV_tunedown_powheg2_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M300_TuneCP5_13TeV_tunedown_powheg2_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/VBF_HToZZTo2L2Nu_M125_TuneCP5_13TeV_tuneup_powheg2_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M300_TuneCP5_13TeV_tuneup_powheg2_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/VBF_HToZZTo2L2Nu_M125_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M160_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M170_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M180_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M190_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M200_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M210_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M230_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M250_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M270_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              // This sample has wrong PU.
              //"/VBF_HToZZTo2L2Nu_M300_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M300_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M350_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M400_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M450_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M500_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M550_TuneCP5_13TeV_powheg2_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M600_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M700_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M800_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M900_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M1000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M1500_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M2000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M2500_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M3000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
            };
          }
        } // End H->ZZ->2l2nu
        else if (sname.Contains("WWTo2L2Nu")){
          // H->WW->2l2nu
          if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/VBFHToWWTo2L2Nu_M125_NNPDF31_TuneCP5Down_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M300_CPS_NNPDF31_TuneCP5Down_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/VBFHToWWTo2L2Nu_M125_NNPDF31_TuneCP5Up_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M300_CPS_NNPDF31_TuneCP5Up_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/VBFHToWWTo2L2Nu_M125_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M160_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M170_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M180_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M190_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M200_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M210_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M230_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M250_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGenV727_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M270_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M300_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M350_13TeV_powheg2_JHUGenV714_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M400_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGenV727_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M450_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M500_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M550_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M600_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M700_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M800_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M900_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M1000_13TeV_powheg2_JHUGenV714_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M1500_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M2000_13TeV_powheg2_JHUGenV714_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M2500_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M3000_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
            };
          }
        } // End H->WW->2l2nu
      } // End VBF
      else if (sname.Contains("ZH")){
        // ZH
        if (sname.Contains("To2Nu2X_2LFilter")){
          // H->ZZ->2nu2x
          if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/ZH_HToZZTo2Nu2X_2LFilter_M125_TuneCP5_13TeV_tunedown_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_M300_TuneCP5_13TeV_tunedown_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/ZH_HToZZTo2Nu2X_2LFilter_M125_TuneCP5_13TeV_tuneup_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_M300_TuneCP5_13TeV_tuneup_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M125_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M160_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M170_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M180_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M190_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M200_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/ZH_HTo2Nu2X_2LFilter_4LVetoFilter_M210_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo2Nu2X_2LFilter_4LVetoFilter_M230_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M250_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/ZH_HTo2Nu2X_2LFilter_4LVetoFilter_M270_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M300_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M350_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M400_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M450_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M500_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M550_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M600_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M700_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M800_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M900_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M1000_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M1500_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M2000_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M2500_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M3000_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
            };
          }
        } // End H->ZZ->2nu2x
        else if (sname.Contains("To2L2Q_2LFilter")){
          // H->ZZ->2l2q
          if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M125_TuneCP5_13TeV_tunedown_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M300_TuneCP5_13TeV_tunedown_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M125_TuneCP5_13TeV_tuneup_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M300_TuneCP5_13TeV_tuneup_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M125_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M160_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M170_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M180_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M190_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M200_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HToZZTo2L2Q_2LFilter_4LVeto_M210_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M230_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M250_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M270_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M300_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M350_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M400_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M450_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M500_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M550_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M600_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M700_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M800_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M900_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M1000_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M1500_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M2000_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M2500_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M3000_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM"
            };
          }
        } // End H->ZZ->2l2q
        else if (sname.Contains("To4Q_2LFilter") && !sname.Contains("WWTo4Q_2LFilter")){
          // H->ZZ->4q
          if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/ZH_HTo4Q_2LFilter_M125_TuneCP5_13TeV_tunedown_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M300_TuneCP5_13TeV_tunedown_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/ZH_HTo4Q_2LFilter_M125_TuneCP5_13TeV_tuneup_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M300_TuneCP5_13TeV_tuneup_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/ZH_HTo4Q_2LFilter_M125_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M160_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M170_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M180_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M190_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M200_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M210_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M230_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M250_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M270_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M350_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M400_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M500_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M550_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M600_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M700_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M800_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M1000_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M2500_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M3000_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM"
            };
          }
        } // End H->ZZ->4q
        else if (sname.Contains("WWTo2L2Nu")){
          // H->WW->2l2nu, better to use the lnuqq sample for the 3l contribution MC
          if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/HZJ_HToWW2L2NuInclusiveZ_M125_NNPDF31_TuneCP5Down_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M300_CPS_NNPDF31_TuneCP5Down_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/HZJ_HToWW2L2NuInclusiveZ_M125_NNPDF31_TuneCP5Up_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M300_CPS_NNPDF31_TuneCP5Up_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/HZJ_HToWW2L2NuInclusiveZ_M125_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M160_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M170_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M180_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M190_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M200_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M210_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M230_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M250_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M270_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M300_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M350_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M400_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M450_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M500_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M550_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M600_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M700_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M800_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M900_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M1000_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M1500_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M2000_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M2500_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M3000_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
            };
          }
        } // End H->WW->2l2nu
        else if (sname.Contains("LNuQQ_2LFilter")){
          // H->WW->lnuqq
          if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/ZH_HToWWToLNuQQ_2LFilter_M125_TuneCP5_13TeV_tunedown_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/ZH_HToWWToLNuQQ_2LFilter_M125_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M160_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M170_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M180_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M190_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M250_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M270_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M1000_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M1500_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M3000_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM"
            };
          }
        } // End H->WW->lnuqq
      } // End ZH
      else if (sname.Contains("WminusH")){
        // WminusH
        if (sname.Contains("ZZTo2L2Nu")){
          // H->ZZ->2l2nu
          if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/WminusH_HToZZTo2L2Nu_M125_TuneCP5_13TeV_tunedown_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M300_TuneCP5_13TeV_tunedown_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/WminusH_HToZZTo2L2Nu_M125_TuneCP5_13TeV_tuneup_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M300_TuneCP5_13TeV_tuneup_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/WminusH_HToZZTo2L2Nu_M125_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M160_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M170_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M180_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M190_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M200_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M210_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M230_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M250_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M270_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M300_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M350_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M400_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M450_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M500_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M550_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M600_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M700_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M800_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M900_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M1000_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M1500_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M2000_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M2500_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M3000_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
            };
          }
        } // End H->ZZ->2l2nu
        else if (sname.Contains("ZZTo2L2Q")){
          // H->ZZ->2l2q
          if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/WminusH_HToZZTo2L2Q_M300_TuneCP5_13TeV_tunedown_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/WminusH_HToZZTo2L2Q_M125_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M160_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M170_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M180_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M190_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M200_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M210_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M230_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M250_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M270_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M300_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M350_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M400_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M450_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M500_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M550_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M600_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M700_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M800_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M900_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M2500_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM"
            };
          }
        } // End H->ZZ->2l2q
        else if (sname.Contains("WW_2LOSFilter")){
          // H->WW->lnu2x with 2l OS filter
          if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/WminusHToWW_2LOSFilter_M125_NNPDF31_TuneCP5Down_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M300_CPS_NNPDF31_TuneCP5Down_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/WminusHToWW_2LOSFilter_M125_NNPDF31_TuneCP5Up_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M300_CPS_NNPDF31_TuneCP5Up_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/WminusHToWW_2LOSFilter_M125_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M160_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M170_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M180_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M190_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M200_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M210_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M230_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M250_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M270_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M300_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M350_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M400_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M450_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M500_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M550_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M600_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M700_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M800_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M900_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M1000_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M1500_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M2000_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M2500_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M3000_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
            };
          }
        } // End H->WW->lnu2x
      } // End WminusH
      else if (sname.Contains("WplusH")){
        // WplusH
        if (sname.Contains("ZZTo2L2Nu")){
          // H->ZZ->2l2nu
          if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/WplusH_HToZZTo2L2Nu_M125_TuneCP5_13TeV_tunedown_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M300_TuneCP5_13TeV_tunedown_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/WplusH_HToZZTo2L2Nu_M125_TuneCP5_13TeV_tuneup_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M300_TuneCP5_13TeV_tuneup_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/WplusH_HToZZTo2L2Nu_M125_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M160_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M170_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M180_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M190_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M200_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M210_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M230_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M250_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M270_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M300_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M350_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M400_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M450_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M500_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M550_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M600_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M700_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M800_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M900_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M1000_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M1500_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M2000_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M2500_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M3000_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
            };
          }
        } // End H->ZZ->2l2nu
        else if (sname.Contains("ZZTo2L2Q")){
          // H->ZZ->2l2q
          if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/WplusH_HToZZTo2L2Q_M300_TuneCP5_13TeV_tunedown_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/WplusH_HToZZTo2L2Q_M125_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M160_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M170_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M180_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M190_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M200_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M210_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M230_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M250_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M270_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M300_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M350_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M400_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M450_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M550_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M600_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M700_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M800_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M900_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M1000_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M1500_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M2000_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_private/MINIAODSIM"
            };
          }
        } // End H->ZZ->2l2q
        else if (sname.Contains("WW_2LOSFilter")){
          // H->WW->lnu2x with 2l OS filter
          if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/WplusHToWW_2LOSFilter_M125_NNPDF31_TuneCP5Down_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M300_CPS_NNPDF31_TuneCP5Down_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/WplusHToWW_2LOSFilter_M125_NNPDF31_TuneCP5Up_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M300_CPS_NNPDF31_TuneCP5Up_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/WplusHToWW_2LOSFilter_M125_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M160_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M170_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M180_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M190_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M200_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M210_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M230_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M250_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M270_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M300_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M350_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M400_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M450_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M500_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M550_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M600_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M700_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M800_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M900_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M1000_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M1500_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M2000_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M2500_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M3000_CPS_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen735_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
            };
          }
        } // End H->WW->lnu2x
      } // End WplusH

      float const mh_req = SampleHelpers::findPoleMass(sname);
      for (unsigned int isample=0; isample<samplelist.size(); isample++){
        auto const& ss = samplelist.at(isample);
        float const mh_sample = SampleHelpers::findPoleMass(ss);
        if (isample>0){
          float const mh_sample_prev = SampleHelpers::findPoleMass(samplelist.at(isample-1));
          if (mh_sample<=mh_sample_prev){
            MELAerr << "SampleHelpers::constructSamplesList: Mass = " << mh_sample << " <= previous mass = " << mh_sample_prev << " during the request for " << sname << endl;
            assert(0);
          }
        }
        if (mh_req<0.f || mh_sample==mh_req) HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ ss });
      }
    }

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
    if (sname == "ZGJets_nunu_nlo_inclusive" || sname == "ZGJets_nunu_nlo") HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ "/ZGTo2NuG_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" });
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
      sname == "EGamma" || sname == "JetMET" || sname == "SingleLepton" // Composite names
      ){
      for (auto const& dp:SampleHelpers::getValidDataPeriods()) constructSamplesList(sname+"_"+dp, syst, samples);
    }
    for (auto const& dp:SampleHelpers::getValidDataPeriods()){
      if (sname == Form("EGamma_%s", dp.Data())){
        constructSamplesList(Form("DoubleEG_%s", dp.Data()), syst, samples);
        constructSamplesList(Form("SingleElectron_%s", dp.Data()), syst, samples);
        constructSamplesList(Form("SinglePhoton_%s", dp.Data()), syst, samples);
      }
      else if (sname == Form("JetMET_%s", dp.Data())){
        constructSamplesList(Form("JetHT_%s", dp.Data()), syst, samples);
        constructSamplesList(Form("HTMHT_%s", dp.Data()), syst, samples);
        constructSamplesList(Form("MET_%s", dp.Data()), syst, samples);
      }
      else if (sname == Form("SingleLepton_%s", dp.Data())){
        constructSamplesList(Form("SingleMuon_%s", dp.Data()), syst, samples);
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
    // POWHEG signal samples
    if (sname.Contains("POWHEG")){
      std::vector<TString> samplelist;

      if (sname.Contains("GGH")){
        // gg->H
        if (sname.Contains("ZZTo2L2Nu")){
          // H->ZZ->2l2nu
          if (syst==tHardJetsDn || syst==tHardJetsUp){
            samplelist = std::vector<TString>{
              "/GluGluHToZZTo2L2Nu_M125_TuneCUETP8M1_13TeV_powheg2_minloHJJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M300_TuneCUETP8M1_13TeV_powheg2_minloHJJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaScaleDn){
            samplelist = std::vector<TString>{
              "/GluGluHToZZTo2L2Nu_M125_TuneCUETP8M1_13TeV_scaledown_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M300_TuneCUETP8M1_13TeV_scaledown_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaScaleUp){
            samplelist = std::vector<TString>{
              "/GluGluHToZZTo2L2Nu_M125_TuneCUETP8M1_13TeV_scaleup_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M300_TuneCUETP8M1_13TeV_scaleup_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/GluGluHToZZTo2L2Nu_M125_TuneCUETP8M1_13TeV_tunedown_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M300_TuneCUETP8M1_13TeV_tunedown_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/GluGluHToZZTo2L2Nu_M125_TuneCUETP8M1_13TeV_tuneup_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M300_TuneCUETP8M1_13TeV_tuneup_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/GluGluHToZZTo2L2Nu_M125_TuneCUETP8M1_13TeV_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M160_TuneCUETP8M1_13TeV_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M170_TuneCUETP8M1_13TeV_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M180_TuneCUETP8M1_13TeV_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M190_TuneCUETP8M1_13TeV_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M200_TuneCUETP8M1_13TeV_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M210_TuneCUETP8M1_13TeV_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M230_TuneCUETP8M1_13TeV_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M250_TuneCUETP8M1_13TeV_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M270_TuneCUETP8M1_13TeV_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M300_13TeV_powheg2_JHUgenV698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M350_TuneCUETP8M1_13TeV_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M400_13TeV_powheg2_JHUgenV698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M450_TuneCUETP8M1_13TeV_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M500_13TeV_powheg2_JHUgenV698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M550_TuneCUETP8M1_13TeV_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M600_13TeV_powheg2_JHUgenV698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M700_13TeV_powheg2_JHUgenV698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M800_13TeV_powheg2_JHUgenV698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M900_13TeV_powheg2_JHUgenV698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M1000_13TeV_powheg2_JHUgenV698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M1500_13TeV_powheg2_JHUgenV698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M2000_13TeV_powheg2_JHUgenV698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M2500_13TeV_powheg2_JHUgenV698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/GluGluHToZZTo2L2Nu_M3000_13TeV_powheg2_JHUgenV698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM"
            };
          }
        } // End H->ZZ->2l2nu
        else if (sname.Contains("WWTo2L2Nu")){
          // H->WW->2l2nu
          if (syst==tHardJetsDn || syst==tHardJetsUp){
            samplelist = std::vector<TString>{
              "/GluGluHToWWTo2L2Nu_M125_TuneCUETP8M1_13TeV_powheg_minloHJJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M300_CPS_TuneCUETP8M1_13TeV_powheg_minlo_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaScaleDn || syst==tPythiaScaleUp){
            samplelist = std::vector<TString>{
              "/GluGluHToWWTo2L2Nu_M125_TuneCUETP8M1_PSweights_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M300_CPS_TuneCUETP8M1_PSweights_13TeV_powheg_JHUGenV698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/GluGluHToWWTo2L2Nu_M125_TuneCUETP8M1Down_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M300_CPS_TuneCUETP8M1Down_13TeV_powheg_JHUGenV698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/GluGluHToWWTo2L2Nu_M125_TuneCUETP8M1Up_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M300_CPS_TuneCUETP8M1Up_13TeV_powheg_JHUGenV698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/GluGluHToWWTo2L2Nu_M125_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M160_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M170_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M180_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M190_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M200_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M210_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M230_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M250_CPS_TuneCUETP8M1_PSweights_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M270_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              //"/GluGluHToWWTo2L2Nu_M300_13TeV_powheg_JHUgenv698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M300_CPS_TuneCUETP8M1_PSweights_13TeV_powheg_JHUGenV698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M350_13TeV_powheg_JHUgenv698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M400_13TeV_powheg_JHUgenv698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              //"/GluGluHToWWTo2L2Nu_M400_13TeV_powheg_JHUgenv727_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M450_13TeV_powheg_JHUgenv698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M500_13TeV_powheg_JHUgenv698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M550_13TeV_powheg_JHUgenv698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M600_13TeV_powheg_JHUgenv698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M700_13TeV_powheg_JHUgenv698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M800_13TeV_powheg_JHUgenv698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M900_13TeV_powheg_JHUgenv698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M1000_13TeV_powheg_JHUgenv698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M1500_13TeV_powheg_JHUgenv698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M2000_13TeV_powheg_JHUgenv698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M2500_13TeV_powheg_JHUgenv698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/GluGluHToWWTo2L2Nu_M3000_13TeV_powheg_JHUgenv698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM"
            };
          }
        } // End H->WW->2l2nu
      } // End gg->H
      else if (sname.Contains("VBF")){
        // VBF
        if (sname.Contains("ZZTo2L2Nu")){
          // H->ZZ->2l2nu
          if (syst==tPythiaScaleDn){
            samplelist = std::vector<TString>{
              "/VBF_HToZZTo2L2Nu_M125_TuneCUETP8M1_13TeV_scaledown_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M300_TuneCUETP8M1_13TeV_scaledown_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaScaleUp){
            samplelist = std::vector<TString>{
              "/VBF_HToZZTo2L2Nu_M125_TuneCUETP8M1_13TeV_scaleup_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M300_TuneCUETP8M1_13TeV_scaleup_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/VBF_HToZZTo2L2Nu_M125_TuneCUETP8M1_13TeV_tunedown_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M300_TuneCUETP8M1_13TeV_tunedown_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/VBF_HToZZTo2L2Nu_M125_TuneCUETP8M1_13TeV_tuneup_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M300_TuneCUETP8M1_13TeV_tuneup_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/VBF_HToZZTo2L2Nu_M125_TuneCUETP8M1_13TeV_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M160_TuneCUETP8M1_13TeV_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M170_TuneCUETP8M1_13TeV_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M180_TuneCUETP8M1_13TeV_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M190_TuneCUETP8M1_13TeV_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M200_TuneCUETP8M1_13TeV_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M210_TuneCUETP8M1_13TeV_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M230_TuneCUETP8M1_13TeV_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M250_TuneCUETP8M1_13TeV_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M300_13TeV_powheg2_JHUgenV698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M350_TuneCUETP8M1_13TeV_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M400_13TeV_powheg2_JHUgenV698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M450_TuneCUETP8M1_13TeV_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M500_13TeV_powheg2_JHUgenV698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M550_TuneCUETP8M1_13TeV_powheg2_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M600_13TeV_powheg2_JHUgenV698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M700_13TeV_powheg2_JHUgenV698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M800_13TeV_powheg2_JHUgenV698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M900_13TeV_powheg2_JHUgenV698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M1000_13TeV_powheg2_JHUgenV698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M1500_13TeV_powheg2_JHUgenV698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M2000_13TeV_powheg2_JHUgenV698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M2500_13TeV_powheg2_JHUgenV698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/VBF_HToZZTo2L2Nu_M3000_13TeV_powheg2_JHUgenV698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM"
            };
          }
        } // End H->ZZ->2l2nu
        else if (sname.Contains("WWTo2L2Nu")){
          // H->WW->2l2nu
          if (syst==tPythiaScaleDn || syst==tPythiaScaleUp){
            samplelist = std::vector<TString>{
              "/VBFHToWWTo2L2Nu_M125_TuneCUETP8M1_PSweights_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M300_CPS_TuneCUETP8M1_PSweights_13TeV_powheg_JHUGenV698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/VBFHToWWTo2L2Nu_M125_TuneCUETP8M1Down_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M300_CPS_TuneCUETP8M1Down_13TeV_powheg_JHUGenV698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/VBFHToWWTo2L2Nu_M125_TuneCUETP8M1Up_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M300_CPS_TuneCUETP8M1Up_13TeV_powheg_JHUGenV698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              //"/VBFHToWWTo2L2Nu_M125_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M125_TuneCUETP8M1_PSweights_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M160_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M170_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M180_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M190_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M200_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M210_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M230_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M250_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M270_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              //"/VBFHToWWTo2L2Nu_M300_13TeV_powheg_JHUgenv698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M300_CPS_TuneCUETP8M1_PSweights_13TeV_powheg_JHUGenV698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M350_13TeV_powheg_JHUgenv698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M400_13TeV_powheg_JHUgenv698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M450_13TeV_powheg_JHUgenv698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M500_13TeV_powheg_JHUgenv698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M550_13TeV_powheg_JHUgenv698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M600_13TeV_powheg_JHUgenv698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M700_13TeV_powheg_JHUgenv698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M800_13TeV_powheg_JHUgenv698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M900_13TeV_powheg_JHUgenv698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M1000_13TeV_powheg_JHUgenv698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M1500_13TeV_powheg_JHUgenv698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M2000_13TeV_powheg_JHUgenv698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M2500_13TeV_powheg_JHUgenv698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
              "/VBFHToWWTo2L2Nu_M3000_13TeV_powheg_JHUgenv698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM"
            };
          }
        } // End H->WW->2l2nu
      } // End VBF
      else if (sname.Contains("ZH")){
        // ZH
        if (sname.Contains("To2Nu2X_2LFilter")){
          // H->ZZ->2nu2x
          if (syst==tPythiaScaleDn){
            samplelist = std::vector<TString>{
              "/ZH_HToZZTo2Nu2X_2LFilter_M125_TuneCUETP8M1_13TeV_scaledown_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_M300_TuneCUETP8M1_13TeV_scaledown_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaScaleUp){
            samplelist = std::vector<TString>{
              "/ZH_HToZZTo2Nu2X_2LFilter_M125_TuneCUETP8M1_13TeV_scaleup_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_M300_TuneCUETP8M1_13TeV_scaleup_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/ZH_HToZZTo2Nu2X_2LFilter_M125_TuneCUETP8M1_13TeV_tunedown_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_M300_TuneCUETP8M1_13TeV_tunedown_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/ZH_HToZZTo2Nu2X_2LFilter_M125_TuneCUETP8M1_13TeV_tuneup_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_M300_TuneCUETP8M1_13TeV_tuneup_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M125_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M160_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M170_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M180_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M190_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M200_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M210_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M230_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M250_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HTo2Nu2X_2LFilter_4LVetoFilter_M270_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M300_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M350_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M400_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M450_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M500_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M550_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M600_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M700_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M800_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M900_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M1000_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M1500_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M2000_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M2500_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2Nu2X_2LFilter_4LVeto_M3000_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
        } // End H->ZZ->2nu2x
        else if (sname.Contains("To2L2Q_2LFilter")){
          // H->ZZ->2l2q
          if (syst==tPythiaScaleDn){
            samplelist = std::vector<TString>{
              "/ZH_HToZZTo2L2Q_2LFilter_M125_TuneCUETP8M1_13TeV_scaledown_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2L2Q_2LFilter_M300_TuneCUETP8M1_13TeV_scaledown_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaScaleUp){
            samplelist = std::vector<TString>{
              "/ZH_HToZZTo2L2Q_2LFilter_M125_TuneCUETP8M1_13TeV_scaleup_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2L2Q_2LFilter_M300_TuneCUETP8M1_13TeV_scaleup_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/ZH_HToZZTo2L2Q_2LFilter_M125_TuneCUETP8M1_13TeV_tunedown_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2L2Q_2LFilter_M300_TuneCUETP8M1_13TeV_tunedown_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/ZH_HToZZTo2L2Q_2LFilter_M125_TuneCUETP8M1_13TeV_tuneup_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2L2Q_2LFilter_M300_TuneCUETP8M1_13TeV_tuneup_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/ZH_HToZZTo2L2Q_2LFilter_4LVeto_M125_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2L2Q_2LFilter_4LVeto_M160_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2L2Q_2LFilter_4LVeto_M170_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2L2Q_2LFilter_4LVeto_M180_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2L2Q_2LFilter_4LVeto_M190_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2L2Q_2LFilter_4LVeto_M200_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2L2Q_2LFilter_4LVeto_M210_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2L2Q_2LFilter_4LVeto_M230_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2L2Q_2LFilter_4LVeto_M250_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HTo2L2Q_2LFilter_4LVetoFilter_M270_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/ZH_HToZZTo2L2Q_2LFilter_4LVeto_M300_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2L2Q_2LFilter_4LVeto_M350_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2L2Q_2LFilter_4LVeto_M400_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2L2Q_2LFilter_4LVeto_M450_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2L2Q_2LFilter_4LVeto_M500_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2L2Q_2LFilter_4LVeto_M550_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2L2Q_2LFilter_4LVeto_M600_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2L2Q_2LFilter_4LVeto_M700_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2L2Q_2LFilter_4LVeto_M800_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2L2Q_2LFilter_4LVeto_M900_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2L2Q_2LFilter_4LVeto_M1000_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2L2Q_2LFilter_4LVeto_M1500_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2L2Q_2LFilter_4LVeto_M2000_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2L2Q_2LFilter_4LVeto_M2500_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo2L2Q_2LFilter_4LVeto_M3000_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
        } // End H->ZZ->2l2q
        else if (sname.Contains("To4Q_2LFilter") && !sname.Contains("WWTo4Q_2LFilter")){
          // H->ZZ->4q
          if (syst==tPythiaScaleDn){
            samplelist = std::vector<TString>{
              "/ZH_HToZZTo4Q_2LFilter_M125_TuneCUETP8M1_13TeV_scaledown_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo4Q_2LFilter_M300_TuneCUETP8M1_13TeV_scaledown_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaScaleUp){
            samplelist = std::vector<TString>{
              "/ZH_HToZZTo4Q_2LFilter_M125_TuneCUETP8M1_13TeV_scaleup_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo4Q_2LFilter_M300_TuneCUETP8M1_13TeV_scaleup_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/ZH_HToZZTo4Q_2LFilter_M125_TuneCUETP8M1_13TeV_tunedown_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo4Q_2LFilter_M300_TuneCUETP8M1_13TeV_tunedown_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/ZH_HToZZTo4Q_2LFilter_M125_TuneCUETP8M1_13TeV_tuneup_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo4Q_2LFilter_M300_TuneCUETP8M1_13TeV_tuneup_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/ZH_HToZZTo4Q_2LFilter_M125_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo4Q_2LFilter_M160_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo4Q_2LFilter_M170_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo4Q_2LFilter_M180_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo4Q_2LFilter_M190_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo4Q_2LFilter_M200_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M210_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/ZH_HToZZTo4Q_2LFilter_M230_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo4Q_2LFilter_M250_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HTo4Q_2LFilter_M270_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/ZH_HToZZTo4Q_2LFilter_M300_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo4Q_2LFilter_M350_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo4Q_2LFilter_M400_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo4Q_2LFilter_M450_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo4Q_2LFilter_M500_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo4Q_2LFilter_M550_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo4Q_2LFilter_M600_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo4Q_2LFilter_M700_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo4Q_2LFilter_M800_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo4Q_2LFilter_M900_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo4Q_2LFilter_M1000_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo4Q_2LFilter_M1500_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo4Q_2LFilter_M2000_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo4Q_2LFilter_M2500_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/ZH_HToZZTo4Q_2LFilter_M3000_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
        } // End H->ZZ->4q
        else if (sname.Contains("WWTo2L2Nu")){
          // H->WW->2l2nu, better to use the lnuqq sample for the 3l contribution MC
          if (syst==tPythiaScaleDn || syst==tPythiaScaleUp){
            samplelist = std::vector<TString>{
              "/HZJ_HToWW2L2NuInclusiveZ_M125_TuneCUETP8M1_PSweights_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M300_CPS_TuneCUETP8M1_PSweights_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/HZJ_HToWW2L2NuInclusiveZ_M125_TuneCUETP8M1Down_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M300_CPS_TuneCUETP8M1Down_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/HZJ_HToWW2L2NuInclusiveZ_M125_TuneCUETP8M1Up_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M300_CPS_TuneCUETP8M1Up_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              //"/HZJ_HToWW2L2NuInclusiveZ_M125_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M125_TuneCUETP8M1_PSweights_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M160_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M170_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M180_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M190_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M200_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M210_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M230_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M250_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M270_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              //"/HZJ_HToWW2L2NuInclusiveZ_M300_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M300_CPS_TuneCUETP8M1_PSweights_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M350_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M400_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M450_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M500_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M550_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M600_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M700_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M800_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M900_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M1000_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M1500_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M2000_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M2500_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/HZJ_HToWW2L2NuInclusiveZ_M3000_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
        } // End H->WW->2l2nu
        else if (sname.Contains("LNuQQ_2LFilter")){
          // H->WW->lnuqq
          if (syst==tPythiaScaleDn){
            samplelist = std::vector<TString>{
              "/ZH_HToWWToLNuQQ_2LFilter_M125_TuneCUETP8M1_13TeV_scaledown_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M300_TuneCUETP8M1_13TeV_scaledown_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM"
            };
          }
          else if (syst==tPythiaScaleUp){
            samplelist = std::vector<TString>{
              "/ZH_HToWWToLNuQQ_2LFilter_M125_TuneCUETP8M1_13TeV_scaleup_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M300_TuneCUETP8M1_13TeV_scaleup_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/ZH_HToWWToLNuQQ_2LFilter_M125_TuneCUETP8M1_13TeV_tunedown_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M300_TuneCUETP8M1_13TeV_tunedown_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/ZH_HToWWToLNuQQ_2LFilter_M125_TuneCUETP8M1_13TeV_tuneup_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M300_TuneCUETP8M1_13TeV_tuneup_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/ZH_HToWWToLNuQQ_2LFilter_M125_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M160_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M170_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M180_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M190_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M200_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M210_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M230_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M250_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M270_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M300_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M350_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M400_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M450_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M500_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M550_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M600_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M700_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M800_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M900_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M1000_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M1500_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M2000_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M2500_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/ZH_HToWWToLNuQQ_2LFilter_M3000_TuneCUETP8M1_13TeV_powheg2-minlo-HZJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM"
            };
          }
        } // End H->WW->lnuqq
      } // End ZH
      else if (sname.Contains("WminusH")){
        // WminusH
        if (sname.Contains("ZZTo2L2Nu")){
          // H->ZZ->2l2nu
          if (syst==tPythiaScaleDn){
            samplelist = std::vector<TString>{
              "/WminusH_HToZZTo2L2Nu_M125_TuneCUETP8M1_13TeV_scaledown_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M300_TuneCUETP8M1_13TeV_scaledown_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaScaleUp){
            samplelist = std::vector<TString>{
              "/WminusH_HToZZTo2L2Nu_M125_TuneCUETP8M1_13TeV_scaleup_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M300_TuneCUETP8M1_13TeV_scaleup_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/WminusH_HToZZTo2L2Nu_M125_TuneCUETP8M1_13TeV_tunedown_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M300_TuneCUETP8M1_13TeV_tunedown_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/WminusH_HToZZTo2L2Nu_M125_TuneCUETP8M1_13TeV_tuneup_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M300_TuneCUETP8M1_13TeV_tuneup_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/WminusH_HToZZTo2L2Nu_M125_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M160_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M170_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M180_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M190_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M200_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M210_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M230_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M250_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M270_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M300_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M350_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M400_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M450_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M500_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M550_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M600_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M700_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M800_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M900_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M1000_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M1500_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M2000_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M2500_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusH_HToZZTo2L2Nu_M3000_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
        } // End H->ZZ->2l2nu
        else if (sname.Contains("ZZTo2L2Q")){
          // H->ZZ->2l2q
          if (syst==tPythiaScaleDn){
            samplelist = std::vector<TString>{
              "/WminusH_HToZZTo2L2Q_M125_TuneCUETP8M1_13TeV_scaledown_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M300_TuneCUETP8M1_13TeV_scaledown_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM"
            };
          }
          else if (syst==tPythiaScaleUp){
            samplelist = std::vector<TString>{
              "/WminusH_HToZZTo2L2Q_M125_TuneCUETP8M1_13TeV_scaleup_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M300_TuneCUETP8M1_13TeV_scaleup_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/WminusH_HToZZTo2L2Q_M125_TuneCUETP8M1_13TeV_tunedown_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M300_TuneCUETP8M1_13TeV_tunedown_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/WminusH_HToZZTo2L2Q_M125_TuneCUETP8M1_13TeV_tuneup_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M300_TuneCUETP8M1_13TeV_tuneup_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/WminusH_HToZZTo2L2Q_M125_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M160_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M170_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M180_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M190_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M200_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M210_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M230_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M250_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M270_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M300_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M350_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M400_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M450_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M500_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M550_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M600_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M700_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M800_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M900_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M1000_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M1500_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M2000_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WminusH_HToZZTo2L2Q_M3000_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM"
            };
          }
        } // End H->ZZ->2l2q
        else if (sname.Contains("WW_2LOSFilter")){
          // H->WW->lnu2x with 2l OS filter
          if (syst==tPythiaScaleDn || syst==tPythiaScaleUp){
            samplelist = std::vector<TString>{
              "/WminusHToWW_2LOSFilter_M125_TuneCUETP8M1_PSweights_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M300_CPS_TuneCUETP8M1_PSweights_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/WminusHToWW_2LOSFilter_M125_TuneCUETP8M1Down_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M300_CPS_TuneCUETP8M1Down_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/WminusHToWW_2LOSFilter_M125_TuneCUETP8M1Up_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M300_CPS_TuneCUETP8M1Up_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              //"/WminusHToWW_2LOSFilter_M125_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M125_TuneCUETP8M1_PSweights_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M160_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M170_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M180_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M190_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M200_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M210_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M230_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M250_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M270_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              //"/WminusHToWW_2LOSFilter_M300_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M300_CPS_TuneCUETP8M1_PSweights_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M350_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M400_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M450_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M500_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M550_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M600_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M700_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M800_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M900_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M1000_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M1500_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M2000_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M2500_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WminusHToWW_2LOSFilter_M3000_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
        } // End H->WW->lnu2x
      } // End WminusH
      else if (sname.Contains("WplusH")){
        // WplusH
        if (sname.Contains("ZZTo2L2Nu")){
          // H->ZZ->2l2nu
          if (syst==tPythiaScaleDn){
            samplelist = std::vector<TString>{
              "/WplusH_HToZZTo2L2Nu_M125_TuneCUETP8M1_13TeV_scaledown_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M300_TuneCUETP8M1_13TeV_scaledown_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaScaleUp){
            samplelist = std::vector<TString>{
              "/WplusH_HToZZTo2L2Nu_M125_TuneCUETP8M1_13TeV_scaleup_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M300_TuneCUETP8M1_13TeV_scaleup_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/WplusH_HToZZTo2L2Nu_M125_TuneCUETP8M1_13TeV_tunedown_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M300_TuneCUETP8M1_13TeV_tunedown_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/WplusH_HToZZTo2L2Nu_M125_TuneCUETP8M1_13TeV_tuneup_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M300_TuneCUETP8M1_13TeV_tuneup_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/WplusH_HToZZTo2L2Nu_M125_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M160_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M170_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M180_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M190_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M200_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M210_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M230_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M250_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M270_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M300_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M350_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M400_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M450_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M500_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M550_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M600_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M700_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M800_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M900_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M1000_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M1500_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M2000_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M2500_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusH_HToZZTo2L2Nu_M3000_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
        } // End H->ZZ->2l2nu
        else if (sname.Contains("ZZTo2L2Q")){
          // H->ZZ->2l2q
          if (syst==tPythiaScaleDn){
            samplelist = std::vector<TString>{
              "/WplusH_HToZZTo2L2Q_M125_TuneCUETP8M1_13TeV_scaledown_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M300_TuneCUETP8M1_13TeV_scaledown_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM"
            };
          }
          else if (syst==tPythiaScaleUp){
            samplelist = std::vector<TString>{
              "/WplusH_HToZZTo2L2Q_M125_TuneCUETP8M1_13TeV_scaleup_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M300_TuneCUETP8M1_13TeV_scaleup_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/WplusH_HToZZTo2L2Q_M125_TuneCUETP8M1_13TeV_tunedown_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M300_TuneCUETP8M1_13TeV_tunedown_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/WplusH_HToZZTo2L2Q_M125_TuneCUETP8M1_13TeV_tuneup_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M300_TuneCUETP8M1_13TeV_tuneup_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              "/WplusH_HToZZTo2L2Q_M125_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M160_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M170_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M180_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M190_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M200_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M210_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M230_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M250_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M270_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M300_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M350_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M400_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M450_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M500_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M550_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M600_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M700_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M800_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M900_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M1000_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M1500_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M2000_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM",
              "/WplusH_HToZZTo2L2Q_M2500_TuneCUETP8M1_13TeV_powheg2-minlo-HWJ_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_private/MINIAODSIM"
            };
          }
        } // End H->ZZ->2l2q
        else if (sname.Contains("WW_2LOSFilter")){
          // H->WW->lnu2x with 2l OS filter
          if (syst==tPythiaScaleDn || syst==tPythiaScaleUp){
            samplelist = std::vector<TString>{
              "/WplusHToWW_2LOSFilter_M125_TuneCUETP8M1_PSweights_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M300_CPS_TuneCUETP8M1_PSweights_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneDn){
            samplelist = std::vector<TString>{
              "/WplusHToWW_2LOSFilter_M125_TuneCUETP8M1Down_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M300_CPS_TuneCUETP8M1Down_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else if (syst==tPythiaTuneUp){
            samplelist = std::vector<TString>{
              "/WplusHToWW_2LOSFilter_M125_TuneCUETP8M1Up_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M300_CPS_TuneCUETP8M1Up_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
          else{
            samplelist = std::vector<TString>{
              //"/WplusHToWW_2LOSFilter_M125_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M125_TuneCUETP8M1_PSweights_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M160_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M170_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M180_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M190_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M200_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M210_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M230_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M250_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M270_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              //"/WplusHToWW_2LOSFilter_M300_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M300_CPS_TuneCUETP8M1_PSweights_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M350_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M400_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M450_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M500_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M550_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M600_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M700_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M800_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M900_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M1000_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M1500_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M2000_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M2500_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
              "/WplusHToWW_2LOSFilter_M3000_CPS_TuneCUETP8M1_13TeV_powheg_JHUGenV735_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM"
            };
          }
        } // End H->WW->lnu2x
      } // End WplusH

      float const mh_req = SampleHelpers::findPoleMass(sname);
      for (unsigned int isample=0; isample<samplelist.size(); isample++){
        auto const& ss = samplelist.at(isample);
        float const mh_sample = SampleHelpers::findPoleMass(ss);
        if (isample>0){
          float const mh_sample_prev = SampleHelpers::findPoleMass(samplelist.at(isample-1));
          if (mh_sample<=mh_sample_prev){
            MELAerr << "SampleHelpers::constructSamplesList: Mass = " << mh_sample << " <= previous mass = " << mh_sample_prev << " during the request for " << sname << endl;
            assert(0);
          }
        }
        if (mh_req<0.f || mh_sample==mh_req) HelperFunctions::appendVector<TString>(samples, std::vector<TString>{ ss });
      }
    }

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

  // Check if any string is accidentally merged
  {
    bool hasInvalidString = false;
    for (auto const& ss:samples){
      if (ss.Contains("MINIAODSIM/") || ss.Contains("MINIAOD/")){
        MELAerr << "SampleHelpers::constructSamplesList: Sample name " << ss << " is invalid." << endl;
        hasInvalidString = true;
      }
    }
    if (hasInvalidString) exit(1);
  }
}
void SampleHelpers::getSamplesList(std::vector<TString> const& s, std::vector<TString>& vs, SystematicsHelpers::SystematicVariationTypes syst, std::vector<size_t>* ns){
  for (auto const& ss:s){
    std::vector<TString> dumappend; constructSamplesList(ss, syst, dumappend);
    HelperFunctions::appendVector<TString>(vs, dumappend);
    if (ns) ns->push_back(dumappend.size());
  }
}

SampleHelpers::HiggsSampleDecayMode SampleHelpers::getHiggsSampleDecayMode(TString const& sname){
  if (sname.Contains("ZZTo4L")) return kZZTo4L;
  else if (sname.Contains("ZZ_4LFilter")) return kZZTo2L2X;

  else if (sname.Contains("ZZTo2L2Nu")) return kZZTo2L2Nu;
  else if (sname.Contains("ZZTo2Nu2X") || sname.Contains("HTo2Nu2X")) return kZZTo2Nu2X;
  else if (sname.Contains("ZZTo2L2Q") || sname.Contains("HTo2L2Q")) return kZZTo2L2Q;
  else if (sname.Contains("ZZTo4Q") || sname.Contains("HTo4Q")) return kZZTo4Q;

  else if (sname.Contains("WWTo2L2N") || sname.Contains("WW2L2N")) return kWWTo2L2Nu;
  else if (sname.Contains("WWToLNuQQ")) return kWWToLNuQQ;
  else if (sname.Contains("WW_2LOSFilter")) return kWWToLNuXX;

  else{
    MELAerr << "SampleHelpers::getHiggsSampleDecayMode: Cannot identify decay mode for " << sname << "." << endl;
    assert(0);
    return nHiggsSampleDecayModes;
  }
}

double SampleHelpers::calculateAdjustedHiggsBREff(TString const& sname, double const& sum_wgts_defaultMemberZero, double const& sum_wgts_defaultLHEEventWeight, bool hasTaus){
  // POWHEG parameters
  constexpr double BR_Z_ll_POWHEG = 0.1004;
  constexpr double BR_W_lnu_POWHEG = 0.3243;

  constexpr double xw = 0.23119;
  constexpr double T3lL = -0.5;
  constexpr double T3lR =  0;
  constexpr double T3nL =  0.5;
  constexpr double T3nR =  0;
  constexpr double T3uL = 0.5;
  constexpr double T3uR = 0;
  constexpr double T3dL = -0.5;
  constexpr double T3dR = 0;
  constexpr double QlL = -1;
  constexpr double QlR = -1;
  constexpr double QnL = 0;
  constexpr double QnR = 0;
  constexpr double QuL = 2./3.;
  constexpr double QuR = 2./3.;
  constexpr double QdL = -1./3.;
  constexpr double QdR = -1./3.;

  // NLO K factors
  constexpr double scale_alpha_Z_qq = 1.03756;
  constexpr double scale_alpha_W_qq = 1.0382;

  constexpr double aR_lep = 2.*(T3lR-QlR*xw);
  constexpr double aL_lep = 2.*(T3lL-QlL*xw);
  constexpr double aR_neu = 2.*(T3nR-QnR*xw);
  constexpr double aL_neu = 2.*(T3nL-QnL*xw);
  const double aR_QUp = 2.*(T3uR-QuR*xw)*std::sqrt(scale_alpha_Z_qq);
  const double aL_QUp = 2.*(T3uL-QuL*xw)*std::sqrt(scale_alpha_Z_qq);
  const double aR_QDn = 2.*(T3dR-QdR*xw)*std::sqrt(scale_alpha_Z_qq);
  const double aL_QDn = 2.*(T3dL-QdL*xw)*std::sqrt(scale_alpha_Z_qq);

  // No color of flavor factors
  const double BR_Z_ll_single = (std::pow(aL_lep, 2)+std::pow(aR_lep, 2));
  const double BR_Z_nunu_single = (std::pow(aL_neu, 2)+std::pow(aR_neu, 2));
  const double BR_Z_uu_single = (std::pow(aL_QUp, 2)+std::pow(aR_QUp, 2));
  const double BR_Z_dd_single = (std::pow(aL_QDn, 2)+std::pow(aR_QDn, 2));

  constexpr double BR_Z_ll_ratio = 1;
  const double BR_Z_nunu_ratio = BR_Z_nunu_single / BR_Z_ll_single;
  const double BR_Z_uu_ratio = BR_Z_uu_single / BR_Z_ll_single;
  const double BR_Z_dd_ratio = BR_Z_dd_single / BR_Z_ll_single;

  constexpr double BR_W_lnu_single = 1;
  constexpr double BR_W_qq_single = scale_alpha_W_qq;
  constexpr double Nflavs_CKM = 2; // 2 flavors from V_CKM without tops

  constexpr double BR_W_lnu_ratio = BR_W_lnu_single/(BR_W_lnu_single*3. + BR_W_qq_single*Nflavs_CKM*3.);
  constexpr double BR_W_qq_ratio = BR_W_qq_single/(BR_W_lnu_single*3. + BR_W_qq_single*Nflavs_CKM*3.);

  double const adj_br = sum_wgts_defaultLHEEventWeight/sum_wgts_defaultMemberZero;

  TString sname_lower = sname; HelperFunctions::lowercase(sname, sname_lower);
  HiggsSampleDecayMode dkmode = getHiggsSampleDecayMode(sname);
  double sampleMH = SampleHelpers::findPoleMass(sname);
  if (!(sampleMH>=0. && sname_lower.Contains("powheg") && sname_lower.Contains("jhugen"))) return 1;

  bool has2LFilter = sname.Contains("2LFilter");
  bool has2LOSFilter = sname.Contains("2LOSFilter");
  bool has4LFilter = sname.Contains("ZZ_4LFilter");
  if (has2LFilter) MELAout << "SampleHelpers::calculateAdjustedHiggsBREff: A 2L filter is detected in " << sname << "." << endl;
  if (has2LOSFilter) MELAout << "SampleHelpers::calculateAdjustedHiggsBREff: A 2L OS filter is detected in " << sname << "." << endl;
  if (has4LFilter) MELAout << "SampleHelpers::calculateAdjustedHiggsBREff: A 4L filter is detected in " << sname << "." << endl;

  std::vector<TString> hypos;
  std::vector<double> BRcorrs;
  std::vector<double> filter_corrs;
  switch (dkmode){
  case kZZTo4L:
    // 4l, different flavors
    hypos.push_back("2e2mu");
    BRcorrs.push_back(std::pow(BR_Z_ll_ratio, 2)*(hasTaus ? 3. : 1.));
    filter_corrs.push_back(1);
    // 4l, same flavors
    hypos.push_back("4e");
    BRcorrs.push_back(std::pow(BR_Z_ll_ratio, 2)*(hasTaus ? 3. : 2.));
    filter_corrs.push_back(1);
    break;
  case kZZTo2L2X:
    // 4l, different flavors
    hypos.push_back("2e2mu");
    BRcorrs.push_back(std::pow(BR_Z_ll_ratio, 2)*(hasTaus ? 3. : 1.));
    filter_corrs.push_back(1);
    // 4l, same flavors
    hypos.push_back("4e");
    BRcorrs.push_back(std::pow(BR_Z_ll_ratio, 2)*(hasTaus ? 3. : 2.));
    filter_corrs.push_back(1);
    // 2l2nu
    hypos.push_back("2e2mu");
    BRcorrs.push_back(BR_Z_nunu_ratio*3. * BR_Z_ll_ratio*(hasTaus ? 3. : 2.));
    filter_corrs.push_back((!has4LFilter ? 1. : BR_Z_ll_POWHEG));
    // 2l2q
    hypos.push_back("2e2mu");
    BRcorrs.push_back(
      (BR_Z_uu_ratio*2.+ BR_Z_dd_ratio*3.)*3. // Reweighting of 2e->2q
      *
      BR_Z_ll_ratio*(hasTaus ? 3. : 2.) // Reweighting of 2mu->2l
    );
    filter_corrs.push_back((!has4LFilter ? 1. : BR_Z_ll_POWHEG));
    break;

  case kZZTo2L2Nu:
    hypos.push_back("2e2mu");
    BRcorrs.push_back(BR_Z_nunu_ratio*3. * BR_Z_ll_ratio*(hasTaus ? 3. : 2.));
    filter_corrs.push_back(1);
    break;
  case kZZTo2Nu2X:
    // 2nu2l
    hypos.push_back("2e2mu");
    BRcorrs.push_back(BR_Z_nunu_ratio*3. * BR_Z_ll_ratio*(hasTaus ? 3. : 2.));
    filter_corrs.push_back((!has2LFilter ? 1. : (1.-BR_Z_ll_POWHEG))); // 4l veto
    // 2nu2q
    hypos.push_back("2e2mu");
    BRcorrs.push_back(
      (BR_Z_uu_ratio*2.+ BR_Z_dd_ratio*3.)*3. // Reweighting of 2e->2q
      *
      BR_Z_nunu_ratio*3. // Reweighting of 2mu->2nu
    );
    filter_corrs.push_back((!has2LFilter ? 1. : BR_Z_ll_POWHEG));
    // 4nu, different flavors
    hypos.push_back("2e2mu");
    BRcorrs.push_back(std::pow(BR_Z_nunu_ratio, 2)*3.);
    filter_corrs.push_back((!has2LFilter ? 1. : BR_Z_ll_POWHEG));
    // 4nu, same flavors
    hypos.push_back("4e");
    BRcorrs.push_back(std::pow(BR_Z_nunu_ratio, 2)*3.);
    filter_corrs.push_back((!has2LFilter ? 1. : BR_Z_ll_POWHEG));
    break;
  case kZZTo2L2Q:
    hypos.push_back("2e2mu");
    BRcorrs.push_back((BR_Z_uu_ratio*2.+ BR_Z_dd_ratio*3.)*3. * BR_Z_ll_ratio*(hasTaus ? 3. : 2.));
    filter_corrs.push_back((!has2LFilter ? 1. : (1.-BR_Z_ll_POWHEG))); // 4l veto
    break;
  case kZZTo4Q:
    // (A) 4q, different flavors or different colors
    hypos.push_back("2e2mu");
    BRcorrs.push_back(
      BR_Z_uu_ratio*BR_Z_uu_ratio * (9.*1. + 6.*2./2.) // (uucc)xNc^2 + (4u+4c)xNcx(Nc-1)(/2 bc. we are using 2e2mu)
      +
      BR_Z_dd_ratio*BR_Z_dd_ratio * (9.*3. + 6.*3./2.) // (ddss + ddbb + ssbb)xNc^2 + (4d+4s+4b)xNcx(Nc-1)(/2 bc. we are using 2e2mu)
      +
      BR_Z_uu_ratio*BR_Z_dd_ratio * (9.*6.) // (uudd + uuss + uubb + ccdd + ccss + ccbb)xNc^2
    );
    filter_corrs.push_back((!has2LFilter ? 1. : BR_Z_ll_POWHEG));
    // (B) 4q, same flavors
    hypos.push_back("4e");
    BRcorrs.push_back(
      BR_Z_uu_ratio*BR_Z_uu_ratio * (3.*2.) // 4u+4c
      +
      BR_Z_dd_ratio*BR_Z_dd_ratio * (3.*3.) // 4d+4s+4b
    );
    filter_corrs.push_back((!has2LFilter ? 1. : BR_Z_ll_POWHEG));
    // Sum of coefficients (ignoring couplings) (A)*2+(B) = 225 = (5*3)^2
    break;

  case kWWTo2L2Nu:
    hypos.push_back("WW");
    BRcorrs.push_back(std::pow(BR_W_lnu_ratio*(hasTaus ? 3. : 2.), 2));
    filter_corrs.push_back(1);
    break;
  case kWWToLNuQQ:
    hypos.push_back("WW");
    BRcorrs.push_back(BR_W_qq_ratio*Nflavs_CKM*3. * BR_W_lnu_ratio*(hasTaus ? 3. : 2.) * 2.); // x2 at the end for l+nuqqb' + l-nubarqbq'
    filter_corrs.push_back((has2LFilter ? BR_Z_ll_POWHEG : 1.));
    break;
  case kWWToLNuXX:
    // lnulnu
    hypos.push_back("WW");
    BRcorrs.push_back(std::pow(BR_W_lnu_ratio*(hasTaus ? 3. : 2.), 2));
    filter_corrs.push_back(1);
    // lnuqq
    hypos.push_back("WW");
    BRcorrs.push_back(BR_W_qq_ratio*Nflavs_CKM*3. * BR_W_lnu_ratio*(hasTaus ? 3. : 2.) * 2.); // x2 at the end for l+nuqqb' + l-nubarqbq'
    filter_corrs.push_back((has2LOSFilter ? 0.5*BR_W_lnu_POWHEG : 1.)); // x0.5 in filter isto cancel the x2 above because the choice of W+ vs W-H in the POWHEG sample determines the sign of the lepton and quarks in H->lnuqq'
    break;

  default:
    MELAerr << "SampleHelpers::calculateAdjustedHiggsBREff: Decay mode " << dkmode << " is not implemented." << endl;
    assert(0);
    return 1;
  }

  double br_sum = 0;
  double br_sum_filtered = 0;
  for (unsigned int ih=0; ih<hypos.size(); ih++){
    HiggsXSBRReader hxsbrReader("${CMSSW_BASE}/src/CMSDataTools/AnalysisTree/data/HiggsXSBR/YR3.csv", hypos.at(ih));
    MELAout << "Evaluating " << hypos.at(ih) << " at " << sampleMH << endl;
    double br_MH = hxsbrReader.eval_br(sampleMH);
    MELAout << "\t- BR raw, BR corr, filter =  " << br_MH << ", " << BRcorrs.at(ih) << ", " << filter_corrs.at(ih) << endl;
    double br_MH_corr = br_MH * BRcorrs.at(ih);
    double br_MH_corr_filtered = br_MH_corr * filter_corrs.at(ih);
    br_sum += br_MH_corr;
    br_sum_filtered += br_MH_corr_filtered;
  }

  MELAout << "SampleHelpers::calculateAdjustedHiggsBREff: Final BR before / after filter: " << adj_br*br_sum << " / " << adj_br*br_sum_filtered << endl;
  MELAout << "\t- Filter efficiency: " << br_sum_filtered/br_sum << endl;

  return adj_br*br_sum_filtered;
}
