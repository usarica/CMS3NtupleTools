#include "common_includes.h"


void plotMET(){
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal;

  TString const cinput_main = "/home/users/usarica/work/Width_AC_Run2/2L2Nu_MiniAOD/CMSSW_10_2_18/src/CMS3/NtupleMaker/test/manualruns/";

  std::unordered_map<int, TString> samples_HWW;
  samples_HWW[125] = "GluGluHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1_0.root";
  samples_HWW[500] = "GluGluHToWWTo2L2Nu_M500_13TeV_powheg2_JHUGenV714_pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1_0.root";
  samples_HWW[3000] = "GluGluHToWWTo2L2Nu_M3000_13TeV_powheg2_JHUGenV714_pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1_0.root";

  std::unordered_map<int, TString> samples_HZZ;
  samples_HZZ[200] = "GluGluHToZZTo2L2Nu_M200_13TeV_powheg2_JHUGenV7011_pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1_0.root";
  samples_HZZ[500] = "GluGluHToZZTo2L2Nu_M500_13TeV_powheg2_JHUGenV7011_pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_0.root";
  samples_HZZ[3000] = "GluGluHToZZTo2L2Nu_M3000_13TeV_powheg2_JHUGenV7011_pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_0.root";

  BaseTree sample_tree(cinput_main + samples_HZZ[200], "cms3ntuple/Events", "", "");

  MELAout << "Getting muon handles" << endl;
  MuonHandler muonHandler;
  muonHandler.bookBranches(&sample_tree);
  muonHandler.wrapTree(&sample_tree);
  
  MELAout << "Getting electron handles" << endl;
  ElectronHandler electronHandler;
  electronHandler.bookBranches(&sample_tree);
  electronHandler.wrapTree(&sample_tree);
  
  MELAout << "Getting photon handles" << endl;
  PhotonHandler photonHandler;
  photonHandler.bookBranches(&sample_tree);
  photonHandler.wrapTree(&sample_tree);
  
  MELAout << "Getting jet and MET handles" << endl;
  JetMETHandler jetHandler;
  jetHandler.bookBranches(&sample_tree);
  jetHandler.wrapTree(&sample_tree);

  MELAout << "Completed getting the handles..." << endl;

  const int nevents = sample_tree.getSelectedNEvents();
  for (int ev=0; ev<nevents; ev++){
    sample_tree.getSelectedEvent(ev);

    muonHandler.constructMuons(theGlobalSyst);
    electronHandler.constructElectrons(theGlobalSyst);
    photonHandler.constructPhotons(theGlobalSyst);
    jetHandler.constructJetMET(theGlobalSyst);

    METObject* pfpuppimet = jetHandler.getPFPUPPIMET();
    METObject* pfchsmet = jetHandler.getPFCHSMET();
    MELAout << "MET values (PFPUPPI, PFCHS) = ( " << pfpuppimet->met() << ", " << pfchsmet->met() << " )" << endl;
  }

}

