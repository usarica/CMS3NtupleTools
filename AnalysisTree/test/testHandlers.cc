#include "common_includes.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TText.h"
#include "TPaveText.h"
#include "TLegend.h"


void testHandlers(int procsel){
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal;

  TString cinput = "/home/users/usarica/work/Width_AC_Run2/Samples/191212/";

  if (procsel == 0) cinput += "DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/allevents.root";
  else if (procsel == 1) cinput += "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/allevents.root";
  else if (procsel == 2) cinput += "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/allevents.root";
  else if (procsel == 3) cinput += "ZZTo2L2Nu_TuneCP5_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v2/allevents.root";
  else if (procsel == 4) cinput += "WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/allevents.root";
  else if (procsel == 5) cinput += "GluGluHToZZTo2L2Nu_M1000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/allevents.root";

  BaseTree sample_tree(cinput, "cms3ntuple/Events", "", "");

  // Get cross section
  sample_tree.bookBranch<float>("xsec", 0.f);

  // Get handlers
  GenInfoHandler genInfoHandler;
  genInfoHandler.bookBranches(&sample_tree);
  genInfoHandler.wrapTree(&sample_tree);

  MuonHandler muonHandler;
  muonHandler.bookBranches(&sample_tree);
  muonHandler.wrapTree(&sample_tree);

  ElectronHandler electronHandler;
  electronHandler.bookBranches(&sample_tree);
  electronHandler.wrapTree(&sample_tree);

  PhotonHandler photonHandler;
  photonHandler.bookBranches(&sample_tree);
  photonHandler.wrapTree(&sample_tree);

  JetMETHandler jetHandler;
  jetHandler.bookBranches(&sample_tree);
  jetHandler.wrapTree(&sample_tree);

  EventFilterHandler eventFilter;
  eventFilter.bookBranches(&sample_tree);
  eventFilter.wrapTree(&sample_tree);

  std::vector<std::string> triggerCheckList{
    "HLT_Ele32_WPTight_Gsf_v*", "HLT_IsoMu24_v*",

    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*", "HLT_DoubleEle25_CaloIdL_MW_v*",
    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*",

    "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*", "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*","HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*",
    "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*", "HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v*", "HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_v*",

    "HLT_TripleMu_10_5_5_DZ_v*", "HLT_TripleMu_12_10_5_v*"
  };

  DileptonHandler dileptonHandler;

  MELAout << "Completed getting the handles..." << endl;
  sample_tree.silenceUnused();

  float xsec=-1;
  double sum_ee=0, sum_mumu=0, sum_emu=0;
  double sum_ee_selected=0, sum_mumu_selected=0, sum_emu_selected=0;
  const int nevents = sample_tree.getSelectedNEvents();
  size_t ev_acc=0, max_ev_acc=10;
  for (int ev=0; ev<nevents; ev++){
    sample_tree.getSelectedEvent(ev);

    eventFilter.constructFilters();
    if (!eventFilter.hasMatchingTriggerPath(triggerCheckList)){
      MELAerr << "No matching trigger paths to " << triggerCheckList << endl;
      MELAerr << "Available trigger paths are as follows:" << endl;
      for (auto const& hltpath:eventFilter.getHLTPaths()) MELAout << "\t- " << hltpath->name << endl;
      break;
    }
    if (eventFilter.getTriggerWeight(triggerCheckList) == 0.){
      /*
      if (ev % 10 == 0){
        for (auto const& hltpath:eventFilter.getHLTPaths()){ if (hltpath->passTrigger) MELAout << "\t- " << hltpath->name << " paased with (L1 prescale, HLT prescale) = (" << hltpath->L1prescale << ", " << hltpath->HLTprescale << ")" << endl; }
      }
      if (ev>=100) break;
      */
      continue;
    }

    MELAout << "================" << endl;
    MELAout << "Event " << ev << ":" << endl;
    MELAout << "================" << endl;

    if (ev_acc==0){
      sample_tree.getVal("xsec", xsec);
      sample_tree.releaseBranch("xsec");
      MELAout << "Sample cross section: " << xsec*1000. << " fb" << endl;
    }

    genInfoHandler.constructGenInfo(theGlobalSyst);
    auto const& genInfo = genInfoHandler.getGenInfo();

    float wgt = genInfo->getGenWeight(true);
    MELAout << "Sample gen. weight: " << wgt << endl;

    muonHandler.constructMuons(theGlobalSyst);
    auto const& muons = muonHandler.getProducts();
    MELAout << "Muons:" << endl;
    for (auto const& part:muons) MELAout << "\t- [" << part->pdgId() << "]: " << part->p4() << ", (L, M, T) = (" << ParticleSelectionHelpers::isLooseParticle(part) << ", " << ParticleSelectionHelpers::isMediumParticle(part) << ", " << ParticleSelectionHelpers::isTightParticle(part) << ")" << endl;

    electronHandler.constructElectrons(theGlobalSyst);
    auto const& electrons = electronHandler.getProducts();
    MELAout << "Electrons:" << endl;
    for (auto const& part:electrons) MELAout << "\t- [" << part->pdgId() << "]: " << part->p4() << ", (L, M, T) = (" << ParticleSelectionHelpers::isLooseParticle(part) << ", " << ParticleSelectionHelpers::isMediumParticle(part) << ", " << ParticleSelectionHelpers::isTightParticle(part) << ")" << endl;

    photonHandler.constructPhotons(theGlobalSyst);
    auto const& photons = photonHandler.getProducts();
    MELAout << "Photons:" << endl;
    for (auto const& part:photons) MELAout << "\t- [" << part->pdgId() << "]: " << part->p4() << ", (L, M, T) = (" << ParticleSelectionHelpers::isLooseParticle(part) << ", " << ParticleSelectionHelpers::isMediumParticle(part) << ", " << ParticleSelectionHelpers::isTightParticle(part) << ")" << endl;

    jetHandler.constructJetMET(theGlobalSyst, &muons, &electrons, &photons);
    auto const& ak4jets = jetHandler.getAK4Jets();
    auto const& ak8jets = jetHandler.getAK8Jets();
    auto const& pfmet = jetHandler.getPFMET();
    auto const& puppimet = jetHandler.getPFPUPPIMET();
    ParticleObject::LorentzVector_t genmet;
    genmet = ParticleObject::PolarLorentzVector_t(genInfo->extras.genmet_met, 0, genInfo->extras.genmet_metPhi, 0);

    MELAout << "ak4 jets:" << endl;
    for (auto const& part:ak4jets) MELAout << "\t- p4 = " << part->p4() << ", (L, T) = (" << ParticleSelectionHelpers::isLooseJet(part) << ", " << ParticleSelectionHelpers::isTightJet(part) << "), btag value = " << part->getBtagValue() << endl;
    MELAout << "ak8 jets:" << endl;
    for (auto const& part:ak8jets) MELAout << "\t- p4 = " << part->p4() << ", (L, T) = (" << ParticleSelectionHelpers::isLooseJet(part) << ", " << ParticleSelectionHelpers::isTightJet(part) << ")" << endl;
    MELAout << "Gen. MET (pT, phi) = (" << genmet.pt() << ", " << genmet.phi() << ")" << endl;
    MELAout << "PF MET (pT, phi) = (" << pfmet->pt() << ", " << pfmet->phi() << ")" << endl;
    MELAout << "PF PUPPI MET (pT, phi) = (" << puppimet->pt() << ", " << puppimet->phi() << ")" << endl;

    dileptonHandler.constructDileptons(&muons, &electrons);
    auto const& dileptons = dileptonHandler.getProducts();

    DileptonObject* theChosenDilepton = nullptr;
    for (auto const& dilepton:dileptons){
      if (dilepton->isValid() && dilepton->isOS() && dilepton->nTightDaughters()==2){
        theChosenDilepton = dilepton;
        break;
      }
    }
    MELAout << "Dileptons:" << endl;
    if (!theChosenDilepton) MELAout << "No valid OS TT dilepton pair is found." << endl;
    else MELAout << "Found " << dileptons.size() << "dileptons. The chosen dilepton pair has p4 = " << theChosenDilepton->p4() << endl;

    MELAout << "Triggers:" << endl;
    for (auto const& strTrigger:triggerCheckList) MELAout << "\t- Trigger weight(" << strTrigger << ") = " << eventFilter.getTriggerWeight({ strTrigger }) << endl;

    ev_acc++;
    if (ev_acc==max_ev_acc) break;
  }
}
