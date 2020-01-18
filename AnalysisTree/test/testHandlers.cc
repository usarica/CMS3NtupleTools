#include "common_includes.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TText.h"
#include "TPaveText.h"
#include "TLegend.h"


void testHandlers(int procsel){
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal;
  SampleHelpers::configure("2018", "hadoop:200101");
  BtagHelpers::setBtagWPType(BtagHelpers::kDeepFlav_Loose);

  TString strSample;
  if (procsel == 0) strSample = "DY_2l_M_10to50";
  else if (procsel == 1) strSample = "DY_2l_M_50";
  else if (procsel == 2) strSample = "TT_2l2nu";
  else if (procsel == 3) strSample = "qqZZ_2l2nu";
  else if (procsel == 4) strSample = "qqWW_2l2nu";
  else if (procsel == 5) strSample = "GGH_M1000_POWHEG";

  std::vector<TString> sampledirs;
  SampleHelpers::constructSamplesList(strSample, theGlobalSyst, sampledirs);

  BaseTree sample_tree(SampleHelpers::getDatasetFileName(sampledirs.front()), EVENTS_TREE_NAME, "", "");
  sample_tree.sampleIdentifier = SampleHelpers::getSampleIdentifier(sampledirs.front());

  // Get cross section
  sample_tree.bookBranch<float>("xsec", 0.f);

  // Get handlers
  SimEventHandler simEventHandler;
  simEventHandler.bookBranches(&sample_tree);
  simEventHandler.wrapTree(&sample_tree);

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

    if (ev==0){
      sample_tree.getVal("xsec", xsec);
      sample_tree.releaseBranch("xsec");
      MELAout << "Sample cross section: " << xsec*1000. << " fb" << endl;
    }

    eventFilter.constructFilters();
    if (ev==0){
      MELAerr << "Available trigger paths are as follows:" << endl;
      for (auto const& hltpath:eventFilter.getHLTPaths()) MELAout << "\t- " << hltpath->name << endl;
    }
    if (!eventFilter.hasMatchingTriggerPath(triggerCheckList)){
      MELAerr << "No matching trigger paths to " << triggerCheckList << endl;
      break;
    }
    if (eventFilter.getTriggerWeight(triggerCheckList) == 0.) continue;

    MELAout << "================" << endl;
    MELAout << "Event " << ev << ":" << endl;
    MELAout << "================" << endl;

    simEventHandler.constructSimEvent(theGlobalSyst);
    MELAout << "Sample chosen data period: " << simEventHandler.getChosenDataPeriod() << " with random number " << simEventHandler.getRandomNumberSeed(SimEventHandler::kDataPeriod) << "." << endl;
    MELAout << "\t- Gen. MET random number: " << simEventHandler.getRandomNumberSeed(SimEventHandler::kGenMETSmear) << endl;
    MELAout << "\t- PU weight: " << simEventHandler.getPileUpWeight() << endl;

    genInfoHandler.constructGenInfo(theGlobalSyst);
    auto const& genInfo = genInfoHandler.getGenInfo();

    float wgt = genInfo->getGenWeight(true);
    MELAout << "Sample gen. weight: " << wgt << endl;

    auto const& lheparticles = genInfoHandler.getLHEParticles();
    MELAout << "LHE particles:" << endl;
    for (auto const& part:lheparticles) MELAout << "\t- [" << part->pdgId() << "]: " << part->p4() << ", status = " << part->status() << endl;

    auto const& genparticles = genInfoHandler.getGenParticles();
    MELAout << "Gen. particles:" << endl;
    for (auto const& part:genparticles){
      MELAout << "\t- [" << part->pdgId() << "]: " << part->p4() << ", status = " << part->status() << endl;
#define GENPARTICLE_VARIABLE(TYPE, NAME, DEFVAL) MELAout << "\t\t- " << #NAME << " = " << part->extras.NAME << endl;
      GENPARTICLE_EXTRA_VARIABLES;
#undef GENPARTICLE_VARIABLE
    }

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

    MELAout << "Event " << (eventFilter.passMETFilters() ? "passes" : "fails") << " MET filters. Available MET filtera are as follows:" << endl;
    for (auto it:eventFilter.getMETFilters()) MELAout << "\t- " << it.first << ": " << it.second << endl;

    ev_acc++;
    if (ev_acc==max_ev_acc) break;
  }
}
