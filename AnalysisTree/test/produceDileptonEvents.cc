#include <cassert>
#include "common_includes.h"
#include "OffshellCutflow.h"
#include <CMS3/MELAHelpers/interface/CMS3MELAHelpers.h>
#include "TStyle.h"


namespace LooperFunctionHelpers{
  using namespace std;
  using namespace MELAStreamHelpers;
  using namespace OffshellCutflow;

  bool looperRule(BaseTreeLooper*, double const&, SimpleEntry&);


  // Helper options for MET
  bool use_MET_Puppi = false;
  bool use_MET_XYCorr = true;
  bool use_MET_JERCorr = true;
  bool use_MET_ParticleMomCorr = true;
  bool use_MET_p4Preservation = true;
  bool use_MET_corrections =true;

  void setMETOptions(bool use_MET_Puppi_, bool use_MET_XYCorr_, bool use_MET_JERCorr_, bool use_MET_ParticleMomCorr_, bool use_MET_p4Preservation_, bool use_MET_corrections_);

  // Helpers for jets
  bool applyPUIdToAK4Jets = true;
  bool applyTightLeptonVetoIdToAK4Jets = false;

  void setAK4JetSelectionOptions(bool applyPUIdToAK4Jets_, bool applyTightLeptonVetoIdToAK4Jets_);


  // Helpers for b-tagging
  float btag_thr_loose = -1;
  float btag_thr_medium = -1;
  float btag_thr_tight = -1;

  void setBtagWPs();

}
bool LooperFunctionHelpers::looperRule(BaseTreeLooper* theLooper, double const& extWgt, SimpleEntry& commonEntry){
  // Define handlers
#define OBJECT_HANDLER_COMMON_DIRECTIVES \
  HANDLER_DIRECTIVE(EventFilterHandler, eventFilter) \
  HANDLER_DIRECTIVE(PFCandidateHandler, pfcandidateHandler) \
  HANDLER_DIRECTIVE(MuonHandler, muonHandler) \
  HANDLER_DIRECTIVE(ElectronHandler, electronHandler) \
  HANDLER_DIRECTIVE(PhotonHandler, photonHandler) \
  /*HANDLER_DIRECTIVE(SuperclusterHandler, superclusterHandler)*/ \
  /*HANDLER_DIRECTIVE(FSRHandler, fsrHandler)*/ \
  HANDLER_DIRECTIVE(JetMETHandler, jetHandler) \
  HANDLER_DIRECTIVE(IsotrackHandler, isotrackHandler) \
  HANDLER_DIRECTIVE(VertexHandler, vertexHandler)
#define OBJECT_HANDLER_SIM_DIRECTIVES \
  HANDLER_DIRECTIVE(SimEventHandler, simEventHandler) \
  HANDLER_DIRECTIVE(GenInfoHandler, genInfoHandler)
#define OBJECT_HANDLER_DIRECTIVES \
  OBJECT_HANDLER_COMMON_DIRECTIVES \
  OBJECT_HANDLER_SIM_DIRECTIVES
#define SCALEFACTOR_HANDLER_COMMON_DIRECTIVES \
  HANDLER_DIRECTIVE(MuonScaleFactorHandler, muonSFHandler) \
  HANDLER_DIRECTIVE(ElectronScaleFactorHandler, electronSFHandler)
#define SCALEFACTOR_HANDLER_SIM_DIRECTIVES \
  HANDLER_DIRECTIVE(PhotonScaleFactorHandler, photonSFHandler) \
  HANDLER_DIRECTIVE(PUJetIdScaleFactorHandler, pujetidSFHandler) \
  HANDLER_DIRECTIVE(BtagScaleFactorHandler, btagSFHandler) \
  HANDLER_DIRECTIVE(METCorrectionHandler, metCorrectionHandler)
#define SCALEFACTOR_HANDLER_DIRECTIVES \
  SCALEFACTOR_HANDLER_COMMON_DIRECTIVES \
  SCALEFACTOR_HANDLER_SIM_DIRECTIVES

  // Get the current tree
  BaseTree* currentTree = theLooper->getWrappedTree();
  if (!currentTree) return false;

  // Acquire global variables
  SystematicsHelpers::SystematicVariationTypes const& theGlobalSyst = theLooper->getSystematic();
  ParticleDisambiguator& particleDisambiguator = theLooper->getParticleDisambiguator();
  DileptonHandler& dileptonHandler = theLooper->getDileptonHandler();
  CMS3MELAHelpers::GMECBlock& MEblock = theLooper->getMEblock();

  // Acquire sample flags
  bool const& isData = theLooper->getCurrentTreeFlag_IsData();
  bool const& isQCD = theLooper->getCurrentTreeFlag_QCDException();
  float pTG_true_exception_range[2]={ -1, -1 };
  bool hasPTGExceptionRange = theLooper->getPTGExceptionRange(pTG_true_exception_range[0], pTG_true_exception_range[1]);
  bool needGenParticleChecks = isQCD || hasPTGExceptionRange;

  // Acquire triggers
  auto const& triggerCheckListMap = theLooper->getHLTMenus();
  auto const& triggerPropsCheckListMap = theLooper->getHLTMenuProperties();
  bool hasSimpleHLTMenus = theLooper->hasSimpleHLTMenus();
  bool hasHLTMenuProperties = theLooper->hasHLTMenuProperties();
  if (hasSimpleHLTMenus && hasHLTMenuProperties){
    MELAerr << "LooperFunctionHelpers::looperRule: Defining both simple HLT menus and menus with properties is not allowed. Choose only one!" << endl;
    assert(0);
  }
  auto it_HLTMenuSimple_SingleLepton = triggerCheckListMap.find("SingleLepton");
  auto it_HLTMenuProps_SingleLepton = triggerPropsCheckListMap.find("SingleLepton");
  if (
    (hasSimpleHLTMenus && it_HLTMenuSimple_SingleLepton == triggerCheckListMap.cend())
    ||
    (hasHLTMenuProperties && it_HLTMenuProps_SingleLepton == triggerPropsCheckListMap.cend())
    ){
    MELAerr << "LooperFunctionHelpers::looperRule: The trigger type 'SingleLepton' has to be defined in this looper rule!" << endl;
    assert(0);
  }
  auto it_HLTMenuSimple_Dilepton_DF = triggerCheckListMap.find("Dilepton_DF");
  auto it_HLTMenuProps_Dilepton_DF = triggerPropsCheckListMap.find("Dilepton_DF");
  if (
    (hasSimpleHLTMenus && it_HLTMenuSimple_Dilepton_DF == triggerCheckListMap.cend())
    ||
    (hasHLTMenuProperties && it_HLTMenuProps_Dilepton_DF == triggerPropsCheckListMap.cend())
    ){
    MELAerr << "LooperFunctionHelpers::looperRule: The trigger type 'Dilepton_DF' has to be defined in this looper rule!" << endl;
    assert(0);
  }
  auto it_HLTMenuSimple_Dilepton_DF_Extra = triggerCheckListMap.find("Dilepton_DF_Extra");
  auto it_HLTMenuProps_Dilepton_DF_Extra = triggerPropsCheckListMap.find("Dilepton_DF_Extra");
  if (
    (hasSimpleHLTMenus && it_HLTMenuSimple_Dilepton_DF_Extra == triggerCheckListMap.cend())
    ||
    (hasHLTMenuProperties && it_HLTMenuProps_Dilepton_DF_Extra == triggerPropsCheckListMap.cend())
    ){
    MELAerr << "LooperFunctionHelpers::looperRule: The trigger type 'Dilepton_DF_Extra' has to be defined in this looper rule!" << endl;
    assert(0);
  }
  auto it_HLTMenuSimple_Dilepton_SF = triggerCheckListMap.find("Dilepton_SF");
  auto it_HLTMenuProps_Dilepton_SF = triggerPropsCheckListMap.find("Dilepton_SF");
  if (
    (hasSimpleHLTMenus && it_HLTMenuSimple_Dilepton_SF == triggerCheckListMap.cend())
    ||
    (hasHLTMenuProperties && it_HLTMenuProps_Dilepton_SF == triggerPropsCheckListMap.cend())
    ){
    MELAerr << "LooperFunctionHelpers::looperRule: The trigger type 'Dilepton_SF' has to be defined in this looper rule!" << endl;
    assert(0);
  }
  auto it_HLTMenuSimple_PFHT_Control = triggerCheckListMap.find("PFHT_Control");
  auto it_HLTMenuProps_PFHT_Control = triggerPropsCheckListMap.find("PFHT_Control");
  if (
    (hasSimpleHLTMenus && it_HLTMenuSimple_PFHT_Control == triggerCheckListMap.cend())
    ||
    (hasHLTMenuProperties && it_HLTMenuProps_PFHT_Control == triggerPropsCheckListMap.cend())
    ){
    MELAerr << "LooperFunctionHelpers::looperRule: The trigger type 'PFHT_Control' has to be defined in this looper rule!" << endl;
    assert(0);
  }
  auto it_HLTMenuSimple_PFMET_MHT_Control = triggerCheckListMap.find("PFMET_MHT_Control");
  auto it_HLTMenuProps_PFMET_MHT_Control = triggerPropsCheckListMap.find("PFMET_MHT_Control");
  if (
    (hasSimpleHLTMenus && it_HLTMenuSimple_PFMET_MHT_Control == triggerCheckListMap.cend())
    ||
    (hasHLTMenuProperties && it_HLTMenuProps_PFMET_MHT_Control == triggerPropsCheckListMap.cend())
    ){
    MELAerr << "LooperFunctionHelpers::looperRule: The trigger type 'PFMET_MHT_Control' has to be defined in this looper rule!" << endl;
    assert(0);
  }

  // Acquire all handlers
#define HANDLER_DIRECTIVE(TYPE, NAME) TYPE* NAME = nullptr;
  OBJECT_HANDLER_DIRECTIVES;
  SCALEFACTOR_HANDLER_DIRECTIVES;
#undef HANDLER_DIRECTIVE
#define HANDLER_DIRECTIVE(TYPE, NAME) TYPE* tmp_##NAME = dynamic_cast<TYPE*>(handler); if (tmp_##NAME){ NAME = tmp_##NAME; continue; }
  for (auto const& handler:theLooper->getObjectHandlers()){
    OBJECT_HANDLER_COMMON_DIRECTIVES;
    if (!isData){
      OBJECT_HANDLER_SIM_DIRECTIVES;
    }
  }
  for (auto const& handler:theLooper->getSFHandlers()){
    SCALEFACTOR_HANDLER_COMMON_DIRECTIVES;
    if (!isData){
      SCALEFACTOR_HANDLER_SIM_DIRECTIVES;
    }
  }
#undef HANDLER_DIRECTIVE
#define HANDLER_DIRECTIVE(TYPE, NAME) \
  if (!NAME){ \
    MELAerr << "LooperFunctionHelpers::looperRule: " << #TYPE << " " << #NAME << " is not registered. Please register and re-run." << endl; \
    assert(0); \
  }
  OBJECT_HANDLER_COMMON_DIRECTIVES;
  SCALEFACTOR_HANDLER_COMMON_DIRECTIVES;
  if (!isData){
    OBJECT_HANDLER_SIM_DIRECTIVES;
    SCALEFACTOR_HANDLER_SIM_DIRECTIVES;
  }
#undef HANDLER_DIRECTIVE

  /************************/
  /* EVENT INTERPRETATION */
  /************************/
  // Recorded variables
#define BRANCH_SCALAR_COMMANDS \
  BRANCH_COMMAND(float, event_wgt) \
  BRANCH_COMMAND(float, event_wgt_adjustment_NNPDF30) \
  BRANCH_COMMAND(float, event_wgt_triggers_SingleLepton) \
  BRANCH_COMMAND(float, event_wgt_triggers_Dilepton) \
  BRANCH_COMMAND(float, event_wgt_triggers_Dilepton_DF_Extra) \
  BRANCH_COMMAND(float, event_wgt_triggers_PFHT_Control) \
  BRANCH_COMMAND(float, event_wgt_triggers_PFMET_MHT_Control) \
  BRANCH_COMMAND(float, event_wgt_SFs) \
  BRANCH_COMMAND(float, event_wgt_SFs_muons) \
  BRANCH_COMMAND(float, event_wgt_SFs_electrons) \
  BRANCH_COMMAND(float, event_wgt_SFs_photons) \
  BRANCH_COMMAND(float, event_wgt_SFs_PUJetId) \
  BRANCH_COMMAND(float, event_wgt_SFs_btagging) \
  BRANCH_COMMAND(float, event_pTmiss) \
  BRANCH_COMMAND(float, event_phimiss) \
  BRANCH_COMMAND(float, event_mTZZ) \
  BRANCH_COMMAND(float, event_mZZ) \
  BRANCH_COMMAND(bool, event_pass_tightMETFilters) \
  BRANCH_COMMAND(float, genmet_pTmiss) \
  BRANCH_COMMAND(float, genmet_phimiss) \
  BRANCH_COMMAND(unsigned int, event_n_vtxs_good) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt30) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt30_btagged_loose) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt30_btagged_medium) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt20) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt20_btagged_loose) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt20_btagged_medium) \
  BRANCH_COMMAND(float, ak4jets_HT) \
  BRANCH_COMMAND(float, ak4jets_MHT) \
  BRANCH_COMMAND(cms3_id_t, dilepton_id) \
  BRANCH_COMMAND(float, dilepton_pt) \
  BRANCH_COMMAND(float, dilepton_eta) \
  BRANCH_COMMAND(float, dilepton_phi) \
  BRANCH_COMMAND(float, dilepton_mass) \
  BRANCH_COMMAND(float, dPhi_pTboson_pTmiss) \
  BRANCH_COMMAND(float, dPhi_pTbosonjets_pTmiss) \
  BRANCH_COMMAND(float, min_abs_dPhi_pTj_pTmiss)
#define BRANCH_VECTOR_COMMANDS \
  BRANCH_COMMAND(bool, leptons_is_genMatched_prompt) \
  BRANCH_COMMAND(bool, leptons_is_TOmatched_SingleLepton) \
  BRANCH_COMMAND(cms3_id_t, leptons_id) \
  BRANCH_COMMAND(float, leptons_pt) \
  BRANCH_COMMAND(float, leptons_eta) \
  BRANCH_COMMAND(float, leptons_phi) \
  BRANCH_COMMAND(float, leptons_mass) \
  BRANCH_COMMAND(float, leptons_eff) \
  BRANCH_COMMAND(float, leptons_eff_DF) \
  BRANCH_COMMAND(float, electrons_full5x5_sigmaIEtaIEta) \
  BRANCH_COMMAND(float, electrons_full5x5_sigmaIPhiIPhi) \
  BRANCH_COMMAND(float, electrons_full5x5_r9) \
  BRANCH_COMMAND(float, electrons_seedTime) \
  BRANCH_COMMAND(bool, ak4jets_is_genMatched) \
  BRANCH_COMMAND(bool, ak4jets_is_genMatched_fullCone) \
  BRANCH_COMMAND(unsigned char, ak4jets_btagWP_Bits) \
  BRANCH_COMMAND(float, ak4jets_NEMF) \
  BRANCH_COMMAND(float, ak4jets_CEMF) \
  BRANCH_COMMAND(float, ak4jets_pt) \
  BRANCH_COMMAND(float, ak4jets_eta) \
  BRANCH_COMMAND(float, ak4jets_phi) \
  BRANCH_COMMAND(float, ak4jets_mass) \
  BRANCH_COMMAND(float, ak8jets_pt) \
  BRANCH_COMMAND(float, ak8jets_eta) \
  BRANCH_COMMAND(float, ak8jets_phi) \
  BRANCH_COMMAND(float, ak8jets_mass)
#define BRANCH_COMMANDS \
  BRANCH_SCALAR_COMMANDS \
  BRANCH_VECTOR_COMMANDS

#define BRANCH_COMMAND(TYPE, NAME) TYPE NAME = 0;
  BRANCH_SCALAR_COMMANDS;
#undef BRANCH_COMMAND
#define BRANCH_COMMAND(TYPE, NAME) std::vector<TYPE> NAME;
  BRANCH_VECTOR_COMMANDS;
#undef BRANCH_COMMAND


  // Always assign the external weight first
  event_wgt = extWgt;
  // Set NNPDF 3.0 adjustment to 1
  event_wgt_adjustment_NNPDF30 = 1;

  if (!isData){
    genInfoHandler->constructGenInfo(theGlobalSyst);
    auto const& genInfo = genInfoHandler->getGenInfo();
    double genwgt_NNPDF30 = genInfo->getGenWeight(false);
    double genwgt_default = genInfo->getGenWeight(true);
    event_wgt_adjustment_NNPDF30 = (genwgt_default!=0. ? genwgt_NNPDF30 / genwgt_default : 0.);
    event_wgt *= genwgt_default;
    genmet_pTmiss = genInfo->extras.genmet_met;
    genmet_phimiss = genInfo->extras.genmet_metPhi;
    auto const& genparticles = genInfoHandler->getGenParticles();

    if (needGenParticleChecks){
      if (isQCD){
        for (auto const& part:genparticles){
          if (PDGHelpers::isAPhoton(part->pdgId()) && part->extras.isPromptFinalState && part->pt()>=25.f){
            event_wgt = 0;
            break;
          }
        }
      }
      if (hasPTGExceptionRange){
        for (auto const& part:genparticles){
          if (PDGHelpers::isAPhoton(part->pdgId()) && part->extras.isHardProcess){
            if ((pTG_true_exception_range[0]>=0.f && part->pt()<pTG_true_exception_range[0]) || (pTG_true_exception_range[1]>=0.f && part->pt()>=pTG_true_exception_range[1])) event_wgt = 0;
            break;
          }
        }
      }
    }

    simEventHandler->constructSimEvent(theGlobalSyst);
    event_wgt *= simEventHandler->getPileUpWeight()*simEventHandler->getL1PrefiringWeight();

    if (event_wgt==0.f) return false;

    // Record LHE MEs and K factors
    for (auto const& it:genInfo->extras.LHE_ME_weights) commonEntry.setNamedVal(it.first, it.second);
    for (auto const& it:genInfo->extras.Kfactors) commonEntry.setNamedVal(it.first, it.second);
  }
  theLooper->incrementSelection("Valid gen. weights");

  vertexHandler->constructVertices();
  if (!vertexHandler->hasGoodPrimaryVertex()) return false;
  event_n_vtxs_good = vertexHandler->getNGoodVertices();
  theLooper->incrementSelection("Good vertices");

  eventFilter->constructFilters(simEventHandler);
  if (isData && !eventFilter->isUniqueDataEvent()) return false;
  theLooper->incrementSelection("Unique data");

  if (!eventFilter->passCommonSkim() || !eventFilter->passMETFilters(EventFilterHandler::kMETFilters_Standard)) return false;
  theLooper->incrementSelection("MET filters");
  event_pass_tightMETFilters = eventFilter->passMETFilters(EventFilterHandler::kMETFilters_Tight);

  pfcandidateHandler->constructPFCandidates(theGlobalSyst);
  auto const& pfcandidates = pfcandidateHandler->getProducts();

  muonHandler->constructMuons(theGlobalSyst, &pfcandidates);
  electronHandler->constructElectrons(theGlobalSyst, &pfcandidates);
  photonHandler->constructPhotons(theGlobalSyst, &pfcandidates);
  particleDisambiguator.disambiguateParticles(muonHandler, electronHandler, photonHandler);

  std::unordered_map<ParticleObject const*, float> lepton_eff_map;

  auto const& muons = muonHandler->getProducts();
  float SF_muons = 1;
  for (auto const& part:muons){
    float theSF = 1;
    float theEff = 1;
    if (!isData){
      muonSFHandler->getIdIsoSFAndEff(theGlobalSyst, part, theSF, &theEff);
      if (theSF == 0.f) continue;
      SF_muons *= theSF;
    }
    else if (ParticleSelectionHelpers::isTightParticle(part)){
      muonSFHandler->getIdIsoSFAndEff(theGlobalSyst, part->pt(), part->eta(), true, true, true, theSF, &theEff);
    }
    lepton_eff_map[part] = theEff;
  }
  event_wgt_SFs_muons = SF_muons;

  auto const& electrons = electronHandler->getProducts();
  float SF_electrons = 1;
  for (auto const& part:electrons){
    float theSF = 1;
    float theEff = 1;
    if (!isData){
      electronSFHandler->getIdIsoSFAndEff(theGlobalSyst, part, theSF, &theEff);
      if (theSF == 0.f) continue;
      SF_electrons *= theSF;
    }
    else if (ParticleSelectionHelpers::isTightParticle(part)){
      electronSFHandler->getIdIsoSFAndEff(theGlobalSyst, part->pt(), part->etaSC(), part->isGap(), true, true, true, theSF, &theEff);
    }
    lepton_eff_map[part] = theEff;
  }
  event_wgt_SFs_electrons = SF_electrons;

  auto const& photons = photonHandler->getProducts();
  unsigned int n_photons_veto = 0;
  float SF_photons = 1;
  for (auto const& part:photons){
    float theSF = 1;
    if (!isData) photonSFHandler->getIdIsoSFAndEff(theGlobalSyst, part, theSF, nullptr);
    if (theSF == 0.f) continue;
    SF_photons *= theSF;

    if (ParticleSelectionHelpers::isVetoParticle(part)) n_photons_veto++;
  }
  event_wgt_SFs_photons = SF_photons;
  if (n_photons_veto!=0) return false;
  theLooper->incrementSelection("Photon veto");

  isotrackHandler->constructIsotracks(&muons, &electrons);
  bool hasVetoIsotrack = false;
  for (auto const& isotrack:isotrackHandler->getProducts()){
    if (isotrack->testSelectionBit(IsotrackSelectionHelpers::kPreselectionVeto)){
      hasVetoIsotrack = true;
      break;
    }
  }
  if (hasVetoIsotrack) return false;
  theLooper->incrementSelection("Isotrack veto");

  dileptonHandler.constructDileptons(&muons, &electrons);
  auto const& dileptons = dileptonHandler.getProducts();
  DileptonObject* theChosenDilepton = nullptr;
  size_t nTightDilep = 0;
  for (auto const& dilepton:dileptons){
    if (dilepton->isValid() && dilepton->isOS() && dilepton->nTightDaughters()==2){
      if (!theChosenDilepton) theChosenDilepton = dilepton;
      nTightDilep++;
    }
  }
  if (!theChosenDilepton || nTightDilep>1) return false;
  theLooper->incrementSelection("Tight, dilepton pair selection");

  bool const dilepton_is_SF = theChosenDilepton->isSF();
  dilepton_id = theChosenDilepton->getDaughter(0)->pdgId()*theChosenDilepton->getDaughter(1)->pdgId();
  dilepton_pt = theChosenDilepton->pt();
  dilepton_eta = theChosenDilepton->eta();
  dilepton_phi = theChosenDilepton->phi();
  dilepton_mass = theChosenDilepton->m();
  // Filter out low-mass resonances
  if (dilepton_mass<12.f) return false;
  theLooper->incrementSelection("Low mass resonance veto");

  jetHandler->constructJetMET(theGlobalSyst, &muons, &electrons, &photons, &pfcandidates);
  auto const& ak4jets = jetHandler->getAK4Jets();
  auto const& ak8jets = jetHandler->getAK8Jets();
  auto const& eventmet = (use_MET_Puppi ? jetHandler->getPFPUPPIMET() : jetHandler->getPFMET());
  if (!isData && use_MET_corrections) metCorrectionHandler->applyCorrections(
    simEventHandler->getChosenDataPeriod(),
    genmet_pTmiss, genmet_phimiss,
    eventmet, !use_MET_Puppi,
    &(simEventHandler->getRandomNumber(SimEventHandler::kGenMETSmear))
  );
  auto event_met_p4 = eventmet->p4(use_MET_XYCorr, use_MET_JERCorr, use_MET_ParticleMomCorr, use_MET_p4Preservation);
  event_pTmiss = event_met_p4.Pt();
  event_phimiss = event_met_p4.Phi();

  std::vector<ParticleObject const*> leptons_TOmatched_SingleLepton;
  if (hasSimpleHLTMenus){
    event_wgt_triggers_SingleLepton = eventFilter->getTriggerWeight(it_HLTMenuSimple_SingleLepton->second);
    event_wgt_triggers_Dilepton = eventFilter->getTriggerWeight((dilepton_is_SF ? it_HLTMenuSimple_Dilepton_SF : it_HLTMenuSimple_Dilepton_DF)->second);
    if (!dilepton_is_SF) event_wgt_triggers_Dilepton_DF_Extra = eventFilter->getTriggerWeight(it_HLTMenuSimple_Dilepton_DF_Extra->second);
    event_wgt_triggers_PFHT_Control = eventFilter->getTriggerWeight(it_HLTMenuSimple_PFHT_Control->second);
    event_wgt_triggers_PFMET_MHT_Control = eventFilter->getTriggerWeight(it_HLTMenuSimple_PFMET_MHT_Control->second);
  }
  else if (hasHLTMenuProperties){
    event_wgt_triggers_SingleLepton = eventFilter->getTriggerWeight(
      it_HLTMenuProps_SingleLepton->second,
      &muons, &electrons, nullptr, nullptr, nullptr, nullptr,
      nullptr, &leptons_TOmatched_SingleLepton
    );
    event_wgt_triggers_Dilepton = eventFilter->getTriggerWeight(
      (dilepton_is_SF ? it_HLTMenuProps_Dilepton_SF : it_HLTMenuProps_Dilepton_DF)->second,
      &muons, &electrons, nullptr, nullptr, nullptr, nullptr
    );
    if (!dilepton_is_SF) event_wgt_triggers_Dilepton_DF_Extra = eventFilter->getTriggerWeight(
      it_HLTMenuProps_Dilepton_DF_Extra->second,
      &muons, &electrons, nullptr, nullptr, nullptr, nullptr
    );
    event_wgt_triggers_PFHT_Control = eventFilter->getTriggerWeight(
      it_HLTMenuProps_PFHT_Control->second,
      nullptr, nullptr, nullptr, &ak4jets, nullptr, nullptr
    );
    event_wgt_triggers_PFMET_MHT_Control = eventFilter->getTriggerWeight(
      it_HLTMenuProps_PFMET_MHT_Control->second,
      nullptr, nullptr, nullptr, &ak4jets, nullptr, eventmet
    );
  }
  if ((event_wgt_triggers_SingleLepton + event_wgt_triggers_Dilepton + event_wgt_triggers_Dilepton_DF_Extra + event_wgt_triggers_PFHT_Control + event_wgt_triggers_PFMET_MHT_Control) == 0.f) return false;
  theLooper->incrementSelection("Trigger");

  // Test HEM filter
  if (!eventFilter->test2018HEMFilter(simEventHandler, nullptr, nullptr, &ak4jets, &ak8jets)) return false;
  theLooper->incrementSelection("HEM15/16 veto");

  // Fill leptons after trigger checks
  for (auto const& dau:theChosenDilepton->getDaughters()){
    MuonObject* dau_muon = dynamic_cast<MuonObject*>(dau);
    ElectronObject* dau_electron = dynamic_cast<ElectronObject*>(dau);

    bool is_genMatched_prompt = (dau_muon ? dau_muon->extras.is_genMatched_prompt : dau_electron->extras.is_genMatched_prompt);
    leptons_is_genMatched_prompt.push_back(is_genMatched_prompt);
    leptons_id.push_back(dau->pdgId());
    leptons_pt.push_back(dau->pt());
    leptons_eta.push_back(dau->eta());
    leptons_phi.push_back(dau->phi());
    leptons_mass.push_back(dau->m());
    leptons_is_TOmatched_SingleLepton.push_back(HelperFunctions::checkListVariable(leptons_TOmatched_SingleLepton, (ParticleObject const*) dau));

    // Fill efficiency for the lepton
    auto it_eff = lepton_eff_map.find(dau);
    leptons_eff.push_back((it_eff==lepton_eff_map.end() ? 1 : it_eff->second));

    // Compute efficiency for the oppsite flavor (mu <-> e)
    {
      float eff_DF = 1, SF_dummy = 1;
      if (dau_muon) electronSFHandler->getIdIsoSFAndEff(theGlobalSyst, dau->pt(), dau->eta(), 2, true, true, true, SF_dummy, &eff_DF);
      else muonSFHandler->getIdIsoSFAndEff(theGlobalSyst, dau->pt(), dau_electron->etaSC(), true, true, true, SF_dummy, &eff_DF);
      leptons_eff_DF.push_back(eff_DF);
    }

    // Extra variables for e/gamma efficiency studies
    if (dau_electron){
      auto const& extras = dau_electron->extras;

      electrons_full5x5_sigmaIEtaIEta.push_back(extras.full5x5_sigmaIEtaIEta);
      electrons_full5x5_sigmaIPhiIPhi.push_back(extras.full5x5_sigmaIPhiIPhi);
      electrons_full5x5_r9.push_back(extras.full5x5_r9);
      electrons_seedTime.push_back(extras.seedTime);
    }
  }

  ParticleObject::LorentzVector_t sump4_ak4jets(0, 0, 0, 0);
  std::vector<AK4JetObject*> ak4jets_tight; ak4jets_tight.reserve(ak4jets.size());
  unsigned int n_ak4jets_tight_pt30_btagged_loose = 0;
  float SF_PUJetId = 1;
  float SF_btagging = 1;
  for (auto const& jet:ak4jets){
    if (!isData){
      float theSF_PUJetId = 1;
      float theSF_btag = 1;
      pujetidSFHandler->getSFAndEff(theGlobalSyst, jet, theSF_PUJetId, nullptr);
      btagSFHandler->getSFAndEff(theGlobalSyst, jet, theSF_btag, nullptr);
      if (theSF_PUJetId != 0.f) SF_PUJetId *= theSF_PUJetId;
      if (theSF_btag != 0.f) SF_btagging *= theSF_btag;
    }

    if (ParticleSelectionHelpers::isTightJet(jet)){
      ak4jets_tight.push_back(jet);
      if (jet->getBtagValue()>=btag_thr_loose) event_n_ak4jets_pt30_btagged_loose++;
      if (jet->getBtagValue()>=btag_thr_medium) event_n_ak4jets_pt30_btagged_medium++;

      ak4jets_HT += jet->pt();
      sump4_ak4jets += jet->p4();
    }

    // Check for pT>20 objects and count them
    if (
      jet->testSelectionBit(AK4JetSelectionHelpers::bit_preselectionTight_id)
      &&
      (!applyPUIdToAK4Jets || jet->testSelectionBit(AK4JetSelectionHelpers::kTightPUJetId))
      &&
      (!applyTightLeptonVetoIdToAK4Jets || jet->testSelectionBit(AK4JetSelectionHelpers::kTightLeptonVetoId))
      &&
      jet->pt()>=20.f && std::abs(jet->eta())<AK4JetSelectionHelpers::etaThr_skim_tight
      ){
      event_n_ak4jets_pt20++;
      if (jet->getBtagValue()>=btag_thr_loose) event_n_ak4jets_pt20_btagged_loose++;
      if (jet->getBtagValue()>=btag_thr_medium) event_n_ak4jets_pt20_btagged_medium++;
    }
  }
  ak4jets_MHT = sump4_ak4jets.Pt();
  event_wgt_SFs_PUJetId = SF_PUJetId;
  event_wgt_SFs_btagging = SF_btagging;
  event_n_ak4jets_pt30 = ak4jets_tight.size();

  // Accumulate ak8 jets as well
  for (auto const& jet:ak8jets){
    if (!ParticleSelectionHelpers::isTightJet(jet)) continue;
    if (jet->mass()<60.f) continue;

    ak8jets_pt.push_back(jet->pt());
    ak8jets_eta.push_back(jet->eta());
    ak8jets_phi.push_back(jet->phi());
    ak8jets_mass.push_back(jet->mass());
  }

  min_abs_dPhi_pTj_pTmiss = TMath::Pi();
  for (auto const& jet:ak4jets_tight){
    ak4jets_pt.push_back(jet->pt());
    ak4jets_eta.push_back(jet->eta());
    ak4jets_phi.push_back(jet->phi());
    ak4jets_mass.push_back(jet->mass());
    ak4jets_is_genMatched.push_back(jet->extras.is_genMatched);
    ak4jets_is_genMatched_fullCone.push_back(jet->extras.is_genMatched_fullCone);
    ak4jets_NEMF.push_back(jet->extras.NEMF);
    ak4jets_CEMF.push_back(jet->extras.CEMF);

    // Determine b-tag WP passing bits
    {
      unsigned char btag_bits=0;
      if (jet->getBtagValue()>=btag_thr_loose) HelperFunctions::set_bit(btag_bits, 0, true);
      if (jet->getBtagValue()>=btag_thr_medium) HelperFunctions::set_bit(btag_bits, 1, true);
      if (jet->getBtagValue()>=btag_thr_tight) HelperFunctions::set_bit(btag_bits, 2, true);
      ak4jets_btagWP_Bits.push_back(btag_bits);
    }

    // Determine min_abs_dPhi_pTj_pTmiss
    float dphi_tmp;
    HelperFunctions::deltaPhi(float(jet->phi()), event_phimiss, dphi_tmp); dphi_tmp = std::abs(dphi_tmp);
    min_abs_dPhi_pTj_pTmiss = std::min(min_abs_dPhi_pTj_pTmiss, dphi_tmp);
  }
  
  // Compute dPhi between the dilepton and pTmiss vector
  dPhi_pTboson_pTmiss = theChosenDilepton->deltaPhi(event_phimiss);
  HelperFunctions::deltaPhi(float((theChosenDilepton->p4()+sump4_ak4jets).Phi()), event_phimiss, dPhi_pTbosonjets_pTmiss);

  // Compute mass variables
  event_mTZZ = std::sqrt(
    std::pow(
    (
      std::sqrt(std::pow(dilepton_pt, 2) + std::pow(dilepton_mass, 2))
      + std::sqrt(std::pow(event_pTmiss, 2) + std::pow(PDGHelpers::Zmass, 2))
      ), 2
    )
    - std::pow((theChosenDilepton->p4() + event_met_p4).Pt(), 2)
  );

  ParticleObject::LorentzVector_t p4_ZZ_approx;
  float const& etamiss_approx = dilepton_eta;
  p4_ZZ_approx = ParticleObject::PolarLorentzVector_t(event_pTmiss, etamiss_approx, event_phimiss, PDGHelpers::Zmass);
  p4_ZZ_approx = p4_ZZ_approx + theChosenDilepton->p4();
  event_mZZ = p4_ZZ_approx.M();

  // Compute MEs
  if (theLooper->hasRecoMEs()){
    SimpleParticleCollection_t daughters;
    daughters.push_back(SimpleParticle_t(25, ParticleObjectHelpers::convertCMSLorentzVectorToTLorentzVector(p4_ZZ_approx)));

    SimpleParticleCollection_t associated;
    for (auto const& jet:ak4jets_tight) associated.push_back(SimpleParticle_t(0, ParticleObjectHelpers::convertCMSLorentzVectorToTLorentzVector(jet->p4())));

    CMS3MELAHelpers::melaHandle->setCandidateDecayMode(TVar::CandidateDecay_Stable);
    CMS3MELAHelpers::melaHandle->setInputEvent(&daughters, &associated, nullptr, false);
    MEblock.computeMELABranches();
    MEblock.pushMELABranches();

    std::unordered_map<std::string, float> ME_values;
    MEblock.getBranchValues(ME_values);
    CMS3MELAHelpers::melaHandle->resetInputEvent();

    // Insert the ME values into commonEntry only when the productTreeList collection is empty.
    // Otherwise, the branches are already made.
    if (!theLooper->hasLinkedOutputTrees()){ for (auto const& it:ME_values) commonEntry.setNamedVal(it.first, it.second); }
  }

  // Set the collection of SFs at the last step
  event_wgt_SFs = SF_muons*SF_electrons*SF_photons*SF_PUJetId*SF_btagging;

  /*********************/
  /* RECORD THE OUTPUT */
  /*********************/
#define BRANCH_COMMAND(TYPE, NAME) commonEntry.setNamedVal(#NAME, NAME);
  BRANCH_COMMANDS;
#undef BRANCH_COMMAND

  return true;

#undef SCALEFACTOR_HANDLER_DIRECTIVES
#undef OBJECT_HANDLER_DIRECTIVES
}

void LooperFunctionHelpers::setMETOptions(bool use_MET_Puppi_, bool use_MET_XYCorr_, bool use_MET_JERCorr_, bool use_MET_ParticleMomCorr_, bool use_MET_p4Preservation_, bool use_MET_corrections_){
  use_MET_Puppi = use_MET_Puppi_;
  use_MET_XYCorr = use_MET_XYCorr_;
  use_MET_JERCorr = use_MET_JERCorr_;
  use_MET_ParticleMomCorr = use_MET_ParticleMomCorr_;
  use_MET_p4Preservation = use_MET_p4Preservation_;
  use_MET_corrections = use_MET_corrections_;
}

void LooperFunctionHelpers::setAK4JetSelectionOptions(bool applyPUIdToAK4Jets_, bool applyTightLeptonVetoIdToAK4Jets_){
  applyPUIdToAK4Jets = applyPUIdToAK4Jets_;
  applyTightLeptonVetoIdToAK4Jets = applyTightLeptonVetoIdToAK4Jets_;
}

void LooperFunctionHelpers::setBtagWPs(){
  std::vector<float> vwps = BtagHelpers::getBtagWPs(false);
  btag_thr_loose = vwps.at(0);
  btag_thr_medium = vwps.at(1);
  btag_thr_tight = vwps.at(2);
}


using namespace SystematicsHelpers;
void getTrees(
  TString strSampleSet, TString period,
  TString prodVersion, TString strdate,
  int ichunk, int nchunks,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal,
  // ME option
  bool computeMEs=false,
  // Jet ID options
  bool applyPUIdToAK4Jets=true, bool applyTightLeptonVetoIdToAK4Jets=false,
  // MET options
  bool use_MET_Puppi=false,
  bool use_MET_XYCorr=true, bool use_MET_JERCorr=true, bool use_MET_ParticleMomCorr=true, bool use_MET_p4Preservation=true, bool use_MET_corrections=true
){
  if (!SampleHelpers::checkRunOnCondor()) std::signal(SIGINT, SampleHelpers::setSignalInterrupt);

  constexpr bool useJetOverlapStripping = false; // Keep overlap removal turned off

  if (nchunks==1){ nchunks = 0; ichunk=0; }
  if (nchunks>0 && (ichunk<0 || ichunk==nchunks)) return;

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  bool useSkims = !(
    strSampleSet.Contains("GluGluH") || strSampleSet.Contains("GGH")
    ||
    strSampleSet.Contains("VBF")
    ||
    strSampleSet.Contains("WminusH")
    ||
    strSampleSet.Contains("WplusH")
    ||
    strSampleSet.Contains("ZH")
    ||
    strSampleSet.Contains("JHUGen") || strSampleSet.Contains("JHUgen") || strSampleSet.Contains("jhugen")
    );
  SampleHelpers::configure(period, Form("%s:%s", (useSkims ? "hadoop_skims" : "hadoop"), prodVersion.Data()));

  const float lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::getDataPeriod());

  std::vector<TString> const validDataPeriods = SampleHelpers::getValidDataPeriods();
  size_t const nValidDataPeriods = validDataPeriods.size();

  BtagHelpers::setBtagWPType(BtagHelpers::kDeepFlav_Loose);
  LooperFunctionHelpers::setBtagWPs();

  std::vector<TriggerHelpers::TriggerType> requiredTriggers_SingleLepton{
    TriggerHelpers::kSingleMu, TriggerHelpers::kSingleMu_HighPt,
    TriggerHelpers::kSingleEle, TriggerHelpers::kSingleEle_HighPt
  };
  std::vector<TriggerHelpers::TriggerType> requiredTriggers_Dilepton_DF{
    TriggerHelpers::kMuEle
  };
  std::vector<TriggerHelpers::TriggerType> requiredTriggers_Dilepton_DF_Extra{
    TriggerHelpers::kMuEle_Extra
  };
  std::vector<TriggerHelpers::TriggerType> requiredTriggers_Dilepton_SF{
    TriggerHelpers::kDoubleMu,
    TriggerHelpers::kDoubleEle, TriggerHelpers::kDoubleEle_HighPt
  };
  std::vector<TriggerHelpers::TriggerType> requiredTriggers_PFHT_Control{ TriggerHelpers::kPFHT_Control };
  std::vector<TriggerHelpers::TriggerType> requiredTriggers_PFMET_MHT_Control{ TriggerHelpers::kPFMET_MHT_Control };
  auto triggerPropsCheckList_SingleLepton = TriggerHelpers::getHLTMenuProperties(requiredTriggers_SingleLepton);
  auto triggerPropsCheckList_Dilepton_DF = TriggerHelpers::getHLTMenuProperties(requiredTriggers_Dilepton_DF);
  auto triggerPropsCheckList_Dilepton_DF_Extra = TriggerHelpers::getHLTMenuProperties(requiredTriggers_Dilepton_DF_Extra);
  auto triggerPropsCheckList_Dilepton_SF = TriggerHelpers::getHLTMenuProperties(requiredTriggers_Dilepton_SF);
  auto triggerPropsCheckList_PFHT_Control = TriggerHelpers::getHLTMenuProperties(requiredTriggers_PFHT_Control);
  auto triggerPropsCheckList_PFMET_MHT_Control = TriggerHelpers::getHLTMenuProperties(requiredTriggers_PFMET_MHT_Control);

  // Get sample specifications
  std::vector<TString> sampledirs;
  SampleHelpers::constructSamplesList(strSampleSet, theGlobalSyst, sampledirs);
  if (sampledirs.empty()) return;
  bool isData = SampleHelpers::checkSampleIsData(sampledirs.front());
  if (isData && (nchunks>0 || theGlobalSyst!=sNominal)) return;

  // Set flags for ak4jet tight id
  AK4JetSelectionHelpers::setPUIdWP(applyPUIdToAK4Jets ? AK4JetSelectionHelpers::kTightPUJetId : AK4JetSelectionHelpers::nSelectionBits); // Default is 'tight'
  AK4JetSelectionHelpers::setApplyTightLeptonVetoIdToJets(applyTightLeptonVetoIdToAK4Jets); // Default is 'false'
  LooperFunctionHelpers::setAK4JetSelectionOptions(applyPUIdToAK4Jets, applyTightLeptonVetoIdToAK4Jets);

  // Set flags for MET
  LooperFunctionHelpers::setMETOptions(use_MET_Puppi, use_MET_XYCorr, use_MET_JERCorr, use_MET_ParticleMomCorr, use_MET_p4Preservation, use_MET_corrections);

  // Set output directory
  TString coutput_main =
    "output/DileptonEvents/SkimTrees/" + strdate
    + "/AK4Jets"
    + "_" + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
    + "_" + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId")
    + "_" + (useJetOverlapStripping ? "ParticleStripped" : "ParticleCleaned")
    + "/" + (use_MET_Puppi ? "PUPPIMET" : "PFMET")
    + "_" + (use_MET_XYCorr ? "WithXY" : "NoXY");
  if (!isData) coutput_main = coutput_main + "_" + (use_MET_JERCorr ? "WithJER" : "NoJER");
  coutput_main = coutput_main
    + "_" + (use_MET_ParticleMomCorr ? "WithPartMomCorr" : "NoPartMomCorr")
    + "_" + (use_MET_p4Preservation ? "P4Preserved" : "P4Default");
  if (!isData) coutput_main = coutput_main + "_" + (use_MET_corrections ? "ResolutionCorrected" : "ResolutionUncorrected");
  coutput_main = coutput_main + "/" + period;

  TDirectory* curdir = gDirectory;
  gSystem->mkdir(coutput_main, true);

  curdir->cd();

  // Create the output file
  TString coutput = SampleHelpers::getSampleIdentifier(strSampleSet);
  HelperFunctions::replaceString(coutput, "_MINIAODSIM", "");
  HelperFunctions::replaceString(coutput, "_MINIAOD", "");
  TString stroutput = Form("%s/%s", coutput_main.Data(), coutput.Data());
  if (nchunks>0) stroutput = stroutput + Form("_%i_of_%i", ichunk, nchunks);
  stroutput += Form("_%s", SystematicsHelpers::getSystName(theGlobalSyst).data());
  stroutput += ".root";
  TFile* foutput = TFile::Open(stroutput, "recreate");
  foutput->cd();
  BaseTree* tout = new BaseTree("SkimTree");
  MELAout << "Created output file " << stroutput << "..." << endl;
  curdir->cd();

  // Declare handlers
  SimEventHandler simEventHandler;
  GenInfoHandler genInfoHandler;
  EventFilterHandler eventFilter;
  PFCandidateHandler pfcandidateHandler;
  MuonHandler muonHandler;
  ElectronHandler electronHandler;
  PhotonHandler photonHandler;
  //SuperclusterHandler superclusterHandler;
  //FSRHandler fsrHandler;
  JetMETHandler jetHandler;
  IsotrackHandler isotrackHandler;
  VertexHandler vertexHandler;

  // Require trigger matching
  eventFilter.setCheckTriggerObjectsForHLTPaths(true);

  OverlapMapHandler<MuonObject, AK4JetObject> overlapMap_muons_ak4jets;
  OverlapMapHandler<MuonObject, AK8JetObject> overlapMap_muons_ak8jets;
  OverlapMapHandler<ElectronObject, AK4JetObject> overlapMap_electrons_ak4jets;
  OverlapMapHandler<ElectronObject, AK8JetObject> overlapMap_electrons_ak8jets;
  OverlapMapHandler<PhotonObject, AK4JetObject> overlapMap_photons_ak4jets;
  OverlapMapHandler<PhotonObject, AK8JetObject> overlapMap_photons_ak8jets;
  if (useJetOverlapStripping){
    muonHandler.registerOverlapMaps(
      overlapMap_muons_ak4jets,
      overlapMap_muons_ak8jets
    );
    electronHandler.registerOverlapMaps(
      overlapMap_electrons_ak4jets,
      overlapMap_electrons_ak8jets
    );
    photonHandler.registerOverlapMaps(
      overlapMap_photons_ak4jets,
      overlapMap_photons_ak8jets
    );
    jetHandler.registerOverlapMaps(
      overlapMap_muons_ak4jets,
      overlapMap_muons_ak8jets,
      overlapMap_electrons_ak4jets,
      overlapMap_electrons_ak8jets,
      overlapMap_photons_ak4jets,
      overlapMap_photons_ak8jets
    );
  }

  MuonScaleFactorHandler muonSFHandler;
  ElectronScaleFactorHandler electronSFHandler;
  PhotonScaleFactorHandler photonSFHandler;
  PUJetIdScaleFactorHandler pujetidSFHandler;
  BtagScaleFactorHandler btagSFHandler;
  METCorrectionHandler metCorrectionHandler;

  genInfoHandler.setAcquireLHEMEWeights(false);
  genInfoHandler.setAcquireLHEParticles(false);
  genInfoHandler.setAcquireGenParticles(true);

  curdir->cd();

  // Configure the looper
  BaseTreeLooper theLooper;
  // Set chunk index
  if (nchunks>0){
    theLooper.setEventIndexRangeBySampleChunks(true);
    theLooper.setEventIndexRange(ichunk, nchunks);
  }
  // Set systematic
  theLooper.setSystematic(theGlobalSyst);
  // Set looper function
  theLooper.setLooperFunction(LooperFunctionHelpers::looperRule);
  // Set object handlers
  theLooper.addObjectHandler(&simEventHandler);
  theLooper.addObjectHandler(&genInfoHandler);
  theLooper.addObjectHandler(&eventFilter);
  theLooper.addObjectHandler(&pfcandidateHandler);
  theLooper.addObjectHandler(&muonHandler);
  theLooper.addObjectHandler(&electronHandler);
  theLooper.addObjectHandler(&photonHandler);
  //theLooper.addObjectHandler(&superclusterHandler);
  //theLooper.addObjectHandler(&fsrHandler);
  theLooper.addObjectHandler(&jetHandler);
  theLooper.addObjectHandler(&isotrackHandler);
  theLooper.addObjectHandler(&vertexHandler);
  if (useJetOverlapStripping){
    theLooper.addObjectHandler(&overlapMap_muons_ak4jets);
    theLooper.addObjectHandler(&overlapMap_muons_ak8jets);
    theLooper.addObjectHandler(&overlapMap_electrons_ak4jets);
    theLooper.addObjectHandler(&overlapMap_electrons_ak8jets);
    theLooper.addObjectHandler(&overlapMap_photons_ak4jets);
    theLooper.addObjectHandler(&overlapMap_photons_ak8jets);
  }

  // Set SF handlers
  theLooper.addSFHandler(&muonSFHandler);
  theLooper.addSFHandler(&electronSFHandler);
  theLooper.addSFHandler(&photonSFHandler);
  theLooper.addSFHandler(&pujetidSFHandler);
  theLooper.addSFHandler(&btagSFHandler);
  theLooper.addSFHandler(&metCorrectionHandler);
  // Set output tree
  theLooper.addOutputTree(tout);
  // Register the HLT menus
  theLooper.addHLTMenu("SingleLepton", triggerPropsCheckList_SingleLepton);
  theLooper.addHLTMenu("Dilepton_DF", triggerPropsCheckList_Dilepton_DF);
  theLooper.addHLTMenu("Dilepton_DF_Extra", triggerPropsCheckList_Dilepton_DF_Extra);
  theLooper.addHLTMenu("Dilepton_SF", triggerPropsCheckList_Dilepton_SF);
  theLooper.addHLTMenu("PFHT_Control", triggerPropsCheckList_PFHT_Control);
  theLooper.addHLTMenu("PFMET_MHT_Control", triggerPropsCheckList_PFMET_MHT_Control);
  // Set the MEs
  if (computeMEs) theLooper.setMatrixElementListFromFile(
    "${CMSSW_BASE}/src/CMS3/AnalysisTree/data/RecoProbabilities/RecoProbabilities.me",
    "AJetsVBFProbabilities_SpinZero_JHUGen,AJetsQCDProbabilities_SpinZero_JHUGen",
    //"AJetsVBFProbabilities_SpinZero_JHUGen,AJetsQCDProbabilities_SpinZero_JHUGen,AJetsVHProbabilities_SpinZero_JHUGen,PMAVJJ_SUPERDIJETMELA",
    false
  );

  curdir->cd();

  std::vector<BaseTree*> sample_trees; sample_trees.reserve(sampledirs.size());
  for (auto const& sname:sampledirs){
    BaseTree* sample_tree = new BaseTree(SampleHelpers::getDatasetFileName(sname), (useSkims ? "cms3ntuple/Dilepton" : "cms3ntuple/Events"), "", ""); sample_trees.push_back(sample_tree);
    sample_tree->sampleIdentifier = SampleHelpers::getSampleIdentifier(sname);

    std::vector<TString> allbranchnames; sample_tree->getValidBranchNamesWithoutAlias(allbranchnames, false);

    const int nEntries = sample_tree->getSelectedNEvents();
    double sum_wgts = (isData ? 1.f : 0.f);
    float xsec = 1;
    float xsec_scale = 1;
    if (!isData){
      // Get cross section
      sample_tree->bookBranch<float>("xsec", 0.f);
      sample_tree->getSelectedEvent(0);
      sample_tree->getVal("xsec", xsec);
      sample_tree->releaseBranch("xsec");
      xsec *= 1000.;

      // Reset gen. and lHE particle settings
      {
        bool has_lheMEweights = false;
        bool has_lheparticles = false;
        bool has_genparticles = false;
        for (auto const& bname:allbranchnames){
          if (bname.Contains("p_Gen") || bname.Contains("LHECandMass")) has_lheMEweights=true;
          else if (bname.Contains(GenInfoHandler::colName_lheparticles)) has_lheparticles = true;
          else if (bname.Contains(GenInfoHandler::colName_genparticles)) has_genparticles = true;
        }
        genInfoHandler.setAcquireLHEMEWeights(has_lheMEweights);
        genInfoHandler.setAcquireLHEParticles(has_lheparticles);
        genInfoHandler.setAcquireGenParticles(has_genparticles);
      }

      // Book branches
      simEventHandler.bookBranches(sample_tree);
      genInfoHandler.bookBranches(sample_tree);

      // Get sum of weights
      std::vector<TString> inputfilenames = SampleHelpers::getDatasetFileNames(sname);
      double sum_wgts_raw_noveto = 0;
      double sum_wgts_raw_withveto = 0;
      bool hasCounters = true;
      {
        int bin_syst = 1 + 1*(theGlobalSyst==SystematicsHelpers::ePUDn) + 2*(theGlobalSyst==SystematicsHelpers::ePUUp);
        int bin_period = 1;
        for (unsigned int iperiod=0; iperiod<nValidDataPeriods; iperiod++){
          if (validDataPeriods.at(iperiod)==SampleHelpers::theDataPeriod){ bin_period += iperiod+1; break; }
        }
        MELAout << "Checking counters histogram bin (" << bin_syst << ", " << bin_period << ") to obtain the sum of weights if the counters histogram exists..." << endl;
        for (auto const& fname:inputfilenames){
          TFile* ftmp = TFile::Open(fname, "read");
          TH2D* hCounters = (TH2D*) ftmp->Get("cms3ntuple/Counters");
          if (!hCounters){
            hasCounters = false;
            sum_wgts = 0;
            break;
          }
          MELAout << "\t- Successfully found the counters histogram in " << fname << endl;
          sum_wgts += hCounters->GetBinContent(bin_syst, bin_period);
          sum_wgts_raw_withveto += hCounters->GetBinContent(0, 0);
          sum_wgts_raw_noveto += hCounters->GetBinContent(0, 0) / (1. - hCounters->GetBinContent(0, 1));
          ftmp->Close();
        }
        if (hasCounters) MELAout << "\t- Obtained the weights from " << inputfilenames.size() << " files..." << endl;
      }
      if (!hasCounters && useSkims){
        MELAerr << "Skims should have contained counters histograms!" << endl;
        assert(0);
      }
      if (!hasCounters){
        MELAout << "No counters histograms are found. Initiation loop over " << nEntries << " events to determine the sample normalization:" << endl;

        simEventHandler.wrapTree(sample_tree);
        genInfoHandler.wrapTree(sample_tree);

        unsigned int n_zero_genwgts=0;
        double frac_zero_genwgts=0;
        for (int ev=0; ev<nEntries; ev++){
          HelperFunctions::progressbar(ev, nEntries);
          sample_tree->getEvent(ev);

          genInfoHandler.constructGenInfo(SystematicsHelpers::sNominal); // Use sNominal here in order to get the weight that corresponds to xsec
          auto const& genInfo = genInfoHandler.getGenInfo();
          double genwgt = genInfo->getGenWeight(true);
          if (genwgt==0.){
            n_zero_genwgts++;
            continue;
          }

          simEventHandler.constructSimEvent(SystematicsHelpers::sNominal);

          sum_wgts_raw_withveto += genwgt;
          sum_wgts += genwgt * simEventHandler.getPileUpWeight();
        }
        if (nEntries>0) frac_zero_genwgts = double(n_zero_genwgts)/double(nEntries);
        sum_wgts_raw_noveto = sum_wgts_raw_withveto / (1. - frac_zero_genwgts);
      }
      xsec_scale = sum_wgts_raw_withveto / sum_wgts_raw_noveto;
    }
    double globalWeight = xsec * xsec_scale * (isData ? 1.f : lumi) / sum_wgts;
    MELAout << "Sample " << sample_tree->sampleIdentifier << " has a gen. weight sum of " << sum_wgts << "." << endl;
    MELAout << "\t- xsec scale = " << xsec_scale << endl;
    MELAout << "\t- Global weight = " << globalWeight << endl;

    // Configure handlers
    pfcandidateHandler.bookBranches(sample_tree);
    muonHandler.bookBranches(sample_tree);
    electronHandler.bookBranches(sample_tree);
    photonHandler.bookBranches(sample_tree);
    jetHandler.bookBranches(sample_tree);
    isotrackHandler.bookBranches(sample_tree);
    vertexHandler.bookBranches(sample_tree);
    eventFilter.bookBranches(sample_tree);

    /*
    bool hasSuperclusters = false;
    for (auto const& bname:allbranchnames){
      if (bname.BeginsWith(SuperclusterHandler::colName.data())){
        hasSuperclusters = true;
        break;
      }
    }
    if (hasSuperclusters) superclusterHandler.bookBranches(sample_tree);
    */
    //fsrHandler.bookBranches(sample_tree);

    if (useJetOverlapStripping){
      overlapMap_muons_ak4jets.bookBranches(sample_tree);
      overlapMap_muons_ak8jets.bookBranches(sample_tree);
      overlapMap_electrons_ak4jets.bookBranches(sample_tree);
      overlapMap_electrons_ak8jets.bookBranches(sample_tree);
      overlapMap_photons_ak4jets.bookBranches(sample_tree);
      overlapMap_photons_ak8jets.bookBranches(sample_tree);
    }

    sample_tree->silenceUnused();

    // Add the input tree to the looper
    theLooper.addTree(sample_tree, globalWeight);
  }

  // Loop over all events
  theLooper.loop(true);

  // No need for the inputs
  for (auto& ss:sample_trees) delete ss;

  // Write output
  foutput->cd();
  tout->writeToFile(foutput);
  delete tout;
  foutput->Close();

  curdir->cd();

  SampleHelpers::addToCondorTransferList(stroutput);
}
