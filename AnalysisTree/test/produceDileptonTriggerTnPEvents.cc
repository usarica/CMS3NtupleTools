#include <cassert>
#include "common_includes.h"
#include "OffshellCutflow.h"
#include <CMS3/MELAHelpers/interface/CMS3MELAHelpers.h>
#include "TStyle.h"


#define CONTROL_TRIGGER_COMMANDS \
CONTROL_TRIGGER_COMMAND(AK8PFJet_Control) \
CONTROL_TRIGGER_COMMAND(VBFJets_Control) \
CONTROL_TRIGGER_COMMAND(PFHT_Control) \
CONTROL_TRIGGER_COMMAND(MET_Control) \
CONTROL_TRIGGER_COMMAND(PFMET_Control) \
CONTROL_TRIGGER_COMMAND(PFHT_PFMET_Control) \
CONTROL_TRIGGER_COMMAND(PFMET_MHT_Control) \
CONTROL_TRIGGER_COMMAND(PFHT_PFMET_MHT_Control)


namespace LooperFunctionHelpers{
  using namespace std;
  using namespace MELAStreamHelpers;
  using namespace OffshellCutflow;

  bool looperRule(BaseTreeLooper*, std::unordered_map<SystematicsHelpers::SystematicVariationTypes, double> const&, SimpleEntry&);


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
bool LooperFunctionHelpers::looperRule(BaseTreeLooper* theLooper, std::unordered_map<SystematicsHelpers::SystematicVariationTypes, double> const& extWgt, SimpleEntry& commonEntry){
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
  bool const& isGJets_HT = theLooper->getCurrentTreeFlag_GJetsHTException();
  float pTG_true_exception_range[2]={ -1, -1 };
  bool hasPTGExceptionRange = theLooper->getPTGExceptionRange(pTG_true_exception_range[0], pTG_true_exception_range[1]);
  bool needGenParticleChecks = isQCD || isGJets_HT || hasPTGExceptionRange;

  // Acquire triggers
  auto const& triggerCheckListMap = theLooper->getHLTMenus();
  auto const& triggerPropsCheckListMap = theLooper->getHLTMenuProperties();
  bool hasSimpleHLTMenus = theLooper->hasSimpleHLTMenus();
  bool hasHLTMenuProperties = theLooper->hasHLTMenuProperties();
  if (hasSimpleHLTMenus || !hasHLTMenuProperties){
    MELAerr << "LooperFunctionHelpers::looperRule: There must be HLT menus with properties." << endl;
    assert(0);
  }
  auto it_HLTMenuProps_SingleLepton = triggerPropsCheckListMap.find("SingleLepton");
  if (it_HLTMenuProps_SingleLepton == triggerPropsCheckListMap.cend()){
    MELAerr << "LooperFunctionHelpers::looperRule: The trigger type 'SingleLepton' has to be defined in this looper rule!" << endl;
    assert(0);
  }
  auto it_HLTMenuProps_Dilepton_DF = triggerPropsCheckListMap.find("Dilepton_DF");
  if (it_HLTMenuProps_Dilepton_DF == triggerPropsCheckListMap.cend()){
    MELAerr << "LooperFunctionHelpers::looperRule: The trigger type 'Dilepton_DF' has to be defined in this looper rule!" << endl;
    assert(0);
  }
  auto it_HLTMenuProps_Dilepton_DF_Extra = triggerPropsCheckListMap.find("Dilepton_DF_Extra");
  if (it_HLTMenuProps_Dilepton_DF_Extra == triggerPropsCheckListMap.cend()){
    MELAerr << "LooperFunctionHelpers::looperRule: The trigger type 'Dilepton_DF_Extra' has to be defined in this looper rule!" << endl;
    assert(0);
  }
  auto it_HLTMenuProps_Dilepton_SF = triggerPropsCheckListMap.find("Dilepton_SF");
  if (it_HLTMenuProps_Dilepton_SF == triggerPropsCheckListMap.cend()){
    MELAerr << "LooperFunctionHelpers::looperRule: The trigger type 'Dilepton_SF' has to be defined in this looper rule!" << endl;
    assert(0);
  }
#define CONTROL_TRIGGER_COMMAND(TYPE) \
  auto it_HLTMenuProps_##TYPE = triggerPropsCheckListMap.find(#TYPE); \
  if (it_HLTMenuProps_##TYPE == triggerPropsCheckListMap.cend()){ \
    MELAerr << "LooperFunctionHelpers::looperRule: The trigger type '" << #TYPE << "' has to be defined in this looper rule!" << endl; \
    assert(0); \
  }
  CONTROL_TRIGGER_COMMANDS;
#undef CONTROL_TRIGGER_COMMAND

  // Trigger leg matches
  std::vector< std::unordered_map<TString, std::vector< std::pair<TriggerHelpers::TriggerType, HLTTriggerPathProperties const*> > >::const_iterator > interestingHLTPropCheckLists{
    it_HLTMenuProps_SingleLepton,
    it_HLTMenuProps_Dilepton_DF, it_HLTMenuProps_Dilepton_DF_Extra,
    it_HLTMenuProps_Dilepton_SF
  };
  std::vector< std::pair<TriggerHelpers::TriggerType, HLTTriggerPathProperties const*> > interestingHLTTypePropPairs;
  std::vector< std::pair<TriggerHelpers::TriggerType, HLTTriggerPathProperties const*> > interestingHLTTypePropPairs_SingleLepton;
  std::unordered_map<HLTTriggerPathProperties const*, bool> trigger_validity_map;
  std::unordered_map<HLTTriggerPathProperties const*, std::vector<float>> dileptons_wgt_triggers_map;
  std::unordered_map<HLTTriggerPathProperties const*, std::vector<float>> leptons_wgt_triggers_SingleLepton_map;
  std::unordered_map<HLTTriggerPathProperties const*, std::vector<float>> tags_wgt_triggers_SingleLepton_map;
  for (auto const& interestingHLTPropCheckList:interestingHLTPropCheckLists){
    for (auto const& enumType_props_pair:interestingHLTPropCheckList->second){
      assert(enumType_props_pair.second != nullptr);
      auto const* hltprop = enumType_props_pair.second;
      trigger_validity_map[hltprop] = false;
      interestingHLTTypePropPairs.push_back(enumType_props_pair);
      dileptons_wgt_triggers_map[hltprop] = std::vector<float>();
      if (enumType_props_pair.first>=TriggerHelpers::kSingleMu && enumType_props_pair.first<TriggerHelpers::kSinglePho){
        interestingHLTTypePropPairs_SingleLepton.push_back(enumType_props_pair);
        leptons_wgt_triggers_SingleLepton_map[hltprop] = std::vector<float>();
        tags_wgt_triggers_SingleLepton_map[hltprop] = std::vector<float>();
      }
    }
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
  BRANCH_COMMAND(float, event_wgt_SFs) \
  BRANCH_COMMAND(float, event_wgt_SFs_muons) \
  BRANCH_COMMAND(float, event_wgt_SFs_electrons) \
  BRANCH_COMMAND(float, event_wgt_SFs_photons) \
  BRANCH_COMMAND(float, event_wgt_SFs_PUJetId) \
  BRANCH_COMMAND(float, event_wgt_SFs_btagging) \
  BRANCH_COMMAND(float, event_pTmiss) \
  BRANCH_COMMAND(float, event_phimiss) \
  BRANCH_COMMAND(bool, event_pass_tightMETFilters) \
  BRANCH_COMMAND(bool, event_pass_isotrackVeto) \
  BRANCH_COMMAND(bool, event_pass_ZWVeto) \
  BRANCH_COMMAND(float, genmet_pTmiss) \
  BRANCH_COMMAND(float, genmet_phimiss) \
  BRANCH_COMMAND(unsigned int, event_n_vtxs_good) \
  BRANCH_COMMAND(unsigned int, event_n_leptons_tight) \
  BRANCH_COMMAND(unsigned int, event_n_photons_veto) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt30) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt30_btagged_loose) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt30_btagged_medium) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt20) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt20_btagged_loose) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt20_btagged_medium) \
  BRANCH_COMMAND(float, ak4jets_HT) \
  BRANCH_COMMAND(float, ak4jets_MHT)
#define BRANCH_VECTOR_COMMANDS \
  BRANCH_COMMAND(cms3_id_t, dileptons_id) \
  BRANCH_COMMAND(bool, dileptons_pass_isotrackVetoWithTags) \
  BRANCH_COMMAND(float, dileptons_pt) \
  BRANCH_COMMAND(float, dileptons_eta) \
  BRANCH_COMMAND(float, dileptons_phi) \
  BRANCH_COMMAND(float, dileptons_mass) \
  BRANCH_COMMAND(float, dileptons_daughters_dR) \
  BRANCH_COMMAND(float, dileptons_dPhi_pTbosonjets_pTmiss) \
  BRANCH_COMMAND(float, dileptons_dPhi_pTboson_pTmiss) \
  BRANCH_COMMAND(std::vector<cms3_listSize_t>, dileptons_daughter_indices) \
  BRANCH_COMMAND(std::vector<cms3_listSize_t>, dileptons_tag_indices) \
  BRANCH_COMMAND(cms3_id_t, leptons_id) \
  BRANCH_COMMAND(bool, leptons_is_genMatched_prompt) \
  BRANCH_COMMAND(float, leptons_pt) \
  BRANCH_COMMAND(float, leptons_eta) \
  BRANCH_COMMAND(float, leptons_etaSC) \
  BRANCH_COMMAND(float, leptons_phi) \
  BRANCH_COMMAND(float, leptons_mass) \
  BRANCH_COMMAND(cms3_id_t, tags_id) \
  BRANCH_COMMAND(bool, tags_is_tight) \
  BRANCH_COMMAND(float, tags_pt) \
  BRANCH_COMMAND(float, tags_eta) \
  BRANCH_COMMAND(float, tags_phi) \
  BRANCH_COMMAND(float, tags_mass)
#define BRANCH_COMMANDS \
  BRANCH_SCALAR_COMMANDS \
  BRANCH_VECTOR_COMMANDS

#define CONTROL_TRIGGER_COMMAND(TYPE) float event_wgt_triggers_##TYPE = 0;
  CONTROL_TRIGGER_COMMANDS;
#undef CONTROL_TRIGGER_COMMAND
#define BRANCH_COMMAND(TYPE, NAME) TYPE NAME = 0;
  BRANCH_SCALAR_COMMANDS;
#undef BRANCH_COMMAND
#define BRANCH_COMMAND(TYPE, NAME) std::vector<TYPE> NAME;
  BRANCH_VECTOR_COMMANDS;
#undef BRANCH_COMMAND


  // Always assign the external weight first
  auto it_extWgt = extWgt.find(theGlobalSyst);
  if (it_extWgt==extWgt.cend()) it_extWgt = extWgt.find(SystematicsHelpers::nSystematicVariations);
  if (it_extWgt==extWgt.cend()){
    MELAerr << "LooperFunctionHelpers::looperRule: External normalization map does not have a proper weight assigned!" << endl;
    assert(0);
  }
  double const& extWgt_central = it_extWgt->second;

  event_wgt = extWgt_central;
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
      if (isGJets_HT){
        for (auto const& part:genparticles){
          if (PDGHelpers::isAPhoton(part->pdgId()) && part->extras.isPromptFinalState){
            double wgt_gjets = std::max(1., 1.71691-0.001221*part->pt());;
            event_wgt *= wgt_gjets;
            break;
          }
        }
      }
      if (hasPTGExceptionRange){
        for (auto const& part:genparticles){
          if (PDGHelpers::isAPhoton(part->pdgId()) && part->extras.isHardProcess){
            if ((pTG_true_exception_range[0]>=0.f && part->pt()<pTG_true_exception_range[0]) || (pTG_true_exception_range[1]>=0.f && part->pt()>=pTG_true_exception_range[1])){
              event_wgt = 0;
            }
            break;
          }
        }
      }
    }

    simEventHandler->constructSimEvent();
    event_wgt *= simEventHandler->getPileUpWeight(theGlobalSyst)*simEventHandler->getL1PrefiringWeight(theGlobalSyst);

    if (event_wgt==0.f) return false;

    // Record LHE MEs and K factors
    for (auto const& it:genInfo->extras.LHE_ME_weights) commonEntry.setNamedVal(it.first, it.second);
    for (auto const& it:genInfo->extras.Kfactors) commonEntry.setNamedVal(it.first, it.second);
  }
  theLooper->incrementSelection("Has valid gen. weights");

  vertexHandler->constructVertices();
  if (!vertexHandler->hasGoodPrimaryVertex()) return false;
  event_n_vtxs_good = vertexHandler->getNGoodVertices();
  theLooper->incrementSelection("Has a good PV");

  eventFilter->constructFilters(simEventHandler);
  auto const& hltpaths = eventFilter->getHLTPaths();
  if (isData && !eventFilter->isUniqueDataEvent()) return false;
  theLooper->incrementSelection("Unique event");

  if (!eventFilter->passCommonSkim() || !eventFilter->passMETFilters(EventFilterHandler::kMETFilters_Standard)) return false;
  theLooper->incrementSelection("Pass MET filters");
  event_pass_tightMETFilters = eventFilter->passMETFilters(EventFilterHandler::kMETFilters_Tight);

  pfcandidateHandler->constructPFCandidates(theGlobalSyst);
  auto const& pfcandidates = pfcandidateHandler->getProducts();

  muonHandler->constructMuons(theGlobalSyst, &pfcandidates);
  electronHandler->constructElectrons(theGlobalSyst, &pfcandidates);
  photonHandler->constructPhotons(theGlobalSyst, &pfcandidates);
  particleDisambiguator.disambiguateParticles(muonHandler, electronHandler, photonHandler);

  auto const& muons = muonHandler->getProducts();
  auto const& electrons = electronHandler->getProducts();
  std::vector<ParticleObject*> leptons_tight; leptons_tight.reserve(muons.size() + electrons.size());
  std::vector<ParticleObject*> tag_candidates;

  float SF_muons = 1;
  for (auto const& part:muons){
    float theSF = 1;
    if (!isData){
      muonSFHandler->getIdIsoSFAndEff(theGlobalSyst, part, theSF, nullptr);
      if (theSF == 0.f) continue;
      SF_muons *= theSF;
    }
    if (ParticleSelectionHelpers::isTightParticle(part)) leptons_tight.push_back(part);
    if (
      (part->testSelectionBit(MuonSelectionHelpers::kProbeId) && part->testSelectionBit(MuonSelectionHelpers::bit_preselectionTight_kin))
      ||
      ParticleSelectionHelpers::isTightParticle(part)
      ) tag_candidates.push_back(part);
  }
  event_wgt_SFs_muons = SF_muons;

  float SF_electrons = 1;
  for (auto const& part:electrons){
    float theSF = 1;
    if (!isData){
      electronSFHandler->getIdIsoSFAndEff(theGlobalSyst, part, theSF, nullptr);
      if (theSF == 0.f) continue;
      SF_electrons *= theSF;
    }
    if (ParticleSelectionHelpers::isTightParticle(part)) leptons_tight.push_back(part);
    if (
      (part->testSelectionBit(ElectronSelectionHelpers::kProbeId) && part->testSelectionBit(ElectronSelectionHelpers::bit_preselectionTight_kin))
      ||
      ParticleSelectionHelpers::isTightParticle(part)
      ) tag_candidates.push_back(part);
  }
  event_wgt_SFs_electrons = SF_electrons;

  event_n_leptons_tight = leptons_tight.size();
  if (event_n_leptons_tight>3 || event_n_leptons_tight<2) return false;
  theLooper->incrementSelection("2<=Nleps<=3");

  auto const& photons = photonHandler->getProducts();
  float SF_photons = 1;
  for (auto const& part:photons){
    float theSF = 1;
    if (!isData){
      photonSFHandler->getIdIsoSFAndEff(theGlobalSyst, part, theSF, nullptr);
      if (theSF == 0.f) continue;
      SF_photons *= theSF;
    }
    if (ParticleSelectionHelpers::isVetoParticle(part)) event_n_photons_veto++;
    if (
      part->testSelectionBit(PhotonSelectionHelpers::bit_preselectionTight_id) && part->testSelectionBit(PhotonSelectionHelpers::bit_preselectionTight_kin)
      &&
      part->testSelectionBit(PhotonSelectionHelpers::kBeamHaloSafe) && part->testSelectionBit(PhotonSelectionHelpers::kSpikeSafe)
      ) tag_candidates.push_back(part);
  }
  event_wgt_SFs_photons = SF_photons;

  // Isotracks
  isotrackHandler->constructIsotracks(&muons, &electrons);
  auto const& isotracks = isotrackHandler->getProducts();
  // Do not clean from isotracks just yet
  event_pass_isotrackVeto = true;
  for (auto const& isotrack:isotracks){
    if (isotrack->testSelectionBit(IsotrackSelectionHelpers::kPreselectionVeto)){
      event_pass_isotrackVeto = false;
      tag_candidates.push_back(isotrack);
    }
  }
  theLooper->incrementSelection("Pass global isotrack veto (no event skip)", event_pass_isotrackVeto);

  // Jets and MET
  jetHandler->constructJetMET(simEventHandler, theGlobalSyst, &muons, &electrons, &photons, &pfcandidates);
  auto const& ak4jets = jetHandler->getAK4Jets();
  auto const& ak8jets = jetHandler->getAK8Jets();

  // Test HEM filter for jets
  if (!eventFilter->test2018HEMFilter(simEventHandler, nullptr, nullptr, &ak4jets, &ak8jets)) return false;
  if (!eventFilter->testNoisyJetFilter(simEventHandler, ak4jets)) return false;
  theLooper->incrementSelection("HEM15/16 and noisy jet vetos");

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

  // Accumulate all tight jets and related quantities
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
      tag_candidates.push_back(jet);
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

  // Accumulate the valid dileptons
  dileptonHandler.constructDileptons(&muons, &electrons);
  auto const& dileptons = dileptonHandler.getProducts();
  std::vector<DileptonObject*> theChosenDileptons; theChosenDileptons.reserve(dileptons.size());
  DileptonObject* dilepton_bestZ = nullptr;
  float min_mll_nosign = -1;
  // Variables to track some counts before rejecting the event.
  unsigned int n_dileptons = dileptons.size();
  unsigned int n_dileptons_valid = 0;
  unsigned int n_dileptons_valid_tight = 0;
  unsigned int n_dileptons_valid_tight_SS = 0;
  unsigned int n_dileptons_valid_tight_OS = 0;
  unsigned int n_dileptons_valid_tight_OSSF = 0;
  unsigned int n_dileptons_valid_tight_OSSF_MW = 0;
  unsigned int n_dileptons_valid_tight_OSDF = 0;
  unsigned int n_dileptons_valid_tight_OSDF_MW = 0;
  for (auto const& dilepton:dileptons){
    if (dilepton->isValid()){
      n_dileptons_valid++;
      if (dilepton->nTightDaughters()==2){
        n_dileptons_valid_tight++;
        float tmp_mll = dilepton->mass();
        if (min_mll_nosign<0.f || min_mll_nosign>tmp_mll) min_mll_nosign = tmp_mll;

        if (dilepton->isOS()){
          n_dileptons_valid_tight_OS++;
          if (dilepton->isSF()){
            n_dileptons_valid_tight_OSSF++;
          }
          else{
            n_dileptons_valid_tight_OSDF++;
          }
        }
        else{
          n_dileptons_valid_tight_SS++;
        }
      }
    }
    if (dilepton->isValid() && dilepton->isOS() && dilepton->nTightDaughters()==2){
      float tmp_mll = dilepton->mass();
      if (dilepton->isSF()){
        //if (std::abs(tmp_mll - 91.2f)>=30.f) continue;
        if (tmp_mll<61.2f) continue;
        theChosenDileptons.push_back(dilepton);
        n_dileptons_valid_tight_OSSF_MW++;
        if (!dilepton_bestZ || std::abs(dilepton_bestZ->mass() - 91.2f)>std::abs(tmp_mll - 91.2f)) dilepton_bestZ = dilepton;
      }
      else{
        if (tmp_mll<76.2) continue;
        theChosenDileptons.push_back(dilepton);
        n_dileptons_valid_tight_OSDF_MW++;
      }
    }
  }
  theLooper->incrementSelection("\t- Total dilepton pairs", n_dileptons);
  theLooper->incrementSelection("\t- Total valid dilepton pairs", n_dileptons_valid);
  theLooper->incrementSelection("\t- Total valid tight dilepton pairs", n_dileptons_valid_tight);
  theLooper->incrementSelection("\t- Total valid tight SS dilepton pairs", n_dileptons_valid_tight_SS);
  theLooper->incrementSelection("\t- Total valid tight OS dilepton pairs", n_dileptons_valid_tight_OS);
  theLooper->incrementSelection("\t- Total valid tight OSSF dilepton pairs", n_dileptons_valid_tight_OSSF);
  theLooper->incrementSelection("\t- Total valid tight OSSF dilepton pairs within the mass window", n_dileptons_valid_tight_OSSF_MW);
  theLooper->incrementSelection("\t- Total valid tight OSDF dilepton pairs", n_dileptons_valid_tight_OSDF);
  theLooper->incrementSelection("\t- Total valid tight OSDF dilepton pairs within the mass window", n_dileptons_valid_tight_OSDF_MW);
  if (theChosenDileptons.empty()) return false;
  theLooper->incrementSelection("Has tight dilepton OS pairs with mass window reqs.");

  // Check main ZW veto
  event_pass_ZWVeto = !(event_n_leptons_tight==3 && event_pass_isotrackVeto && dilepton_bestZ && std::abs(dilepton_bestZ->mass() - 91.2f)<30.f && min_mll_nosign>4.f && event_pTmiss>=20.f && event_n_ak4jets_pt30_btagged_loose==0);
  theLooper->incrementSelection("Has 3 tight dileptons", event_n_leptons_tight==3);
  theLooper->incrementSelection("Has 3 tight dileptons passing ZW vetoes", event_n_leptons_tight==3 && event_pass_ZWVeto);

  bool passOrthogonalControlTriggers = false;
#define CONTROL_TRIGGER_COMMAND(TYPE) \
  event_wgt_triggers_##TYPE = eventFilter->getTriggerWeight( \
    it_HLTMenuProps_##TYPE->second, \
    nullptr, nullptr, nullptr, &ak4jets, &ak8jets, eventmet \
  ); \
  passOrthogonalControlTriggers |= (event_wgt_triggers_##TYPE != 0.f);
  CONTROL_TRIGGER_COMMANDS;
#undef CONTROL_TRIGGER_COMMAND

  // Record validity
  for (auto& it:trigger_validity_map){
    auto const& hltprop = it.first;
    for (auto const& hltpath:hltpaths){
      if (hltprop->isSameTrigger(hltpath->name)){
        it.second = hltpath->isValid();
        break;
      }
    }
  }
  // Acquire valid single lepton paths
  std::unordered_map<HLTTriggerPathProperties const*, HLTTriggerPathObject*> hltprop_hltpath_map_SingleLepton;
  for (auto const& enumType_props_pair:interestingHLTTypePropPairs_SingleLepton){
    auto const& trigger_type = enumType_props_pair.first;
    auto const& hltprop = enumType_props_pair.second;
    for (auto const& hltpath:hltpaths){
      if (hltpath->isValid() && hltprop->isSameTrigger(hltpath->name)){
        hltprop_hltpath_map_SingleLepton[hltprop] = hltpath;
        break;
      }
    }
  }

  bool passSingleLeptonTriggers = false;
  std::vector<ParticleObject*> tags_selected;
  for (auto const& part:tag_candidates){
    std::unordered_map<HLTTriggerPathProperties const*, float> hltprop_wgt_map;
    for (auto const& enumType_props_pair:interestingHLTTypePropPairs_SingleLepton){
      auto const& trigger_type = enumType_props_pair.first;
      auto const& hltprop = enumType_props_pair.second;

      HLTTriggerPathObject* hltpath = nullptr;
      std::unordered_map<HLTTriggerPathProperties const*, HLTTriggerPathObject*>::const_iterator it_hltprop_hltpath_map_SingleLepton;
      float wgt_trigger = 0;
      bool doTest = (
        !(
          (std::abs(part->pdgId())==13 && !(trigger_type>=TriggerHelpers::kSingleMu && trigger_type<TriggerHelpers::kSingleEle))
          ||
          ((std::abs(part->pdgId())==11 || std::abs(part->pdgId())==22) && !(trigger_type>=TriggerHelpers::kSingleEle && trigger_type<TriggerHelpers::kSinglePho))
          )
        &&
        HelperFunctions::getUnorderedMapIterator(hltprop, hltprop_hltpath_map_SingleLepton, it_hltprop_hltpath_map_SingleLepton)
        &&
        (hltpath=it_hltprop_hltpath_map_SingleLepton->second) && hltpath->passTrigger
        );
      if (doTest){
        HLTTriggerPathProperties::TriggerObjectExceptionType const& TOexception = hltprop->getTOException();
        auto const& passedTriggerObjects = hltpath->getPassedTriggerObjects();

        bool hasTORecovery = false;
        std::vector<TriggerObject const*> passedTriggerObjectsWithRecovery;
        if (TOexception == HLTTriggerPathProperties::toRecoverObjectsFromFailing){
          // Determine what to recover
          unsigned short n_TOmuons = 0;
          unsigned short n_TOelectrons = 0;
          unsigned short n_TOphotons = 0;
          for (auto const& TOobj:passedTriggerObjects){
            if (TOobj->isTriggerObjectType(trigger::TriggerMuon)) n_TOmuons++;
            else if (TOobj->isTriggerObjectType(trigger::TriggerElectron)) n_TOelectrons++;
            else if (TOobj->isTriggerObjectType(trigger::TriggerPhoton) || TOobj->isTriggerObjectType(trigger::TriggerCluster)){
              n_TOelectrons++;
              n_TOphotons++;
            }
          }
          auto const& props_map = hltprop->getObjectProperties();
          unsigned short n_TOmuons_req = (props_map.find(HLTObjectProperties::kMuon)!=props_map.end() ? props_map.find(HLTObjectProperties::kMuon)->second.size() : 0);
          unsigned short n_TOelectrons_req = (props_map.find(HLTObjectProperties::kElectron)!=props_map.end() ? props_map.find(HLTObjectProperties::kElectron)->second.size() : 0);
          unsigned short n_TOphotons_req = (props_map.find(HLTObjectProperties::kPhoton)!=props_map.end() ? props_map.find(HLTObjectProperties::kPhoton)->second.size() : 0);
          bool needMuonRecovery = (n_TOmuons<n_TOmuons_req);
          bool needElectronRecovery = (n_TOelectrons<n_TOelectrons_req);
          bool needPhotonRecovery = (n_TOphotons<n_TOphotons_req);
          // Add existing passing objects
          for (auto const& TOobj:passedTriggerObjects) passedTriggerObjectsWithRecovery.push_back(TOobj);
          // Add from failing objects
          for (auto const& TOobj:hltpath->getFailedTriggerObjects()){
            if (
              (needMuonRecovery && TOobj->isTriggerObjectType(trigger::TriggerMuon))
              ||
              (needElectronRecovery && (TOobj->isTriggerObjectType(trigger::TriggerElectron) || TOobj->isTriggerObjectType(trigger::TriggerPhoton) || TOobj->isTriggerObjectType(trigger::TriggerCluster)))
              ||
              (needPhotonRecovery && (TOobj->isTriggerObjectType(trigger::TriggerPhoton) || TOobj->isTriggerObjectType(trigger::TriggerCluster)))
              ) passedTriggerObjectsWithRecovery.push_back(TOobj);
          }
          hasTORecovery = true;
        }

        std::vector<TriggerObject const*> const& passedTriggerObjects_final = (!hasTORecovery ? passedTriggerObjects : passedTriggerObjectsWithRecovery);
        std::vector<ParticleObject*> tmpvec{ part };
        std::vector<ParticleObject*> tmpvec_matched;

        if (std::abs(part->pdgId())==13) TriggerObject::getMatchedPhysicsObjects(
          passedTriggerObjects_final, { trigger::TriggerMuon }, 0.2,
          tmpvec, tmpvec_matched
        );
        else if (std::abs(part->pdgId())==11 || std::abs(part->pdgId())==22) TriggerObject::getMatchedPhysicsObjects(
          passedTriggerObjects_final, { trigger::TriggerElectron, trigger::TriggerPhoton, trigger::TriggerCluster }, 0.2,
          tmpvec, tmpvec_matched
        );
        else TriggerObject::getMatchedPhysicsObjects(
          passedTriggerObjects_final, { trigger::TriggerMuon, trigger::TriggerElectron, trigger::TriggerPhoton, trigger::TriggerCluster }, 0.2,
          tmpvec, tmpvec_matched
        );

        if (!tmpvec_matched.empty()) wgt_trigger = hltpath->L1prescale * hltpath->HLTprescale;
      }
      if (wgt_trigger!=0.f) hltprop_wgt_map[hltprop] = wgt_trigger;
    }
    if (!hltprop_wgt_map.empty()){
      tags_selected.push_back(part);
      for (auto& it:tags_wgt_triggers_SingleLepton_map){
        auto const& hltprop = it.first;
        auto& vec_wgts = it.second;
        std::unordered_map<HLTTriggerPathProperties const*, float>::const_iterator it_hltprop_wgt_map;
        vec_wgts.push_back((HelperFunctions::getUnorderedMapIterator(hltprop, hltprop_wgt_map, it_hltprop_wgt_map) ? it_hltprop_wgt_map->second : 0.f));
        passSingleLeptonTriggers |= (vec_wgts.back() != 0.f);
      }
    }
    // All tight leptons are tag candidates, and the order of tight leptons within the tag candidates collection is the same.
    if (HelperFunctions::checkListVariable(leptons_tight, part)){
      for (auto& it:leptons_wgt_triggers_SingleLepton_map){
        auto const& hltprop = it.first;
        auto& vec_wgts = it.second;
        std::unordered_map<HLTTriggerPathProperties const*, float>::const_iterator it_hltprop_wgt_map;
        vec_wgts.push_back((HelperFunctions::getUnorderedMapIterator(hltprop, hltprop_wgt_map, it_hltprop_wgt_map) ? it_hltprop_wgt_map->second : 0.f));
        passSingleLeptonTriggers |= (vec_wgts.back() != 0.f);
      }
    }
  }
  // Must pass at least one single lepton or orthogonal trigger
  if (!passSingleLeptonTriggers && !passOrthogonalControlTriggers) return false;
  theLooper->incrementSelection("Pass single lepton or control triggers");
  theLooper->incrementSelection("\t- Pass single lepton triggers", passSingleLeptonTriggers);
  theLooper->incrementSelection("\t- Pass control triggers", passOrthogonalControlTriggers);

  // Fill leptons
  for (auto const& part:leptons_tight){
    MuonObject* theMuon = dynamic_cast<MuonObject*>(part);
    ElectronObject* theElectron = dynamic_cast<ElectronObject*>(part);
    bool is_genMatched_prompt = (theMuon ? theMuon->extras.is_genMatched_prompt : theElectron->extras.is_genMatched_prompt);
    float etaSC = (theElectron ? theElectron->etaSC() : theMuon->eta());

    leptons_id.push_back(part->pdgId());
    leptons_is_genMatched_prompt.push_back(part->pdgId());
    leptons_pt.push_back(part->pt());
    leptons_eta.push_back(part->eta());
    leptons_etaSC.push_back(etaSC);
    leptons_phi.push_back(part->phi());
    leptons_mass.push_back(part->mass());
  }

  // Fill selected tags
  for (auto const& part:tags_selected){
    MuonObject* theMuon = dynamic_cast<MuonObject*>(part);
    ElectronObject* theElectron = dynamic_cast<ElectronObject*>(part);
    PhotonObject* thePhoton = dynamic_cast<PhotonObject*>(part);
    IsotrackObject* theIsotrack = dynamic_cast<IsotrackObject*>(part);
    AK4JetObject* theAK4Jet = dynamic_cast<AK4JetObject*>(part);

    cms3_id_t tag_id = -9000;
    bool tag_is_tight = true;
    if (theMuon || theElectron || thePhoton){
      tag_id = std::abs(part->pdgId());
      tag_is_tight = ParticleSelectionHelpers::isTightParticle(part);
    }
    else if (theIsotrack) tag_id = -std::abs(part->pdgId());
    else if (theAK4Jet) tag_id = 1;
    else MELAerr << "Tag type is not defined!" << endl;

    tags_id.push_back(tag_id);
    tags_is_tight.push_back(tag_is_tight);
    tags_pt.push_back(part->pt());
    tags_eta.push_back(part->eta());
    tags_phi.push_back(part->phi());
    tags_mass.push_back(part->mass());
  }

  for (auto const& theChosenDilepton:theChosenDileptons){
    auto const& daughters = theChosenDilepton->getDaughters();
    ParticleObject* dau_first = theChosenDilepton->getDaughter_leadingPt();
    ParticleObject* dau_second = theChosenDilepton->getDaughter_subleadingPt();
    if (std::abs(dau_first->pdgId())==11 && std::abs(dau_second->pdgId())==13) std::swap(dau_first, dau_second);
    std::vector<MuonObject*> daughters_muon;
    std::vector<ElectronObject*> daughters_electron;
    for (auto const& dau:daughters){
      MuonObject* dau_muon = dynamic_cast<MuonObject*>(dau);
      ElectronObject* dau_electron = dynamic_cast<ElectronObject*>(dau);
      if (dau_muon) daughters_muon.push_back(dau_muon);
      else if (dau_electron) daughters_electron.push_back(dau_electron);
    }

    cms3_id_t const dilepton_id = daughters.front()->pdgId()*daughters.back()->pdgId();
    bool const dilepton_is_SF = dilepton_id!=-143;
    float const dPhi_pTboson_pTmiss = theChosenDilepton->deltaPhi(event_phimiss);
    float dPhi_pTbosonjets_pTmiss = TMath::Pi();
    HelperFunctions::deltaPhi(float((theChosenDilepton->p4()+sump4_ak4jets).Phi()), event_phimiss, dPhi_pTbosonjets_pTmiss);

    theLooper->incrementSelection("\t- Total number of dilepton pairs selected");
    theLooper->incrementSelection("\t- Total number of dilepton pairs selected (ee)", dilepton_id==-121);
    theLooper->incrementSelection("\t- Total number of dilepton pairs selected (mue)", dilepton_id==-143);
    theLooper->incrementSelection("\t- Total number of dilepton pairs selected (mumu)", dilepton_id==-169);

    // Apply veto from SR-like selection
    OffshellCutflow::setActiveFinalState(OffshellCutflow::fs_ZZ_2l2nu);
    bool passZZ2l2nuSRlikeSelection = (
      event_n_leptons_tight==2
      &&
      OffshellCutflow::check_pTboson(theChosenDilepton->pt())
      &&
      OffshellCutflow::check_pTmiss(event_pTmiss)
      &&
      OffshellCutflow::check_dPhi_pTll_pTmiss(dPhi_pTboson_pTmiss)
      &&
      OffshellCutflow::check_dPhi_pTlljets_pTmiss(dPhi_pTbosonjets_pTmiss)
      );
    OffshellCutflow::setActiveFinalState(OffshellCutflow::fs_WW_2l2nu);
    bool passWW2l2nuSRlikeSelection = (
      event_n_leptons_tight==2
      &&
      OffshellCutflow::check_pTboson(theChosenDilepton->pt())
      &&
      OffshellCutflow::check_pTmiss(event_pTmiss)
      &&
      OffshellCutflow::check_dPhi_pTll_pTmiss(dPhi_pTboson_pTmiss)
      &&
      OffshellCutflow::check_dPhi_pTlljets_pTmiss(dPhi_pTbosonjets_pTmiss)
      );
    bool passBtaggedJets = (event_n_ak4jets_pt30_btagged_medium==0);
    bool passSRVeto = !((passZZ2l2nuSRlikeSelection || passWW2l2nuSRlikeSelection) && passBtaggedJets);
    theLooper->incrementSelection("\t- Number of dilepton pairs passing SR veto", passSRVeto);
    theLooper->incrementSelection("\t- Number of dilepton pairs passing SR veto (ee)", passSRVeto && dilepton_id==-121);
    theLooper->incrementSelection("\t- Number of dilepton pairs passing SR veto (mue)", passSRVeto && dilepton_id==-143);
    theLooper->incrementSelection("\t- Number of dilepton pairs passing SR veto (mumu)", passSRVeto && dilepton_id==-169);

    // Do the iteration twice to preserve order
    std::vector<cms3_listSize_t> daughter_indices;
    for (unsigned int ilep=0; ilep<leptons_tight.size(); ilep++){
      auto const& part = leptons_tight.at(ilep);
      if (part == dau_first){ daughter_indices.push_back(ilep); break; }
    }
    for (unsigned int ilep=0; ilep<leptons_tight.size(); ilep++){
      auto const& part = leptons_tight.at(ilep);
      if (part == dau_second){ daughter_indices.push_back(ilep); break; }
    }
    assert(daughter_indices.size()==2);

    std::vector<cms3_listSize_t> tag_indices;
    std::vector<ParticleObject*> tags_dilepton;
    for (unsigned int itag=0; itag<tags_selected.size(); itag++){
      auto const& part = tags_selected.at(itag);
      if (theChosenDilepton->hasDaughter(part)) continue;
      if (dau_first->deltaR(part)<0.4f || dau_second->deltaR(part)<0.4f) continue;
      tag_indices.push_back(itag);
      tags_dilepton.push_back(part);
    }
    theLooper->incrementSelection("\t- Dilepton pairs with a matched extra taggable object", !tag_indices.empty());
    theLooper->incrementSelection("\t- Dilepton pairs with a matched extra taggable object (ee)", !tag_indices.empty() && dilepton_id==-121);
    theLooper->incrementSelection("\t- Dilepton pairs with a matched extra taggable object (mue)", !tag_indices.empty() && dilepton_id==-143);
    theLooper->incrementSelection("\t- Dilepton pairs with a matched extra taggable object (mumu)", !tag_indices.empty() && dilepton_id==-169);

    theLooper->incrementSelection("\t- Dilepton pairs passing orthogonal trigger selections", passOrthogonalControlTriggers);
    theLooper->incrementSelection("\t- Dilepton pairs passing orthogonal trigger selections (ee)", passOrthogonalControlTriggers && dilepton_id==-121);
    theLooper->incrementSelection("\t- Dilepton pairs passing orthogonal trigger selections (mue)", passOrthogonalControlTriggers && dilepton_id==-143);
    theLooper->incrementSelection("\t- Dilepton pairs passing orthogonal trigger selections (mumu)", passOrthogonalControlTriggers && dilepton_id==-169);

    bool hasVetoIsotrack = false;
    for (auto const& isotrack:isotracks){
      if (isotrack->testSelectionBit(IsotrackSelectionHelpers::kPreselectionVeto)){
        bool doSkip = false;
        for (auto const& tag:tags_dilepton){
          if (tag==isotrack || tag->deltaR(isotrack)<0.4f){
            doSkip = true;
            break;
          }
        }
        if (doSkip) continue;
        hasVetoIsotrack = true;
        break;
      }
    }
    theLooper->incrementSelection("\t- Number of dilepton (+ taggable object) collections passing isotrack veto", !hasVetoIsotrack);

    // Apply all vetoes before filling the vectors
    if (!passSRVeto) continue;

    dileptons_id.push_back(dilepton_id);
    dileptons_pass_isotrackVetoWithTags.push_back(!hasVetoIsotrack);
    dileptons_pt.push_back(theChosenDilepton->pt());
    dileptons_eta.push_back(theChosenDilepton->eta());
    dileptons_phi.push_back(theChosenDilepton->phi());
    dileptons_mass.push_back(theChosenDilepton->m());
    dileptons_daughters_dR.push_back(dau_first->deltaR(dau_second));
    dileptons_dPhi_pTbosonjets_pTmiss.push_back(dPhi_pTbosonjets_pTmiss);
    dileptons_dPhi_pTboson_pTmiss.push_back(dPhi_pTboson_pTmiss);
    dileptons_daughter_indices.push_back(daughter_indices);
    dileptons_tag_indices.push_back(tag_indices);

    for (auto const& enumType_props_pair:interestingHLTTypePropPairs){
      auto const& trigger_type = enumType_props_pair.first;
      auto const& hltprop = enumType_props_pair.second;
      float tmp_wgt_trigger = eventFilter->getTriggerWeight(
        { enumType_props_pair },
        &daughters_muon, &daughters_electron, nullptr, &ak4jets, nullptr, nullptr
      );
      dileptons_wgt_triggers_map.find(hltprop)->second.push_back(tmp_wgt_trigger);
    }
  }
  if (dileptons_id.empty()) return false;
  theLooper->incrementSelection("Has valid tight dileptons");

  // Set the collection of SFs at the last step
  event_wgt_SFs = SF_muons*SF_electrons*SF_photons*SF_PUJetId*SF_btagging;

  /*********************/
  /* RECORD THE OUTPUT */
  /*********************/
#define BRANCH_COMMAND(TYPE, NAME) commonEntry.setNamedVal(#NAME, NAME);
  BRANCH_COMMANDS;
#undef BRANCH_COMMAND
#define CONTROL_TRIGGER_COMMAND(TYPE) commonEntry.setNamedVal(Form("event_wgt_triggers_%s", #TYPE), event_wgt_triggers_##TYPE);
  CONTROL_TRIGGER_COMMANDS;
#undef CONTROL_TRIGGER_COMMAND

  for (auto const& it:trigger_validity_map){
    TString bname = "is_valid_trigger_"; bname = bname + it.first->getName().data();
    commonEntry.setNamedVal(bname, it.second);
  }
  for (auto const& it:dileptons_wgt_triggers_map){
    TString bname = "dileptons_wgt_triggers_"; bname = bname + it.first->getName().data();
    commonEntry.setNamedVal(bname, it.second);
  }
  for (auto const& it:leptons_wgt_triggers_SingleLepton_map){
    TString bname = "leptons_wgt_triggers_"; bname = bname + it.first->getName().data();
    commonEntry.setNamedVal(bname, it.second);
  }
  for (auto const& it:tags_wgt_triggers_SingleLepton_map){
    TString bname = "tags_wgt_triggers_"; bname = bname + it.first->getName().data();
    commonEntry.setNamedVal(bname, it.second);
  }

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
  //useSkims = false;
  SampleHelpers::configure(period, Form("%s:%s", (useSkims ? "store_skims" : "store"), prodVersion.Data()));

  const float lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::getDataPeriod());

  std::vector<TString> const validDataPeriods = SampleHelpers::getValidDataPeriods();
  size_t const nValidDataPeriods = validDataPeriods.size();

  BtagHelpers::setBtagWPType(BtagHelpers::kDeepFlav_Loose);
  LooperFunctionHelpers::setBtagWPs();

  std::vector<TriggerHelpers::TriggerType> requiredTriggers_SingleLepton{
    TriggerHelpers::kSingleMu, TriggerHelpers::kSingleMu_HighPt, TriggerHelpers::kSingleMu_Control,
    TriggerHelpers::kSingleEle, TriggerHelpers::kSingleEle_HighPt, TriggerHelpers::kSingleEle_Control
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
#define CONTROL_TRIGGER_COMMAND(TYPE) std::vector<TriggerHelpers::TriggerType> requiredTriggers_##TYPE{ TriggerHelpers::k##TYPE };
  CONTROL_TRIGGER_COMMANDS;
#undef CONTROL_TRIGGER_COMMAND
  TriggerHelpers::dropSelectionCuts(TriggerHelpers::kMET_Control);
  TriggerHelpers::dropSelectionCuts(TriggerHelpers::kPFMET_Control);
  TriggerHelpers::dropSelectionCuts(TriggerHelpers::kPFHT_PFMET_Control);
  TriggerHelpers::dropSelectionCuts(TriggerHelpers::kPFMET_MHT_Control);
  TriggerHelpers::dropSelectionCuts(TriggerHelpers::kPFHT_PFMET_MHT_Control);
  auto triggerPropsCheckList_SingleLepton = TriggerHelpers::getHLTMenuProperties(requiredTriggers_SingleLepton);
  auto triggerPropsCheckList_Dilepton_DF = TriggerHelpers::getHLTMenuProperties(requiredTriggers_Dilepton_DF);
  auto triggerPropsCheckList_Dilepton_DF_Extra = TriggerHelpers::getHLTMenuProperties(requiredTriggers_Dilepton_DF_Extra);
  auto triggerPropsCheckList_Dilepton_SF = TriggerHelpers::getHLTMenuProperties(requiredTriggers_Dilepton_SF);
#define CONTROL_TRIGGER_COMMAND(TYPE) auto triggerPropsCheckList_##TYPE = TriggerHelpers::getHLTMenuProperties(requiredTriggers_##TYPE);
  CONTROL_TRIGGER_COMMANDS;
#undef CONTROL_TRIGGER_COMMAND

  // Get sample specifications
  std::vector<TString> sampledirs;
  SampleHelpers::constructSamplesList(strSampleSet, theGlobalSyst, sampledirs);
  if (sampledirs.empty()) return;
  bool isData = SampleHelpers::checkSampleIsData(sampledirs.front());
  if (isData && theGlobalSyst!=sNominal) return;

  // Set flags for ak4jet tight id
  AK4JetSelectionHelpers::setPUIdWP(applyPUIdToAK4Jets ? AK4JetSelectionHelpers::kTightPUJetId : AK4JetSelectionHelpers::nSelectionBits); // Default is 'tight'
  AK4JetSelectionHelpers::setApplyTightLeptonVetoIdToJets(applyTightLeptonVetoIdToAK4Jets); // Default is 'false'
  LooperFunctionHelpers::setAK4JetSelectionOptions(applyPUIdToAK4Jets, applyTightLeptonVetoIdToAK4Jets);

  // Set flags for MET
  LooperFunctionHelpers::setMETOptions(use_MET_Puppi, use_MET_XYCorr, use_MET_JERCorr, use_MET_ParticleMomCorr, use_MET_p4Preservation, use_MET_corrections);

  // Set output directory
  TString coutput_main =
    "output/DileptonTriggerTnPEvents/SkimTrees/" + strdate
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
  stroutput += Form("_%s", SystematicsHelpers::getSystName(theGlobalSyst).data());
  if (nchunks>0) stroutput = stroutput + Form("_%i_of_%i", ichunk, nchunks);
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
#define CONTROL_TRIGGER_COMMAND(TYPE) theLooper.addHLTMenu(#TYPE, triggerPropsCheckList_##TYPE);
  CONTROL_TRIGGER_COMMANDS;
#undef CONTROL_TRIGGER_COMMAND

  curdir->cd();

  std::vector<BaseTree*> sample_trees; sample_trees.reserve(sampledirs.size());
  for (auto const& sname:sampledirs){
    BaseTree* sample_tree = new BaseTree(SampleHelpers::getDatasetFileName(sname), (useSkims ? "cms3ntuple/Dilepton" : "cms3ntuple/Events"), "", ""); sample_trees.push_back(sample_tree);
    sample_tree->sampleIdentifier = SampleHelpers::getSampleIdentifier(sname);

    std::vector<TString> allbranchnames; sample_tree->getValidBranchNamesWithoutAlias(allbranchnames, false);

    const int nEntries = sample_tree->getSelectedNEvents();
    double sum_wgts = (isData ? 1.f : 0.f);
    double sum_wgts_PUDn = sum_wgts;
    double sum_wgts_PUUp = sum_wgts;
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
            sum_wgts = sum_wgts_PUDn = sum_wgts_PUUp = 0;
            break;
          }
          MELAout << "\t- Successfully found the counters histogram in " << fname << endl;
          sum_wgts += hCounters->GetBinContent(bin_syst, bin_period);
          sum_wgts_PUDn += hCounters->GetBinContent(2, bin_period);
          sum_wgts_PUUp += hCounters->GetBinContent(3, bin_period);
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

          simEventHandler.constructSimEvent();

          sum_wgts_raw_withveto += genwgt;
          sum_wgts += genwgt * simEventHandler.getPileUpWeight(theGlobalSyst);
          sum_wgts_PUDn += genwgt * simEventHandler.getPileUpWeight(SystematicsHelpers::ePUDn);
          sum_wgts_PUUp += genwgt * simEventHandler.getPileUpWeight(SystematicsHelpers::ePUUp);
        }
        if (nEntries>0) frac_zero_genwgts = double(n_zero_genwgts)/double(nEntries);
        sum_wgts_raw_noveto = sum_wgts_raw_withveto / (1. - frac_zero_genwgts);
      }
      xsec_scale = sum_wgts_raw_withveto / sum_wgts_raw_noveto;
    }
    std::unordered_map<SystematicsHelpers::SystematicVariationTypes, double> globalWeights;
    double globalWeight = xsec * xsec_scale * (isData ? 1.f : lumi) / sum_wgts; globalWeights[theGlobalSyst] = globalWeight;
    double globalWeight_PUDn = xsec * xsec_scale * (isData ? 1.f : lumi) / sum_wgts_PUDn; globalWeights[SystematicsHelpers::ePUDn] = globalWeight_PUDn;
    double globalWeight_PUUp = xsec * xsec_scale * (isData ? 1.f : lumi) / sum_wgts_PUUp; globalWeights[SystematicsHelpers::ePUUp] = globalWeight_PUUp;
    MELAout << "Sample " << sample_tree->sampleIdentifier << " has a gen. weight sum of " << sum_wgts << " (PU dn: " << sum_wgts_PUDn << ", PU up: " << sum_wgts_PUUp << ")." << endl;
    MELAout << "\t- xsec scale = " << xsec_scale << endl;
    MELAout << "\t- Global weight = " << globalWeight << endl;
    MELAout << "\t- Global weight (PU dn) = " << globalWeight_PUDn << endl;
    MELAout << "\t- Global weight (PU up) = " << globalWeight_PUUp << endl;

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
    theLooper.addTree(sample_tree, globalWeights);
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

  splitFileAndAddForTransfer(stroutput);
}


#undef CONTROL_TRIGGER_COMMANDS
