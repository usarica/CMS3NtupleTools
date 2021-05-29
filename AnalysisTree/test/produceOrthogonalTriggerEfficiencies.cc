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
#define CONTROL_TRIGGER_COMMAND(TYPE) \
  auto it_HLTMenuProps_##TYPE = triggerPropsCheckListMap.find(#TYPE); \
  if (it_HLTMenuProps_##TYPE == triggerPropsCheckListMap.cend()){ \
    MELAerr << "LooperFunctionHelpers::looperRule: The trigger type '" << #TYPE << "' has to be defined in this looper rule!" << endl; \
    assert(0); \
  }
  CONTROL_TRIGGER_COMMANDS;
#undef CONTROL_TRIGGER_COMMAND

  // Trigger leg matches
  std::vector< std::unordered_map<TString, std::vector< std::pair<TriggerHelpers::TriggerType, HLTTriggerPathProperties const*> > >::const_iterator > interestingHLTPropCheckLists;
#define CONTROL_TRIGGER_COMMAND(TYPE) interestingHLTPropCheckLists.push_back(it_HLTMenuProps_##TYPE);
  CONTROL_TRIGGER_COMMANDS;
#undef CONTROL_TRIGGER_COMMAND
  std::unordered_map<HLTTriggerPathProperties const*, std::pair<bool, float>> trigger_validity_wgt_map;
  for (auto const& interestingHLTPropCheckList:interestingHLTPropCheckLists){
    for (auto const& enumType_props_pair:interestingHLTPropCheckList->second){
      assert(enumType_props_pair.second != nullptr);
      auto const* hltprop = enumType_props_pair.second;
      trigger_validity_wgt_map[hltprop] = std::pair<bool, float>(false, 0);
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
  BRANCH_COMMAND(float, event_wgt_triggers) \
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
  BRANCH_COMMAND(float, genmet_pTmiss) \
  BRANCH_COMMAND(float, genmet_phimiss) \
  BRANCH_COMMAND(unsigned int, event_n_vtxs_good) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt30) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt30_btagged_loose) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt30_btagged_medium) \
  BRANCH_COMMAND(float, ak4jets_HT) \
  BRANCH_COMMAND(float, ak4jets_MHT) \
  BRANCH_COMMAND(float, dPhi_pTboson_pTmiss) \
  BRANCH_COMMAND(float, dPhi_pTbosonjets_pTmiss) \
  BRANCH_COMMAND(float, min_abs_dPhi_pTj_pTmiss) \
  BRANCH_COMMAND(cms3_id_t, lepton_id) \
  BRANCH_COMMAND(bool, lepton_is_genMatched_prompt) \
  BRANCH_COMMAND(float, lepton_pt) \
  BRANCH_COMMAND(float, lepton_eta) \
  BRANCH_COMMAND(float, lepton_etaSC) \
  BRANCH_COMMAND(float, lepton_phi) \
  BRANCH_COMMAND(float, lepton_mass)
#define BRANCH_VECTOR_COMMANDS

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
  unsigned int n_leptons_veto = 0;

  float SF_muons = 1;
  for (auto const& part:muons){
    float theSF = 1;
    if (!isData){
      muonSFHandler->getIdIsoSFAndEff(theGlobalSyst, part, theSF, nullptr);
      if (theSF == 0.f) continue;
      SF_muons *= theSF;
    }
    if (ParticleSelectionHelpers::isTightParticle(part)) leptons_tight.push_back(part);
    else if (ParticleSelectionHelpers::isVetoParticle(part)) n_leptons_veto++;
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
    else if (ParticleSelectionHelpers::isVetoParticle(part)) n_leptons_veto++;
  }
  event_wgt_SFs_electrons = SF_electrons;

  if (leptons_tight.size()!=1 || n_leptons_veto>0) return false;
  theLooper->incrementSelection("Nleps==1 and no other veto leptons");

  auto const& photons = photonHandler->getProducts();
  float SF_photons = 1;
  for (auto const& part:photons){
    float theSF = 1;
    if (!isData){
      photonSFHandler->getIdIsoSFAndEff(theGlobalSyst, part, theSF, nullptr);
      if (theSF == 0.f) continue;
      SF_photons *= theSF;
    }
    if (ParticleSelectionHelpers::isVetoParticle(part)) return false;
  }
  event_wgt_SFs_photons = SF_photons;

  // Isotracks
  isotrackHandler->constructIsotracks(&muons, &electrons);
  auto const& isotracks = isotrackHandler->getProducts();
  // Do not clean from isotracks just yet
  for (auto const& isotrack:isotracks){
    if (isotrack->testSelectionBit(IsotrackSelectionHelpers::kPreselectionVeto)) return false;
  }
  theLooper->incrementSelection("Pass global isotrack veto (no event skip)");

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
  min_abs_dPhi_pTj_pTmiss = TMath::Pi();
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

    // Determine min_abs_dPhi_pTj_pTmiss
    float dphi_tmp;
    HelperFunctions::deltaPhi(float(jet->phi()), event_phimiss, dphi_tmp); dphi_tmp = std::abs(dphi_tmp);
    min_abs_dPhi_pTj_pTmiss = std::min(min_abs_dPhi_pTj_pTmiss, dphi_tmp);
  }
  ak4jets_MHT = sump4_ak4jets.Pt();
  event_wgt_SFs_PUJetId = SF_PUJetId;
  event_wgt_SFs_btagging = SF_btagging;
  event_n_ak4jets_pt30 = ak4jets_tight.size();

  bool passOrthogonalControlTriggers = false;
#define CONTROL_TRIGGER_COMMAND(TYPE) \
  event_wgt_triggers_##TYPE = eventFilter->getTriggerWeight( \
    it_HLTMenuProps_##TYPE->second, \
    nullptr, nullptr, nullptr, &ak4jets, &ak8jets, eventmet \
  ); \
  passOrthogonalControlTriggers |= (event_wgt_triggers_##TYPE != 0.f);
  CONTROL_TRIGGER_COMMANDS;
#undef CONTROL_TRIGGER_COMMAND

  event_wgt_triggers = eventFilter->getTriggerWeight(
    it_HLTMenuProps_SingleLepton->second,
    &muons, &electrons, nullptr, nullptr, nullptr, nullptr
  );
  if (event_wgt_triggers==0.) return false;
  theLooper->incrementSelection("Pass single lepton triggers");

  // Record validity
  for (auto& it:trigger_validity_wgt_map){
    auto const& hltprop = it.first;
    for (auto const& hltpath:hltpaths){
      if (hltprop->isSameTrigger(hltpath->name)){
        it.second = std::pair<bool, float>(hltpath->isValid(), (hltpath->passTrigger ? hltpath->L1prescale * hltpath->HLTprescale : 0.));
        break;
      }
    }
  }

  // Fill leptons
  ParticleObject* theChosenLepton = leptons_tight.front();
  MuonObject* theMuon = dynamic_cast<MuonObject*>(theChosenLepton);
  ElectronObject* theElectron = dynamic_cast<ElectronObject*>(theChosenLepton);
  lepton_id = theChosenLepton->pdgId();
  lepton_is_genMatched_prompt = (theMuon ? theMuon->extras.is_genMatched_prompt : theElectron->extras.is_genMatched_prompt);
  lepton_pt = theChosenLepton->pt();
  lepton_eta = theChosenLepton->eta();
  lepton_etaSC = (theElectron ? theElectron->etaSC() : theMuon->eta());
  lepton_phi = theChosenLepton->phi();
  lepton_mass = theChosenLepton->mass();

  // Compute dPhi between the dilepton and pTmiss vector
  dPhi_pTboson_pTmiss = theChosenLepton->deltaPhi(event_phimiss);
  HelperFunctions::deltaPhi(float((theChosenLepton->p4()+sump4_ak4jets).Phi()), event_phimiss, dPhi_pTbosonjets_pTmiss);

  // Set the collection of SFs at the last step
  event_wgt_SFs = SF_muons*SF_electrons*SF_photons*SF_PUJetId*SF_btagging;

  /*********************/
  /* RECORD THE OUTPUT */
  /*********************/
#define BRANCH_COMMAND(TYPE, NAME) commonEntry.setNamedVal(#NAME, NAME);
  BRANCH_COMMANDS;
#undef BRANCH_COMMAND

  for (auto const& it:trigger_validity_wgt_map){
    TString bname = "is_valid_trigger_"; bname = bname + it.first->getName().data();
    commonEntry.setNamedVal(bname, it.second.first);
    bname = "event_wgt_trigger_"; bname = bname + it.first->getName().data();
    commonEntry.setNamedVal(bname, it.second.second);
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
  bool use_MET_XYCorr=true, bool use_MET_JERCorr=false, bool use_MET_ParticleMomCorr=true, bool use_MET_p4Preservation=true, bool use_MET_corrections=true
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
    //TriggerHelpers::kSingleMu, TriggerHelpers::kSingleMu_HighPt,
    TriggerHelpers::kSingleEle, TriggerHelpers::kSingleEle_HighPt
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
    "output/OrthogonalTriggerEfficiencies/SkimTrees/" + strdate
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
#define CONTROL_TRIGGER_COMMAND(TYPE) theLooper.addHLTMenu(#TYPE, triggerPropsCheckList_##TYPE);
  CONTROL_TRIGGER_COMMANDS;
#undef CONTROL_TRIGGER_COMMAND

  curdir->cd();

  BaseTree::setRobustInputCheck(false);
  std::vector<BaseTree*> sample_trees; sample_trees.reserve(sampledirs.size());
  for (auto const& sname:sampledirs){
    TString strdsetfname = SampleHelpers::getDatasetFileName(sname);
    MELAout << "=> Accessing the input trees from " << strdsetfname << "..." << endl;
    BaseTree* sample_tree = new BaseTree(strdsetfname, (useSkims ? "cms3ntuple/SingleLepton" : "cms3ntuple/Events"), "", ""); sample_trees.push_back(sample_tree);
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
