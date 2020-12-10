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

  // Fakeable ID usage
  bool applyFakeables = false;
  void setApplyFakeableId(bool applyFakeables_);

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
  auto const& triggerPropsCheckListMap = theLooper->getHLTMenuProperties();
  bool hasHLTMenuProperties = theLooper->hasHLTMenuProperties();
  if (!hasHLTMenuProperties){
    MELAerr << "LooperFunctionHelpers::looperRule: There must be HLT menus with properties." << endl;
    assert(0);
  }
  auto it_HLTMenuProps = triggerPropsCheckListMap.find("SingleLepton");
  if (it_HLTMenuProps == triggerPropsCheckListMap.cend()){
    MELAerr << "LooperFunctionHelpers::looperRule: The trigger type 'SingleLepton' has to be defined in this looper rule!" << endl;
    assert(0);
  }
#define CONTROL_TRIGGER_COMMAND(TYPE) auto it_HLTMenuProps_##TYPE = triggerPropsCheckListMap.find(#TYPE);
  CONTROL_TRIGGER_COMMANDS;
#undef CONTROL_TRIGGER_COMMAND
#define CONTROL_TRIGGER_COMMAND(TYPE) \
  if (it_HLTMenuProps_##TYPE == triggerPropsCheckListMap.cend()){ \
    MELAerr << "LooperFunctionHelpers::looperRule: The trigger type '" << #TYPE << "' has to be defined in this looper rule!" << endl; \
    assert(0); \
  }
  if (applyFakeables){
    CONTROL_TRIGGER_COMMANDS;
  }
#undef CONTROL_TRIGGER_COMMAND

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
  if (!isData){
    OBJECT_HANDLER_SIM_DIRECTIVES;
    SCALEFACTOR_HANDLER_DIRECTIVES;
  }
#undef HANDLER_DIRECTIVE

  /************************/
  /* EVENT INTERPRETATION */
  /************************/
  // Recorded variables
#define BRANCH_SCALAR_COMMANDS \
  BRANCH_COMMAND(float, event_wgt) \
  BRANCH_COMMAND(float, event_wgt_L1PrefiringDn) \
  BRANCH_COMMAND(float, event_wgt_L1PrefiringUp) \
  BRANCH_COMMAND(float, event_wgt_PUDn) \
  BRANCH_COMMAND(float, event_wgt_PUUp) \
  BRANCH_COMMAND(float, event_wgt_adjustment_NNPDF30) \
  BRANCH_COMMAND(float, event_wgt_triggers) \
  BRANCH_COMMAND(float, event_wgt_SFs_muons) \
  BRANCH_COMMAND(float, event_wgt_SFs_muons_StatDn) \
  BRANCH_COMMAND(float, event_wgt_SFs_muons_StatUp) \
  BRANCH_COMMAND(float, event_wgt_SFs_muons_SystDn) \
  BRANCH_COMMAND(float, event_wgt_SFs_muons_SystUp) \
  BRANCH_COMMAND(float, event_wgt_SFs_muons_AltMCDn) \
  BRANCH_COMMAND(float, event_wgt_SFs_muons_AltMCUp) \
  BRANCH_COMMAND(float, event_wgt_SFs_electrons) \
  BRANCH_COMMAND(float, event_wgt_SFs_electrons_StatDn) \
  BRANCH_COMMAND(float, event_wgt_SFs_electrons_StatUp) \
  BRANCH_COMMAND(float, event_wgt_SFs_electrons_SystDn) \
  BRANCH_COMMAND(float, event_wgt_SFs_electrons_SystUp) \
  BRANCH_COMMAND(float, event_wgt_SFs_electrons_AltMCDn) \
  BRANCH_COMMAND(float, event_wgt_SFs_electrons_AltMCUp) \
  BRANCH_COMMAND(float, event_wgt_SFs_photons) \
  BRANCH_COMMAND(float, event_wgt_SFs_photons_EffDn) \
  BRANCH_COMMAND(float, event_wgt_SFs_photons_EffUp) \
  BRANCH_COMMAND(float, event_wgt_SFs_PUJetId) \
  BRANCH_COMMAND(float, event_wgt_SFs_PUJetId_EffDn) \
  BRANCH_COMMAND(float, event_wgt_SFs_PUJetId_EffUp) \
  BRANCH_COMMAND(float, event_wgt_SFs_btagging) \
  BRANCH_COMMAND(float, event_wgt_SFs_btagging_EffDn) \
  BRANCH_COMMAND(float, event_wgt_SFs_btagging_EffUp) \
  BRANCH_COMMAND(float, event_pTmiss) \
  BRANCH_COMMAND(float, event_phimiss) \
  BRANCH_COMMAND(float, event_mTZZ) \
  BRANCH_COMMAND(float, event_mZZ) \
  BRANCH_COMMAND(bool, event_pass_tightMETFilters) \
  BRANCH_COMMAND(float, genmet_pTmiss) \
  BRANCH_COMMAND(float, genmet_phimiss) \
  BRANCH_COMMAND(unsigned int, event_n_vtxs_good) \
  BRANCH_COMMAND(unsigned int, event_n_leptons_fakeableBase) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt30) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt30_btagged_loose) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt30_btagged_medium) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt20) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt20_btagged_loose) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt20_btagged_medium) \
  BRANCH_COMMAND(float, ak4jets_HT) \
  BRANCH_COMMAND(float, ak4jets_MHT) \
  BRANCH_COMMAND(cms3_id_t, lepton_id) \
  BRANCH_COMMAND(float, lepton_pt) \
  BRANCH_COMMAND(float, lepton_eta) \
  BRANCH_COMMAND(float, lepton_phi) \
  BRANCH_COMMAND(float, lepton_mass) \
  BRANCH_COMMAND(float, lepton_relIso) \
  BRANCH_COMMAND(bool, lepton_is_genMatched_prompt) \
  BRANCH_COMMAND(float, electron_full5x5_sigmaIEtaIEta) \
  BRANCH_COMMAND(float, electron_full5x5_sigmaIPhiIPhi) \
  BRANCH_COMMAND(float, electron_full5x5_r9) \
  BRANCH_COMMAND(float, electron_seedTime) \
  BRANCH_COMMAND(float, dPhi_pTboson_pTmiss) \
  BRANCH_COMMAND(float, dPhi_pTbosonjets_pTmiss) \
  BRANCH_COMMAND(float, min_abs_dPhi_pTj_pTmiss)
#define BRANCH_VECTOR_COMMANDS \
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

#define CONTROL_TRIGGER_COMMAND(TYPE) float event_wgt_triggers_##TYPE = 0;
  CONTROL_TRIGGER_COMMANDS;
#undef CONTROL_TRIGGER_COMMAND
#define BRANCH_COMMAND(TYPE, NAME) TYPE NAME = 0;
  BRANCH_SCALAR_COMMANDS;
#undef BRANCH_COMMAND
#define BRANCH_COMMAND(TYPE, NAME) std::vector<TYPE> NAME;
  BRANCH_VECTOR_COMMANDS;
#undef BRANCH_COMMAND
  std::vector<TString> bnames_exclude;
#define BRANCH_COMMAND(TYPE, NAME) if (theGlobalSyst!=SystematicsHelpers::sNominal){ TString strtmp = #NAME; if (strtmp.EndsWith("Up") || strtmp.EndsWith("Dn")) bnames_exclude.emplace_back(strtmp); }
  BRANCH_SCALAR_COMMANDS;
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

  auto it_extWgt_PUDn = extWgt.find(SystematicsHelpers::ePUDn);
  if (it_extWgt_PUDn==extWgt.cend()) it_extWgt_PUDn = it_extWgt;
  double const& extWgt_PUDn = it_extWgt_PUDn->second;

  auto it_extWgt_PUUp = extWgt.find(SystematicsHelpers::ePUUp);
  if (it_extWgt_PUUp==extWgt.cend()) it_extWgt_PUUp = it_extWgt;
  double const& extWgt_PUUp = it_extWgt_PUUp->second;

  event_wgt = event_wgt_L1PrefiringDn = event_wgt_L1PrefiringUp = extWgt_central;
  event_wgt_PUDn = extWgt_PUDn; event_wgt_PUUp = extWgt_PUUp;
  // Set NNPDF 3.0 adjustment to 1
  event_wgt_adjustment_NNPDF30 = 1;

  if (!isData){
    genInfoHandler->constructGenInfo(theGlobalSyst);
    auto const& genInfo = genInfoHandler->getGenInfo();
    double genwgt_NNPDF30 = genInfo->getGenWeight(false);
    double genwgt_default = genInfo->getGenWeight(true);
    event_wgt_adjustment_NNPDF30 = (genwgt_default!=0. ? genwgt_NNPDF30 / genwgt_default : 0.);
    event_wgt *= genwgt_default;
    event_wgt_PUDn *= genwgt_default;
    event_wgt_PUUp *= genwgt_default;
    genmet_pTmiss = genInfo->extras.genmet_met;
    genmet_phimiss = genInfo->extras.genmet_metPhi;
    auto const& genparticles = genInfoHandler->getGenParticles();

    if (needGenParticleChecks){
      if (isQCD){
        for (auto const& part:genparticles){
          if (PDGHelpers::isAPhoton(part->pdgId()) && part->extras.isPromptFinalState && part->pt()>=25.f){
            event_wgt = event_wgt_PUDn = event_wgt_PUUp = 0;
            break;
          }
        }
      }
      if (isGJets_HT){
        for (auto const& part:genparticles){
          if (PDGHelpers::isAPhoton(part->pdgId()) && part->extras.isPromptFinalState){
            double wgt_gjets = std::max(1., 1.71691-0.001221*part->pt());;
            event_wgt *= wgt_gjets;
            event_wgt_PUDn *= wgt_gjets;
            event_wgt_PUUp *= wgt_gjets;
            break;
          }
        }
      }
      if (hasPTGExceptionRange){
        for (auto const& part:genparticles){
          if (PDGHelpers::isAPhoton(part->pdgId()) && part->extras.isHardProcess){
            if ((pTG_true_exception_range[0]>=0.f && part->pt()<pTG_true_exception_range[0]) || (pTG_true_exception_range[1]>=0.f && part->pt()>=pTG_true_exception_range[1])){
              event_wgt = event_wgt_PUDn = event_wgt_PUUp = 0;
            }
            break;
          }
        }
      }
    }
    event_wgt_L1PrefiringDn = event_wgt_L1PrefiringUp = event_wgt;

    simEventHandler->constructSimEvent();
    event_wgt *= simEventHandler->getPileUpWeight(theGlobalSyst)*simEventHandler->getL1PrefiringWeight(theGlobalSyst);
    event_wgt_L1PrefiringDn *= simEventHandler->getPileUpWeight(theGlobalSyst)*simEventHandler->getL1PrefiringWeight(SystematicsHelpers::eL1PrefiringDn);
    event_wgt_L1PrefiringUp *= simEventHandler->getPileUpWeight(theGlobalSyst)*simEventHandler->getL1PrefiringWeight(SystematicsHelpers::eL1PrefiringUp);
    event_wgt_PUDn *= simEventHandler->getPileUpWeight(SystematicsHelpers::ePUDn)*simEventHandler->getL1PrefiringWeight(theGlobalSyst);
    event_wgt_PUUp *= simEventHandler->getPileUpWeight(SystematicsHelpers::ePUUp)*simEventHandler->getL1PrefiringWeight(theGlobalSyst);

    if (
      event_wgt==0.f
      &&
      event_wgt_L1PrefiringDn==0.f
      &&
      event_wgt_L1PrefiringUp==0.f
      &&
      event_wgt_PUDn==0.f
      &&
      event_wgt_PUUp==0.f
      ) return false;

    // Record LHE MEs and K factors
    for (auto const& it:genInfo->extras.LHE_ME_weights) commonEntry.setNamedVal(it.first, it.second);
    for (auto const& it:genInfo->extras.Kfactors) commonEntry.setNamedVal(it.first, it.second);
  }
  theLooper->incrementSelection("Valid gen. weights");

  vertexHandler->constructVertices();
  if (!vertexHandler->hasGoodPrimaryVertex()) return false;
  event_n_vtxs_good = vertexHandler->getNGoodVertices();
  theLooper->incrementSelection("Good vertices");

  pfcandidateHandler->constructPFCandidates(theGlobalSyst);
  auto const& pfcandidates = pfcandidateHandler->getProducts();

  muonHandler->constructMuons(theGlobalSyst, &pfcandidates);
  electronHandler->constructElectrons(theGlobalSyst, &pfcandidates);
  photonHandler->constructPhotons(theGlobalSyst, &pfcandidates);
  particleDisambiguator.disambiguateParticles(muonHandler, electronHandler, photonHandler);

  ParticleObject* theChosenLepton = nullptr;
  unsigned short n_leptons_tight = 0;
  unsigned short n_leptons_veto = 0;

  auto const& muons = muonHandler->getProducts();
  float SF_muons = 1;
  float SF_muons_StatDn = 1;
  float SF_muons_StatUp = 1;
  float SF_muons_SystDn = 1;
  float SF_muons_SystUp = 1;
  float SF_muons_AltMCDn = 1;
  float SF_muons_AltMCUp = 1;
  for (auto const& part:muons){
    float theSF = 1;
    if (!isData){
      theSF = 1; muonSFHandler->getIdIsoSFAndEff(theGlobalSyst, part, theSF, nullptr); theSF = std::max(1e-5f, theSF); SF_muons *= theSF;
      theSF = 1; muonSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eMuEffStatDn, part, theSF, nullptr); theSF = std::max(1e-5f, theSF); SF_muons_StatDn *= theSF;
      theSF = 1; muonSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eMuEffStatUp, part, theSF, nullptr); theSF = std::max(1e-5f, theSF); SF_muons_StatUp *= theSF;
      theSF = 1; muonSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eMuEffSystDn, part, theSF, nullptr); theSF = std::max(1e-5f, theSF); SF_muons_SystDn *= theSF;
      theSF = 1; muonSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eMuEffSystUp, part, theSF, nullptr); theSF = std::max(1e-5f, theSF); SF_muons_SystUp *= theSF;
      theSF = 1; muonSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eMuEffAltMCDn, part, theSF, nullptr); theSF = std::max(1e-5f, theSF); SF_muons_AltMCDn *= theSF;
      theSF = 1; muonSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eMuEffAltMCUp, part, theSF, nullptr); theSF = std::max(1e-5f, theSF); SF_muons_AltMCUp *= theSF;
    }

    if ((!applyFakeables && ParticleSelectionHelpers::isTightParticle(part)) || (applyFakeables && ParticleSelectionHelpers::isLooseParticle(part))){
      theChosenLepton = part;
      n_leptons_tight++;
    }
    else if (ParticleSelectionHelpers::isVetoParticle(part)) n_leptons_veto++;
    else if (!ParticleSelectionHelpers::isTightParticle(part) && part->testSelectionBit(MuonSelectionHelpers::kFakeableBase)) event_n_leptons_fakeableBase++;
  }
  event_wgt_SFs_muons = SF_muons;
  event_wgt_SFs_muons_StatDn = SF_muons_StatDn;
  event_wgt_SFs_muons_StatUp = SF_muons_StatUp;
  event_wgt_SFs_muons_SystDn = SF_muons_SystDn;
  event_wgt_SFs_muons_SystUp = SF_muons_SystUp;
  event_wgt_SFs_muons_AltMCDn = SF_muons_AltMCDn;
  event_wgt_SFs_muons_AltMCUp = SF_muons_AltMCUp;

  auto const& electrons = electronHandler->getProducts();
  float SF_electrons = 1;
  float SF_electrons_StatDn = 1;
  float SF_electrons_StatUp = 1;
  float SF_electrons_SystDn = 1;
  float SF_electrons_SystUp = 1;
  float SF_electrons_AltMCDn = 1;
  float SF_electrons_AltMCUp = 1;
  for (auto const& part:electrons){
    float theSF = 1;
    if (!isData){
      theSF = 1; electronSFHandler->getIdIsoSFAndEff(theGlobalSyst, part, theSF, nullptr); theSF = std::max(1e-5f, theSF); SF_electrons *= theSF;
      theSF = 1; electronSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eEleEffStatDn, part, theSF, nullptr); theSF = std::max(1e-5f, theSF); SF_electrons_StatDn *= theSF;
      theSF = 1; electronSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eEleEffStatUp, part, theSF, nullptr); theSF = std::max(1e-5f, theSF); SF_electrons_StatUp *= theSF;
      theSF = 1; electronSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eEleEffSystDn, part, theSF, nullptr); theSF = std::max(1e-5f, theSF); SF_electrons_SystDn *= theSF;
      theSF = 1; electronSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eEleEffSystUp, part, theSF, nullptr); theSF = std::max(1e-5f, theSF); SF_electrons_SystUp *= theSF;
      theSF = 1; electronSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eEleEffAltMCDn, part, theSF, nullptr); theSF = std::max(1e-5f, theSF); SF_electrons_AltMCDn *= theSF;
      theSF = 1; electronSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eEleEffAltMCUp, part, theSF, nullptr); theSF = std::max(1e-5f, theSF); SF_electrons_AltMCUp *= theSF;
    }

    if ((!applyFakeables && ParticleSelectionHelpers::isTightParticle(part)) || (applyFakeables && ParticleSelectionHelpers::isLooseParticle(part))){
      theChosenLepton = part;
      n_leptons_tight++;
    }
    else if (ParticleSelectionHelpers::isVetoParticle(part)) n_leptons_veto++;
    else if (!ParticleSelectionHelpers::isTightParticle(part) && part->testSelectionBit(ElectronSelectionHelpers::kFakeableBase)) event_n_leptons_fakeableBase++;
  }
  event_wgt_SFs_electrons = SF_electrons;
  event_wgt_SFs_electrons_StatDn = SF_electrons_StatDn;
  event_wgt_SFs_electrons_StatUp = SF_electrons_StatUp;
  event_wgt_SFs_electrons_SystDn = SF_electrons_SystDn;
  event_wgt_SFs_electrons_SystUp = SF_electrons_SystUp;
  event_wgt_SFs_electrons_AltMCDn = SF_electrons_AltMCDn;
  event_wgt_SFs_electrons_AltMCUp = SF_electrons_AltMCUp;

  if (!(n_leptons_tight==1 && n_leptons_veto==0)) return false;
  theLooper->incrementSelection("Exactly one lepton");

  auto const& photons = photonHandler->getProducts();
  unsigned int n_photons_veto = 0;
  float SF_photons = 1;
  float SF_photons_EffDn = 1;
  float SF_photons_EffUp = 1;
  for (auto const& part:photons){
    if (!isData){
      float theSF = 1;
      photonSFHandler->getIdIsoSFAndEff(theGlobalSyst, part, theSF, nullptr); theSF = std::max(1e-5f, theSF); SF_photons *= theSF;
      theSF = 1; photonSFHandler->getIdIsoSFAndEff(SystematicsHelpers::ePhoEffDn, part, theSF, nullptr); theSF = std::max(1e-5f, theSF); SF_photons_EffDn *= theSF;
      theSF = 1; photonSFHandler->getIdIsoSFAndEff(SystematicsHelpers::ePhoEffUp, part, theSF, nullptr); theSF = std::max(1e-5f, theSF); SF_photons_EffUp *= theSF;
    }

    if (ParticleSelectionHelpers::isVetoParticle(part)) n_photons_veto++;
  }
  event_wgt_SFs_photons = SF_photons;
  event_wgt_SFs_photons_EffDn = SF_photons_EffDn;
  event_wgt_SFs_photons_EffUp = SF_photons_EffUp;
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

  // Record lepton variables
  MuonObject* theMuon = dynamic_cast<MuonObject*>(theChosenLepton);
  ElectronObject* theElectron = dynamic_cast<ElectronObject*>(theChosenLepton);
  lepton_id = theChosenLepton->pdgId();
  lepton_pt = theChosenLepton->pt();
  lepton_eta = theChosenLepton->eta();
  lepton_phi = theChosenLepton->phi();
  lepton_mass = theChosenLepton->m();
  lepton_is_genMatched_prompt = (theMuon ? theMuon->extras.is_genMatched_prompt : theElectron->extras.is_genMatched_prompt);
  if (theElectron){
    electron_full5x5_sigmaIEtaIEta = theElectron->extras.full5x5_sigmaIEtaIEta;
    electron_full5x5_sigmaIPhiIPhi = theElectron->extras.full5x5_sigmaIPhiIPhi;
    electron_full5x5_r9 = theElectron->extras.full5x5_r9;
    electron_seedTime = theElectron->extras.seedTime;
    lepton_relIso = ElectronSelectionHelpers::computeIso(*theElectron);
  }
  else{
    lepton_relIso = MuonSelectionHelpers::computeIso(*theMuon);
  }

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

  ParticleObject::LorentzVector_t sump4_ak4jets(0, 0, 0, 0);
  std::vector<AK4JetObject*> ak4jets_tight; ak4jets_tight.reserve(ak4jets.size());
  unsigned int n_ak4jets_tight_pt30_btagged_loose = 0;
  float SF_PUJetId = 1;
  float SF_PUJetId_EffDn = 1;
  float SF_PUJetId_EffUp = 1;
  float SF_btagging = 1;
  float SF_btagging_EffDn = 1;
  float SF_btagging_EffUp = 1;
  for (auto const& jet:ak4jets){
    if (!isData){
      float theSF_PUJetId = 1;
      float theSF_btag = 1;
      pujetidSFHandler->getSFAndEff(theGlobalSyst, jet, theSF_PUJetId, nullptr); theSF_PUJetId = std::max(theSF_PUJetId, 1e-5f); SF_PUJetId *= theSF_PUJetId;
      btagSFHandler->getSFAndEff(theGlobalSyst, jet, theSF_btag, nullptr); theSF_btag = std::max(theSF_btag, 1e-5f); SF_btagging *= theSF_btag;

      theSF_PUJetId = 1;
      theSF_btag = 1;
      pujetidSFHandler->getSFAndEff(SystematicsHelpers::ePUJetIdEffDn, jet, theSF_PUJetId, nullptr); theSF_PUJetId = std::max(theSF_PUJetId, 1e-5f); SF_PUJetId_EffDn *= theSF_PUJetId;
      btagSFHandler->getSFAndEff(SystematicsHelpers::eBTagSFDn, jet, theSF_btag, nullptr); theSF_btag = std::max(theSF_btag, 1e-5f); SF_btagging_EffDn *= theSF_btag;

      theSF_PUJetId = 1;
      theSF_btag = 1;
      pujetidSFHandler->getSFAndEff(SystematicsHelpers::ePUJetIdEffUp, jet, theSF_PUJetId, nullptr); theSF_PUJetId = std::max(theSF_PUJetId, 1e-5f); SF_PUJetId_EffUp *= theSF_PUJetId;
      btagSFHandler->getSFAndEff(SystematicsHelpers::eBTagSFUp, jet, theSF_btag, nullptr); theSF_btag = std::max(theSF_btag, 1e-5f); SF_btagging_EffUp *= theSF_btag;
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
  event_wgt_SFs_PUJetId = std::min(SF_PUJetId, 3.f);
  event_wgt_SFs_PUJetId_EffDn = std::min(SF_PUJetId_EffDn, 3.f);
  event_wgt_SFs_PUJetId_EffUp = std::min(SF_PUJetId_EffUp, 3.f);
  event_wgt_SFs_btagging = SF_btagging;
  event_wgt_SFs_btagging_EffDn = SF_btagging_EffDn;
  event_wgt_SFs_btagging_EffUp = SF_btagging_EffUp;
  event_n_ak4jets_pt30 = ak4jets_tight.size();
  ak4jets_MHT = sump4_ak4jets.Pt();

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
  dPhi_pTboson_pTmiss = theChosenLepton->deltaPhi(event_phimiss);
  HelperFunctions::deltaPhi(float((theChosenLepton->p4()+sump4_ak4jets).Phi()), event_phimiss, dPhi_pTbosonjets_pTmiss);

  // Compute mass variables
  float const& etamiss_approx = lepton_eta;
  ParticleObject::LorentzVector_t p4_lepton_Z_approx; p4_lepton_Z_approx = ParticleObject::PolarLorentzVector_t(lepton_pt, lepton_eta, lepton_phi, PDGHelpers::Zmass);
  ParticleObject::LorentzVector_t p4_ZZ_approx; p4_ZZ_approx = ParticleObject::PolarLorentzVector_t(event_pTmiss, etamiss_approx, event_phimiss, PDGHelpers::Zmass);
  p4_ZZ_approx = p4_ZZ_approx + p4_lepton_Z_approx;

  event_mTZZ = std::sqrt(
    std::pow(
    (
      std::sqrt(std::pow(lepton_pt, 2) + std::pow(p4_lepton_Z_approx.M(), 2))
      + std::sqrt(std::pow(event_pTmiss, 2) + std::pow(PDGHelpers::Zmass, 2))
      ), 2
    )
    - std::pow((p4_lepton_Z_approx + event_met_p4).Pt(), 2)
  );
  event_mZZ = p4_ZZ_approx.M();

  // Apply veto from SR-like selection
  OffshellCutflow::setActiveFinalState(OffshellCutflow::fs_ZZ_2l2nu);
  std::vector<bool> const v_passZZ2l2nuSRlikeSelection={
    OffshellCutflow::check_pTboson(lepton_pt),
    OffshellCutflow::check_pTmiss(event_pTmiss),
    OffshellCutflow::check_dPhi_pTll_pTmiss(dPhi_pTboson_pTmiss),
    OffshellCutflow::check_dPhi_pTlljets_pTmiss(dPhi_pTbosonjets_pTmiss),
    OffshellCutflow::check_Nb_veto(event_n_ak4jets_pt30_btagged_loose)
  };
  OffshellCutflow::setActiveFinalState(OffshellCutflow::fs_WW_2l2nu);
  std::vector<bool> const v_passWW2l2nuSRlikeSelection={
    OffshellCutflow::check_pTboson(lepton_pt),
    OffshellCutflow::check_pTmiss(event_pTmiss),
    OffshellCutflow::check_dPhi_pTll_pTmiss(dPhi_pTboson_pTmiss),
    OffshellCutflow::check_dPhi_pTlljets_pTmiss(dPhi_pTbosonjets_pTmiss),
    OffshellCutflow::check_Nb_veto(event_n_ak4jets_pt30_btagged_loose)
  };
  unsigned short n_passZZ2l2nuSRlikeSelection = 0;
  unsigned short n_passWW2l2nuSRlikeSelection = 0;
  for (auto const& selreq:v_passZZ2l2nuSRlikeSelection){ if (selreq) n_passZZ2l2nuSRlikeSelection++; }
  for (auto const& selreq:v_passWW2l2nuSRlikeSelection){ if (selreq) n_passWW2l2nuSRlikeSelection++; }
  bool const pass_SRSel_Nminus1 = (n_passZZ2l2nuSRlikeSelection >= v_passZZ2l2nuSRlikeSelection.size()-1) || (n_passWW2l2nuSRlikeSelection >= v_passWW2l2nuSRlikeSelection.size()-1);

  // Do not record sigle lepton CR events which do not pass N-1 selection requirements
  if (!applyFakeables){
    if (!pass_SRSel_Nminus1) return false;
    theLooper->incrementSelection("N-1 veto");
  }

  // Accumulate ak8 jets as well
  for (auto const& jet:ak8jets){
    if (!ParticleSelectionHelpers::isTightJet(jet)) continue;
    if (jet->mass()<60.f) continue;

    ak8jets_pt.push_back(jet->pt());
    ak8jets_eta.push_back(jet->eta());
    ak8jets_phi.push_back(jet->phi());
    ak8jets_mass.push_back(jet->mass());
  }

  eventFilter->constructFilters(simEventHandler);
  if (isData && !eventFilter->isUniqueDataEvent()) return false;
  theLooper->incrementSelection("Unique data");

  if (!eventFilter->passCommonSkim() || !eventFilter->passMETFilters(EventFilterHandler::kMETFilters_Standard)) return false;
  theLooper->incrementSelection("MET filters");
  event_pass_tightMETFilters = eventFilter->passMETFilters(EventFilterHandler::kMETFilters_Tight);

  // Test HEM filter
  if (!eventFilter->test2018HEMFilter(simEventHandler, nullptr, nullptr, &ak4jets, &ak8jets)) return false;
  theLooper->incrementSelection("HEM15/16 veto");

  if (!applyFakeables) event_wgt_triggers = eventFilter->getTriggerWeight(
    it_HLTMenuProps->second,
    &muons, &electrons, nullptr, nullptr, nullptr, nullptr
  );
  else{
#define CONTROL_TRIGGER_COMMAND(TYPE) \
    event_wgt_triggers_##TYPE = eventFilter->getTriggerWeight( \
      it_HLTMenuProps_##TYPE->second, \
      nullptr, nullptr, nullptr, &ak4jets, &ak8jets, eventmet \
    ); \
    if (event_wgt_triggers_##TYPE != 0.f) event_wgt_triggers = (event_wgt_triggers==0.f ? event_wgt_triggers_##TYPE : std::min(event_wgt_triggers_##TYPE, event_wgt_triggers));
    CONTROL_TRIGGER_COMMANDS;
#undef CONTROL_TRIGGER_COMMAND

    // Check matching of single lepton triggers to jets
    float event_wgt_triggers_SingleLepton_jetMatched = 0;
    auto const& hltpaths = eventFilter->getHLTPaths();
    std::vector<ParticleObject*> triggerjet_candidates; // In order to match to single lepton triggers manually
    for (auto const& jet:ak4jets_tight) triggerjet_candidates.push_back(jet);
    for (auto const& enumType_props_pair:it_HLTMenuProps->second){
      if (event_wgt_triggers_SingleLepton_jetMatched==1.f) break;

      auto const& trigger_type = enumType_props_pair.first;
      auto const& hltprop = enumType_props_pair.second;
      for (auto const& hltpath:hltpaths){
        if (event_wgt_triggers_SingleLepton_jetMatched==1.f) break;

        if (hltpath->isValid() && hltpath->passTrigger && hltprop->isSameTrigger(hltpath->name)){
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
          std::vector<ParticleObject*> tmpvec_matched;
          TriggerObject::getMatchedPhysicsObjects(
            passedTriggerObjects_final, { trigger::TriggerMuon, trigger::TriggerElectron, trigger::TriggerPhoton, trigger::TriggerCluster }, 0.2,
            triggerjet_candidates, tmpvec_matched
          );
          if (!tmpvec_matched.empty()){
            float wgt_trigger = hltpath->L1prescale * hltpath->HLTprescale;
            if (wgt_trigger!=0.f) event_wgt_triggers_SingleLepton_jetMatched = (event_wgt_triggers_SingleLepton_jetMatched==0.f ? wgt_trigger : std::min(wgt_trigger, event_wgt_triggers_SingleLepton_jetMatched));
          }
        }
      }
    }
    commonEntry.setNamedVal("event_wgt_triggers_SingleLepton_jetMatched", event_wgt_triggers_SingleLepton_jetMatched);
    event_wgt_triggers = (event_wgt_triggers==0.f ? event_wgt_triggers_SingleLepton_jetMatched : std::min(event_wgt_triggers_SingleLepton_jetMatched, event_wgt_triggers));
  }
  if (event_wgt_triggers == 0.f) return false;
  theLooper->incrementSelection("Trigger");

  // Compute MEs
  bool computeMEs = theLooper->hasRecoMEs() && pass_SRSel_Nminus1;
  if (computeMEs){
    SimpleParticleCollection_t daughters;
    daughters.push_back(SimpleParticle_t(25, ParticleObjectHelpers::convertCMSLorentzVectorToTLorentzVector(p4_ZZ_approx)));

    SimpleParticleCollection_t associated;
    for (auto const& jet:ak4jets_tight) associated.push_back(SimpleParticle_t(0, ParticleObjectHelpers::convertCMSLorentzVectorToTLorentzVector(jet->p4())));

    CMS3MELAHelpers::melaHandle->setCandidateDecayMode(TVar::CandidateDecay_Stable);
    CMS3MELAHelpers::melaHandle->setInputEvent(&daughters, &associated, nullptr, false);
    MEblock.computeMELABranches();
    MEblock.pushMELABranches();

    CMS3MELAHelpers::melaHandle->resetInputEvent();
  }
  else{
    MEblock.computeMELABranches();
    MEblock.pushMELABranches();
  }
  // Insert the ME values into commonEntry only when the productTreeList collection is empty.
  // Otherwise, the branches are already made.
  if (!theLooper->hasLinkedOutputTrees()){
    std::unordered_map<std::string, float> ME_values;
    MEblock.getBranchValues(ME_values);
    for (auto const& it:ME_values) commonEntry.setNamedVal(it.first, it.second);
  }

  // For the case of fakeable ids, add extra bracnhes for fakeable flavor
  if (applyFakeables){
    // All leptons are loose, so record which ones are tight.
    bool lepton_is_tight = ParticleSelectionHelpers::isTightParticle(theChosenLepton); commonEntry.setNamedVal("lepton_is_tight", lepton_is_tight);
    float event_MT = std::sqrt(2.*lepton_pt*event_pTmiss*(1.f - std::cos(dPhi_pTboson_pTmiss))); commonEntry.setNamedVal("event_MT", event_MT);

    std::vector<bool> lepton_pass_fakeable_ids; // These are the fakeable id variations for trigger
    if (theElectron){
      using namespace ElectronTriggerCutEnums;

      cms3_electron_cutbasedbits_triggeremulation_t const& trigbits_v1 = theElectron->extras.id_cutBased_triggerEmulationV1_Bits;
      cms3_electron_cutbasedbits_triggeremulation_t const& trigbits_v2 = theElectron->extras.id_cutBased_triggerEmulationV2_Bits;
      std::vector<ElectronTriggerCutTypes> testbits_v1;
      std::vector<ElectronTriggerCutTypes> testbits_v2;
      switch (SampleHelpers::getDataYear()){
      case 2016:
        // In 2016, EAs for ECAL and HCAL PF isos differ between the two version, so test all valid bits in both.
        testbits_v1 = std::vector<ElectronTriggerCutTypes>{
          kCuts_CaloIdL_TrackIdL,
          kCuts_CaloIdL_TrackIdL_IsoVL_v1,
          kCuts_CaloIdL_TrackIdL_IsoVL_v2,
          kCuts_CaloIdL_TrackIdL_IsoVL_v3,
          kCuts_CaloIdL_GsfTrkIdVL,
          kCuts_CaloIdL_MW,
          kCuts_WPLoose,
          kCuts_WPTight_v1,
          kCuts_DoublePhoton
        };
        testbits_v2 = testbits_v1;
        break;
      case 2017:
        // In 2017, EAs for ECAL and HCAL PF isos are the same between the two version, so test only the bits that are have different selection properties.
        testbits_v1 = std::vector<ElectronTriggerCutTypes>{
          kCuts_CaloIdL_TrackIdL,
          kCuts_CaloIdL_TrackIdL_IsoVL_v1,
          kCuts_CaloIdL_TrackIdL_IsoVL_v2,
          kCuts_CaloIdL_MW,
          kCuts_WPTight_v1,
          kCuts_WPTight_v2,
          kCuts_DoublePhoton
        };
        testbits_v2 = std::vector<ElectronTriggerCutTypes>{
          kCuts_CaloIdL_TrackIdL,
          kCuts_CaloIdL_TrackIdL_IsoVL_v1,
          kCuts_WPTight_v1,
          kCuts_WPTight_v2
        };
        break;
      case 2018:
        // In 2018, EAs for ECAL and HCAL PF isos are the same between the two version, so test only the bits that are have different selection properties.
        testbits_v1 = std::vector<ElectronTriggerCutTypes>{
          kCuts_CaloIdL_TrackIdL,
          kCuts_CaloIdL_TrackIdL_IsoVL_v1,
          kCuts_CaloIdL_MW,
          kCuts_WPTight_v1,
          kCuts_DoublePhoton
        };
        testbits_v2 = std::vector<ElectronTriggerCutTypes>{
          kCuts_CaloIdL_TrackIdL,
          kCuts_CaloIdL_TrackIdL_IsoVL_v1,
          kCuts_WPTight_v1
        };
        break;
      default:
        MELAerr << "LooperFunctionHelpers::looperRule: Year " << SampleHelpers::getDataYear() << " is not defined for the electron trigger emulation bit selections." << endl;
        assert(0);
      }

      lepton_pass_fakeable_ids.reserve(testbits_v1.size() + testbits_v2.size());
      for (auto const& testbit:testbits_v1) lepton_pass_fakeable_ids.push_back(HelperFunctions::test_bit(trigbits_v1, testbit));
      for (auto const& testbit:testbits_v2) lepton_pass_fakeable_ids.push_back(HelperFunctions::test_bit(trigbits_v2, testbit));
    }
    else{
      bool isFakeableBaseWithTrkIso03 = theMuon->extras.trkIso03_trackerSumPt<0.4f*lepton_pt; // i.e. TrkIsoVVL
      lepton_pass_fakeable_ids.reserve(1);
      lepton_pass_fakeable_ids.push_back(isFakeableBaseWithTrkIso03);
    }

    commonEntry.setNamedVal("lepton_pass_fakeable_ids", lepton_pass_fakeable_ids);
  }


  /*********************/
  /* RECORD THE OUTPUT */
  /*********************/
#define BRANCH_COMMAND(TYPE, NAME) commonEntry.setNamedVal(#NAME, NAME);
  BRANCH_COMMANDS;
#undef BRANCH_COMMAND
#define CONTROL_TRIGGER_COMMAND(TYPE) commonEntry.setNamedVal(Form("event_wgt_triggers_%s", #TYPE), event_wgt_triggers_##TYPE);
  if (applyFakeables){
    CONTROL_TRIGGER_COMMANDS;
  }
#undef CONTROL_TRIGGER_COMMAND

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

void LooperFunctionHelpers::setApplyFakeableId(bool applyFakeables_){
  applyFakeables = applyFakeables_;
  if (applyFakeables){
    MuonSelectionHelpers::setAllowFakeableInLooseSelection(true);
    MuonSelectionHelpers::doRequireTrackerIsolationInFakeable(-1); // Set rel. trk. iso. to -1 so that we can test it separately.
    ElectronSelectionHelpers::setAllowFakeableInLooseSelection(true);
  }
}


using namespace SystematicsHelpers;
void getTrees(
  TString strSampleSet, TString period,
  TString prodVersion, TString strdate,
  int ichunk, int nchunks,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal,
  // ME option
  bool computeMEs=false,
  // Additional option to include fakeable ids
  bool useFakeables=false,
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

  std::string systName = SystematicsHelpers::getSystName(theGlobalSyst);
  std::vector<SystematicsHelpers::SystematicVariationTypes> const disallowedSysts{
    eEleEffDn, eEleEffUp,
    eEleEffStatDn, eEleEffStatUp,
    eEleEffSystDn, eEleEffSystUp,
    eEleEffAltMCDn, eEleEffAltMCUp,

    eMuEffDn, eMuEffUp,
    eMuEffStatDn, eMuEffStatUp,
    eMuEffSystDn, eMuEffSystUp,
    eMuEffAltMCDn, eMuEffAltMCUp,

    ePhoEffDn, ePhoEffUp,

    ePUJetIdEffDn, ePUJetIdEffUp,
    eBTagSFDn, eBTagSFUp,

    ePUDn, ePUUp,
    eL1PrefiringDn, eL1PrefiringUp
  };
  if (HelperFunctions::checkListVariable(disallowedSysts, theGlobalSyst)){
    MELAout << "Systematic type " << systName << " is not allowed because the set of weights already cover it." << endl;
    return;
  }

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  bool const useSkims = !(
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

  std::vector<TriggerHelpers::TriggerType> requiredTriggers{
    TriggerHelpers::kSingleMu, TriggerHelpers::kSingleMu_HighPt,
    TriggerHelpers::kSingleEle, TriggerHelpers::kSingleEle_HighPt
  };
  auto triggerPropsCheckList = TriggerHelpers::getHLTMenuProperties(requiredTriggers);

#define CONTROL_TRIGGER_COMMAND(TYPE) std::vector<TriggerHelpers::TriggerType> requiredTriggers_##TYPE{ TriggerHelpers::k##TYPE };
  CONTROL_TRIGGER_COMMANDS;
#undef CONTROL_TRIGGER_COMMAND
  TriggerHelpers::dropSelectionCuts(TriggerHelpers::kMET_Control);
  TriggerHelpers::dropSelectionCuts(TriggerHelpers::kPFMET_Control);
  TriggerHelpers::dropSelectionCuts(TriggerHelpers::kPFHT_PFMET_Control);
  TriggerHelpers::dropSelectionCuts(TriggerHelpers::kPFMET_MHT_Control);
  TriggerHelpers::dropSelectionCuts(TriggerHelpers::kPFHT_PFMET_MHT_Control);
#define CONTROL_TRIGGER_COMMAND(TYPE) auto triggerPropsCheckList_##TYPE = TriggerHelpers::getHLTMenuProperties(requiredTriggers_##TYPE);
  CONTROL_TRIGGER_COMMANDS;
#undef CONTROL_TRIGGER_COMMAND

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

  // Set fakeable ID option
  LooperFunctionHelpers::setApplyFakeableId(useFakeables);
  if (useFakeables) computeMEs=false;

  // Set output directory
  TString coutput_main =
    "output/SingleLeptonEvents/SkimTrees/" + strdate
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
  coutput_main = coutput_main + "/" + (useFakeables ? "WithFakeableLeptons" : "DefaultLeptons");
  coutput_main = coutput_main + "/" + period;

  TDirectory* curdir = gDirectory;
  gSystem->mkdir(coutput_main, true);

  curdir->cd();

  // Create the output file
  TString coutput = SampleHelpers::getSampleIdentifier(strSampleSet);
  HelperFunctions::replaceString(coutput, "_MINIAODSIM", "");
  HelperFunctions::replaceString(coutput, "_MINIAOD", "");
  TString stroutput = Form("%s/%s", coutput_main.Data(), coutput.Data());
  stroutput += Form("_%s", systName.data());
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
  theLooper.addHLTMenu("SingleLepton", triggerPropsCheckList);
#define CONTROL_TRIGGER_COMMAND(TYPE) theLooper.addHLTMenu(#TYPE, triggerPropsCheckList_##TYPE);
  if (useFakeables){
    CONTROL_TRIGGER_COMMANDS;
  }
#undef CONTROL_TRIGGER_COMMAND
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
    TString strdsetfname = SampleHelpers::getDatasetFileName(sname);
    MELAout << "=> Accessing the input trees from " << strdsetfname << "..." << endl;
    BaseTree* sample_tree;
    if (useSkims) sample_tree = new BaseTree(
      strdsetfname,
      {
        "cms3ntuple/SingleLepton", // Includes base fakeable objects, so no need to add an exception for useFakeables=true
        "cms3ntuple/Dilepton_Control", // The other two trees are included in cases systematic variations move a few percent of events around
        "cms3ntuple/Dilepton"
      },
      ""
    );
    else sample_tree = new BaseTree(strdsetfname, "cms3ntuple/Events", "", "");
    sample_trees.push_back(sample_tree);
    sample_tree->sampleIdentifier = SampleHelpers::getSampleIdentifier(sname);
    MELAout << "\t- Sample identifier (is data ? " << isData << "): " << sample_tree->sampleIdentifier << endl;

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
          if (bname.Contains("p_Gen") || bname.Contains("LHECandMass")) has_lheMEweights = true;
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
