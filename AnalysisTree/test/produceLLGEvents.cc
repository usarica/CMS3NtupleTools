#include <cassert>
#include "common_includes.h"
#include "OffshellCutflow.h"
#include <IvyFramework/IvyAutoMELA/interface/IvyMELAHelpers.h>
#include "TStyle.h"


namespace LooperFunctionHelpers{
  using namespace std;
  using namespace IvyStreamHelpers;
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

  // Flag to apply thte low dilepton mass requirement
  bool applyLowDileptonMassReq = true;
  void setApplyLowDileptonMassReq(bool applyLowDileptonMassReq_);

  // Helpers for b-tagging
  float btag_thr_loose = -1;
  float btag_thr_medium = -1;
  float btag_thr_tight = -1;

  void setBtagWPs();

  // Helpers to keep or discard gen. ak4jet quantities
  bool keepGenAK4JetInfo = false;
  void setKeepGenAK4JetInfo(bool keepGenAK4JetInfo_);

  // Helpers to record the LHE Higgs information
  bool keepLHEGenPartInfo = false;
  void setKeepLHEGenPartInfo(bool keepLHEGenPartInfo_);

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
  IvyMELAHelpers::GMECBlock& MEblock = theLooper->getMEblock();

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
  if (hasSimpleHLTMenus && hasHLTMenuProperties){
    IVYerr << "LooperFunctionHelpers::looperRule: Defining both simple HLT menus and menus with properties is not allowed. Choose only one!" << endl;
    assert(0);
  }
  auto it_HLTMenuSimple_SingleLepton = triggerCheckListMap.find("SingleLepton");
  auto it_HLTMenuProps_SingleLepton = triggerPropsCheckListMap.find("SingleLepton");
  if (
    (hasSimpleHLTMenus && it_HLTMenuSimple_SingleLepton == triggerCheckListMap.cend())
    ||
    (hasHLTMenuProperties && it_HLTMenuProps_SingleLepton == triggerPropsCheckListMap.cend())
    ){
    IVYerr << "LooperFunctionHelpers::looperRule: The trigger type 'SingleLepton' has to be defined in this looper rule!" << endl;
    assert(0);
  }
  auto it_HLTMenuSimple_Dilepton = triggerCheckListMap.find("Dilepton");
  auto it_HLTMenuProps_Dilepton = triggerPropsCheckListMap.find("Dilepton");
  if (
    (hasSimpleHLTMenus && it_HLTMenuSimple_Dilepton == triggerCheckListMap.cend())
    ||
    (hasHLTMenuProperties && it_HLTMenuProps_Dilepton == triggerPropsCheckListMap.cend())
    ){
    IVYerr << "LooperFunctionHelpers::looperRule: The trigger type 'Dilepton' has to be defined in this looper rule!" << endl;
    assert(0);
  }
  auto it_HLTMenuSimple_SinglePhoton = triggerCheckListMap.find("SinglePhoton");
  auto it_HLTMenuProps_SinglePhoton = triggerPropsCheckListMap.find("SinglePhoton");
  if (
    (hasSimpleHLTMenus && it_HLTMenuSimple_SinglePhoton == triggerCheckListMap.cend())
    ||
    (hasHLTMenuProperties && it_HLTMenuProps_SinglePhoton == triggerPropsCheckListMap.cend())
    ){
    IVYerr << "LooperFunctionHelpers::looperRule: The trigger type 'Dilepton_SinglePhoton' has to be defined in this looper rule!" << endl;
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
    IVYerr << "LooperFunctionHelpers::looperRule: " << #TYPE << " " << #NAME << " is not registered. Please register and re-run." << endl; \
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
  BRANCH_COMMAND(float, event_wgt_L1PrefiringDn) \
  BRANCH_COMMAND(float, event_wgt_L1PrefiringUp) \
  BRANCH_COMMAND(float, event_wgt_PUDn) \
  BRANCH_COMMAND(float, event_wgt_PUUp) \
  /* Pythia weight adjustments are independent of PDF choice */ \
  BRANCH_COMMAND(float, event_wgt_adjustment_PythiaScaleDn) \
  BRANCH_COMMAND(float, event_wgt_adjustment_PythiaScaleUp) \
  /* Factorization and renormalization scale weight adjustments are independent of PDF choice (because they are only done for the default PDF set) */ \
  BRANCH_COMMAND(float, event_wgt_adjustment_PDFScaleDn) \
  BRANCH_COMMAND(float, event_wgt_adjustment_PDFScaleUp) \
  BRANCH_COMMAND(float, event_wgt_adjustment_QCDScaleDn) \
  BRANCH_COMMAND(float, event_wgt_adjustment_QCDScaleUp) \
  /* a_s(mZ) and PDF replica weight adjustments come from the specific PDF set, so they are split between 'default' vs 'NNPDF3.0'. */ \
  BRANCH_COMMAND(float, event_wgt_adjustment_AsMZDn) \
  BRANCH_COMMAND(float, event_wgt_adjustment_AsMZUp) \
  BRANCH_COMMAND(float, event_wgt_adjustment_PDFReplicaDn) \
  BRANCH_COMMAND(float, event_wgt_adjustment_PDFReplicaUp) \
  BRANCH_COMMAND(float, event_wgt_adjustment_NNPDF30) \
  BRANCH_COMMAND(float, event_wgt_adjustment_NNPDF30_AsMZDn) \
  BRANCH_COMMAND(float, event_wgt_adjustment_NNPDF30_AsMZUp) \
  BRANCH_COMMAND(float, event_wgt_adjustment_NNPDF30_PDFReplicaDn) \
  BRANCH_COMMAND(float, event_wgt_adjustment_NNPDF30_PDFReplicaUp) \
  BRANCH_COMMAND(float, event_wgt_triggers_SingleLepton) \
  BRANCH_COMMAND(float, event_wgt_triggers_Dilepton) \
  BRANCH_COMMAND(float, event_wgt_triggers_SinglePhoton) \
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
  BRANCH_COMMAND(float, event_mllg) \
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
  BRANCH_COMMAND(cms3_id_t, dilepton_id) \
  BRANCH_COMMAND(float, dilepton_pt) \
  BRANCH_COMMAND(float, dilepton_eta) \
  BRANCH_COMMAND(float, dilepton_phi) \
  BRANCH_COMMAND(float, dilepton_mass) \
  BRANCH_COMMAND(float, photon_pt) \
  BRANCH_COMMAND(float, photon_eta) \
  BRANCH_COMMAND(float, photon_phi) \
  BRANCH_COMMAND(float, photon_mass) \
  BRANCH_COMMAND(float, photon_full5x5_sigmaIEtaIEta) \
  BRANCH_COMMAND(float, photon_full5x5_sigmaIPhiIPhi) \
  BRANCH_COMMAND(float, photon_full5x5_r9) \
  BRANCH_COMMAND(float, photon_seedTime) \
  BRANCH_COMMAND(float, photon_MIPTotalEnergy) \
  BRANCH_COMMAND(bool, photon_is_genMatched_prompt) \
  BRANCH_COMMAND(bool, photon_is_conversionSafe) \
  BRANCH_COMMAND(bool, photon_is_inTime) \
  BRANCH_COMMAND(bool, photon_is_beamHaloSafe) \
  BRANCH_COMMAND(bool, photon_is_spikeSafe) \
  BRANCH_COMMAND(bool, photon_is_PFID) \
  BRANCH_COMMAND(bool, photon_is_METSafe) \
  BRANCH_COMMAND(bool, photon_pass_HGGSelection) \
  BRANCH_COMMAND(bool, photon_isGap) \
  BRANCH_COMMAND(bool, photon_isEB) \
  BRANCH_COMMAND(bool, photon_isEE) \
  BRANCH_COMMAND(bool, photon_isEBEEGap) \
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
  BRANCH_COMMAND(float, leptons_eff_StatDn) \
  BRANCH_COMMAND(float, leptons_eff_StatUp) \
  BRANCH_COMMAND(float, leptons_eff_SystDn) \
  BRANCH_COMMAND(float, leptons_eff_SystUp) \
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
  std::vector<TString> bnames_exclude;
#define BRANCH_COMMAND(TYPE, NAME) if (theGlobalSyst!=SystematicsHelpers::sNominal){ TString strtmp = #NAME; if (strtmp.EndsWith("Up") || strtmp.EndsWith("Dn")) bnames_exclude.emplace_back(strtmp); }
  BRANCH_SCALAR_COMMANDS;
  BRANCH_VECTOR_COMMANDS;
#undef BRANCH_COMMAND


  // Always assign the external weight first
  auto it_extWgt = extWgt.find(theGlobalSyst);
  if (it_extWgt==extWgt.cend()) it_extWgt = extWgt.find(SystematicsHelpers::nSystematicVariations);
  if (it_extWgt==extWgt.cend()){
    IVYerr << "LooperFunctionHelpers::looperRule: External normalization map does not have a proper weight assigned!" << endl;
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

    auto const& genInfoExtras = genInfo->extras;
    event_wgt_adjustment_PythiaScaleDn = genInfoExtras.PythiaWeight_isr_muR0p25 * genInfoExtras.PythiaWeight_fsr_muR0p25;
    event_wgt_adjustment_PythiaScaleUp = genInfoExtras.PythiaWeight_isr_muR4 * genInfoExtras.PythiaWeight_fsr_muR4;
    event_wgt_adjustment_PDFScaleDn = genInfoExtras.LHEweight_QCDscale_muR1_muF0p5;
    event_wgt_adjustment_PDFScaleUp = genInfoExtras.LHEweight_QCDscale_muR1_muF2;
    event_wgt_adjustment_QCDScaleDn = genInfoExtras.LHEweight_QCDscale_muR0p5_muF1;
    event_wgt_adjustment_QCDScaleUp = genInfoExtras.LHEweight_QCDscale_muR2_muF1;
    event_wgt_adjustment_AsMZDn = genInfoExtras.LHEweight_AsMZ_Dn_default;
    event_wgt_adjustment_AsMZUp = genInfoExtras.LHEweight_AsMZ_Up_default;
    event_wgt_adjustment_NNPDF30_AsMZDn = genInfoExtras.LHEweight_AsMZ_Dn_NNPDF30;
    event_wgt_adjustment_NNPDF30_AsMZUp = genInfoExtras.LHEweight_AsMZ_Up_NNPDF30;
    event_wgt_adjustment_PDFReplicaDn = genInfoExtras.LHEweight_PDFVariation_Dn_default;
    event_wgt_adjustment_PDFReplicaUp = genInfoExtras.LHEweight_PDFVariation_Up_default;
    event_wgt_adjustment_NNPDF30_PDFReplicaDn = genInfoExtras.LHEweight_PDFVariation_Dn_NNPDF30;
    event_wgt_adjustment_NNPDF30_PDFReplicaUp = genInfoExtras.LHEweight_PDFVariation_Up_NNPDF30;
    // Adjust for cases where NNPDF 3.0 does not exist.
    if (
      event_wgt_adjustment_NNPDF30==1.f
      &&
      event_wgt_adjustment_NNPDF30_AsMZDn==1.f && event_wgt_adjustment_NNPDF30_AsMZUp==1.f
      &&
      event_wgt_adjustment_NNPDF30_PDFReplicaDn==1.f && event_wgt_adjustment_NNPDF30_PDFReplicaUp==1.f
      ){
      event_wgt_adjustment_NNPDF30_AsMZDn = event_wgt_adjustment_AsMZDn;
      event_wgt_adjustment_NNPDF30_AsMZUp = event_wgt_adjustment_AsMZUp;
      event_wgt_adjustment_NNPDF30_PDFReplicaDn = event_wgt_adjustment_PDFReplicaDn;
      event_wgt_adjustment_NNPDF30_PDFReplicaUp = event_wgt_adjustment_PDFReplicaUp;
    }

    genmet_pTmiss = genInfoExtras.genmet_met;
    genmet_phimiss = genInfoExtras.genmet_metPhi;
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

  vertexHandler->constructVertices();
  if (!vertexHandler->hasGoodPrimaryVertex()) return false;
  event_n_vtxs_good = vertexHandler->getNGoodVertices();

  eventFilter->constructFilters(simEventHandler);
  if (isData && !eventFilter->isUniqueDataEvent()) return false;

  if (!eventFilter->passCommonSkim() || !eventFilter->passMETFilters(EventFilterHandler::kMETFilters_Standard)) return false;
  event_pass_tightMETFilters = eventFilter->passMETFilters(EventFilterHandler::kMETFilters_Tight);

  pfcandidateHandler->constructPFCandidates(theGlobalSyst);
  auto const& pfcandidates = pfcandidateHandler->getProducts();

  muonHandler->constructMuons(theGlobalSyst, &pfcandidates);
  electronHandler->constructElectrons(theGlobalSyst, &pfcandidates);
  photonHandler->constructPhotons(theGlobalSyst, &pfcandidates);
  particleDisambiguator.disambiguateParticles(muonHandler, electronHandler, photonHandler);

  std::unordered_map<ParticleObject const*, float> lepton_eff_map;
  std::unordered_map<ParticleObject const*, float> lepton_eff_StatDn_map;
  std::unordered_map<ParticleObject const*, float> lepton_eff_StatUp_map;
  std::unordered_map<ParticleObject const*, float> lepton_eff_SystDn_map;
  std::unordered_map<ParticleObject const*, float> lepton_eff_SystUp_map;

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
    float theEff = 1;
    if (!isData){
      theSF = 1; theEff = 1; muonSFHandler->getIdIsoSFAndEff(theGlobalSyst, part, theSF, &theEff); theSF = std::max(1e-5f, theSF); SF_muons *= theSF; lepton_eff_map[part] = theEff;
      theSF = 1; theEff = 1; muonSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eMuEffStatDn, part, theSF, &theEff); theSF = std::max(1e-5f, theSF); SF_muons_StatDn *= theSF; lepton_eff_StatDn_map[part] = theEff;
      theSF = 1; theEff = 1; muonSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eMuEffStatUp, part, theSF, &theEff); theSF = std::max(1e-5f, theSF); SF_muons_StatUp *= theSF; lepton_eff_StatUp_map[part] = theEff;
      theSF = 1; theEff = 1; muonSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eMuEffSystDn, part, theSF, &theEff); theSF = std::max(1e-5f, theSF); SF_muons_SystDn *= theSF; lepton_eff_SystDn_map[part] = theEff;
      theSF = 1; theEff = 1; muonSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eMuEffSystUp, part, theSF, &theEff); theSF = std::max(1e-5f, theSF); SF_muons_SystUp *= theSF; lepton_eff_SystUp_map[part] = theEff;
      theSF = 1; theEff = 1; muonSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eMuEffAltMCDn, part, theSF, &theEff); theSF = std::max(1e-5f, theSF); SF_muons_AltMCDn *= theSF;
      theSF = 1; theEff = 1; muonSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eMuEffAltMCUp, part, theSF, &theEff); theSF = std::max(1e-5f, theSF); SF_muons_AltMCUp *= theSF;
    }
    else if (ParticleSelectionHelpers::isTightParticle(part)){
      theEff = 1; muonSFHandler->getIdIsoSFAndEff(theGlobalSyst, part->pt(), part->eta(), true, true, true, theSF, &theEff); lepton_eff_map[part] = theEff;
      theEff = 1; muonSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eMuEffStatDn, part->pt(), part->eta(), true, true, true, theSF, &theEff); lepton_eff_StatDn_map[part] = theEff;
      theEff = 1; muonSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eMuEffStatUp, part->pt(), part->eta(), true, true, true, theSF, &theEff); lepton_eff_StatUp_map[part] = theEff;
      theEff = 1; muonSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eMuEffSystDn, part->pt(), part->eta(), true, true, true, theSF, &theEff); lepton_eff_SystDn_map[part] = theEff;
      theEff = 1; muonSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eMuEffSystUp, part->pt(), part->eta(), true, true, true, theSF, &theEff); lepton_eff_SystUp_map[part] = theEff;
    }

    if (!ParticleSelectionHelpers::isTightParticle(part) && part->testSelectionBit(MuonSelectionHelpers::kFakeableBase)) event_n_leptons_fakeableBase++;
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
    float theEff = 1;
    if (!isData){
      theSF = 1; theEff = 1; electronSFHandler->getIdIsoSFAndEff(theGlobalSyst, part, theSF, &theEff); theSF = std::max(1e-5f, theSF); SF_electrons *= theSF; lepton_eff_map[part] = theEff;
      theSF = 1; theEff = 1; electronSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eEleEffStatDn, part, theSF, &theEff); theSF = std::max(1e-5f, theSF); SF_electrons_StatDn *= theSF; lepton_eff_StatDn_map[part] = theEff;
      theSF = 1; theEff = 1; electronSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eEleEffStatUp, part, theSF, &theEff); theSF = std::max(1e-5f, theSF); SF_electrons_StatUp *= theSF; lepton_eff_StatUp_map[part] = theEff;
      theSF = 1; theEff = 1; electronSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eEleEffSystDn, part, theSF, &theEff); theSF = std::max(1e-5f, theSF); SF_electrons_SystDn *= theSF; lepton_eff_SystDn_map[part] = theEff;
      theSF = 1; theEff = 1; electronSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eEleEffSystUp, part, theSF, &theEff); theSF = std::max(1e-5f, theSF); SF_electrons_SystUp *= theSF; lepton_eff_SystUp_map[part] = theEff;
      theSF = 1; theEff = 1; electronSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eEleEffAltMCDn, part, theSF, &theEff); theSF = std::max(1e-5f, theSF); SF_electrons_AltMCDn *= theSF;
      theSF = 1; theEff = 1; electronSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eEleEffAltMCUp, part, theSF, &theEff); theSF = std::max(1e-5f, theSF); SF_electrons_AltMCUp *= theSF;
    }
    else if (ParticleSelectionHelpers::isTightParticle(part)){
      theEff = 1; electronSFHandler->getIdIsoSFAndEff(theGlobalSyst, part->pt(), part->etaSC(), part->isGap(), true, true, true, theSF, &theEff); lepton_eff_map[part] = theEff;
      theEff = 1; electronSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eEleEffStatDn, part->pt(), part->etaSC(), part->isGap(), true, true, true, theSF, &theEff); lepton_eff_StatDn_map[part] = theEff;
      theEff = 1; electronSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eEleEffStatUp, part->pt(), part->etaSC(), part->isGap(), true, true, true, theSF, &theEff); lepton_eff_StatUp_map[part] = theEff;
      theEff = 1; electronSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eEleEffSystDn, part->pt(), part->etaSC(), part->isGap(), true, true, true, theSF, &theEff); lepton_eff_SystDn_map[part] = theEff;
      theEff = 1; electronSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eEleEffSystUp, part->pt(), part->etaSC(), part->isGap(), true, true, true, theSF, &theEff); lepton_eff_SystUp_map[part] = theEff;
    }

    if (!ParticleSelectionHelpers::isTightParticle(part) && part->testSelectionBit(ElectronSelectionHelpers::kFakeableBase)) event_n_leptons_fakeableBase++;
  }
  event_wgt_SFs_electrons = SF_electrons;
  event_wgt_SFs_electrons_StatDn = SF_electrons_StatDn;
  event_wgt_SFs_electrons_StatUp = SF_electrons_StatUp;
  event_wgt_SFs_electrons_SystDn = SF_electrons_SystDn;
  event_wgt_SFs_electrons_SystUp = SF_electrons_SystUp;
  event_wgt_SFs_electrons_AltMCDn = SF_electrons_AltMCDn;
  event_wgt_SFs_electrons_AltMCUp = SF_electrons_AltMCUp;

  auto const& photons = photonHandler->getProducts();
  unsigned int n_photons_veto = 0;
  float SF_photons = 1;
  float SF_photons_EffDn = 1;
  float SF_photons_EffUp = 1;
  PhotonObject const* theChosenPhoton = nullptr;
  for (auto const& part:photons){
    if (!isData){
      float theSF = 1;
      photonSFHandler->getIdIsoSFAndEff(theGlobalSyst, part, theSF, nullptr); theSF = std::max(1e-5f, theSF); SF_photons *= theSF;
      theSF = 1; photonSFHandler->getIdIsoSFAndEff(SystematicsHelpers::ePhoEffDn, part, theSF, nullptr); theSF = std::max(1e-5f, theSF); SF_photons_EffDn *= theSF;
      theSF = 1; photonSFHandler->getIdIsoSFAndEff(SystematicsHelpers::ePhoEffUp, part, theSF, nullptr); theSF = std::max(1e-5f, theSF); SF_photons_EffUp *= theSF;
    }

    if (!theChosenPhoton && ParticleSelectionHelpers::isTightParticle(part)) theChosenPhoton = part;
    else if (ParticleSelectionHelpers::isVetoParticle(part)) n_photons_veto++;
  }
  event_wgt_SFs_photons = SF_photons;
  event_wgt_SFs_photons_EffDn = SF_photons_EffDn;
  event_wgt_SFs_photons_EffUp = SF_photons_EffUp;
  if (!theChosenPhoton || n_photons_veto!=0) return false;

  photon_pt = theChosenPhoton->pt();
  photon_eta = theChosenPhoton->eta();
  photon_phi = theChosenPhoton->phi();
  photon_mass = theChosenPhoton->m();
  photon_is_genMatched_prompt = theChosenPhoton->extras.is_genMatched_prompt;
  photon_is_conversionSafe = theChosenPhoton->testSelection(PhotonSelectionHelpers::kConversionSafe);
  photon_is_inTime = theChosenPhoton->testSelection(PhotonSelectionHelpers::kInTimeSeed);
  photon_is_beamHaloSafe = theChosenPhoton->testSelection(PhotonSelectionHelpers::kBeamHaloSafe);
  photon_is_spikeSafe = theChosenPhoton->testSelection(PhotonSelectionHelpers::kSpikeSafe);
  photon_is_PFID = theChosenPhoton->testSelection(PhotonSelectionHelpers::kPFPhotonId);
  photon_is_METSafe = theChosenPhoton->testSelection(PhotonSelectionHelpers::kPFMETSafe);
  photon_pass_HGGSelection = theChosenPhoton->extras.id_cutBased_HGG_Bits;
  photon_isGap = theChosenPhoton->isGap();
  photon_isEBEEGap = theChosenPhoton->isEBEEGap();
  photon_isEB = theChosenPhoton->isEB();
  photon_isEE = theChosenPhoton->isEE();
  photon_full5x5_sigmaIEtaIEta = theChosenPhoton->extras.full5x5_sigmaIEtaIEta;
  photon_full5x5_sigmaIPhiIPhi = theChosenPhoton->extras.full5x5_sigmaIPhiIPhi;
  photon_full5x5_r9 = theChosenPhoton->extras.full5x5_r9;
  photon_seedTime = theChosenPhoton->extras.seedTime;
  photon_MIPTotalEnergy = theChosenPhoton->extras.MIPTotalEnergy;

  isotrackHandler->constructIsotracks(&muons, &electrons);
  bool hasVetoIsotrack = false;
  for (auto const& isotrack:isotrackHandler->getProducts()){
    if (isotrack->testSelectionBit(IsotrackSelectionHelpers::kPreselectionVeto)){
      hasVetoIsotrack = true;
      break;
    }
  }
  if (hasVetoIsotrack) return false;

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
  // Veto also emu final states
  if (!theChosenDilepton || !theChosenDilepton->isSF() || nTightDilep>1) return false;

  dilepton_id = theChosenDilepton->getDaughter(0)->pdgId()*theChosenDilepton->getDaughter(1)->pdgId();
  dilepton_pt = theChosenDilepton->pt();
  dilepton_eta = theChosenDilepton->eta();
  dilepton_phi = theChosenDilepton->phi();
  dilepton_mass = theChosenDilepton->m();
  // Filter out low-mass resonances
  if (applyLowDileptonMassReq && dilepton_mass<12.f) return false;

  jetHandler->constructJetMET(simEventHandler, theGlobalSyst, &muons, &electrons, &photons, &pfcandidates);
  auto const& ak4jets = jetHandler->getAK4Jets();
  auto const& ak8jets = jetHandler->getAK8Jets();

  std::vector<ParticleObject const*> leptons_TOmatched_SingleLepton;
  if (hasSimpleHLTMenus){
    event_wgt_triggers_SingleLepton = eventFilter->getTriggerWeight(it_HLTMenuSimple_SingleLepton->second);
    event_wgt_triggers_Dilepton = eventFilter->getTriggerWeight(it_HLTMenuSimple_Dilepton->second);
    event_wgt_triggers_SinglePhoton = eventFilter->getTriggerWeight(it_HLTMenuSimple_SinglePhoton->second);
  }
  else if (hasHLTMenuProperties){
    event_wgt_triggers_SingleLepton = eventFilter->getTriggerWeight(
      it_HLTMenuProps_SingleLepton->second,
      &muons, &electrons, nullptr, nullptr, nullptr, nullptr,
      nullptr, &leptons_TOmatched_SingleLepton
    );
    event_wgt_triggers_Dilepton = eventFilter->getTriggerWeight(
      it_HLTMenuProps_Dilepton->second,
      &muons, &electrons, nullptr, nullptr, nullptr, nullptr
    );
    event_wgt_triggers_SinglePhoton = eventFilter->getTriggerWeight(
      it_HLTMenuProps_SinglePhoton->second,
      nullptr, nullptr, &photons, nullptr, nullptr, nullptr
    );
  }
  if ((event_wgt_triggers_SingleLepton + event_wgt_triggers_Dilepton + event_wgt_triggers_SinglePhoton) == 0.f) return false;

  // Test HEM filter
  if (!eventFilter->test2018HEMFilter(simEventHandler, nullptr, nullptr, &ak4jets, &ak8jets)) return false;
  if (!eventFilter->testNoisyJetFilter(simEventHandler, ak4jets)) return false;

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
    it_eff = lepton_eff_StatDn_map.find(dau);
    leptons_eff_StatDn.push_back((it_eff==lepton_eff_StatDn_map.end() ? 1 : it_eff->second));
    it_eff = lepton_eff_StatUp_map.find(dau);
    leptons_eff_StatUp.push_back((it_eff==lepton_eff_StatUp_map.end() ? 1 : it_eff->second));
    it_eff = lepton_eff_SystDn_map.find(dau);
    leptons_eff_SystDn.push_back((it_eff==lepton_eff_SystDn_map.end() ? 1 : it_eff->second));
    it_eff = lepton_eff_SystUp_map.find(dau);
    leptons_eff_SystUp.push_back((it_eff==lepton_eff_SystUp_map.end() ? 1 : it_eff->second));

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

  // Accumulate ak8 jets as well
  for (auto const& jet:ak8jets){
    if (!ParticleSelectionHelpers::isTightJet(jet)) continue;
    if (jet->mass()<60.f) continue;

    ak8jets_pt.push_back(jet->pt());
    ak8jets_eta.push_back(jet->eta());
    ak8jets_phi.push_back(jet->phi());
    ak8jets_mass.push_back(jet->mass());
  }

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

  // dPhi variables between pTmiss and boson are corrected in this looper by adding the momentum of the dilepton object
  ParticleObject::LorentzVector_t p4_llmet = event_met_p4 + theChosenDilepton->p4();
  float pt_llmet = p4_llmet.Pt();
  float phi_llmet = p4_llmet.Phi();

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
    HelperFunctions::deltaPhi(float(jet->phi()), phi_llmet, dphi_tmp); dphi_tmp = std::abs(dphi_tmp);
    min_abs_dPhi_pTj_pTmiss = std::min(min_abs_dPhi_pTj_pTmiss, dphi_tmp);
  }
  
  // Compute dPhi between the dilepton and pTmiss vector
  dPhi_pTboson_pTmiss = theChosenPhoton->deltaPhi(phi_llmet);
  HelperFunctions::deltaPhi(float((theChosenPhoton->p4()+sump4_ak4jets).Phi()), phi_llmet, dPhi_pTbosonjets_pTmiss);

  // Compute mass variables
  event_mllg = (theChosenPhoton->p4() + theChosenDilepton->p4()).M();

  float const& etamiss_approx = photon_eta;
  ParticleObject::LorentzVector_t p4_photon_Z_approx; p4_photon_Z_approx = ParticleObject::PolarLorentzVector_t(photon_pt, photon_eta, photon_phi, PDGHelpers::Zmass);
  ParticleObject::LorentzVector_t p4_ZZ_approx; p4_ZZ_approx = ParticleObject::PolarLorentzVector_t(pt_llmet, etamiss_approx, phi_llmet, PDGHelpers::Zmass);
  p4_ZZ_approx = p4_ZZ_approx + p4_photon_Z_approx;

  event_mTZZ = std::sqrt(
    std::pow(
    (
      std::sqrt(std::pow(photon_pt, 2) + std::pow(p4_photon_Z_approx.M(), 2))
      + std::sqrt(std::pow(pt_llmet, 2) + std::pow(PDGHelpers::Zmass, 2)) // Do not use the invariant mass of the ll+MET p4 because ll+MET is now the MET proxy, and the formula uses Z mass for this part.
      ), 2
    )
    - std::pow((p4_photon_Z_approx + p4_llmet).Pt(), 2)
  );
  event_mZZ = p4_ZZ_approx.M();

  // Apply veto from SR-like selection
  OffshellCutflow::setActiveFinalState(OffshellCutflow::fs_ZZ_2l2nu);
  std::vector<bool> const v_passZZ2l2nuSRlikeSelection={
    OffshellCutflow::check_pTboson(photon_pt),
    OffshellCutflow::check_pTmiss(pt_llmet),
    OffshellCutflow::check_dPhi_pTll_pTmiss(dPhi_pTboson_pTmiss),
    OffshellCutflow::check_dPhi_pTlljets_pTmiss(dPhi_pTbosonjets_pTmiss),
    OffshellCutflow::check_Nb_veto(event_n_ak4jets_pt30_btagged_loose)
  };
  OffshellCutflow::setActiveFinalState(OffshellCutflow::fs_WW_2l2nu);
  std::vector<bool> const v_passWW2l2nuSRlikeSelection={
    OffshellCutflow::check_pTboson(photon_pt),
    OffshellCutflow::check_pTmiss(pt_llmet),
    OffshellCutflow::check_dPhi_pTll_pTmiss(dPhi_pTboson_pTmiss),
    OffshellCutflow::check_dPhi_pTlljets_pTmiss(dPhi_pTbosonjets_pTmiss),
    OffshellCutflow::check_Nb_veto(event_n_ak4jets_pt30_btagged_loose)
  };
  unsigned short n_passZZ2l2nuSRlikeSelection = 0;
  unsigned short n_passWW2l2nuSRlikeSelection = 0;
  for (auto const& selreq:v_passZZ2l2nuSRlikeSelection){ if (selreq) n_passZZ2l2nuSRlikeSelection++; }
  for (auto const& selreq:v_passWW2l2nuSRlikeSelection){ if (selreq) n_passWW2l2nuSRlikeSelection++; }
  bool const pass_SRSel_Nminus1 = (n_passZZ2l2nuSRlikeSelection >= v_passZZ2l2nuSRlikeSelection.size()-1) || (n_passWW2l2nuSRlikeSelection >= v_passWW2l2nuSRlikeSelection.size()-1);

  // Compute MEs
  bool computeMEs = theLooper->hasRecoMEs() && pass_SRSel_Nminus1;
  if (computeMEs){
    SimpleParticleCollection_t daughters;
    daughters.push_back(SimpleParticle_t(25, ParticleObjectHelpers::convertCMSLorentzVectorToTLorentzVector(p4_ZZ_approx)));

    SimpleParticleCollection_t associated;
    for (auto const& jet:ak4jets_tight) associated.push_back(SimpleParticle_t(0, ParticleObjectHelpers::convertCMSLorentzVectorToTLorentzVector(jet->p4())));

    IvyMELAHelpers::melaHandle->setCandidateDecayMode(TVar::CandidateDecay_Stable);
    IvyMELAHelpers::melaHandle->setInputEvent(&daughters, &associated, nullptr, false);
    MEblock.computeMELABranches();
    MEblock.pushMELABranches();

    IvyMELAHelpers::melaHandle->resetInputEvent();
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

  if (!isData && keepGenAK4JetInfo){
    auto const& genak4jets = genInfoHandler->getGenAK4Jets();
    std::vector<float> genak4jets_pt; genak4jets_pt.reserve(genak4jets.size());
    std::vector<float> genak4jets_eta; genak4jets_eta.reserve(genak4jets.size());
    std::vector<float> genak4jets_phi; genak4jets_phi.reserve(genak4jets.size());
    std::vector<float> genak4jets_mass; genak4jets_mass.reserve(genak4jets.size());
    for (auto const& jet:genak4jets){
      float const jet_pt = jet->pt();
      float const jet_eta = jet->eta();
      float const jet_phi = jet->phi();
      float const jet_mass = jet->mass();
      genak4jets_pt.push_back(jet_pt);
      genak4jets_eta.push_back(jet_eta);
      genak4jets_phi.push_back(jet_phi);
      genak4jets_mass.push_back(jet_mass);
    }
    commonEntry.setNamedVal("genak4jets_pt", genak4jets_pt);
    commonEntry.setNamedVal("genak4jets_eta", genak4jets_eta);
    commonEntry.setNamedVal("genak4jets_phi", genak4jets_phi);
    commonEntry.setNamedVal("genak4jets_mass", genak4jets_mass);
  }
  if (!isData && keepLHEGenPartInfo){
    auto const& lheparticles = genInfoHandler->getLHEParticles();
    for (auto const& part:lheparticles){
      if (PDGHelpers::isAHiggs(part->pdgId())){
        commonEntry.setNamedVal<float>("lheHiggs_pt", part->pt());
        commonEntry.setNamedVal<float>("lheHiggs_rapidity", part->rapidity());
        commonEntry.setNamedVal<float>("lheHiggs_mass", part->mass());
        break;
      }
    }

    ParticleObject::LorentzVector_t genpromptparticles_sump4;
    auto const& genparticles = genInfoHandler->getGenParticles();
    for (auto const& part:genparticles){
      if (
        part->testSelectionBit(GenParticleSelectionHelpers::kHardPromptFinalVisibleParticle)
        ||
        (part->extras.isPromptFinalState && PDGHelpers::isANeutrino(part->pdgId()))
        ) genpromptparticles_sump4 += part->p4();
    }

    commonEntry.setNamedVal<float>("genpromptparticles_sump4_pt", genpromptparticles_sump4.Pt());
    commonEntry.setNamedVal<float>("genpromptparticles_sump4_mass", genpromptparticles_sump4.M());
    commonEntry.setNamedVal<float>("genpromptparticles_sump4_rapidity", genpromptparticles_sump4.Rapidity());
  }

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

void LooperFunctionHelpers::setApplyLowDileptonMassReq(bool applyLowDileptonMassReq_){ applyLowDileptonMassReq = applyLowDileptonMassReq_; }

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

void LooperFunctionHelpers::setKeepGenAK4JetInfo(bool keepGenAK4JetInfo_){ keepGenAK4JetInfo = keepGenAK4JetInfo_; }

void LooperFunctionHelpers::setKeepLHEGenPartInfo(bool keepLHEGenPartInfo_){ keepLHEGenPartInfo = keepLHEGenPartInfo_; }


using namespace SystematicsHelpers;
void getTrees(
  TString strSampleSet, TString period,
  TString prodVersion, TString strdate,
  int ichunk, int nchunks,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal,
  // ME option
  bool computeMEs=false,
  // Apply low mass dilepton cut
  bool applyLowDileptonMassReq=true,
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
    IVYout << "Systematic type " << systName << " is not allowed because the set of weights already cover it." << endl;
    return;
  }

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  TString strSampleSet_lower; HelperFunctions::lowercase(strSampleSet, strSampleSet_lower);
  bool const isHighMassPOWHEG =
    strSampleSet_lower.Contains("powheg")
    &&
    (
      strSampleSet_lower.Contains("jhugen")
      ||
      strSampleSet.BeginsWith("GGH")
      ||
      strSampleSet.BeginsWith("VBF")
      ||
      strSampleSet.BeginsWith("WplusH")
      ||
      strSampleSet.BeginsWith("WminusH")
      ||
      strSampleSet.BeginsWith("ZH")
      );
  bool const useSkims = !isHighMassPOWHEG;

  SampleHelpers::configure(period, Form("%s:%s", (useSkims ? "store_skims" : "store"), prodVersion.Data()));

  const float lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::getDataPeriod());

  std::vector<TString> const validDataPeriods = SampleHelpers::getValidDataPeriods();
  size_t const nValidDataPeriods = validDataPeriods.size();

  BtagHelpers::setBtagWPType(BtagHelpers::kDeepFlav_Loose);
  LooperFunctionHelpers::setBtagWPs();

  std::vector<TriggerHelpers::TriggerType> requiredTriggers_SingleLepton{
    TriggerHelpers::kSingleMu, TriggerHelpers::kSingleMu_HighPt,
    TriggerHelpers::kSingleEle, TriggerHelpers::kSingleEle_HighPt
  };
  std::vector<TriggerHelpers::TriggerType> requiredTriggers_Dilepton{
    TriggerHelpers::kDoubleMu,
    TriggerHelpers::kDoubleEle, TriggerHelpers::kDoubleEle_HighPt
  };
  std::vector<TriggerHelpers::TriggerType> requiredTriggers_SinglePhoton{ TriggerHelpers::kSinglePho };
  auto triggerPropsCheckList_SingleLepton = TriggerHelpers::getHLTMenuProperties(requiredTriggers_SingleLepton);
  auto triggerPropsCheckList_Dilepton = TriggerHelpers::getHLTMenuProperties(requiredTriggers_Dilepton);
  auto triggerPropsCheckList_SinglePhoton = TriggerHelpers::getHLTMenuProperties(requiredTriggers_SinglePhoton);

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

  // Set flag to apply the low dilepton mass requirement
  LooperFunctionHelpers::setApplyLowDileptonMassReq(applyLowDileptonMassReq);

  // Set output directory
  TString coutput_main =
    "output/LLGEvents/SkimTrees/" + strdate
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
  stroutput += Form("_%s", systName.data());
  if (nchunks>0) stroutput = stroutput + Form("_%i_of_%i", ichunk, nchunks);
  stroutput += ".root";
  TFile* foutput = TFile::Open(stroutput, "recreate");
  foutput->cd();
  BaseTree* tout = new BaseTree("SkimTree");
  IVYout << "Created output file " << stroutput << "..." << endl;
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
  theLooper.addHLTMenu("Dilepton", triggerPropsCheckList_Dilepton);
  theLooper.addHLTMenu("SinglePhoton", triggerPropsCheckList_SinglePhoton);
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
    IVYout << "=> Accessing the input trees from " << strdsetfname << "..." << endl;
    BaseTree* sample_tree = new BaseTree(strdsetfname, (useSkims ? "cms3ntuple/Dilepton" : "cms3ntuple/Events"), "", ""); sample_trees.push_back(sample_tree);
    sample_tree->sampleIdentifier = SampleHelpers::getSampleIdentifier(sname);
    float const sampleMH = SampleHelpers::findPoleMass(sample_tree->sampleIdentifier);
    IVYout << "\t- Sample identifier (is data ? " << isData << "): " << sample_tree->sampleIdentifier << endl;

    std::vector<TString> allbranchnames; sample_tree->getValidBranchNamesWithoutAlias(allbranchnames, false);

    const int nEntries = sample_tree->getNEvents();
    bool hasTaus = false;
    double sum_wgts = (isData ? 1.f : 0.f);
    double sum_wgts_PUDn = sum_wgts;
    double sum_wgts_PUUp = sum_wgts;
    double sum_wgts_raw_noveto = sum_wgts;
    double sum_wgts_raw_withveto = sum_wgts;
    double sum_wgts_raw_withveto_defaultMemberZero = sum_wgts;
    float xsec = 1;
    float xsec_scale = 1;
    float BR_scale = 1;
    if (!isData){
      // Get cross section
      sample_tree->bookBranch<float>("xsec", 0.f);
      sample_tree->getSelectedEvent(0);
      sample_tree->getVal("xsec", xsec);
      sample_tree->releaseBranch("xsec");
      xsec *= 1000.;

      bool has_lheMEweights = false;
      bool has_lheparticles = false;
      bool has_genparticles = false;
      bool has_genak4jets = false;
      for (auto const& bname:allbranchnames){
        if (bname.Contains("p_Gen") || bname.Contains("LHECandMass")) has_lheMEweights = true;
        else if (bname.Contains(GenInfoHandler::colName_lheparticles)) has_lheparticles = true;
        else if (bname.Contains(GenInfoHandler::colName_genparticles)) has_genparticles = true;
        else if (bname.Contains(GenInfoHandler::colName_genak4jets)) has_genak4jets = true;
      }

      // Book branches
      simEventHandler.bookBranches(sample_tree);

      bool const isHighMassPOWHEGSample = (isHighMassPOWHEG && sampleMH>0.f);
      bool const checkForTaus = has_lheparticles && isHighMassPOWHEGSample;

      genInfoHandler.setAcquireLHEMEWeights(false);
      genInfoHandler.setAcquireLHEParticles(checkForTaus);
      genInfoHandler.setAcquireGenParticles(false);
      genInfoHandler.bookBranches(sample_tree);

      // Get sum of weights
      std::vector<TString> inputfilenames = SampleHelpers::getDatasetFileNames(sname);
      bool hasCounters = true;
      {
        int bin_syst = 1 + 1*(theGlobalSyst==SystematicsHelpers::ePUDn) + 2*(theGlobalSyst==SystematicsHelpers::ePUUp);
        int bin_period = 1;
        for (unsigned int iperiod=0; iperiod<nValidDataPeriods; iperiod++){
          if (validDataPeriods.at(iperiod)==SampleHelpers::theDataPeriod){ bin_period += iperiod+1; break; }
        }
        IVYout << "Checking counters histogram bin (" << bin_syst << ", " << bin_period << ") to obtain the sum of weights if the counters histogram exists..." << endl;
        for (auto const& fname:inputfilenames){
          TFile* ftmp = TFile::Open(fname, "read");
          TH2D* hCounters = (TH2D*) ftmp->Get("cms3ntuple/Counters");
          if (!hCounters){
            hasCounters = false;
            sum_wgts = sum_wgts_PUDn = sum_wgts_PUUp = 0;
            break;
          }
          IVYout << "\t- Successfully found the counters histogram in " << fname << endl;
          sum_wgts += hCounters->GetBinContent(bin_syst, bin_period);
          sum_wgts_PUDn += hCounters->GetBinContent(2, bin_period);
          sum_wgts_PUUp += hCounters->GetBinContent(3, bin_period);
          sum_wgts_raw_withveto += hCounters->GetBinContent(0, 0);
          sum_wgts_raw_noveto += hCounters->GetBinContent(0, 0) / (1. - hCounters->GetBinContent(0, 1));
          ftmp->Close();
        }
        if (hasCounters) IVYout << "\t- Obtained the weights from " << inputfilenames.size() << " files..." << endl;
      }
      if (!hasCounters && useSkims){
        IVYerr << "Skims should have contained counters histograms!" << endl;
        assert(0);
      }
      if (!hasCounters){
        IVYout << "No counters histograms are found. Initiation loop over " << nEntries << " events to determine the sample normalization:" << endl;

        simEventHandler.wrapTree(sample_tree);
        genInfoHandler.wrapTree(sample_tree);

        unsigned int n_zero_genwgts=0;
        double frac_zero_genwgts=0;
        for (int ev=0; ev<nEntries; ev++){
          HelperFunctions::progressbar(ev, nEntries);
          sample_tree->getEvent(ev);

          genInfoHandler.constructGenInfo(SystematicsHelpers::sNominal); // Use sNominal here in order to get the weight that corresponds to xsec
          auto const& genInfo = genInfoHandler.getGenInfo();

          if (checkForTaus){
            auto const& lheparticles = genInfoHandler.getLHEParticles();
            if (!hasTaus){
              for (auto const& part:lheparticles){
                if (part->status()==1 && std::abs(part->pdgId())==15){
                  hasTaus = true;
                  genInfoHandler.setAcquireLHEParticles(false);
                  break;
                }
              }
            }
          }

          double genwgt = genInfo->getGenWeight(true);
          double genwgt_defaultMemberZero = genInfo->extras.LHEweight_defaultMemberZero;
          if (genwgt==0.){
            n_zero_genwgts++;
            continue;
          }

          simEventHandler.constructSimEvent();

          sum_wgts_raw_withveto_defaultMemberZero += genwgt_defaultMemberZero;
          sum_wgts_raw_withveto += genwgt;
          sum_wgts += genwgt * simEventHandler.getPileUpWeight(theGlobalSyst);
          sum_wgts_PUDn += genwgt * simEventHandler.getPileUpWeight(SystematicsHelpers::ePUDn);
          sum_wgts_PUUp += genwgt * simEventHandler.getPileUpWeight(SystematicsHelpers::ePUUp);
        }
        if (nEntries>0) frac_zero_genwgts = double(n_zero_genwgts)/double(nEntries);
        sum_wgts_raw_noveto = sum_wgts_raw_withveto / (1. - frac_zero_genwgts);
      }
      xsec_scale = sum_wgts_raw_withveto / sum_wgts_raw_noveto;
      if (isHighMassPOWHEGSample) BR_scale = SampleHelpers::calculateAdjustedHiggsBREff(sname, sum_wgts_raw_withveto_defaultMemberZero, sum_wgts_raw_withveto, hasTaus);

      // Reset gen. and LHE particle settings
      genInfoHandler.setAcquireLHEMEWeights(has_lheMEweights);
      genInfoHandler.setAcquireLHEParticles(has_lheparticles);
      genInfoHandler.setAcquireGenParticles(has_genparticles);
      genInfoHandler.setAcquireGenAK4Jets(has_genak4jets && theGlobalSyst==sNominal);
      genInfoHandler.setDoGenJetsVDecayCleaning(has_genak4jets && theGlobalSyst==sNominal);
      genInfoHandler.bookBranches(sample_tree);

      LooperFunctionHelpers::setKeepGenAK4JetInfo(theGlobalSyst==sNominal);
      LooperFunctionHelpers::setKeepLHEGenPartInfo(theGlobalSyst==sNominal);
    }
    std::unordered_map<SystematicsHelpers::SystematicVariationTypes, double> globalWeights;
    double globalWeight = xsec * xsec_scale * BR_scale * (isData ? 1.f : lumi) / sum_wgts; globalWeights[theGlobalSyst] = globalWeight;
    double globalWeight_PUDn = xsec * xsec_scale * BR_scale * (isData ? 1.f : lumi) / sum_wgts_PUDn; globalWeights[SystematicsHelpers::ePUDn] = globalWeight_PUDn;
    double globalWeight_PUUp = xsec * xsec_scale * BR_scale * (isData ? 1.f : lumi) / sum_wgts_PUUp; globalWeights[SystematicsHelpers::ePUUp] = globalWeight_PUUp;
    IVYout << "Sample " << sample_tree->sampleIdentifier << " has a gen. weight sum of " << sum_wgts << " (PU dn: " << sum_wgts_PUDn << ", PU up: " << sum_wgts_PUUp << ")." << endl;
    IVYout << "\t- Raw xsec = " << xsec << endl;
    IVYout << "\t- xsec scale * BR scale = " << xsec_scale * BR_scale << endl;
    IVYout << "\t- xsec * BR * lumi = " << xsec * xsec_scale * BR_scale * (isData ? 1.f : lumi) << endl;
    IVYout << "\t- Global weight = " << globalWeight << endl;
    IVYout << "\t- Global weight (PU dn) = " << globalWeight_PUDn << endl;
    IVYout << "\t- Global weight (PU up) = " << globalWeight_PUUp << endl;

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
