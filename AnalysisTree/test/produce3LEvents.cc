#include <cassert>
#include <chrono>
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
  bool use_MET_corrections = true;

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

  void setApplyFakeableId(bool applyFakeables);

  // Helpers to keep or discard gen. ak4jet quantities
  bool keepGenAK4JetInfo = false;
  void setKeepGenAK4JetInfo(bool keepGenAK4JetInfo_);

  // Helpers to record the LHE Higgs information
  bool keepLHEGenPartInfo = false;
  void setKeepLHEGenPartInfo(bool keepLHEGenPartInfo_);

  std::vector<std::pair<TString, std::chrono::microseconds>> type_accTime_pairs;
  void addTimeDuration(TString const& strname, std::chrono::microseconds const& dur);

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
  bool hasHLTMenuProperties = theLooper->hasHLTMenuProperties();
  if (!hasHLTMenuProperties){
    IVYerr << "LooperFunctionHelpers::looperRule: Has to have HLT menus with properties..." << endl;
    assert(0);
  }
  auto it_HLTMenuSimple_SingleLepton = triggerCheckListMap.find("SingleLepton");
  auto it_HLTMenuProps_SingleLepton = triggerPropsCheckListMap.find("SingleLepton");
  if (it_HLTMenuProps_SingleLepton == triggerPropsCheckListMap.cend()){
    IVYerr << "LooperFunctionHelpers::looperRule: The trigger type 'SingleLepton' has to be defined in this looper rule!" << endl;
    assert(0);
  }
  auto it_HLTMenuSimple_Dilepton_DF = triggerCheckListMap.find("Dilepton_DF");
  auto it_HLTMenuProps_Dilepton_DF = triggerPropsCheckListMap.find("Dilepton_DF");
  if (it_HLTMenuProps_Dilepton_DF == triggerPropsCheckListMap.cend()){
    IVYerr << "LooperFunctionHelpers::looperRule: The trigger type 'Dilepton_DF' has to be defined in this looper rule!" << endl;
    assert(0);
  }
  auto it_HLTMenuSimple_Dilepton_DF_Extra = triggerCheckListMap.find("Dilepton_DF_Extra");
  auto it_HLTMenuProps_Dilepton_DF_Extra = triggerPropsCheckListMap.find("Dilepton_DF_Extra");
  if (it_HLTMenuProps_Dilepton_DF_Extra == triggerPropsCheckListMap.cend()){
    IVYerr << "LooperFunctionHelpers::looperRule: The trigger type 'Dilepton_DF_Extra' has to be defined in this looper rule!" << endl;
    assert(0);
  }
  auto it_HLTMenuSimple_Dilepton_SF = triggerCheckListMap.find("Dilepton_SF");
  auto it_HLTMenuProps_Dilepton_SF = triggerPropsCheckListMap.find("Dilepton_SF");
  if (it_HLTMenuProps_Dilepton_SF == triggerPropsCheckListMap.cend()){
    IVYerr << "LooperFunctionHelpers::looperRule: The trigger type 'Dilepton_SF' has to be defined in this looper rule!" << endl;
    assert(0);
  }
  auto it_HLTMenuSimple_Trilepton = triggerCheckListMap.find("Trilepton");
  auto it_HLTMenuProps_Trilepton = triggerPropsCheckListMap.find("Trilepton");
  if (it_HLTMenuProps_Trilepton == triggerPropsCheckListMap.cend()){
    IVYerr << "LooperFunctionHelpers::looperRule: The trigger type 'Trilepton' has to be defined in this looper rule!" << endl;
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
  BRANCH_COMMAND(float, event_wgt_triggers_Trilepton) \
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
  BRANCH_COMMAND(cms3_id_t, dilepton_id) \
  BRANCH_COMMAND(float, dilepton_pt) \
  BRANCH_COMMAND(float, dilepton_eta) \
  BRANCH_COMMAND(float, dilepton_phi) \
  BRANCH_COMMAND(float, dilepton_mass) \
  BRANCH_COMMAND(float, ak4jets_HT) \
  BRANCH_COMMAND(float, ak4jets_MHT) \
  BRANCH_COMMAND(float, event_m3l) \
  BRANCH_COMMAND(float, dPhi_Z_W) \
  BRANCH_COMMAND(float, dPhi_lepW_pTmiss) \
  BRANCH_COMMAND(float, event_mTl) \
  BRANCH_COMMAND(float, event_mWVis) \
  BRANCH_COMMAND(float, dPhi_pTleptonsjets_pTmiss) \
  BRANCH_COMMAND(float, min_abs_dPhi_pTj_pTmiss)
#define BRANCH_VECTOR_COMMANDS \
  BRANCH_COMMAND(float, event_wgt_triggers_Dilepton_SF) \
  BRANCH_COMMAND(float, event_wgt_triggers_Dilepton_DF) \
  BRANCH_COMMAND(float, event_wgt_triggers_Dilepton_DF_Extra) \
  BRANCH_COMMAND(cms3_listSize_t, dilepton_daughter_indices) \
  BRANCH_COMMAND(float, event_wgt_triggers_SingleLepton) \
  BRANCH_COMMAND(bool, leptons_is_genMatched_prompt) \
  BRANCH_COMMAND(bool, leptons_is_fakeableBase) \
  BRANCH_COMMAND(std::vector<bool>, leptons_pass_fakeable_ids) \
  BRANCH_COMMAND(cms3_id_t, leptons_id) \
  BRANCH_COMMAND(float, leptons_pt) \
  BRANCH_COMMAND(float, leptons_eta) \
  BRANCH_COMMAND(float, leptons_phi) \
  BRANCH_COMMAND(float, leptons_mass) \
  BRANCH_COMMAND(float, leptons_relIso) \
  BRANCH_COMMAND(float, leptons_eff) \
  BRANCH_COMMAND(float, leptons_eff_StatDn) \
  BRANCH_COMMAND(float, leptons_eff_StatUp) \
  BRANCH_COMMAND(float, leptons_eff_SystDn) \
  BRANCH_COMMAND(float, leptons_eff_SystUp) \
  BRANCH_COMMAND(float, leptons_eff_DF) \
  BRANCH_COMMAND(float, leptons_eff_DF_StatDn) \
  BRANCH_COMMAND(float, leptons_eff_DF_StatUp) \
  BRANCH_COMMAND(float, leptons_eff_DF_SystDn) \
  BRANCH_COMMAND(float, leptons_eff_DF_SystUp) \
  BRANCH_COMMAND(bool, electrons_conv_vtx_flag) \
  BRANCH_COMMAND(bool, electrons_pass_tightCharge) \
  BRANCH_COMMAND(cms3_electron_missinghits_t, electrons_n_missing_inner_hits) \
  BRANCH_COMMAND(cms3_electron_missinghits_t, electrons_n_all_missing_inner_hits) \
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
  std::unordered_map<ParticleObject const*, float> lepton_eff_StatDn_map;
  std::unordered_map<ParticleObject const*, float> lepton_eff_StatUp_map;
  std::unordered_map<ParticleObject const*, float> lepton_eff_SystDn_map;
  std::unordered_map<ParticleObject const*, float> lepton_eff_SystUp_map;

  auto const& muons = muonHandler->getProducts();
  auto const& electrons = electronHandler->getProducts();
  std::vector<MuonObject*> muons_selected; muons_selected.reserve(muons.size());
  std::vector<ElectronObject*> electrons_selected; electrons_selected.reserve(electrons.size());
  unsigned int nleptons_tight = 0;

  float SF_muons = 1;
  float SF_muons_StatDn = 1;
  float SF_muons_StatUp = 1;
  float SF_muons_SystDn = 1;
  float SF_muons_SystUp = 1;
  float SF_muons_AltMCDn = 1;
  float SF_muons_AltMCUp = 1;
  for (auto const& part:muons){
    bool is_tight = ParticleSelectionHelpers::isTightParticle(part);
    bool is_fakeable = ParticleSelectionHelpers::isLooseParticle(part);

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
    else if (is_tight){
      theEff = 1; muonSFHandler->getIdIsoSFAndEff(theGlobalSyst, part->pt(), part->eta(), true, true, true, theSF, &theEff); lepton_eff_map[part] = theEff;
      theEff = 1; muonSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eMuEffStatDn, part->pt(), part->eta(), true, true, true, theSF, &theEff); lepton_eff_StatDn_map[part] = theEff;
      theEff = 1; muonSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eMuEffStatUp, part->pt(), part->eta(), true, true, true, theSF, &theEff); lepton_eff_StatUp_map[part] = theEff;
      theEff = 1; muonSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eMuEffSystDn, part->pt(), part->eta(), true, true, true, theSF, &theEff); lepton_eff_SystDn_map[part] = theEff;
      theEff = 1; muonSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eMuEffSystUp, part->pt(), part->eta(), true, true, true, theSF, &theEff); lepton_eff_SystUp_map[part] = theEff;
    }

    if (is_tight || is_fakeable){
      if (is_tight) nleptons_tight++;
      else event_n_leptons_fakeableBase++;
      muons_selected.push_back(part);
    }
  }
  event_wgt_SFs_muons = SF_muons;
  event_wgt_SFs_muons_StatDn = SF_muons_StatDn;
  event_wgt_SFs_muons_StatUp = SF_muons_StatUp;
  event_wgt_SFs_muons_SystDn = SF_muons_SystDn;
  event_wgt_SFs_muons_SystUp = SF_muons_SystUp;
  event_wgt_SFs_muons_AltMCDn = SF_muons_AltMCDn;
  event_wgt_SFs_muons_AltMCUp = SF_muons_AltMCUp;

  float SF_electrons = 1;
  float SF_electrons_StatDn = 1;
  float SF_electrons_StatUp = 1;
  float SF_electrons_SystDn = 1;
  float SF_electrons_SystUp = 1;
  float SF_electrons_AltMCDn = 1;
  float SF_electrons_AltMCUp = 1;
  for (auto const& part:electrons){
    bool is_tight = ParticleSelectionHelpers::isTightParticle(part);
    bool is_fakeable = ParticleSelectionHelpers::isLooseParticle(part);

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
    else if (is_tight){
      theEff = 1; electronSFHandler->getIdIsoSFAndEff(theGlobalSyst, part->pt(), part->etaSC(), part->isGap(), true, true, true, theSF, &theEff); lepton_eff_map[part] = theEff;
      theEff = 1; electronSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eEleEffStatDn, part->pt(), part->etaSC(), part->isGap(), true, true, true, theSF, &theEff); lepton_eff_StatDn_map[part] = theEff;
      theEff = 1; electronSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eEleEffStatUp, part->pt(), part->etaSC(), part->isGap(), true, true, true, theSF, &theEff); lepton_eff_StatUp_map[part] = theEff;
      theEff = 1; electronSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eEleEffSystDn, part->pt(), part->etaSC(), part->isGap(), true, true, true, theSF, &theEff); lepton_eff_SystDn_map[part] = theEff;
      theEff = 1; electronSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eEleEffSystUp, part->pt(), part->etaSC(), part->isGap(), true, true, true, theSF, &theEff); lepton_eff_SystUp_map[part] = theEff;
    }

    if (is_tight || is_fakeable){
      if (is_tight) nleptons_tight++;
      else event_n_leptons_fakeableBase++;
      electrons_selected.push_back(part);
    }
  }
  event_wgt_SFs_electrons = SF_electrons;
  event_wgt_SFs_electrons_StatDn = SF_electrons_StatDn;
  event_wgt_SFs_electrons_StatUp = SF_electrons_StatUp;
  event_wgt_SFs_electrons_SystDn = SF_electrons_SystDn;
  event_wgt_SFs_electrons_SystUp = SF_electrons_SystUp;
  event_wgt_SFs_electrons_AltMCDn = SF_electrons_AltMCDn;
  event_wgt_SFs_electrons_AltMCUp = SF_electrons_AltMCUp;
  if (nleptons_tight<1 || nleptons_tight>3) return false;
  if ((nleptons_tight + event_n_leptons_fakeableBase)!=3) return false;
  theLooper->incrementSelection("3P+2P1F+1P2F leptons");

  std::vector<ParticleObject*> leptons_selected; leptons_selected.reserve(muons_selected.size() + electrons_selected.size());
  for (auto const& part:muons_selected) leptons_selected.push_back(part);
  for (auto const& part:electrons_selected) leptons_selected.push_back(part);
  ParticleObjectHelpers::sortByGreaterPt(leptons_selected);
  assert(leptons_selected.size()==3);
  if (!OffshellCutflow::check_pTl1(leptons_selected.front()->pt()) || !OffshellCutflow::check_pTl2(leptons_selected.at(1)->pt()) || !OffshellCutflow::check_pTl3(leptons_selected.back()->pt())) return false;
  theLooper->incrementSelection("Trigger efficiency plateau veto");

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

  dileptonHandler.constructDileptons(&muons_selected, &electrons_selected);
  auto const& dileptons = dileptonHandler.getProducts();

  DileptonObject* theChosenDilepton = nullptr;
  bool passMLLQCDSuppression = true;
  for (auto const& dilepton:dileptons){
    if (dilepton->isValid()){
      float tmp_mass = dilepton->m();
      float diff_mass = std::abs(tmp_mass - OffshellCutflow::MZ_VAL_CUTS);

      // Check if all ll pairs, regardless of flavor or charge, satisfy QCD suppression.
      passMLLQCDSuppression &= OffshellCutflow::check_mll_QCDsuppression(tmp_mass);

      // Check if the dilepton is the best Z->ll by mass
      if (OffshellCutflow::check_mll(tmp_mass, true)){
        if (!theChosenDilepton || std::abs(theChosenDilepton->m() - OffshellCutflow::MZ_VAL_CUTS)>diff_mass){
          if (!theChosenDilepton || (dilepton->isOS() && dilepton->isSF()) || !(theChosenDilepton && theChosenDilepton->isOS() && theChosenDilepton->isSF())) theChosenDilepton = dilepton;
        }
      }
    }
  }
  if (!passMLLQCDSuppression) return false;
  theLooper->incrementSelection("QCD veto");

  if (!theChosenDilepton) return false;
  theLooper->incrementSelection("Dilepton pair close to the Z mass");

  dilepton_id = theChosenDilepton->getDaughter(0)->pdgId()*theChosenDilepton->getDaughter(1)->pdgId();
  dilepton_pt = theChosenDilepton->pt();
  dilepton_eta = theChosenDilepton->eta();
  dilepton_phi = theChosenDilepton->phi();
  dilepton_mass = theChosenDilepton->m();
  dilepton_daughter_indices.reserve(2);
  for (unsigned short idau=0; idau<2; idau++){
    for (cms3_listSize_t idx=0; idx<leptons_selected.size(); idx++){
      if (leptons_selected.at(idx) == theChosenDilepton->getDaughter(idau)){
        dilepton_daughter_indices.push_back(idx);
        break;
      }
    }
  }
  assert(dilepton_daughter_indices.size()==2);
  ParticleObject* lepton_W = leptons_selected.at(3-dilepton_daughter_indices.front()-dilepton_daughter_indices.back());

  jetHandler->constructJetMET(simEventHandler, theGlobalSyst, &muons, &electrons, &photons, &pfcandidates);
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

  float event_wgt_Dilepton_SingleLepton_Combined = 0.f;
  for (auto it_part_i=leptons_selected.begin(); it_part_i!=leptons_selected.end(); it_part_i++){
    ParticleObject* part_i = *it_part_i;
    MuonObject* theMuon_i = dynamic_cast<MuonObject*>(part_i);
    ElectronObject* theElectron_i = dynamic_cast<ElectronObject*>(part_i);

    {
      std::vector<MuonObject*> muons_trigger; muons_trigger.reserve(1);
      std::vector<ElectronObject*> electrons_trigger; electrons_trigger.reserve(1);
      if (theMuon_i) muons_trigger.push_back(theMuon_i);
      if (theElectron_i) electrons_trigger.push_back(theElectron_i);

      event_wgt_triggers_SingleLepton.push_back(
        eventFilter->getTriggerWeight(
          it_HLTMenuProps_SingleLepton->second,
          &muons_trigger, &electrons_trigger, nullptr, nullptr, nullptr, nullptr
        )
      ); event_wgt_Dilepton_SingleLepton_Combined += event_wgt_triggers_SingleLepton.back();
    }

    for (auto it_part_j=it_part_i; it_part_j!=leptons_selected.end(); it_part_j++){
      if (it_part_i == it_part_j) continue;

      ParticleObject* part_j = *it_part_j;
      MuonObject* theMuon_j = dynamic_cast<MuonObject*>(part_j);
      ElectronObject* theElectron_j = dynamic_cast<ElectronObject*>(part_j);

      std::vector<MuonObject*> muons_trigger; muons_trigger.reserve(2);
      std::vector<ElectronObject*> electrons_trigger; electrons_trigger.reserve(2);
      if (theMuon_i) muons_trigger.push_back(theMuon_i);
      if (theMuon_j) muons_trigger.push_back(theMuon_j);
      if (theElectron_i) electrons_trigger.push_back(theElectron_i);
      if (theElectron_j) electrons_trigger.push_back(theElectron_j);

      event_wgt_triggers_Dilepton_SF.push_back(
        eventFilter->getTriggerWeight(
          it_HLTMenuProps_Dilepton_SF->second,
          &muons_trigger, &electrons_trigger, nullptr, nullptr, nullptr, nullptr
        )
      ); event_wgt_Dilepton_SingleLepton_Combined += event_wgt_triggers_Dilepton_SF.back();
      event_wgt_triggers_Dilepton_DF.push_back(
        eventFilter->getTriggerWeight(
          it_HLTMenuProps_Dilepton_DF->second,
          &muons_trigger, &electrons_trigger, nullptr, nullptr, nullptr, nullptr
        )
      ); event_wgt_Dilepton_SingleLepton_Combined += event_wgt_triggers_Dilepton_DF.back();
      event_wgt_triggers_Dilepton_DF_Extra.push_back(
        eventFilter->getTriggerWeight(
          it_HLTMenuProps_Dilepton_DF_Extra->second,
          &muons_trigger, &electrons_trigger, nullptr, nullptr, nullptr, nullptr
        )
      ); event_wgt_Dilepton_SingleLepton_Combined += event_wgt_triggers_Dilepton_DF_Extra.back();
    }
  }
  event_wgt_triggers_Trilepton = eventFilter->getTriggerWeight(
    it_HLTMenuProps_Trilepton->second,
    &muons, &electrons, nullptr, nullptr, nullptr, nullptr
  );
  if ((event_wgt_Dilepton_SingleLepton_Combined + event_wgt_triggers_Trilepton) == 0.f) return false;
  theLooper->incrementSelection("Trigger");

  // Test HEM filter
  if (!eventFilter->test2018HEMFilter(simEventHandler, nullptr, nullptr, &ak4jets, &ak8jets)) return false;
  if (!eventFilter->testNoisyJetFilter(simEventHandler, ak4jets)) return false;
  theLooper->incrementSelection("HEM15/16 and noisy jet vetos");

  // Fill leptons after trigger checks
  ParticleObject::LorentzVector_t sump4_leptons;
  for (auto const& part:leptons_selected){
    sump4_leptons = sump4_leptons + part->p4();

    float lepton_pt = part->pt();

    MuonObject* theMuon = dynamic_cast<MuonObject*>(part);
    ElectronObject* theElectron = dynamic_cast<ElectronObject*>(part);

    bool is_genMatched_prompt = (theMuon ? theMuon->extras.is_genMatched_prompt : theElectron->extras.is_genMatched_prompt);
    bool is_fakeableBase = !ParticleSelectionHelpers::isTightParticle(part);
    leptons_is_genMatched_prompt.push_back(is_genMatched_prompt);
    leptons_is_fakeableBase.push_back(is_fakeableBase);
    leptons_id.push_back(part->pdgId());
    leptons_pt.push_back(lepton_pt);
    leptons_eta.push_back(part->eta());
    leptons_phi.push_back(part->phi());
    leptons_mass.push_back(part->m());

    // Fill efficiency for the lepton
    auto it_eff = lepton_eff_map.find(part);
    leptons_eff.push_back((it_eff==lepton_eff_map.end() ? 1 : it_eff->second));
    it_eff = lepton_eff_StatDn_map.find(part);
    leptons_eff_StatDn.push_back((it_eff==lepton_eff_StatDn_map.end() ? 1 : it_eff->second));
    it_eff = lepton_eff_StatUp_map.find(part);
    leptons_eff_StatUp.push_back((it_eff==lepton_eff_StatUp_map.end() ? 1 : it_eff->second));
    it_eff = lepton_eff_SystDn_map.find(part);
    leptons_eff_SystDn.push_back((it_eff==lepton_eff_SystDn_map.end() ? 1 : it_eff->second));
    it_eff = lepton_eff_SystUp_map.find(part);
    leptons_eff_SystUp.push_back((it_eff==lepton_eff_SystUp_map.end() ? 1 : it_eff->second));

    {
      float lepton_relIso = 0;
      std::vector<bool> lepton_pass_fakeable_ids;
      float eff_DF = 1, eff_DF_StatDn = 1, eff_DF_StatUp = 1, eff_DF_SystDn = 1, eff_DF_SystUp = 1, SF_dummy = 1;
      if (theMuon){
        lepton_relIso = MuonSelectionHelpers::computeIso(*theMuon);

        electronSFHandler->getIdIsoSFAndEff(theGlobalSyst, part->pt(), part->eta(), 2, true, true, true, SF_dummy, &eff_DF);
        electronSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eEleEffStatDn, part->pt(), part->eta(), 2, true, true, true, SF_dummy, &eff_DF_StatDn);
        electronSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eEleEffStatUp, part->pt(), part->eta(), 2, true, true, true, SF_dummy, &eff_DF_StatUp);
        electronSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eEleEffSystDn, part->pt(), part->eta(), 2, true, true, true, SF_dummy, &eff_DF_SystDn);
        electronSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eEleEffSystUp, part->pt(), part->eta(), 2, true, true, true, SF_dummy, &eff_DF_SystUp);

        // Fakeable ids
        bool const isFakeableBaseWithTrkIso03 = theMuon->extras.trkIso03_trackerSumPt<0.4f*lepton_pt; // i.e. TrkIsoVVL
        lepton_pass_fakeable_ids.reserve(1);
        lepton_pass_fakeable_ids.push_back(isFakeableBaseWithTrkIso03);
      }
      else{
        lepton_relIso = ElectronSelectionHelpers::computeIso(*theElectron);

        using namespace ElectronTriggerCutEnums;

        muonSFHandler->getIdIsoSFAndEff(theGlobalSyst, part->pt(), theElectron->etaSC(), true, true, true, SF_dummy, &eff_DF);
        muonSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eMuEffStatDn, part->pt(), theElectron->etaSC(), true, true, true, SF_dummy, &eff_DF_StatDn);
        muonSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eMuEffStatUp, part->pt(), theElectron->etaSC(), true, true, true, SF_dummy, &eff_DF_StatUp);
        muonSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eMuEffSystDn, part->pt(), theElectron->etaSC(), true, true, true, SF_dummy, &eff_DF_SystDn);
        muonSFHandler->getIdIsoSFAndEff(SystematicsHelpers::eMuEffSystUp, part->pt(), theElectron->etaSC(), true, true, true, SF_dummy, &eff_DF_SystUp);

        // Fakeable ids
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
          IVYerr << "LooperFunctionHelpers::looperRule: Year " << SampleHelpers::getDataYear() << " is not defined for the electron trigger emulation bit selections." << endl;
          assert(0);
        }

        lepton_pass_fakeable_ids.reserve(testbits_v1.size() + testbits_v2.size());
        for (auto const& testbit:testbits_v1) lepton_pass_fakeable_ids.push_back(HelperFunctions::test_bit(trigbits_v1, testbit));
        for (auto const& testbit:testbits_v2) lepton_pass_fakeable_ids.push_back(HelperFunctions::test_bit(trigbits_v2, testbit));
      }
      leptons_relIso.push_back(lepton_relIso);

      leptons_eff_DF.push_back(eff_DF);
      leptons_eff_DF_StatDn.push_back(eff_DF_StatDn);
      leptons_eff_DF_StatUp.push_back(eff_DF_StatUp);
      leptons_eff_DF_SystDn.push_back(eff_DF_SystDn);
      leptons_eff_DF_SystUp.push_back(eff_DF_SystUp);

      leptons_pass_fakeable_ids.push_back(lepton_pass_fakeable_ids);
    }

    // Extra variables for e/gamma efficiency studies
    if (theElectron){
      auto const& extras = theElectron->extras;

      electrons_conv_vtx_flag.push_back(extras.conv_vtx_flag);
      electrons_pass_tightCharge.push_back(HelperFunctions::test_bit(extras.charge_consistency_bits, 2));
      electrons_n_missing_inner_hits.push_back(extras.n_missing_inner_hits);
      electrons_n_all_missing_inner_hits.push_back(extras.n_all_missing_inner_hits);
      electrons_full5x5_sigmaIEtaIEta.push_back(extras.full5x5_sigmaIEtaIEta);
      electrons_full5x5_sigmaIPhiIPhi.push_back(extras.full5x5_sigmaIPhiIPhi);
      electrons_full5x5_r9.push_back(extras.full5x5_r9);
      electrons_seedTime.push_back(extras.seedTime);
    }
  }

  ParticleObject::LorentzVector_t sump4_ak4jets(0, 0, 0, 0);
  std::vector<AK4JetObject*> ak4jets_tight; ak4jets_tight.reserve(ak4jets.size());
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
  
  ParticleObject::LorentzVector_t sump4_leptons_ak4jets = sump4_leptons + sump4_ak4jets;
  ParticleObject::LorentzVector_t p4_W = lepton_W->p4() + event_met_p4;
  HelperFunctions::deltaPhi(float(sump4_leptons_ak4jets.Phi()), event_phimiss, dPhi_pTleptonsjets_pTmiss);
  HelperFunctions::deltaPhi(float(theChosenDilepton->phi()), float(p4_W.Phi()), dPhi_Z_W);
  HelperFunctions::deltaPhi(float(lepton_W->phi()), event_phimiss, dPhi_lepW_pTmiss);
  event_m3l = sump4_leptons.M();
  event_mTl = std::sqrt(std::pow(lepton_W->pt() + event_pTmiss, 2) - std::pow(p4_W.Pt(), 2));
  event_mWVis = p4_W.M();

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
#define BRANCH_COMMAND(TYPE, NAME) if (!HelperFunctions::checkListVariable<TString>(bnames_exclude, #NAME)) commonEntry.setNamedVal(#NAME, NAME);
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

void LooperFunctionHelpers::setKeepGenAK4JetInfo(bool keepGenAK4JetInfo_){ keepGenAK4JetInfo = keepGenAK4JetInfo_; }

void LooperFunctionHelpers::setKeepLHEGenPartInfo(bool keepLHEGenPartInfo_){ keepLHEGenPartInfo = keepLHEGenPartInfo_; }

void LooperFunctionHelpers::addTimeDuration(TString const& strname, std::chrono::microseconds const& dur){
  for (auto& pp:type_accTime_pairs){
    if (pp.first==strname){
      pp.second = std::chrono::duration_cast<std::chrono::microseconds>(pp.second + dur);
      return;
    }
  }
  type_accTime_pairs.emplace_back(strname, dur);
}

void LooperFunctionHelpers::setApplyFakeableId(bool applyFakeables){
  MuonSelectionHelpers::setAllowFakeableInLooseSelection(applyFakeables);
  MuonSelectionHelpers::doRequireTrackerIsolationInFakeable(-1); // Set rel. trk. iso. to -1 so that we can test it separately.
  ElectronSelectionHelpers::setAllowFakeableInLooseSelection(applyFakeables);
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

  LooperFunctionHelpers::setApplyFakeableId(true); // Always apply fakeable id in loose selection/object cleaning
  OffshellCutflow::setActiveFinalState(OffshellCutflow::fs_ZW_3l1nu); // Set selection scheme

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
  std::vector<TriggerHelpers::TriggerType> requiredTriggers_Trilepton{
    TriggerHelpers::kTripleLep
  };
  auto triggerPropsCheckList_SingleLepton = TriggerHelpers::getHLTMenuProperties(requiredTriggers_SingleLepton);
  auto triggerPropsCheckList_Dilepton_DF = TriggerHelpers::getHLTMenuProperties(requiredTriggers_Dilepton_DF);
  auto triggerPropsCheckList_Dilepton_DF_Extra = TriggerHelpers::getHLTMenuProperties(requiredTriggers_Dilepton_DF_Extra);
  auto triggerPropsCheckList_Dilepton_SF = TriggerHelpers::getHLTMenuProperties(requiredTriggers_Dilepton_SF);
  auto triggerPropsCheckList_Trilepton = TriggerHelpers::getHLTMenuProperties(requiredTriggers_Trilepton);

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
    "output/3LEvents/SkimTrees/" + strdate
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
  theLooper.addHLTMenu("Dilepton_DF", triggerPropsCheckList_Dilepton_DF);
  theLooper.addHLTMenu("Dilepton_DF_Extra", triggerPropsCheckList_Dilepton_DF_Extra);
  theLooper.addHLTMenu("Dilepton_SF", triggerPropsCheckList_Dilepton_SF);
  theLooper.addHLTMenu("Trilepton", triggerPropsCheckList_Trilepton);

  curdir->cd();

  std::vector<BaseTree*> sample_trees; sample_trees.reserve(sampledirs.size());
  for (auto const& sname:sampledirs){
    bool const include_SingleLepton = !sname.Contains("ZZTo4L");
    std::vector<TString> inskimtrees;
    if (include_SingleLepton) inskimtrees = std::vector<TString>{
      "cms3ntuple/Dilepton", // 3P + 2P1F
      "cms3ntuple/SingleLepton", // Needed for 1P2F
      "cms3ntuple/Dilepton_Control" // Needed for residual events from systematics migrations
    };
    else inskimtrees = std::vector<TString>{
      "cms3ntuple/Dilepton",
      "cms3ntuple/Dilepton_Control"
    };

    TString strdsetfname = SampleHelpers::getDatasetFileName(sname);
    IVYout << "=> Accessing the input trees from " << strdsetfname << "..." << endl;
    BaseTree* sample_tree;
    if (useSkims) sample_tree = new BaseTree(strdsetfname, inskimtrees, "");
    else sample_tree = new BaseTree(strdsetfname, "cms3ntuple/Events", "", "");
    sample_trees.push_back(sample_tree);
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

  for (auto const& pp:LooperFunctionHelpers::type_accTime_pairs) IVYout << pp.first << " duration: " << pp.second.count() << endl;

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
