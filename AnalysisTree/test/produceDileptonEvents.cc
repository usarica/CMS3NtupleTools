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
#define SCALEFACTOR_HANDLER_DIRECTIVES \
  HANDLER_DIRECTIVE(MuonScaleFactorHandler, muonSFHandler) \
  HANDLER_DIRECTIVE(ElectronScaleFactorHandler, electronSFHandler) \
  HANDLER_DIRECTIVE(PhotonScaleFactorHandler, photonSFHandler) \
  HANDLER_DIRECTIVE(BtagScaleFactorHandler, btagSFHandler) \
  HANDLER_DIRECTIVE(METCorrectionHandler, metCorrectionHandler)

  // Get the current tree
  BaseTree* currentTree = theLooper->getWrappedTree();
  if (!currentTree) return false;

  // Acquire global variables
  SystematicsHelpers::SystematicVariationTypes const& theGlobalSyst = theLooper->getSystematic();
  ParticleDisambiguator& particleDisambiguator = theLooper->getParticleDisambiguator();
  DileptonHandler& dileptonHandler = theLooper->getDileptonHandler();

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
  auto it_HLTMenuSimple_OSDF = triggerCheckListMap.find("Dilepton_OSDF");
  auto it_HLTMenuProps_OSDF = triggerPropsCheckListMap.find("Dilepton_OSDF");
  if (
    (hasSimpleHLTMenus && it_HLTMenuSimple_OSDF == triggerCheckListMap.cend())
    ||
    (hasHLTMenuProperties && it_HLTMenuProps_OSDF == triggerPropsCheckListMap.cend())
    ){
    MELAerr << "LooperFunctionHelpers::looperRule: The trigger type 'Dilepton_OSDF' has to be defined in this looper rule!" << endl;
    assert(0);
  }
  auto it_HLTMenuSimple_OSSF = triggerCheckListMap.find("Dilepton_OSSF");
  auto it_HLTMenuProps_OSSF = triggerPropsCheckListMap.find("Dilepton_OSSF");
  if (
    (hasSimpleHLTMenus && it_HLTMenuSimple_OSSF == triggerCheckListMap.cend())
    ||
    (hasHLTMenuProperties && it_HLTMenuProps_OSSF == triggerPropsCheckListMap.cend())
    ){
    MELAerr << "LooperFunctionHelpers::looperRule: The trigger type 'Dilepton_OSSF' has to be defined in this looper rule!" << endl;
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
  if (!isData){
    for (auto const& handler:theLooper->getSFHandlers()){
      SCALEFACTOR_HANDLER_DIRECTIVES;
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
  BRANCH_COMMAND(float, event_wgt_triggers) \
  BRANCH_COMMAND(float, event_wgt_SFs) \
  BRANCH_COMMAND(float, event_pTmiss) \
  BRANCH_COMMAND(float, event_phimiss) \
  BRANCH_COMMAND(bool, event_passTightMETFilters) \
  BRANCH_COMMAND(float, genmet_pTmiss) \
  BRANCH_COMMAND(float, genmet_phimiss) \
  BRANCH_COMMAND(unsigned int, event_n_vtxs_good) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt30) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt20) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt20_btagged_loose) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt20_btagged_medium) \
  BRANCH_COMMAND(float, ak4jets_HT) \
  BRANCH_COMMAND(cms3_id_t, dilepton_id) \
  BRANCH_COMMAND(float, dilepton_pt) \
  BRANCH_COMMAND(float, dilepton_eta) \
  BRANCH_COMMAND(float, dilepton_phi) \
  BRANCH_COMMAND(float, dilepton_mass) \
  BRANCH_COMMAND(float, dPhi_pTll_pTmiss) \
  BRANCH_COMMAND(float, dPhi_pTlljets_pTmiss) \
  BRANCH_COMMAND(float, min_abs_dPhi_pTj_pTmiss)
#define BRANCH_VECTOR_COMMANDS \
  BRANCH_COMMAND(cms3_id_t, leptons_id) \
  BRANCH_COMMAND(float, leptons_pt) \
  BRANCH_COMMAND(float, leptons_eta) \
  BRANCH_COMMAND(float, leptons_phi) \
  BRANCH_COMMAND(float, leptons_mass) \
  BRANCH_COMMAND(float, ak4jets_pt) \
  BRANCH_COMMAND(float, ak4jets_eta) \
  BRANCH_COMMAND(float, ak4jets_phi) \
  BRANCH_COMMAND(float, ak4jets_mass)
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

  if (!isData){
    genInfoHandler->constructGenInfo(theGlobalSyst);
    auto const& genInfo = genInfoHandler->getGenInfo();
    event_wgt *= genInfo->getGenWeight(true);
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
  }

  vertexHandler->constructVertices();
  if (!vertexHandler->hasGoodPrimaryVertex()) return false;
  event_n_vtxs_good = vertexHandler->getNGoodVertices();

  eventFilter->constructFilters(simEventHandler);
  if (isData && !eventFilter->isUniqueDataEvent()) return false;

  if (!eventFilter->passCommonSkim() || !eventFilter->passMETFilters(EventFilterHandler::kMETFilters_Standard)) return false;
  event_passTightMETFilters = eventFilter->passMETFilters(EventFilterHandler::kMETFilters_Tight);

  pfcandidateHandler->constructPFCandidates(theGlobalSyst);
  auto const& pfcandidates = pfcandidateHandler->getProducts();

  muonHandler->constructMuons(theGlobalSyst, &pfcandidates);
  electronHandler->constructElectrons(theGlobalSyst, &pfcandidates);
  photonHandler->constructPhotons(theGlobalSyst, &pfcandidates);
  particleDisambiguator.disambiguateParticles(muonHandler, electronHandler, photonHandler);

  auto const& muons = muonHandler->getProducts();
  float SF_muons = 1;
  for (auto const& part:muons){
    float theSF = 1;
    if (!isData) muonSFHandler->getIdIsoSFAndEff(theGlobalSyst, part, theSF, nullptr);
    if (theSF == 0.f) continue;
    SF_muons *= theSF;
  }

  auto const& electrons = electronHandler->getProducts();
  float SF_electrons = 1;
  for (auto const& part:electrons){
    float theSF = 1;
    //if (!isData) electronSFHandler.getIdIsoSFAndEff(theGlobalSyst, part, theSF, nullptr);
    if (theSF == 0.f) continue;
    SF_electrons *= theSF;
  }

  auto const& photons = photonHandler->getProducts();
  unsigned int n_photons_tight = 0;
  float SF_photons = 1;
  for (auto const& part:photons){
    float theSF = 1;
    if (!isData) photonSFHandler->getIdIsoSFAndEff(theGlobalSyst, part, theSF, nullptr);
    if (theSF == 0.f) continue;
    SF_photons *= theSF;

    if (ParticleSelectionHelpers::isTightParticle(part)) n_photons_tight++;
  }
  if (n_photons_tight!=0) return false;

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
  if (!theChosenDilepton || nTightDilep>1) return false;

  bool const dilepton_is_SF = theChosenDilepton->isSF();
  dilepton_id = theChosenDilepton->getDaughter(0)->pdgId()*theChosenDilepton->getDaughter(1)->pdgId();
  dilepton_pt = theChosenDilepton->pt();
  dilepton_eta = theChosenDilepton->eta();
  dilepton_phi = theChosenDilepton->phi();
  dilepton_mass = theChosenDilepton->m();
  for (auto const& dau:theChosenDilepton->getDaughters()){
    leptons_id.push_back(dau->pdgId());
    leptons_pt.push_back(dau->pt());
    leptons_eta.push_back(dau->eta());
    leptons_phi.push_back(dau->phi());
  }

  jetHandler->constructJetMET(theGlobalSyst, &muons, &electrons, &photons, &pfcandidates);
  auto const& ak4jets = jetHandler->getAK4Jets();
  auto const& ak8jets = jetHandler->getAK8Jets();

  if (hasSimpleHLTMenus) event_wgt_triggers = eventFilter->getTriggerWeight((dilepton_is_SF ? it_HLTMenuSimple_OSSF : it_HLTMenuSimple_OSDF)->second);
  else if (hasHLTMenuProperties) event_wgt_triggers = eventFilter->getTriggerWeight(
    (dilepton_is_SF ? it_HLTMenuProps_OSSF : it_HLTMenuProps_OSDF)->second,
    &muons, &electrons, nullptr, nullptr, nullptr, nullptr
  );
  if (event_wgt_triggers == 0.f) return false;

  // Test HEM filter
  if (!eventFilter->test2018HEMFilter(simEventHandler, nullptr, nullptr, &ak4jets, &ak8jets)) return false;

  ParticleObject::LorentzVector_t sump4_ak4jets(0, 0, 0, 0);
  std::vector<AK4JetObject*> ak4jets_tight; ak4jets_tight.reserve(ak4jets.size());
  unsigned int n_ak4jets_tight_pt30_btagged_loose = 0;
  float SF_btagging = 1;
  for (auto const& jet:ak4jets){
    float theSF = 1;
    if (!isData) btagSFHandler->getSFAndEff(theGlobalSyst, jet, theSF, nullptr);
    if (theSF != 0.f) SF_btagging *= theSF;

    if (ParticleSelectionHelpers::isTightJet(jet)){
      ak4jets_tight.push_back(jet);
      if (jet->getBtagValue()>=btag_thr_loose) n_ak4jets_tight_pt30_btagged_loose++;

      ak4jets_HT += jet->pt();
      sump4_ak4jets += jet->p4();
    }

    // Check for pT>20 objects and count them
    if (
      jet->testSelectionBit(AK4JetSelectionHelpers::kTightId)
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
  if (n_ak4jets_tight_pt30_btagged_loose>0) return false;
  event_n_ak4jets_pt30 = ak4jets_tight.size();

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

  min_abs_dPhi_pTj_pTmiss = TMath::Pi();
  for (auto const& jet:ak4jets_tight){
    ak4jets_pt.push_back(jet->pt());
    ak4jets_eta.push_back(jet->eta());
    ak4jets_phi.push_back(jet->phi());
    ak4jets_mass.push_back(jet->mass());

    float dphi_tmp;
    HelperFunctions::deltaPhi(float(jet->phi()), event_phimiss, dphi_tmp); dphi_tmp = std::abs(dphi_tmp);
    min_abs_dPhi_pTj_pTmiss = std::min(min_abs_dPhi_pTj_pTmiss, dphi_tmp);
  }

  dPhi_pTll_pTmiss = theChosenDilepton->deltaPhi(event_phimiss);
  HelperFunctions::deltaPhi(float((theChosenDilepton->p4()+sump4_ak4jets).Phi()), event_phimiss, dPhi_pTlljets_pTmiss);

  // Set the collection of SFs at the last step
  event_wgt_SFs = SF_muons*SF_electrons*SF_photons*SF_btagging;

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
}


using namespace SystematicsHelpers;
void getTrees(
  TString strSampleSet, TString period,
  TString prodVersion, TString strdate,
  int ichunk, int nchunks,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal,
  // Jet ID options
  bool applyPUIdToAK4Jets=true, bool applyTightLeptonVetoIdToAK4Jets=false, bool useJetOverlapStripping=true,
  // MET options
  bool use_MET_Puppi=false,
  bool use_MET_XYCorr=true, bool use_MET_JERCorr=true, bool use_MET_ParticleMomCorr=true, bool use_MET_p4Preservation=true, bool use_MET_corrections=true
){
  if (nchunks==1){ nchunks = 0; ichunk=0; }
  if (nchunks>0 && (ichunk<0 || ichunk==nchunks)) return;

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SampleHelpers::configure(period, "hadoop_skims:"+prodVersion);

  const float lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::theDataPeriod);

  std::vector<TString> const validDataPeriods = SampleHelpers::getValidDataPeriods();
  size_t const nValidDataPeriods = validDataPeriods.size();

  BtagHelpers::setBtagWPType(BtagHelpers::kDeepFlav_Loose);
  LooperFunctionHelpers::setBtagWPs();

  std::vector<TriggerHelpers::TriggerType> requiredTriggers_OSDF{
    TriggerHelpers::kMuEle,
    TriggerHelpers::kSingleMu, TriggerHelpers::kSingleMu_HighPt,
    TriggerHelpers::kSingleEle, TriggerHelpers::kSingleEle_HighPt
  };
  std::vector<TriggerHelpers::TriggerType> requiredTriggers_OSSF{
    TriggerHelpers::kDoubleMu,
    TriggerHelpers::kDoubleEle, TriggerHelpers::kDoubleEle_HighPt,
    TriggerHelpers::kSingleMu, TriggerHelpers::kSingleMu_HighPt,
    TriggerHelpers::kSingleEle, TriggerHelpers::kSingleEle_HighPt
  };
  auto triggerPropsCheckList_OSDF = TriggerHelpers::getHLTMenuProperties(requiredTriggers_OSDF);
  auto triggerPropsCheckList_OSSF = TriggerHelpers::getHLTMenuProperties(requiredTriggers_OSSF);

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
  theLooper.addSFHandler(&btagSFHandler);
  theLooper.addSFHandler(&metCorrectionHandler);
  // Set output tree
  theLooper.addOutputTree(tout);
  // Register the HLT menus
  theLooper.addHLTMenu("Dilepton_OSDF", triggerPropsCheckList_OSDF);
  theLooper.addHLTMenu("Dilepton_OSSF", triggerPropsCheckList_OSSF);

  curdir->cd();

  std::vector<BaseTree*> sample_trees; sample_trees.reserve(sampledirs.size());
  for (auto const& sname:sampledirs){
    BaseTree* sample_tree = new BaseTree(SampleHelpers::getDatasetFileName(sname), "cms3ntuple/Dilepton", "", ""); sample_trees.push_back(sample_tree);
    sample_tree->sampleIdentifier = SampleHelpers::getSampleIdentifier(sname);

    std::vector<TString> allbranchnames; sample_tree->getValidBranchNamesWithoutAlias(allbranchnames, false);

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
      if (!hasCounters){
        MELAerr << "Please use skim ntuples!" << endl;
        assert(0);
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
