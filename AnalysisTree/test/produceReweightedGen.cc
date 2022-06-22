#include <cassert>
#include "common_includes.h"
#include "OffshellCutflow.h"
#include <IvyFramework/IvyAutoMELA/interface/IvyMELAHelpers.h>
#include <MelaAnalytics/CandidateLOCaster/interface/MELACandidateRecaster.h>
#include <MelaAnalytics/EventContainer/interface/HiggsComparators.h>
#include <MelaAnalytics/EventContainer/interface/TopComparators.h>
#include "TStyle.h"


// Recorded variables
#define BRANCH_SCALAR_COMMANDS \
  BRANCH_COMMAND(float, event_wgt) \
  BRANCH_COMMAND(float, event_wgt_adjustment_NNPDF30) \
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
  BRANCH_COMMAND(float, event_wgt_adjustment_NNPDF30_AsMZDn) \
  BRANCH_COMMAND(float, event_wgt_adjustment_NNPDF30_AsMZUp) \
  BRANCH_COMMAND(float, event_wgt_adjustment_NNPDF30_PDFReplicaDn) \
  BRANCH_COMMAND(float, event_wgt_adjustment_NNPDF30_PDFReplicaUp) \
  /* External adjustment factors from the registered histograms */ \
  BRANCH_COMMAND(float, event_wgt_adjustment_ext_HardJetsDn) \
  BRANCH_COMMAND(float, event_wgt_adjustment_ext_HardJetsUp) \
  BRANCH_COMMAND(float, event_wgt_adjustment_ext_PythiaScaleDn) \
  BRANCH_COMMAND(float, event_wgt_adjustment_ext_PythiaScaleUp) \
  BRANCH_COMMAND(float, event_wgt_adjustment_ext_PDFScaleDn) \
  BRANCH_COMMAND(float, event_wgt_adjustment_ext_PDFScaleUp) \
  BRANCH_COMMAND(float, event_wgt_adjustment_ext_QCDScaleDn) \
  BRANCH_COMMAND(float, event_wgt_adjustment_ext_QCDScaleUp) \
  BRANCH_COMMAND(float, event_wgt_adjustment_ext_AsMZDn) \
  BRANCH_COMMAND(float, event_wgt_adjustment_ext_AsMZUp) \
  BRANCH_COMMAND(float, event_wgt_adjustment_ext_PDFReplicaDn) \
  BRANCH_COMMAND(float, event_wgt_adjustment_ext_PDFReplicaUp) \
  /* ME adjustment weights */ \
  BRANCH_COMMAND(bool, invalidReweightingWgts) \
  BRANCH_COMMAND(float, sample_wgt) \
  BRANCH_COMMAND(float, sample_wgt_pairwiseComponent) \
  /* Other variables */ \
  BRANCH_COMMAND(bool, passGenCompSelection) \
  BRANCH_COMMAND(float, lheHiggs_mass) \
  BRANCH_COMMAND(float, lheHiggs_pt) \
  BRANCH_COMMAND(float, lheLeptonicDecay_pt) \
  BRANCH_COMMAND(float, lheLeptonicDecay_mass) \
  BRANCH_COMMAND(float, genLeptonicDecay_pt) \
  BRANCH_COMMAND(float, genLeptonicDecay_mass)
#define BRANCH_VECTOR_COMMANDS \
  BRANCH_COMMAND(cms3_id_t, lheparticles_id) \
  BRANCH_COMMAND(cms3_id_t, lheparticles_status) \
  BRANCH_COMMAND(cms3_id_t, lheQCDLOparticles_id) \
  BRANCH_COMMAND(cms3_id_t, lheQCDLOparticles_status) \
  BRANCH_COMMAND(cms3_id_t, genparticles_id) \
  BRANCH_COMMAND(float, lheQCDLOparticles_px) \
  BRANCH_COMMAND(float, lheQCDLOparticles_py) \
  BRANCH_COMMAND(float, lheQCDLOparticles_pz) \
  BRANCH_COMMAND(float, lheQCDLOparticles_E) \
  BRANCH_COMMAND(float, lheparticles_px) \
  BRANCH_COMMAND(float, lheparticles_py) \
  BRANCH_COMMAND(float, lheparticles_pz) \
  BRANCH_COMMAND(float, lheparticles_E) \
  BRANCH_COMMAND(float, genparticles_pt) \
  BRANCH_COMMAND(float, genparticles_eta) \
  BRANCH_COMMAND(float, genparticles_phi) \
  BRANCH_COMMAND(float, genparticles_mass) \
  BRANCH_COMMAND(float, genak4jets_pt) \
  BRANCH_COMMAND(float, genak4jets_eta) \
  BRANCH_COMMAND(float, genak4jets_phi) \
  BRANCH_COMMAND(float, genak4jets_mass)
#define BRANCH_COMMANDS \
  BRANCH_SCALAR_COMMANDS \
  BRANCH_VECTOR_COMMANDS


namespace LooperFunctionHelpers{
  using namespace std;
  using namespace IvyStreamHelpers;
  using namespace OffshellCutflow;

  bool recastLHETopology = false;
  TVar::Production candScheme = TVar::nProductions; // Recasting scheme

  // Candidate interpretation
  MELAEvent::CandidateVVMode VVMode = MELAEvent::nCandidateVVModes;
  int VVDecayMode = -1;

  std::vector<TString> selectedMEs;

  bool applyPythiaScaleExternally = false;
  std::vector<SystematicsHelpers::SystematicVariationTypes> registeredSysts;
  std::unordered_map<SystematicsHelpers::SystematicVariationTypes, TH1F*> syst_h1d_map;
  std::unordered_map<SystematicsHelpers::SystematicVariationTypes, TH2F*> syst_h2d_map;
  std::unordered_map<SystematicsHelpers::SystematicVariationTypes, TH3F*> syst_h3d_map;
  void registerSystHistogram(SystematicsHelpers::SystematicVariationTypes const& syst, TH1F* hist);
  void registerSystHistogram(SystematicsHelpers::SystematicVariationTypes const& syst, TH2F* hist);
  void registerSystHistogram(SystematicsHelpers::SystematicVariationTypes const& syst, TH3F* hist);

  double EvalSystHistogram(TH1F* hist, float const& xvar);
  double EvalSystHistogram(TH2F* hist, float const& xvar, float const& yvar);
  double EvalSystHistogram(TH3F* hist, float const& xvar, float const& yvar, float const& zvar);

  bool looperRule(BaseTreeLooper*, std::unordered_map<SystematicsHelpers::SystematicVariationTypes, double> const&, SimpleEntry&);

  void cleanup();
}
bool LooperFunctionHelpers::looperRule(BaseTreeLooper* theLooper, std::unordered_map<SystematicsHelpers::SystematicVariationTypes, double> const& extWgt, SimpleEntry& commonEntry){
  // Define handlers
#define OBJECT_HANDLER_DIRECTIVES \
  HANDLER_DIRECTIVE(GenInfoHandler, genInfoHandler)

  // Get the current tree
  BaseTree* currentTree = theLooper->getWrappedTree();
  if (!currentTree) return false;

  // Acquire global variables
  SystematicsHelpers::SystematicVariationTypes const& theGlobalSyst = theLooper->getSystematic();

  // Acquire all handlers
#define HANDLER_DIRECTIVE(TYPE, NAME) TYPE* NAME = nullptr;
  OBJECT_HANDLER_DIRECTIVES;
#undef HANDLER_DIRECTIVE
#define HANDLER_DIRECTIVE(TYPE, NAME) TYPE* tmp_##NAME = dynamic_cast<TYPE*>(handler); if (tmp_##NAME){ NAME = tmp_##NAME; continue; }
  for (auto const& handler:theLooper->getObjectHandlers()){
    OBJECT_HANDLER_DIRECTIVES;
  }
#undef HANDLER_DIRECTIVE

  BulkReweightingBuilder* rewgtBuilder = nullptr;
  {
    auto const& rewgt_map = theLooper->getRegisteredRewgtBuilders();
    auto it_rewgt_map = rewgt_map.find("MERewgt");
    if (it_rewgt_map!=rewgt_map.end()) rewgtBuilder = it_rewgt_map->second;
  }

  /************************/
  /* EVENT INTERPRETATION */
  /************************/
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
    IVYerr << "LooperFunctionHelpers::looperRule: External normalization map does not have a proper weight assigned!" << endl;
    assert(0);
  }
  double const& extWgt_central = it_extWgt->second;

  event_wgt = extWgt_central;
  // Set NNPDF 3.0 adjustment to 1
  event_wgt_adjustment_NNPDF30 = 1;

  genInfoHandler->constructGenInfo(theGlobalSyst);
  auto const& genInfo = genInfoHandler->getGenInfo();
  double genwgt_NNPDF30 = genInfo->getGenWeight(false);
  double genwgt_default = genInfo->getGenWeight(true);
  event_wgt_adjustment_NNPDF30 = (genwgt_default!=0. ? genwgt_NNPDF30 / genwgt_default : 0.);
  event_wgt *= genwgt_default;
  if (event_wgt==0.f) return false;


  event_wgt_adjustment_ext_PDFScaleDn = event_wgt_adjustment_ext_PDFScaleUp
  = event_wgt_adjustment_ext_QCDScaleDn = event_wgt_adjustment_ext_QCDScaleUp
  = event_wgt_adjustment_ext_AsMZDn = event_wgt_adjustment_ext_AsMZUp
  = event_wgt_adjustment_ext_PDFReplicaDn = event_wgt_adjustment_ext_PDFReplicaUp
  = event_wgt_adjustment_ext_PythiaScaleDn = event_wgt_adjustment_ext_PythiaScaleUp
  = event_wgt_adjustment_ext_HardJetsDn = event_wgt_adjustment_ext_HardJetsUp = 1;

  auto const& genInfoExtras = genInfo->extras;
  event_wgt_adjustment_PythiaScaleDn = (!applyPythiaScaleExternally ? genInfoExtras.PythiaWeight_isr_muR0p25 * genInfoExtras.PythiaWeight_fsr_muR0p25 : 1.f);
  event_wgt_adjustment_PythiaScaleUp = (!applyPythiaScaleExternally ? genInfoExtras.PythiaWeight_isr_muR4 * genInfoExtras.PythiaWeight_fsr_muR4 : 1.f);
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

  passGenCompSelection = true;

  // Variables for systematics reweighting
  float genpromptparticles_sump4_pt = -1;
  unsigned int n_genak4jets = 0;

  bool foundLHEHiggs = false;
  bool hasTaus = false;
  unsigned int n_leps_nus=0;
  ParticleObject::LorentzVector_t p4_lheHiggs;
  ParticleObject::LorentzVector_t p4_lheLeptonicDecay;
  auto const& lheparticles = genInfoHandler->getLHEParticles();
  auto const& genparticles = genInfoHandler->getGenParticles();
  if (!lheparticles.empty()){
    for (auto const& part:lheparticles){
      if (part->status()==2 && PDGHelpers::isAHiggs(part->pdgId())){
        p4_lheHiggs = part->p4();
        foundLHEHiggs = true;
      }
      else if (part->status()==1 && (PDGHelpers::isALepton(part->pdgId()) || PDGHelpers::isANeutrino(part->pdgId()))){
        p4_lheLeptonicDecay = p4_lheLeptonicDecay + part->p4();
        if (std::abs(part->pdgId())==15) hasTaus = true;
        n_leps_nus++;
      }
    }
  }
  else if (!genparticles.empty()){
    for (auto const& part:genparticles){
      if (!part->extras.isHardProcess) continue;
      if (PDGHelpers::isAHiggs(part->pdgId())){
        p4_lheHiggs = part->p4();
        foundLHEHiggs = true;
      }
      else if (part->status()!=21 && (PDGHelpers::isALepton(part->pdgId()) || PDGHelpers::isANeutrino(part->pdgId()))){
        p4_lheLeptonicDecay = p4_lheLeptonicDecay + part->p4();
        if (std::abs(part->pdgId())==15) hasTaus = true;
        n_leps_nus++;
      }
    }
  }
  else{
    IVYerr << "LooperFunctionHelpers::looperRule: Both gen. and LHE particle collections are empty!" << endl;
    assert(0);
  }
  passGenCompSelection &= (n_leps_nus==4);
  passGenCompSelection &= (!hasTaus);
  if (foundLHEHiggs){
    lheHiggs_mass = p4_lheHiggs.M();
    lheHiggs_pt = p4_lheHiggs.Pt();
  }
  lheLeptonicDecay_pt = p4_lheLeptonicDecay.Pt();
  lheLeptonicDecay_mass = p4_lheLeptonicDecay.M();

  {
    ParticleObject::LorentzVector_t p4_genLeptonicDecay;
    ParticleObject::LorentzVector_t genpromptparticles_sump4;
    auto const& genparticles = genInfoHandler->getGenParticles();
    for (auto const& part:genparticles){
      if (
        part->testSelectionBit(GenParticleSelectionHelpers::kHardPromptFinalVisibleParticle)
        ||
        (part->extras.isPromptFinalState && PDGHelpers::isANeutrino(part->pdgId()))
        ) genpromptparticles_sump4 += part->p4();
      if (part->extras.isPromptFinalState && (PDGHelpers::isAPhoton(part->pdgId()) || PDGHelpers::isANeutrino(part->pdgId()) || PDGHelpers::isALepton(part->pdgId()))){
        p4_genLeptonicDecay += part->p4();
        genparticles_id.push_back(part->pdgId());
        genparticles_pt.push_back(part->pt());
        genparticles_eta.push_back(part->eta());
        genparticles_phi.push_back(part->phi());
        genparticles_mass.push_back(part->mass());
      }
    }
    genLeptonicDecay_pt = p4_genLeptonicDecay.Pt();
    genLeptonicDecay_mass = p4_genLeptonicDecay.M();
    genpromptparticles_sump4_pt = genpromptparticles_sump4.Pt();
  }

  auto const& genak4jets = genInfoHandler->getGenAK4Jets();
  for (auto const& jet:genak4jets){
    n_genak4jets++;

    if (jet->pt()<20.f) continue;

    genak4jets_pt.push_back(jet->pt());
    genak4jets_eta.push_back(jet->eta());
    genak4jets_phi.push_back(jet->phi());
    genak4jets_mass.push_back(jet->mass());
  }

  {
    std::vector<MELAParticle*> particleList;
    {
      std::unordered_map<IvyParticle*, MELAParticle*> particle_translation_map;
      if (!lheparticles.empty()){
        for (auto const& part:lheparticles){
          MELAParticle* onePart = new MELAParticle(part->pdgId(), part->p4_TLV());
          onePart->setGenStatus(part->status());
          particleList.push_back(onePart);
          particle_translation_map[part] = onePart;
          if (lheHiggs_mass<0. && PDGHelpers::isAHiggs(part->pdgId())){
            lheHiggs_mass = part->mass();
            lheHiggs_pt = part->pt();
          }
        }
        // Loop to also assign mothers
        for (auto const& part:lheparticles){
          MELAParticle* const& thePart = particle_translation_map.find(part)->second;
          for (auto const& mother:part->getMothers()){
            if (!mother) continue;
            MELAParticle* const& theMother = particle_translation_map.find(mother)->second;
            thePart->addMother(theMother);
          }
        }
      }
      else{
        for (auto const& part:genparticles){
          if (!part->extras.isHardProcess) continue;
          MELAParticle* onePart = new MELAParticle(part->pdgId(), part->p4_TLV());
          int lhe_status = 1;
          if (part->status()!=2) lhe_status = PDGHelpers::convertPythiaStatus(part->status());
          onePart->setGenStatus(lhe_status);
          particleList.push_back(onePart);
          particle_translation_map[part] = onePart;
          if (lheHiggs_mass<0. && PDGHelpers::isAHiggs(part->pdgId())){
            lheHiggs_mass = part->mass();
            lheHiggs_pt = part->pt();
          }
        }
        // Loop to also assign mothers
        for (auto const& part:genparticles){
          if (!part->extras.isHardProcess) continue;
          MELAParticle* const& thePart = particle_translation_map.find(part)->second;
          for (auto const& mother:part->getMothers()){
            if (!mother) continue;
            MELAParticle* const& theMother = particle_translation_map.find(mother)->second;
            thePart->addMother(theMother);
          }
        }
      }
    }

    MELAEvent* genEvent = new MELAEvent();
    std::vector<MELAParticle*> writtenGenCands;
    std::vector<MELAParticle*> writtenGenTopCands;
    {
      using namespace PDGHelpers;

      for (MELAParticle* genPart:particleList){
        if (isATopQuark(genPart->id)){
          writtenGenTopCands.push_back(genPart);
          if (genPart->genStatus==1) genEvent->addIntermediate(genPart);
        }
        if (isAHiggs(genPart->id)){
          writtenGenCands.push_back(genPart);
          if (VVMode==MELAEvent::UndecayedMode && (genPart->genStatus==1 || genPart->genStatus==2)) genEvent->addIntermediate(genPart);
        }
        if (genPart->genStatus==1){
          if (isALepton(genPart->id)) genEvent->addLepton(genPart);
          else if (isANeutrino(genPart->id)) genEvent->addNeutrino(genPart);
          else if (isAPhoton(genPart->id)) genEvent->addPhoton(genPart);
          else if (isAKnownJet(genPart->id) && !isATopQuark(genPart->id)) genEvent->addJet(genPart);
        }
        else if (genPart->genStatus==-1) genEvent->addMother(genPart);
      }
    }

    genEvent->constructTopCandidates();
    {
      std::vector<MELATopCandidate_t*> matchedTops;
      for (auto* writtenGenTopCand:writtenGenTopCands){
        MELATopCandidate_t* tmpCand = TopComparators::matchATopToParticle(*genEvent, writtenGenTopCand);
        if (tmpCand) matchedTops.push_back(tmpCand);
      }
      for (MELATopCandidate_t* tmpCand:genEvent->getTopCandidates()){
        if (std::find(matchedTops.begin(), matchedTops.end(), tmpCand)==matchedTops.end()) tmpCand->setSelected(false);
      }
    }

    MELACandidate* genCand = nullptr;
    MELACandidate* candModified = nullptr;

    genEvent->constructVVCandidates(VVMode, VVDecayMode);
    genEvent->addVVCandidateAppendages();
    for (auto* writtenGenCand:writtenGenCands){
      MELACandidate* tmpCand = HiggsComparators::matchAHiggsToParticle(*genEvent, writtenGenCand);
      if (tmpCand){
        if (!genCand) genCand = tmpCand;
        else genCand = HiggsComparators::candComparator(genCand, tmpCand, HiggsComparators::BestZ1ThenZ2ScSumPt, VVMode);
      }
    }
    if (!genCand) genCand = HiggsComparators::candidateSelector(*genEvent, HiggsComparators::BestZ1ThenZ2ScSumPt, VVMode);

    if (lheHiggs_mass<0. && genCand){
      lheHiggs_mass = genCand->m();
      lheHiggs_pt = genCand->pt();
    }

    // Write the default LHE particles
    {
      std::vector<MELAParticle*> particles_to_write;
      for (auto const& part:genCand->getMothers()) particles_to_write.push_back(part);
      for (auto const& part:genCand->getAssociatedPhotons()) particles_to_write.push_back(part);
      // No need to return neutrinos. They are already part of the leptons collection.
      //for (auto const& part:genCand->getAssociatedNeutrinos()) particles_to_write.push_back(part);
      for (auto const& part:genCand->getAssociatedLeptons()) particles_to_write.push_back(part);
      for (auto const& part:genCand->getAssociatedJets()) particles_to_write.push_back(part);
      for (auto const& part:genCand->getSortedDaughters()) particles_to_write.push_back(part);

      for (auto const& part:particles_to_write){
        lheparticles_id.push_back(part->id);
        lheparticles_status.push_back(part->genStatus);
        lheparticles_px.push_back(part->x());
        lheparticles_py.push_back(part->y());
        lheparticles_pz.push_back(part->z());
        lheparticles_E.push_back(part->t());
      }
    }

    MELACandidate* candToWrite = genCand;
    MELACandidateRecaster* recaster = nullptr;
    if (recastLHETopology){
      recaster = new MELACandidateRecaster(candScheme);
      recaster->copyCandidate(genCand, candModified);

      if (candScheme==TVar::JJVBF) recaster->reduceJJtoQuarks(candModified);
      else if (candScheme==TVar::Had_ZH || candScheme==TVar::Had_WH){
        MELAParticle* bestAV = MELACandidateRecaster::getBestAssociatedV(genCand, candScheme);
        if (bestAV) recaster->deduceLOVHTopology(candModified);
        else{
          IVYerr << "ERROR: No associated V can be found in the VH recasting scheme." << endl;
          exit(1);
        }
      }

      candToWrite = candModified;
    }

    // Write the recasted LHE particles
    {
      std::vector<MELAParticle*> particles_to_write;
      for (auto const& part:candToWrite->getMothers()) particles_to_write.push_back(part);
      for (auto const& part:candToWrite->getAssociatedPhotons()) particles_to_write.push_back(part);
      // No need to return neutrinos. They are already part of the leptons collection.
      //for (auto const& part:candToWrite->getAssociatedNeutrinos()) particles_to_write.push_back(part);
      for (auto const& part:candToWrite->getAssociatedLeptons()) particles_to_write.push_back(part);
      for (auto const& part:candToWrite->getAssociatedJets()) particles_to_write.push_back(part);
      for (auto const& part:candToWrite->getSortedDaughters()) particles_to_write.push_back(part);

      for (auto const& part:particles_to_write){
        lheQCDLOparticles_id.push_back(part->id);
        lheQCDLOparticles_status.push_back(part->genStatus);
        lheQCDLOparticles_px.push_back(part->x());
        lheQCDLOparticles_py.push_back(part->y());
        lheQCDLOparticles_pz.push_back(part->z());
        lheQCDLOparticles_E.push_back(part->t());
      }
    }

    delete candModified;
    delete recaster;
    delete genEvent;
    for (auto& part:particleList) delete part;
  }

  sample_wgt = (rewgtBuilder ? rewgtBuilder->getOverallReweightingNormalization(currentTree) : 1.);
  invalidReweightingWgts = (rewgtBuilder ? !rewgtBuilder->checkWeightsBelowThreshold(currentTree) : false);
  // The weight below is only for bookkeeping purposes. It should not be multiplied.
  sample_wgt_pairwiseComponent = (rewgtBuilder ? rewgtBuilder->getSamplePairwiseNormalization(currentTree) : 1.);

  // Systematics reweighting
  for (auto const& syst:registeredSysts){
    using namespace SystematicsHelpers;
    double syst_corr = 1;
    if (
      syst==tPDFScaleDn || syst==tPDFScaleUp
      || syst==tQCDScaleDn || syst==tQCDScaleUp
      || syst==tAsMZDn || syst==tAsMZUp
      || syst==tPDFReplicaDn || syst==tPDFReplicaUp
      ) syst_corr = EvalSystHistogram(syst_h1d_map[syst], lheHiggs_mass);
    else if (
      syst==tPythiaTuneDn || syst==tPythiaTuneUp || syst==tPythiaScaleDn || syst==tPythiaScaleUp
      ) syst_corr = EvalSystHistogram(syst_h3d_map[syst], lheHiggs_mass, n_genak4jets, genpromptparticles_sump4_pt / lheHiggs_mass);
    else if (
      syst==tHardJetsDn || syst==tHardJetsUp
      ) syst_corr = EvalSystHistogram(syst_h3d_map[syst], lheHiggs_mass, n_genak4jets, lheHiggs_pt / lheHiggs_mass);

    switch (syst){
    case tPDFScaleDn:
      event_wgt_adjustment_ext_PDFScaleDn *= syst_corr;
      break;
    case tPDFScaleUp:
      event_wgt_adjustment_ext_PDFScaleUp *= syst_corr;
      break;
    case tQCDScaleDn:
      event_wgt_adjustment_ext_QCDScaleDn *= syst_corr;
      break;
    case tQCDScaleUp:
      event_wgt_adjustment_ext_QCDScaleUp *= syst_corr;
      break;
    case tAsMZDn:
      event_wgt_adjustment_ext_AsMZDn *= syst_corr;
      break;
    case tAsMZUp:
      event_wgt_adjustment_ext_AsMZUp *= syst_corr;
      break;
    case tPDFReplicaDn:
      event_wgt_adjustment_ext_PDFReplicaDn *= syst_corr;
      break;
    case tPDFReplicaUp:
      event_wgt_adjustment_ext_PDFReplicaUp *= syst_corr;
      break;
    case tPythiaScaleDn:
      event_wgt_adjustment_ext_PythiaScaleDn *= syst_corr;
      break;
    case tPythiaScaleUp:
      event_wgt_adjustment_ext_PythiaScaleUp *= syst_corr;
      break;
    case tHardJetsDn:
      event_wgt_adjustment_ext_HardJetsDn *= syst_corr;
      break;
    case tHardJetsUp:
      event_wgt_adjustment_ext_HardJetsUp *= syst_corr;
      break;
    default:
      break;
    }
  }

  // Record LHE MEs and K factors
  for (auto const& it:genInfo->extras.LHE_ME_weights) commonEntry.setNamedVal(it.first, it.second);
  for (auto const& it:genInfo->extras.Kfactors) commonEntry.setNamedVal(it.first, it.second);

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

void LooperFunctionHelpers::registerSystHistogram(SystematicsHelpers::SystematicVariationTypes const& syst, TH1F* hist){
  if (!hist) return;
  syst_h1d_map[syst] = hist;
  registeredSysts.push_back(syst);
  IVYout << "LooperFunctionHelpers::registerSystHistogram: Registered 1D histogram for " << SystematicsHelpers::getSystName(syst) << endl;
  if (syst==SystematicsHelpers::tPythiaScaleDn || syst==SystematicsHelpers::tPythiaScaleUp) applyPythiaScaleExternally = true;
}
void LooperFunctionHelpers::registerSystHistogram(SystematicsHelpers::SystematicVariationTypes const& syst, TH2F* hist){
  if (!hist) return;
  syst_h2d_map[syst] = hist;
  registeredSysts.push_back(syst);
  IVYout << "LooperFunctionHelpers::registerSystHistogram: Registered 2D histogram for " << SystematicsHelpers::getSystName(syst) << endl;
  if (syst==SystematicsHelpers::tPythiaScaleDn || syst==SystematicsHelpers::tPythiaScaleUp) applyPythiaScaleExternally = true;
}
void LooperFunctionHelpers::registerSystHistogram(SystematicsHelpers::SystematicVariationTypes const& syst, TH3F* hist){
  if (!hist) return;
  syst_h3d_map[syst] = hist;
  registeredSysts.push_back(syst);
  IVYout << "LooperFunctionHelpers::registerSystHistogram: Registered 3D histogram for " << SystematicsHelpers::getSystName(syst) << endl;
  if (syst==SystematicsHelpers::tPythiaScaleDn || syst==SystematicsHelpers::tPythiaScaleUp) applyPythiaScaleExternally = true;
}

double LooperFunctionHelpers::EvalSystHistogram(TH1F* hist, float const& xvar){
  if (!hist) return 1;

  const int nbinsx = hist->GetNbinsX();
  int ibin = hist->GetXaxis()->FindBin(xvar);
  if (ibin<1) ibin = 1;
  else if (ibin>nbinsx) ibin = nbinsx;

  return hist->GetBinContent(ibin);
}
double LooperFunctionHelpers::EvalSystHistogram(TH2F* hist, float const& xvar, float const& yvar){
  if (!hist) return 1;

  const int nbinsx = hist->GetNbinsX();
  int ibin = hist->GetXaxis()->FindBin(xvar);
  if (ibin<1) ibin = 1;
  else if (ibin>nbinsx) ibin = nbinsx;

  const int nbinsy = hist->GetNbinsY();
  int jbin = hist->GetYaxis()->FindBin(yvar);
  if (jbin<1) jbin = 1;
  else if (jbin>nbinsy) jbin = nbinsy;

  return hist->GetBinContent(ibin, jbin);
}
double LooperFunctionHelpers::EvalSystHistogram(TH3F* hist, float const& xvar, float const& yvar, float const& zvar){
  if (!hist) return 1;

  const int nbinsx = hist->GetNbinsX();
  int ibin = hist->GetXaxis()->FindBin(xvar);
  if (ibin<1) ibin = 1;
  else if (ibin>nbinsx) ibin = nbinsx;

  const int nbinsy = hist->GetNbinsY();
  int jbin = hist->GetYaxis()->FindBin(yvar);
  if (jbin<1) jbin = 1;
  else if (jbin>nbinsy) jbin = nbinsy;

  const int nbinz = hist->GetNbinsZ();
  int kbin = hist->GetZaxis()->FindBin(zvar);
  if (kbin<1) kbin = 1;
  else if (kbin>nbinz) kbin = nbinz;

  return hist->GetBinContent(ibin, jbin, kbin);
}

void LooperFunctionHelpers::cleanup(){
  for (auto& hh:syst_h1d_map) delete hh.second;
  for (auto& hh:syst_h2d_map) delete hh.second;
  for (auto& hh:syst_h3d_map) delete hh.second;
}

using namespace PhysicsProcessHelpers;


PhysicsProcessHandler* getPhysicsProcessHandler(TString strSampleSet, ACHypothesisHelpers::DecayType dktype){
  PhysicsProcessHandler* res = nullptr;
  if (strSampleSet.Contains("GGH") || strSampleSet.Contains("GluGluH")) res = new GGProcessHandler(dktype);
  else if (strSampleSet.Contains("VBF")) res = new VVProcessHandler(dktype, kProcess_VBF);
  else if (strSampleSet.Contains("ZH")) res = new VVProcessHandler(dktype, kProcess_ZH);
  else if (strSampleSet.Contains("WminusH") || strSampleSet.Contains("WplusH")) res = new VVProcessHandler(dktype, kProcess_WH);
  else{
    IVYerr << "getPhysicsProcessHandler: Cannot identify process " << strSampleSet;
    assert(0);
  }
  return res;
}

// Acquire special systematics reweighting histograms
void registerSystematics(TString strSampleSet){
  using namespace SystematicsHelpers;

  TDirectory* curdir = gDirectory;

  TString strinput_customSyst_main = ANALYSISTREEPKGDATAPATH + "ScaleFactors/SystematicsCustomReweighting/";
  HostHelpers::ExpandEnvironmentVariables(strinput_customSyst_main);
  std::vector<SystematicVariationTypes> const allowedSysts_1D{
    tPDFScaleDn, tPDFScaleUp,
    tQCDScaleDn, tQCDScaleUp,
    tAsMZDn, tAsMZUp,
    tPDFReplicaDn, tPDFReplicaUp
  };
  std::vector<SystematicVariationTypes> const allowedSysts_2D;
  std::vector<SystematicVariationTypes> const allowedSysts_3D{
    tPythiaScaleDn, tPythiaScaleUp,
    //tPythiaTuneDn, tPythiaTuneUp,
    tHardJetsDn, tHardJetsUp
  };
  std::vector<SystematicVariationTypes> allowedSysts = allowedSysts_1D;
  HelperFunctions::appendVector(allowedSysts, allowedSysts_2D);
  HelperFunctions::appendVector(allowedSysts, allowedSysts_3D);
  for (auto const& syst:allowedSysts){
    TString const strSyst = getSystName(syst).data();
    TString strinput_customSyst;
    TH1F* h_ratio_syst_1D = nullptr;
    TH2F* h_ratio_syst_2D = nullptr;
    TH3F* h_ratio_syst_3D = nullptr;

    if (SampleHelpers::getDataYear()==2016){
      strinput_customSyst = strinput_customSyst_main + "2016/" + strSampleSet + "_" + strSyst + ".root";
      if (!HostHelpers::FileReadable(strinput_customSyst)) strinput_customSyst = strinput_customSyst_main + "2016_2017_2018/" + strSampleSet + "_" + strSyst + ".root";
      if (!HostHelpers::FileReadable(strinput_customSyst)) strinput_customSyst = "";
    }
    else if (SampleHelpers::getDataYear()==2017 || SampleHelpers::getDataYear()==2018){
      strinput_customSyst = strinput_customSyst_main + "2017_2018/" + strSampleSet + "_" + strSyst + ".root";
      if (!HostHelpers::FileReadable(strinput_customSyst)) strinput_customSyst = strinput_customSyst_main + "2016_2017_2018/" + strSampleSet + "_" + strSyst + ".root";
      if (!HostHelpers::FileReadable(strinput_customSyst)) strinput_customSyst = "";
    }
    if (strinput_customSyst!=""){
      IVYout << "Acquiring the systematics file " << strinput_customSyst << "..." << endl;
      TFile* finput_syst = TFile::Open(strinput_customSyst, "read");
      if (HelperFunctions::checkListVariable(allowedSysts_1D, syst)){
        TH1F* htmp = (TH1F*) finput_syst->Get("h_ratio");
        curdir->cd();
        h_ratio_syst_1D = (TH1F*) htmp->Clone(htmp->GetName());
      }
      else if (HelperFunctions::checkListVariable(allowedSysts_2D, syst)){
        TH2F* htmp = (TH2F*) finput_syst->Get("h_ratio");
        curdir->cd();
        h_ratio_syst_2D = (TH2F*) htmp->Clone(htmp->GetName());
      }
      else if (HelperFunctions::checkListVariable(allowedSysts_3D, syst)){
        TH3F* htmp = (TH3F*) finput_syst->Get("h_ratio");
        curdir->cd();
        h_ratio_syst_3D = (TH3F*) htmp->Clone(htmp->GetName());
      }
      finput_syst->Close();
      if (!h_ratio_syst_1D && !h_ratio_syst_2D && !h_ratio_syst_3D){
        IVYerr << "\t- No 1, 2, or 3D reweighting could be acquired. Aborting..." << endl;
        assert(0);
      }
    }

    curdir->cd();

    LooperFunctionHelpers::registerSystHistogram(syst, h_ratio_syst_1D);
    LooperFunctionHelpers::registerSystHistogram(syst, h_ratio_syst_2D);
    LooperFunctionHelpers::registerSystHistogram(syst, h_ratio_syst_3D);
  }
}

void produceReweightedGen(
  TString strSampleSet, TString period,
  TString prodVersion, TString strdate
){
  SystematicsHelpers::SystematicVariationTypes const theGlobalSyst = SystematicsHelpers::sNominal;

  if (!SampleHelpers::checkRunOnCondor()) std::signal(SIGINT, SampleHelpers::setSignalInterrupt);

  if (strSampleSet.Contains("/MINIAOD")){
    IVYerr << "Processing single samples is not the design goal of produceReweightingRecords." << endl;
    return;
  }

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SampleHelpers::configure(period, Form("%s:%s", "store", prodVersion.Data()));

  constexpr float lumi = 1;

  std::vector<TString> const validDataPeriods = SampleHelpers::getValidDataPeriods();
  size_t const nValidDataPeriods = validDataPeriods.size();

  bool const isGG = strSampleSet.Contains("GluGluH") || strSampleSet.Contains("GGH");
  bool const isVBF = strSampleSet.Contains("VBF");
  bool const isWH = strSampleSet.Contains("WminusH") || strSampleSet.Contains("WplusH");
  //bool const isVH = strSampleSet.Contains("WminusH") || strSampleSet.Contains("WplusH") || strSampleSet.Contains("ZH") || strSampleSet.Contains("HZJ");
  SampleHelpers::HiggsSampleDecayMode const hdecaymode = SampleHelpers::getHiggsSampleDecayMode(strSampleSet);
  bool const hasHWWDecay = SampleHelpers::isHiggsToWWDecay(hdecaymode);
  bool const isPowheg = strSampleSet.Contains("POWHEG");
  bool const hasDirectHWW = isWH || hasHWWDecay;
  bool const isWHWW = isWH && hasDirectHWW;

  // Get sample specifications
  std::vector<TString> sampledirs;
  SampleHelpers::constructSamplesList(strSampleSet, theGlobalSyst, sampledirs);

  if (sampledirs.empty()) return;
  bool isData = SampleHelpers::checkSampleIsData(sampledirs.front());
  if (isData) return;

  // Set output directory
  TString coutput_main = "output/ReweightedGenTrees/" + strdate + "/" + period;

  TDirectory* curdir = gDirectory;
  gSystem->mkdir(coutput_main, true);

  curdir->cd();

  // Register external systematics shape reweighting
  if (isPowheg) registerSystematics(strSampleSet);

  // Set variables to allow recasting of topology
  LooperFunctionHelpers::recastLHETopology = isPowheg && !isGG;
  if (LooperFunctionHelpers::recastLHETopology) LooperFunctionHelpers::candScheme = (isVBF ? TVar::JJVBF : (isWH ? TVar::Had_WH : TVar::Had_ZH));

  // Set decay mode interpretation
  // Leave VVDecayMode as -1 in POWHEG samples
  LooperFunctionHelpers::VVMode = (hasHWWDecay ? MELAEvent::WWMode : MELAEvent::ZZMode);
  if (!isPowheg){
    switch (hdecaymode){
    case SampleHelpers::kZZTo4L:
      LooperFunctionHelpers::VVDecayMode = 0;
      break;
    default:
      break;
    }
  }

  // Create the output file
  TString coutput = SampleHelpers::getSampleIdentifier(strSampleSet);
  HelperFunctions::replaceString(coutput, "_MINIAODSIM", "");
  HelperFunctions::replaceString(coutput, "_MINIAOD", "");
  TString stroutput = Form("%s/%s", coutput_main.Data(), coutput.Data());
  stroutput += Form("_%s", SystematicsHelpers::getSystName(theGlobalSyst).data());
  stroutput += ".root";
  TString stroutput_txt = stroutput;
  HelperFunctions::replaceString<TString, TString const>(stroutput_txt, ".root", ".txt");
  TString stroutput_rewgtRcd = Form("%s/rewgtRcd_%s%s", coutput_main.Data(), coutput.Data(), ".root");
  TString stroutput_rewgtRcd_txt = Form("%s/rewgtRcd_%s%s", coutput_main.Data(), coutput.Data(), ".txt");
  TFile* foutput = TFile::Open(stroutput, "recreate");
  foutput->cd();
  BaseTree* tout = new BaseTree("SkimTree");
  IVYout << "Created output file " << stroutput << "..." << endl;
  curdir->cd();

  // Declare handlers
  GenInfoHandler genInfoHandler;
  genInfoHandler.setAcquireLHEMEWeights(true);
  genInfoHandler.setAcquireLHEParticles(true);
  genInfoHandler.setAcquireGenParticles(true);

  curdir->cd();

  // Configure the looper
  BaseTreeLooper theLooper;
  // Set systematic
  theLooper.setSystematic(theGlobalSyst);
  // Set looper function
  theLooper.setLooperFunction(LooperFunctionHelpers::looperRule);
  // Set object handlers
  theLooper.addObjectHandler(&genInfoHandler);

  // Set output tree
  theLooper.addOutputTree(tout);

  ExtendedBinning binning_rewgt;
  binning_rewgt.addBinBoundary(70);
  binning_rewgt.addBinBoundary(13000);
  {
    std::vector<double> sample_masses;
    for (auto const& sname:sampledirs){
      TString sid = SampleHelpers::getSampleIdentifier(sname);
      double sample_mass = SampleHelpers::findPoleMass(sid);
      if (sample_mass>0.) HelperFunctions::addByLowest(sample_masses, sample_mass, true);
    }
    if (sample_masses.size()>1){
      for (unsigned int im=0; im<sample_masses.size()-1; im++){
        double const& thisMass = sample_masses.at(im);
        double const& nextMass = sample_masses.at(im+1);

        // Special treatment for mH=125 GeV
        // Need to account for the relatively large gap with the next mass point
        if (std::abs(thisMass-125.)<0.8 && std::abs(nextMass-125.)>0.8){
          binning_rewgt.addBinBoundary(thisMass + (nextMass - thisMass)/3.);
          binning_rewgt.addBinBoundary(thisMass + (nextMass - thisMass)*2./3.);
        }
        else binning_rewgt.addBinBoundary((thisMass + nextMass)/2.);
      }
    }
  }
  BulkReweightingBuilder* rewgtBuilder = nullptr;
  if (isPowheg){
    IVYout << "Reweighting bin boundaries: " << binning_rewgt.getBinningVector() << endl;
    rewgtBuilder = new BulkReweightingBuilder(
      binning_rewgt,
      { "LHECandMass" },
      { "genHEPMCweight_default" },
      { "xsec" },
      ReweightingFunctions::getSimpleVariableBin,
      ReweightingFunctions::getSimpleWeight,
      ReweightingFunctions::getSimpleWeight
    );
    theLooper.addReweightingBuilder("MERewgt", rewgtBuilder);
  }

  // Acquire the MEs
  std::vector<std::pair<TString, bool>> strMEs;
  if (rewgtBuilder){
    PhysicsProcessHandler* proc_handler = getPhysicsProcessHandler(strSampleSet, (hdecaymode==SampleHelpers::kZZTo4L ? ACHypothesisHelpers::kZZ4l_offshell : ACHypothesisHelpers::kZZ2l2nu_offshell));
    for (unsigned int iac=0; iac<(unsigned int) ACHypothesisHelpers::nACHypotheses; iac++){
      ACHypothesisHelpers::ACHypothesis hypo = (ACHypothesisHelpers::ACHypothesis) iac;
      std::vector<TString> strMEs_hypo = proc_handler->getMELAHypothesisWeights(hypo, false);
      std::vector<bool> excludeFromNormRewgt_list(strMEs_hypo.size(), (hasDirectHWW && hypo==ACHypothesisHelpers::kL1ZGs));
      // WW-only interaction can also happen in VBF, but that concerns only one of the L1ZGs ME types.
      if (isVBF && hypo==ACHypothesisHelpers::kL1ZGs){
        for (unsigned int ibsi=0; ibsi<strMEs_hypo.size(); ibsi++){
          if (ibsi==(static_cast<unsigned int>(VVProcessHandler::VVTplSigBSM)-static_cast<unsigned int>(VVProcessHandler::nVVTplSMTypes))) excludeFromNormRewgt_list.at(ibsi)=true;
        }
      }
      std::vector<std::pair<TString, bool>> strMEs_tmp;
      HelperFunctions::zipVectors(strMEs_hypo, excludeFromNormRewgt_list, strMEs_tmp);
      HelperFunctions::appendVector(strMEs, strMEs_tmp);
    }
    delete proc_handler;
  }

  //int iRefHypo = -1;
  int iRefHypo = 0;
  constexpr double tol_wgt = 5;
  float const thr_frac_Neff = (isGG ? 0.005 : 0.01);
  {
    unsigned int ihypo = 0;
    for (auto const& strME_pp:strMEs){
      auto const& strME = strME_pp.first;
      auto const& excludeFromNormRewgt = strME_pp.second;
      double thr_wgt = 0.9995;
      if (isGG){
        //if (strME == "p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_MCFM") iRefHypo = ihypo;
        /* Do nothing, 0.9995 is good enough. */
      }
      else if (strME == "p_Gen_JJEW_SIG_ghv1_1_MCFM"){
        //iRefHypo = ihypo;
        /* Do nothing to weight threshold, 0.9995 is good enough. */
      }
      else if (isVBF) thr_wgt = 0.999;
      else if (isWHWW) thr_wgt = 0.995;

      std::vector<TString> wgtcoll{ strME };
      if (isPowheg) wgtcoll.push_back("p_Gen_CPStoBWPropRewgt");
      if (rewgtBuilder) rewgtBuilder->addReweightingWeights(
        wgtcoll,
        ReweightingFunctions::getSimpleWeight,
        thr_wgt, tol_wgt,
        excludeFromNormRewgt
      );

      ihypo++;
    }
    if (iRefHypo>=0) IVYout << "Reweighting builder reference hypothesis index is " << iRefHypo << " (" << strMEs.at(iRefHypo).first << ")." << endl;
    else{
      if (rewgtBuilder){
        IVYerr << "Reweighting builder reference hypothesis is not found. Aborting... " << endl;
        delete rewgtBuilder;
        delete tout;
        foutput->Close();
        return;
      }
    }
  }

  curdir->cd();

  bool allTreesValid = true;
  BaseTree* tree_MH125 = nullptr;
  BaseTree* tree_MHLowestOffshell = nullptr;
  std::vector<BaseTree*> sample_trees; sample_trees.reserve(sampledirs.size());
  for (auto const& sname:sampledirs){
    TString strinput = SampleHelpers::getDatasetFileName(sname);
    IVYout << "Acquiring " << sname << " from input file(s) " << strinput << "..." << endl;
    BaseTree* sample_tree = new BaseTree(strinput, "cms3ntuple/Events", "", ""); sample_trees.push_back(sample_tree);
    if (!sample_tree->isValid()){
      IVYerr << "\t- Tree is invalid. Aborting..." << endl;
      delete sample_tree;
      for (auto& ss:sample_trees) delete ss;
      allTreesValid = false;
      break;
    }
    sample_tree->sampleIdentifier = SampleHelpers::getSampleIdentifier(sname);
    float const sampleMH = (isPowheg ? SampleHelpers::findPoleMass(sample_tree->sampleIdentifier) : -1.f);
    if (std::abs(sampleMH-125.f)<0.8f) tree_MH125 = sample_tree;
    else if (!tree_MHLowestOffshell && sampleMH>=200.f) tree_MHLowestOffshell = sample_tree;

    std::vector<TString> allbranchnames; sample_tree->getValidBranchNamesWithoutAlias(allbranchnames, false);

    const int nEntries = sample_tree->getNEvents();
    bool hasTaus = false;
    double sum_wgts_raw_noveto = 0;
    double sum_wgts_raw_withveto = 0;
    double sum_wgts_raw_withveto_defaultMemberZero = 0;
    float xsec = 1;
    float xsec_scale = 1;
    float BR_scale = 1;
    if (!isData){
      // Get cross section
      sample_tree->bookBranch<float>("xsec", 0.f);
      sample_tree->getSelectedEvent(0);
      sample_tree->getVal("xsec", xsec);
      xsec *= 1000.;

      // Book branches
      genInfoHandler.setAcquireLHEMEWeights(false);
      genInfoHandler.setAcquireLHEParticles(true);
      genInfoHandler.setAcquireGenParticles(true);
      genInfoHandler.bookBranches(sample_tree);

      sample_tree->silenceUnused();

      // Get sum of weights
      {
        IVYout << "No counters histograms are found. Initiation loop over " << nEntries << " events to determine the sample normalization:" << endl;

        genInfoHandler.wrapTree(sample_tree);

        unsigned int n_zero_genwgts=0;
        double frac_zero_genwgts=0;
        for (int ev=0; ev<nEntries; ev++){
          HelperFunctions::progressbar(ev, nEntries);
          sample_tree->getEvent(ev);

          genInfoHandler.constructGenInfo(SystematicsHelpers::sNominal); // Use sNominal here in order to get the weight that corresponds to xsec
          auto const& genInfo = genInfoHandler.getGenInfo();

          auto const& lheparticles = genInfoHandler.getLHEParticles();
          auto const& genparticles = genInfoHandler.getGenParticles();
          if (!hasTaus){
            if (!lheparticles.empty()){
              genInfoHandler.setAcquireGenParticles(false);
              for (auto const& part:lheparticles){
                if (part->status()==1 && std::abs(part->pdgId())==15){
                  hasTaus = true;
                  genInfoHandler.setAcquireLHEParticles(false);
                  break;
                }
              }
            }
            else{
              genInfoHandler.setAcquireLHEParticles(false);
              if (!genparticles.empty()){
                for (auto const& part:genparticles){
                  if (!part->extras.isHardProcess) continue;
                  if ((part->status()==1 || part->status()==2) && std::abs(part->pdgId())==15){
                    hasTaus = true;
                    genInfoHandler.setAcquireGenParticles(false);
                    break;
                  }
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

          sum_wgts_raw_withveto_defaultMemberZero += genwgt_defaultMemberZero;
          sum_wgts_raw_withveto += genwgt;
        }
        if (nEntries>0) frac_zero_genwgts = double(n_zero_genwgts)/double(nEntries);
        sum_wgts_raw_noveto = sum_wgts_raw_withveto / (1. - frac_zero_genwgts);
      }
      xsec_scale = sum_wgts_raw_withveto / sum_wgts_raw_noveto;
      if (sampleMH>0.f) BR_scale = SampleHelpers::calculateAdjustedHiggsBREff(sname, sum_wgts_raw_withveto_defaultMemberZero, sum_wgts_raw_withveto, hasTaus);
    }
    std::unordered_map<SystematicsHelpers::SystematicVariationTypes, double> globalWeights;
    double globalWeight = xsec * xsec_scale * BR_scale / sum_wgts_raw_withveto; globalWeights[theGlobalSyst] = globalWeight;
    IVYout << "Sample " << sample_tree->sampleIdentifier << " has a gen. weight sum of " << sum_wgts_raw_withveto << "." << endl;
    IVYout << "\t- Raw xsec = " << xsec << endl;
    IVYout << "\t- xsec scale * BR scale = " << xsec_scale * BR_scale << endl;
    IVYout << "\t- xsec * BR = " << xsec * xsec_scale * BR_scale << endl;
    IVYout << "\t- Global weight = " << globalWeight << endl;

    // Reset gen. and LHE particle settings, and book those branches as well
    {
      genInfoHandler.setAcquireLHEMEWeights(true);
      genInfoHandler.setAcquireLHEParticles(true);
      genInfoHandler.setAcquireGenParticles(true);
      genInfoHandler.setAcquireGenAK4Jets(true);
      // Specify this flag to omit V->qq->jj
      genInfoHandler.setDoGenJetsVDecayCleaning(false);
      // Do not clean jets, allow user to do this on their own
      genInfoHandler.setDoGenJetsCleaning(false);
      genInfoHandler.bookBranches(sample_tree);
    }

    // Register tree
    if (rewgtBuilder){
      IVYout << "\t- Registering the sample for reweighting..." << endl;
      rewgtBuilder->registerTree(sample_tree, xsec_scale * BR_scale / sum_wgts_raw_withveto);
    }

    sample_tree->silenceUnused();

    IVYout << "\t- Registering the sample to the looper..." << endl;
    // Add the input tree to the looper
    theLooper.addTree(sample_tree, globalWeights);
    //if (sample_trees.size()==2) break;
  }

  if (rewgtBuilder){
    std::vector< std::pair<BaseTree*, BaseTree*> > tree_normTree_pairs; tree_normTree_pairs.reserve(sample_trees.size()-1);
    for (unsigned int itree=0; itree<sample_trees.size(); itree++){
      BaseTree* const& sample_tree = sample_trees.at(itree);
      if (sample_tree==tree_MH125 || sample_tree==tree_MHLowestOffshell) continue;
      float const sampleMH = SampleHelpers::findPoleMass(sample_tree->sampleIdentifier);
      if ((tree_MH125 && sampleMH>SampleHelpers::findPoleMass(tree_MH125->sampleIdentifier) && sampleMH<160.f) || (tree_MHLowestOffshell && sampleMH>SampleHelpers::findPoleMass(tree_MHLowestOffshell->sampleIdentifier))){
        tree_normTree_pairs.emplace_back(sample_trees.at(itree), sample_trees.at(itree-1));
        IVYout << "Normalizing mass " << sampleMH << " to mass " << SampleHelpers::findPoleMass(sample_trees.at(itree-1)->sampleIdentifier) << endl;
      }
    }

    // Record the reweighting weights or use the existing record
    bool const rwgtrcd_exists = HostHelpers::FileReadable(stroutput_rewgtRcd);
    if (rwgtrcd_exists) rewgtBuilder->setupFromFile(stroutput_rewgtRcd);
    else{
      IVYout.open(stroutput_rewgtRcd_txt);
      rewgtBuilder->setup(iRefHypo, &tree_normTree_pairs, thr_frac_Neff);
      IVYout.close();

      TFile* foutput_rewgtRcd = TFile::Open(stroutput_rewgtRcd, "recreate");
      rewgtBuilder->writeToFile(foutput_rewgtRcd);
      foutput_rewgtRcd->Close();
    }
  }

  curdir->cd();

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

  LooperFunctionHelpers::cleanup();

  SampleHelpers::addToCondorTransferList(stroutput);
}

void produceHistograms(TString strSampleSet, TString period, TString strdate){
  SystematicsHelpers::SystematicVariationTypes const theGlobalSyst = SystematicsHelpers::sNominal;

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  bool const isGG = strSampleSet.Contains("GluGluH") || strSampleSet.Contains("GGH");
  //bool isVVVV = strSampleSet.Contains("VBF") || strSampleSet.Contains("ZH") || strSampleSet.Contains("WminusH") || strSampleSet.Contains("WplusH");
  IVYout << "Sample is " << (isGG ? "a gg" : "an EW") << " sample." << endl;

  // Set output directory
  TString cinput_main = "output/ReweightedGenTrees/" + strdate + "/" + period;
  TString coutput_main = "output/ReweightedGenTrees/" + strdate + "/" + period + "/Plots";

  TDirectory* curdir = gDirectory;
  gSystem->mkdir(coutput_main, true);

  curdir->cd();

  // Create the output file
  TString coutput = SampleHelpers::getSampleIdentifier(strSampleSet);
  HelperFunctions::replaceString(coutput, "_MINIAODSIM", "");
  HelperFunctions::replaceString(coutput, "_MINIAOD", "");
  TString stroutput = Form("%s/%s", coutput_main.Data(), coutput.Data());
  stroutput += Form("_%s", SystematicsHelpers::getSystName(theGlobalSyst).data());
  stroutput += ".root";
  TFile* foutput = TFile::Open(stroutput, "recreate");
  IVYout << "Created output file " << stroutput << "..." << endl;
  foutput->cd();
  
  TString strinput = cinput_main + "/" + coutput + "_" + SystematicsHelpers::getSystName(theGlobalSyst).data() + ".root";
  IVYout << "Acquiring input " << strinput << "..." << endl;
  BaseTree* tin = new BaseTree(strinput, "SkimTree", "", "");
  if (!tin->isValid()){
    IVYerr << "\t- Failed to acquire." << endl;
    delete tin;
    foutput->Close();
    exit(1);
  }

#define BRANCH_COMMAND(TYPE, NAME) TYPE* NAME = nullptr; tin->bookBranch<TYPE>(#NAME, 0);
  BRANCH_SCALAR_COMMANDS;
#undef BRANCH_COMMAND
#define BRANCH_COMMAND(TYPE, NAME) std::vector<TYPE>** NAME=nullptr; tin->bookBranch<std::vector<TYPE>*>(#NAME, nullptr);
  BRANCH_VECTOR_COMMANDS;
#undef BRANCH_COMMAND
#define BRANCH_COMMAND(TYPE, NAME) tin->getValRef(#NAME, NAME);
  BRANCH_COMMANDS;
#undef BRANCH_COMMAND

  std::unordered_map<TString, float*> ME_Kfactor_values;
  std::vector<TString> allbranchnames;
  tin->getValidBranchNamesWithoutAlias(allbranchnames, false);
  for (auto const& bname:allbranchnames){
    if (
      (bname.BeginsWith("p_") && (bname.Contains("JHUGen") || bname.Contains("MCFM")))
      ||
      bname.BeginsWith("p_Gen")
      ||
      bname.Contains("LHECandMass")
      ||
      bname.BeginsWith("KFactor")
      ){
      tin->bookBranch<float>(bname, -1.f);
      ME_Kfactor_values[bname] = nullptr;
      tin->getValRef(bname, ME_Kfactor_values[bname]);
    }
  }

  tin->silenceUnused();

  float* val_Kfactor_QCD = nullptr;
  float* val_ME_SIG = nullptr;
  if (isGG){
    val_Kfactor_QCD = ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_Nominal")->second;
    val_ME_SIG = ME_Kfactor_values.find("p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_MCFM")->second;
  }
  else{
    val_ME_SIG = ME_Kfactor_values.find("p_Gen_JJEW_SIG_ghv1_1_MCFM")->second;
  }
  float* val_ME_CPS = ME_Kfactor_values.find("p_Gen_CPStoBWPropRewgt")->second;
  float* LHECandMass = ME_Kfactor_values.find("LHECandMass")->second;

  bool hasError = false;
  if (!val_ME_SIG){
    IVYerr << "val_ME_SIG is null!" << endl;
    hasError = true;
  }
  if (!val_ME_CPS){
    IVYerr << "val_ME_CPS is null!" << endl;
    hasError = true;
  }
  if (isGG && !val_Kfactor_QCD){
    IVYerr << "val_Kfactor_QCD is null!" << endl;
    hasError = true;
  }
  if (hasError){
    delete tin;
    foutput->Close();
    exit(1);
  }

  std::vector<TString> const strCatNames{ "Nj_eq_0", "Nj_eq_1", "Nj_geq_2" };
  std::vector<TString> const strCatLabels{ "N_{j}=0", "N_{j}=1", "N_{j} #geq 2" };
  size_t const nCats = strCatNames.size();

  ExtendedBinning binning_mass("genmass", "m_{ZZ} (GeV)");
  constexpr bool useConstBinning = false;
  if (!useConstBinning){
    binning_mass.addBinBoundary(100);
    binning_mass.addBinBoundary(124);
    binning_mass.addBinBoundary(126);
    binning_mass.addBinBoundary(150);
    {
      double bin_low_edge = 180;
      while (bin_low_edge<=500.){
        binning_mass.addBinBoundary(bin_low_edge);
        bin_low_edge += 20.;
      }
      while (bin_low_edge<=1000.){
        binning_mass.addBinBoundary(bin_low_edge);
        bin_low_edge += 50.;
      }
      while (bin_low_edge<=2000.){
        binning_mass.addBinBoundary(bin_low_edge);
        bin_low_edge += 100.;
      }
      while (bin_low_edge<=3000.){
        binning_mass.addBinBoundary(bin_low_edge);
        bin_low_edge += 500.;
      }
    }
  }
  else{
    double bin_low_edge = 100;
    while (bin_low_edge<=3000.){
      binning_mass.addBinBoundary(bin_low_edge);
      bin_low_edge += 25.;
    }
  }

  std::vector<TH1F> hlist_genmass; hlist_genmass.reserve(strCatNames.size());
  std::vector<TH1F> hlist_genak4jetpt; hlist_genak4jetpt.reserve(strCatNames.size());
  std::vector<TH1F> hlist_genak4jetselectedpt; hlist_genak4jetselectedpt.reserve(strCatNames.size());
  std::vector<TH1F> hlist_genak4jetselected_leadingpt; hlist_genak4jetselected_leadingpt.reserve(strCatNames.size());
  std::vector<TH1F> hlist_genak4jetselected_subleadingpt; hlist_genak4jetselected_subleadingpt.reserve(strCatNames.size());
  for (unsigned short icat=0; icat<nCats; icat++){
    hlist_genmass.emplace_back(
      Form("h_genmass_%s", strCatNames.at(icat).Data()), strCatLabels.at(icat),
      binning_mass.getNbins(), binning_mass.getBinning()
    );
    hlist_genak4jetpt.emplace_back(
      Form("h_genak4jetpt_%s", strCatNames.at(icat).Data()), strCatLabels.at(icat),
      100, 20, 1020
    );
    hlist_genak4jetselectedpt.emplace_back(
      Form("h_genak4jetselectedpt_%s", strCatNames.at(icat).Data()), strCatLabels.at(icat),
      100, 20, 1020
    );
    hlist_genak4jetselected_leadingpt.emplace_back(
      Form("h_genak4jetselected_leadingpt_%s", strCatNames.at(icat).Data()), strCatLabels.at(icat),
      100, 20, 1020
    );
    hlist_genak4jetselected_subleadingpt.emplace_back(
      Form("h_genak4jetselected_subleadingpt_%s", strCatNames.at(icat).Data()), strCatLabels.at(icat),
      100, 20, 1020
    );
  }

  float sum_wgts = 0;
  float sum_wgts_accepted = 0;
  float sum_wgts_rejected_Nleptons_lt_4 = 0;
  float sum_wgts_rejected_Nleptons_gt_4 = 0;
  float sum_wgts_rejected_noZZ = 0;
  float sum_wgts_hasHardPhotons = 0;
  int n_printouts = -1;
  int nEntries = tin->getNEvents();
  IVYout << "Looping over " << nEntries << " events:" << endl;
  for (int ev=0; ev<nEntries; ev++){
    tin->getEvent(ev);
    HelperFunctions::progressbar(ev, nEntries);

    float wgt = (*event_wgt) * (*sample_wgt) * (*invalidReweightingWgts ? 0.f : 1.f) * (val_Kfactor_QCD ? *val_Kfactor_QCD : 1.f) * (*val_ME_SIG) * (*val_ME_CPS);
    sum_wgts += wgt;

    unsigned int n_genparticles = (*genparticles_id)->size();
    std::vector<MELAParticle> genparticles; genparticles.reserve(n_genparticles);
    for (unsigned int ipart=0; ipart<n_genparticles; ipart++){
      TLorentzVector tmp_p4;
      tmp_p4.SetPtEtaPhiM((*genparticles_pt)->at(ipart), (*genparticles_eta)->at(ipart), (*genparticles_phi)->at(ipart), (*genparticles_mass)->at(ipart));
      genparticles.emplace_back((*genparticles_id)->at(ipart), tmp_p4);
    }

    int sumid_genleptons_selected = 0;
    TLorentzVector sump4_genleptons_selected;
    std::vector<MELAParticle const*> genleptons_selected; genleptons_selected.reserve(n_genparticles);
    std::vector<MELAParticle const*> genphotons_selected; genphotons_selected.reserve(n_genparticles);
    for (auto const& part:genparticles){
      int pid = part.id;
      if (PDGHelpers::isALepton(pid) && part.pt()>=5. && std::abs(part.eta())<(std::abs(pid)==11 ? 2.5 : 2.4)){
        sumid_genleptons_selected += part.id;
        sump4_genleptons_selected += part.p4;
        genleptons_selected.push_back(&part);
      }
      else if (PDGHelpers::isAPhoton(pid) && part.pt()>=20. && std::abs(part.eta())<2.5) genphotons_selected.push_back(&part);
    }

    if (genleptons_selected.size()!=4){
      if (n_printouts>=0 && n_printouts<100){
        IVYout << "Event " << ev << " is rejected because it has " << genleptons_selected.size() << " leptons." << endl;
        n_printouts++;
      }
      if (genleptons_selected.size()<4) sum_wgts_rejected_Nleptons_lt_4 += wgt;
      if (genleptons_selected.size()>4) sum_wgts_rejected_Nleptons_gt_4 += wgt;
      continue;
    }
    if (sumid_genleptons_selected!=0){
      if (n_printouts>=0 && n_printouts<100){
        IVYout << "Event " << ev << " is rejected because it doesn't have a proper ZZ candidate." << endl;
        n_printouts++;
      }
      sum_wgts_rejected_noZZ += wgt;
      continue;
    }
    sum_wgts_accepted += wgt;

    if (!genphotons_selected.empty()){
      bool hasHardPhotons = false;
      for (auto const& photon:genphotons_selected){
        bool isSeparated = true;
        for (auto const& lepton:genleptons_selected){
          if (photon->deltaR(lepton)<0.4){
            isSeparated = false;
            break;
          }
        }
        if (isSeparated){
          hasHardPhotons = true;
          break;
        }
      }
      if (hasHardPhotons) sum_wgts_hasHardPhotons += wgt;
    }

    unsigned int n_genak4jets = (*genak4jets_pt)->size();
    std::vector<MELAParticle> genak4jets; genak4jets.reserve(n_genak4jets);
    for (unsigned int ipart=0; ipart<n_genak4jets; ipart++){
      TLorentzVector tmp_p4;
      tmp_p4.SetPtEtaPhiM((*genak4jets_pt)->at(ipart), (*genak4jets_eta)->at(ipart), (*genak4jets_phi)->at(ipart), (*genak4jets_mass)->at(ipart));
      genak4jets.emplace_back(0, tmp_p4);
    }

    std::vector<MELAParticle const*> genak4jets_selected; genak4jets_selected.reserve(genak4jets.size());
    for (auto const& jet:genak4jets){
      if (jet.pt()<30. || std::abs(jet.eta())>=4.7) continue;
      bool doSkip = false;
      for (auto const& part:genleptons_selected){
        if (jet.deltaR(part)<0.4){ doSkip = true; break; }
      }
      if (!doSkip) genak4jets_selected.push_back(&jet);
    }

    unsigned int icat = std::min(genak4jets_selected.size(), nCats-1);

    hlist_genmass.at(icat).Fill(
      //*LHECandMass,
      sump4_genleptons_selected.M(),
      wgt
    );
    if (/*(*LHECandMass)*/sump4_genleptons_selected.M()>200.){
      for (auto const& jet:genak4jets) hlist_genak4jetpt.at(icat).Fill(jet.pt(), wgt);
      for (auto const& jet:genak4jets_selected) hlist_genak4jetselectedpt.at(icat).Fill(jet->pt(), wgt);
      if (genak4jets_selected.size()>=1){
        hlist_genak4jetselected_leadingpt.at(icat).Fill(genak4jets_selected.at(0)->pt(), wgt);
        if (genak4jets_selected.size()>=2){
          hlist_genak4jetselected_subleadingpt.at(icat).Fill(genak4jets_selected.at(1)->pt(), wgt);
        }
      }
    }
  }

  IVYout << "Sum of weights: " << sum_wgts << endl;
  IVYout << "\t- Fraction of events accepted: " << sum_wgts_accepted / sum_wgts << endl;
  IVYout << "\t- Fraction of events rejected because selected number of leptons < 4: " << sum_wgts_rejected_Nleptons_lt_4 / sum_wgts << endl;
  IVYout << "\t- Fraction of events rejected because selected number of leptons > 4: " << sum_wgts_rejected_Nleptons_gt_4 / sum_wgts << endl;
  IVYout << "\t- Fraction of events rejected because no ZZ candidate can be constructed: " << sum_wgts_rejected_noZZ / sum_wgts << endl;
  IVYout << "\t- Fraction of accepted events with at least one hard photon: " << sum_wgts_hasHardPhotons / sum_wgts_accepted << endl;

  for (auto& hh:hlist_genmass){
    HelperFunctions::divideBinWidth(&hh);
    foutput->WriteTObject(&hh);
  }
  for (auto& hh:hlist_genak4jetpt){
    HelperFunctions::divideBinWidth(&hh);
    foutput->WriteTObject(&hh);
  }
  for (auto& hh:hlist_genak4jetselectedpt){
    HelperFunctions::divideBinWidth(&hh);
    foutput->WriteTObject(&hh);
  }
  for (auto& hh:hlist_genak4jetselected_leadingpt){
    HelperFunctions::divideBinWidth(&hh);
    foutput->WriteTObject(&hh);
  }
  for (auto& hh:hlist_genak4jetselected_subleadingpt){
    HelperFunctions::divideBinWidth(&hh);
    foutput->WriteTObject(&hh);
  }

  delete tin;
  foutput->Close();
}

void produceBWHistograms(TString strSampleSet, TString period, TString prodVersion, TString strdate){
  SystematicsHelpers::SystematicVariationTypes const theGlobalSyst = SystematicsHelpers::sNominal;

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  TDirectory* curdir = gDirectory;
  //gSystem->mkdir(coutput_main, true);
  //curdir->cd();

  SampleHelpers::configure(period, Form("%s:%s", "store", prodVersion.Data()));

  std::vector<TString> sampledirs;
  SampleHelpers::constructSamplesList(strSampleSet, theGlobalSyst, sampledirs);
  if (sampledirs.size()!=1) return;

  TString infiles = SampleHelpers::getDatasetFileName(sampledirs.front());
  //HelperFunctions::replaceString<TString, TString const>(infiles, "/ceph/", "/hadoop/");
  BaseTree sample_tree(infiles, "cms3ntuple/Events", "", "");
  sample_tree.sampleIdentifier = SampleHelpers::getSampleIdentifier(sampledirs.front());
  bool const isData = SampleHelpers::checkSampleIsData(sample_tree.sampleIdentifier);
  if (isData) return;
  float const sampleMH = SampleHelpers::findPoleMass(sample_tree.sampleIdentifier);

  curdir->cd();

  TString coutput_main = "output/BWRewgtPlots/" + strdate;
  gSystem->mkdir(coutput_main, true);
  TString coutput = coutput_main + "/hists_" + strSampleSet + ".root";
  TFile* foutput = TFile::Open(coutput, "recreate");

#define HIST_COMMAND(NAME) \
  TH1F* h_##NAME = new TH1F(Form("h_%s", #NAME), "", 70, 200, 1600); \
  h_##NAME->Sumw2(); \
  h_##NAME->GetXaxis()->SetTitle("m_{ZZ} (GeV)"); \
  h_##NAME->GetYaxis()->SetTitle("Events / GeV"); \
  h_##NAME->SetLineWidth(2);
  HIST_COMMAND(norewgt);
  HIST_COMMAND(BRrewgt);
  HIST_COMMAND(fullproprewgt);
#undef HIST_COMMAND

  curdir->cd();

  bool hasTaus = false;
  double sum_wgts_raw_noveto = 0;
  double sum_wgts_raw_withveto = 0;
  double sum_wgts_raw_withveto_defaultMemberZero = 0;
  float xsec = 1;
  float xsec_scale = 1;
  float BR_scale = 1;
  float* LHECandMass = nullptr;
  float* p_Gen_CPStoBWPropRewgt = nullptr;

  // Get cross section
  sample_tree.bookBranch<float>("LHECandMass", -1.f);
  sample_tree.bookBranch<float>("p_Gen_CPStoBWPropRewgt", -1.f);
  sample_tree.getValRef("LHECandMass", LHECandMass);
  sample_tree.getValRef("p_Gen_CPStoBWPropRewgt", p_Gen_CPStoBWPropRewgt);

  sample_tree.silenceUnused();

  GenInfoHandler genInfoHandler;
  genInfoHandler.setAcquireLHEMEWeights(false);
  genInfoHandler.setAcquireLHEParticles(true);
  genInfoHandler.setAcquireGenParticles(false);
  genInfoHandler.bookBranches(&sample_tree);
  genInfoHandler.wrapTree(&sample_tree);

  unsigned int n_zero_genwgts=0;
  double frac_zero_genwgts=0;
  const int nEntries = sample_tree.getNEvents();
  for (int ev=0; ev<nEntries; ev++){
    HelperFunctions::progressbar(ev, nEntries);
    sample_tree.getEvent(ev);

    genInfoHandler.constructGenInfo(SystematicsHelpers::sNominal); // Use sNominal here in order to get the weight that corresponds to xsec
    auto const& genInfo = genInfoHandler.getGenInfo();
    if (ev==0) xsec = genInfo->extras.xsec*1000.;

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

    double genwgt = genInfo->getGenWeight(true);
    double genwgt_defaultMemberZero = genInfo->extras.LHEweight_defaultMemberZero;
    if (genwgt==0.){
      n_zero_genwgts++;
      continue;
    }
    sum_wgts_raw_withveto_defaultMemberZero += genwgt_defaultMemberZero;
    sum_wgts_raw_withveto += genwgt;

    h_norewgt->Fill(*LHECandMass, genwgt_defaultMemberZero);
    h_BRrewgt->Fill(*LHECandMass, genwgt);
    h_fullproprewgt->Fill(*LHECandMass, genwgt*(*p_Gen_CPStoBWPropRewgt));
  }
  if (nEntries>0) frac_zero_genwgts = double(n_zero_genwgts)/double(nEntries);
  sum_wgts_raw_noveto = sum_wgts_raw_withveto / (1. - frac_zero_genwgts);
  xsec_scale = sum_wgts_raw_withveto / sum_wgts_raw_noveto;
  if (sampleMH>0.f){
    IVYout << "Calculating BR scale..." << endl;
    BR_scale = SampleHelpers::calculateAdjustedHiggsBREff(sampledirs.front(), sum_wgts_raw_withveto_defaultMemberZero, sum_wgts_raw_withveto, hasTaus);
  }
  IVYout << "xsec scale = " << xsec_scale << endl;
  IVYout << "BR scale = " << BR_scale << endl;
  double globalWeight_fullproprewgt = xsec * xsec_scale * BR_scale / sum_wgts_raw_withveto;
  double globalWeight_BRrewgt = globalWeight_fullproprewgt;
  double globalWeight_norewgt = globalWeight_fullproprewgt * sum_wgts_raw_withveto_defaultMemberZero/sum_wgts_raw_withveto;
  IVYout << "Sample " << sample_tree.sampleIdentifier << " has a gen. weight sum of " << sum_wgts_raw_withveto << "." << endl;
  IVYout << "\t- Raw xsec = " << xsec << endl;
  IVYout << "\t- xsec scale * BR scale = " << xsec_scale * BR_scale << endl;
  IVYout << "\t- xsec * BR = " << xsec * xsec_scale * BR_scale << endl;
  IVYout << "\t- Global weights = " << globalWeight_norewgt << ", " << globalWeight_BRrewgt << ", " << globalWeight_fullproprewgt << endl;
  h_norewgt->Scale(globalWeight_norewgt); HelperFunctions::divideBinWidth(h_norewgt);
  h_BRrewgt->Scale(globalWeight_BRrewgt); HelperFunctions::divideBinWidth(h_BRrewgt);
  h_fullproprewgt->Scale(globalWeight_fullproprewgt); HelperFunctions::divideBinWidth(h_fullproprewgt);

  curdir->cd();
  foutput->WriteTObject(h_norewgt); delete h_norewgt;
  foutput->WriteTObject(h_BRrewgt); delete h_BRrewgt;
  foutput->WriteTObject(h_fullproprewgt); delete h_fullproprewgt;
  foutput->Close();
}


#include "PlottingHelpers.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TText.h"
#include "TPaveText.h"
#include "TLegend.h"

void makePlot(
  TString const& coutput_main,
  TString canvasname,
  std::vector<TH1F*> hlist,
  std::vector<TH1F*> hratios,
  std::vector<TString> hlabels,
  TString selectionLabels,
  TString drawopts,
  bool useLogY,
  bool adjustYLow,
  float factorYHigh
){
  using namespace PlottingHelpers;

  bool const addRatioPanel = !hratios.empty();

  size_t nplottables = hlist.size();
  if (hlabels.size()!=nplottables) return;
  for (auto const& hlabel:hlabels){ if (hlabel=="") nplottables--; }

  std::vector<TString> selectionList;
  if (selectionLabels!="") HelperFunctions::splitOptionRecursive(selectionLabels, selectionList, '|');

  bool hasData = false;
  if (useLogY) adjustYLow = true;

  std::vector<bool> hHasErrors;

  int nbins = -1;
  double ymin = 0;
  if (adjustYLow) ymin=9e9;
  double ymax = -9e9;
  for (size_t is=0; is<hlist.size(); is++){
    TH1F* hist = hlist.at(is);
    TString hname = hist->GetName();
    if (hname.Contains("Data") || hname.Contains("data")) hasData = true;
    bool hasErrors=false;
    if (nbins<0) nbins = hist->GetNbinsX();
    for (int ix=1; ix<=nbins; ix++){
      double bc = hist->GetBinContent(ix);
      double be = hist->GetBinError(ix);
      if (be>0.2*std::abs(bc)) be = 0.2*std::abs(bc);
      if (be!=0.f) hasErrors = true;
      ymax = std::max(ymax, bc+be);
      double bclow=bc; if (be<=bclow) bclow -= be;
      if (adjustYLow && bc>=1e-2) ymin = std::min(ymin, bclow);
    }
    hHasErrors.push_back(hasErrors);
    //IVYout << "ymin, ymax after " << hname << ": " << ymin << ", " << ymax << endl;
  }
  if (ymax>=0.) ymax *= (factorYHigh>0.f ? factorYHigh : (useLogY ? 150. : 1.5));
  else ymax /= (factorYHigh>0.f ? factorYHigh : 1.5);
  ymin *= (ymin>=0. ? 0.95 : 1.05);
  for (TH1F* const& hist:hlist) hist->GetYaxis()->SetRangeUser(ymin, ymax);

  TString varlabel = "";
  TString quantlabel = "";
  std::vector<TH1F*> hnum_MC_list;
  std::unordered_map<TH1F*, TGraphAsymmErrors*> hist_tg_map;
  for (size_t is=0; is<hlist.size(); is++){
    TH1F* hist = hlist.at(is);
    TString hname = hist->GetName();

    if (varlabel=="") varlabel = hist->GetXaxis()->GetTitle();
    if (quantlabel=="") quantlabel = hist->GetYaxis()->GetTitle();

    hist->GetXaxis()->SetTitle("");
    hist->GetYaxis()->SetTitle("");

    if (hname.Contains("Data")){
      TGraphAsymmErrors* tg = nullptr;
      HelperFunctions::convertTH1FToTGraphAsymmErrors(hist, tg, false, true);
      tg->SetName(Form("%s_noZeros", tg->GetName()));
      tg->SetTitle("");
      tg->GetYaxis()->SetRangeUser(ymin, ymax);

      tg->GetXaxis()->SetTitle("");
      tg->GetYaxis()->SetTitle("");

      hist_tg_map[hist] = tg;
    }
  }

  constexpr double npixels_pad_xy = 800;
  CMSLogoStep cmslogotype = kSimulation;
  PlotCanvas plot(canvasname, npixels_pad_xy, npixels_pad_xy, 1, (addRatioPanel ? 2 : 1), 0.2, 0.05, 0.15, 0.07, 0., 0.1, 0.2);
  plot.addCMSLogo(cmslogotype, 13, -1, 0, "13 TeV");

  for (size_t is=0; is<hlist.size(); is++){
    TH1F* hist = hlist.at(is);

    hist->GetXaxis()->SetNdivisions(505);
    hist->GetXaxis()->SetLabelFont(43);
    hist->GetXaxis()->SetLabelOffset(plot.getStdOffset_XLabel());
    hist->GetXaxis()->SetLabelSize(plot.getStdPixelSize_XYLabel());
    hist->GetYaxis()->SetNdivisions(505);
    hist->GetYaxis()->SetLabelFont(43);
    hist->GetYaxis()->SetLabelOffset(plot.getStdOffset_YLabel());
    hist->GetYaxis()->SetLabelSize(plot.getStdPixelSize_XYLabel());

    if (addRatioPanel) hist->GetXaxis()->SetLabelSize(0);
  }
  for (auto& pp:hist_tg_map){
    TGraphAsymmErrors* tg = pp.second;

    tg->GetXaxis()->SetNdivisions(505);
    tg->GetXaxis()->SetLabelFont(43);
    tg->GetXaxis()->SetLabelOffset(plot.getStdOffset_XLabel());
    tg->GetXaxis()->SetLabelSize(plot.getStdPixelSize_XYLabel());
    tg->GetYaxis()->SetNdivisions(505);
    tg->GetYaxis()->SetLabelFont(43);
    tg->GetYaxis()->SetLabelOffset(plot.getStdOffset_YLabel());
    tg->GetYaxis()->SetLabelSize(plot.getStdPixelSize_XYLabel());

    if (addRatioPanel) tg->GetXaxis()->SetLabelSize(0);
  }

  TH1F* hdummy_ratio = nullptr;
  std::vector<TGraphAsymmErrors*> tgratios;
  std::vector<TString> stropts_ratios;
  if (addRatioPanel){
    double ymin_ratio = 0;
    double ymax_ratio = 2;
    hdummy_ratio = dynamic_cast<TH1F*>(hratios.front()->Clone("hdummy_ratio")); hdummy_ratio->Reset("ICESM");
    hdummy_ratio->GetYaxis()->SetRangeUser(ymin_ratio, ymax_ratio);

    hdummy_ratio->GetXaxis()->SetNdivisions(505);
    hdummy_ratio->GetXaxis()->SetLabelFont(43);
    hdummy_ratio->GetXaxis()->SetLabelOffset(plot.getStdOffset_XLabel());
    hdummy_ratio->GetXaxis()->SetLabelSize(plot.getStdPixelSize_XYLabel());
    hdummy_ratio->GetYaxis()->SetNdivisions(505);
    hdummy_ratio->GetYaxis()->SetLabelFont(43);
    hdummy_ratio->GetYaxis()->SetLabelOffset(plot.getStdOffset_YLabel());
    hdummy_ratio->GetYaxis()->SetLabelSize(plot.getStdPixelSize_XYLabel());
  }

  TPad* pad_hists = plot.getInsidePanels().front().back();
  TPad* pad_ratios = (addRatioPanel ? plot.getInsidePanels().front().front() : nullptr);

  // Add x and y titles
  TPad* pad_xtitle = plot.getBorderPanels().at(0); pad_xtitle->cd();
  TLatex* xtitle = new TLatex(); plot.addText(xtitle);
  xtitle->SetTextAlign(22);
  xtitle->SetTextFont(PlotCanvas::getStdFont_XYTitle());
  xtitle->SetTextSize(plot.getStdPixelSize_XYTitle());
  plot.addText(xtitle->DrawLatexNDC(0.5, 0.5, varlabel));

  TPad* pad_ytitle = plot.getBorderPanels().at(1); pad_ytitle->cd();
  TLatex* ytitle = new TLatex(); plot.addText(ytitle);
  ytitle->SetTextAlign(22);
  ytitle->SetTextFont(PlotCanvas::getStdFont_XYTitle());
  ytitle->SetTextSize(plot.getStdPixelSize_XYTitle());
  ytitle->SetTextAngle(90);
  plot.addText(ytitle->DrawLatexNDC(0.5, (addRatioPanel ? 1.-0.5/1.4 : 0.5), quantlabel));
  if (addRatioPanel) plot.addText(ytitle->DrawLatexNDC(0.5, 0.13/1.4, "Ratio"));

  pad_hists->cd();
  if (useLogY) pad_hists->SetLogy(true);

  constexpr double legend_ymax = 0.90;
  double legend_pixelsize = plot.getStdPixelSize_XYTitle();
  double legend_reldy = legend_pixelsize/npixels_pad_xy*1.3;
  TLegend* legend = new TLegend(
    0.55,
    legend_ymax-legend_reldy*float(nplottables),
    0.90,
    legend_ymax,
    "", "NDC"
  );
  legend->SetBorderSize(0);
  legend->SetTextFont(43);
  legend->SetTextAlign(12);
  legend->SetTextSize(legend_pixelsize);
  legend->SetLineColor(1);
  legend->SetLineStyle(1);
  legend->SetLineWidth(1);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  plot.addLegend(legend);
  TText* text;

  pad_hists->cd();

  bool firstHist = true;
  for (size_t is=0; is<hlist.size(); is++){
    TH1F* hist = hlist.at(is);
    TString const& hlabel = hlabels.at(is);

    bool hasGraph = hist_tg_map.find(hist)!=hist_tg_map.end();
    bool hasErrors = hHasErrors.at(is);

    TString stropt = drawopts;
    if (!hasErrors) stropt = "hist";

    hist->SetTitle("");
    if (hlabel!=""){
      if (!hasGraph){
        if (stropt=="hist") legend->AddEntry(hist, hlabel, "f");
        else legend->AddEntry(hist, hlabel, "elp");
      }
      else legend->AddEntry(hist_tg_map[hist], hlabel, "e1p");
    }

    if (hasGraph) continue;

    if (firstHist){
      if (!hasGraph) hist->Draw(stropt);
      //else hist_tg_map[hist]->Draw("ae1p");
      firstHist = false;
    }
    else{
      if (!hasGraph) hist->Draw(stropt+"same");
      //else hist_tg_map[hist]->Draw("e1psame");
    }
  }

  pad_hists->cd();

  // Re-draw data.
  // Draw in reverse in order to make sure real data is drawn the last.
  for (int is=hlist.size()-1; is>=0; is--){
    TH1F* hist = hlist.at(is);
    TString hname = hist->GetName();
    if (!hname.Contains("Data")) continue;
    bool hasGraph = hist_tg_map.find(hist)!=hist_tg_map.end();
    bool hasErrors = hHasErrors.at(is);
    TString stropt = drawopts;
    if (!hasErrors) stropt = "hist";
    if (!hasGraph) hist->Draw(stropt+"same");
    //else hist_tg_map[hist]->Draw("e1psame");
  }

  pad_hists->cd();
  legend->Draw();

  pad_hists->cd();
  TLatex* selectionstitle = new TLatex(); plot.addText(selectionstitle);
  selectionstitle->SetTextAlign(12);
  selectionstitle->SetTextFont(43);
  selectionstitle->SetTextSize(legend_pixelsize);
  {
    double pt_ymax = legend_ymax;
    double pt_dy = legend_reldy;
    for (auto const& strSel:selectionList){
      plot.addText(selectionstitle->DrawLatexNDC(0.22/(1.+0.25+0.0625)+0.05, pt_ymax-pt_dy/2., strSel));
      pt_ymax -= pt_dy;
    }
  }

  if (pad_ratios){
    pad_ratios->cd();
    hdummy_ratio->SetTitle("");
    hdummy_ratio->Draw("hist");
    for (auto& hh:hratios){
      hh->SetTitle("");
      hh->Draw("histsame");
    }
  }

  plot.update();
  plot.save(coutput_main, "png");
  plot.save(coutput_main, "pdf");

  delete hdummy_ratio;
  for (auto& pp:hist_tg_map) delete pp.second;
}


void makeBWRewgtPlots(TString strSampleSet, TString strdate){
  bool isGG = strSampleSet.Contains("GGH");
  bool isVBF = strSampleSet.Contains("VBFH");

  TString cinput_main = "output/BWRewgtPlots/" + strdate;
  TString coutput_main = cinput_main;
  TString cinput = cinput_main + "/hists_" + strSampleSet + ".root";
  TFile* finput = TFile::Open(cinput, "read");

  TH1F* h_norewgt = (TH1F*) finput->Get("h_norewgt"); h_norewgt->SetLineColor(kBlack); h_norewgt->SetMarkerColor(kBlack); h_norewgt->GetYaxis()->SetTitle("d#sigma / dm_{ZZ} (fb/GeV)");
  TH1F* h_BRrewgt = (TH1F*) finput->Get("h_BRrewgt"); h_BRrewgt->SetLineColor(kBlue); h_BRrewgt->SetMarkerColor(kBlue); h_BRrewgt->GetYaxis()->SetTitle("d#sigma / dm_{ZZ} (fb/GeV)");
  TH1F* h_fullproprewgt = (TH1F*) finput->Get("h_fullproprewgt"); h_fullproprewgt->SetLineColor(kViolet); h_fullproprewgt->SetMarkerColor(kViolet); h_fullproprewgt->GetYaxis()->SetTitle("d#sigma / dm_{ZZ} (fb/GeV)");

  TH1F* hratio_BRrewgt = (TH1F*) h_BRrewgt->Clone("hratio_BRrewgt"); hratio_BRrewgt->Divide(h_norewgt); hratio_BRrewgt->GetXaxis()->SetTitle(""); hratio_BRrewgt->GetYaxis()->SetTitle("");
  TH1F* hratio_fullproprewgt = (TH1F*) h_fullproprewgt->Clone("hratio_fullproprewgt"); hratio_fullproprewgt->Divide(h_BRrewgt); hratio_fullproprewgt->GetXaxis()->SetTitle(""); hratio_fullproprewgt->GetYaxis()->SetTitle("");

  {
    std::vector<TH1F*> hlist{
      h_norewgt,
      h_BRrewgt,
      h_fullproprewgt
    };
    std::vector<TH1F*> hratios{
      hratio_BRrewgt,
      hratio_fullproprewgt
    };
    std::vector<TString> hlabels{
      "Default POWHEG",
      "+ BR evolution",
      "+ CPS#rightarrowBW"
    };
    makePlot(
      coutput_main,
      Form("c_%s", strSampleSet.Data()),
      hlist,
      hratios,
      hlabels,
      (isGG ? "gg#rightarrowH#rightarrow2l2#nu" : "VBF, H#rightarrow2l2#nu"),
      "hist",
      false, false, -1.
    );
  }

  finput->Close();
}


void computeIntegratedXsecAfterSelection(
  TString strSampleSet, TString period, TString strdate,
  double minMass, double maxMass,
  bool use4lSel = false, int nleps_req = -1,
  bool useJetCuts = false, bool useGenJetExtraCuts = false,
  bool useOnlyWrongComb = false
){
  SystematicsHelpers::SystematicVariationTypes const theGlobalSyst = SystematicsHelpers::sNominal;

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  bool const isGG = strSampleSet.Contains("GluGluH") || strSampleSet.Contains("GGH");
  bool const isVH = strSampleSet.Contains("WminusH") || strSampleSet.Contains("WplusH") || strSampleSet.Contains("ZH") || strSampleSet.Contains("HZJ");
  bool const isPowheg = strSampleSet.Contains("POWHEG");
  bool const acquireMEs = isPowheg;

  // Set output directory
  TString cinput_main = "output/ReweightedGenTrees/" + strdate + "/" + period;

  TDirectory* curdir = gDirectory;
  curdir->cd();

  TFile* foutput = TFile::Open("bla.root", "recreate");
  TH1F htmp("hh", "", 900, 100, 1000);
  curdir->cd();

  // Create the output file
  TString cinput = SampleHelpers::getSampleIdentifier(strSampleSet);
  HelperFunctions::replaceString(cinput, "_MINIAODSIM", "");
  HelperFunctions::replaceString(cinput, "_MINIAOD", "");
  TString strinput = cinput_main + "/" + cinput + "_" + SystematicsHelpers::getSystName(theGlobalSyst).data() + ".root";
  IVYout << "Acquiring input " << strinput << "..." << endl;
  BaseTree* tin = new BaseTree(strinput, "SkimTree", "", "");
  if (!tin->isValid()){
    IVYerr << "\t- Failed to acquire." << endl;
    delete tin;
    exit(1);
  }

#define BRANCH_COMMAND(TYPE, NAME) TYPE* NAME = nullptr; tin->bookBranch<TYPE>(#NAME, 0);
  BRANCH_SCALAR_COMMANDS;
#undef BRANCH_COMMAND
#define BRANCH_COMMAND(TYPE, NAME) std::vector<TYPE>** NAME=nullptr; tin->bookBranch<std::vector<TYPE>*>(#NAME, nullptr);
  BRANCH_VECTOR_COMMANDS;
#undef BRANCH_COMMAND
#define BRANCH_COMMAND(TYPE, NAME) tin->getValRef(#NAME, NAME);
  BRANCH_COMMANDS;
#undef BRANCH_COMMAND

  std::unordered_map<TString, float*> ME_Kfactor_values;
  std::vector<TString> allbranchnames;
  tin->getValidBranchNamesWithoutAlias(allbranchnames, false);
  for (auto const& bname:allbranchnames){
    if (
      (bname.BeginsWith("p_") && (bname.Contains("JHUGen") || bname.Contains("MCFM")))
      ||
      bname.BeginsWith("p_Gen")
      ||
      bname.Contains("LHECandMass")
      ||
      bname.BeginsWith("KFactor")
      ){
      tin->bookBranch<float>(bname, -1.f);
      ME_Kfactor_values[bname] = nullptr;
      tin->getValRef(bname, ME_Kfactor_values[bname]);
    }
  }

  tin->silenceUnused();

  float* val_Kfactor_QCD = nullptr;
  float* val_ME_SIG = nullptr;
  float* val_ME_CPS = nullptr;
  if (isGG){
    val_Kfactor_QCD = ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_Nominal")->second;
    if (acquireMEs){
      val_ME_SIG = ME_Kfactor_values.find("p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_MCFM")->second;
      val_ME_CPS = ME_Kfactor_values.find("p_Gen_CPStoBWPropRewgt")->second;
    }
  }
  else{
    if (acquireMEs){
      val_ME_SIG = ME_Kfactor_values.find("p_Gen_JJEW_SIG_ghv1_1_MCFM")->second;
      val_ME_CPS = ME_Kfactor_values.find("p_Gen_CPStoBWPropRewgt")->second;
    }
  }

  double sum_wgts = 0;
  double sum_wgts_passGenCompSelection = 0;
  double sum_wgts_accepted = 0;
  int nEntries = tin->getNEvents();
  IVYout << "Looping over " << nEntries << " events:" << endl;
  for (int ev=0; ev<nEntries; ev++){
    tin->getEvent(ev);
    HelperFunctions::progressbar(ev, nEntries);

    double wgt =
      (*event_wgt) * (*event_wgt_adjustment_NNPDF30)
      * (*sample_wgt) * (*invalidReweightingWgts ? 0.f : 1.f)
      * (val_Kfactor_QCD ? *val_Kfactor_QCD : 1.f)
      * (val_ME_SIG ? *val_ME_SIG : 1.f) * (val_ME_CPS ? *val_ME_CPS : 1.f);
    sum_wgts += wgt;

    if (!(use4lSel && isPowheg && isVH)){
      if (!(*passGenCompSelection)) continue;
    }
    sum_wgts_passGenCompSelection += wgt;

    bool pass_LHEFinalState = true;

    std::vector<ParticleObject> mothers; mothers.reserve(2);
    std::vector<ParticleObject> leptons; leptons.reserve((*lheQCDLOparticles_id)->size());
    std::vector<ParticleObject> jets; jets.reserve((*lheQCDLOparticles_id)->size());
    ParticleObject::LorentzVector_t leptons_sump4;
    for (unsigned int ip=0; ip<(*lheQCDLOparticles_id)->size(); ip++){
      auto const& st = (*lheQCDLOparticles_status)->at(ip);
      auto const& pid = (*lheQCDLOparticles_id)->at(ip);
      ParticleObject::LorentzVector_t pp4((*lheQCDLOparticles_px)->at(ip), (*lheQCDLOparticles_py)->at(ip), (*lheQCDLOparticles_pz)->at(ip), (*lheQCDLOparticles_E)->at(ip));

      if (st<0){
        mothers.emplace_back(pid, pp4);
        continue;
      }

      if ((PDGHelpers::isALepton(pid) && std::abs(pid)!=15) || PDGHelpers::isANeutrino(pid)){
        leptons.emplace_back(pid, pp4);
        leptons_sump4 += pp4;
      }
      else if (PDGHelpers::isAJet(pid)) jets.emplace_back(pid, pp4);
    }

    if (useOnlyWrongComb){
      std::vector<std::pair<ParticleObject*, ParticleObject*>> dilepton_pairs;
      std::vector<std::pair<ParticleObject*, ParticleObject*>> dinu_pairs;
      std::vector<std::pair<ParticleObject*, ParticleObject*>> dijet_pairs;
      for (auto it_l1=leptons.begin(); it_l1!=leptons.end();it_l1++){
        for (auto it_l2=it_l1+1; it_l2!=leptons.end(); it_l2++){
          if (it_l1->pdgId()+it_l2->pdgId()==0){
            if (PDGHelpers::isANeutrino(it_l1->pdgId())) dinu_pairs.emplace_back(&(*it_l1), &(*it_l2));
            else dilepton_pairs.emplace_back(&(*it_l1), &(*it_l2));
          }
        }
      }
      for (auto it_j1=jets.begin(); it_j1!=jets.end(); it_j1++){
        for (auto it_j2=it_j1+1; it_j2!=jets.end(); it_j2++){
          if (it_j1->pdgId()+it_j2->pdgId()==0) dijet_pairs.emplace_back(&(*it_j1), &(*it_j2));
        }
      }
      bool has_llqq = false;
      bool has_llnn = false;
      for (auto const& pl:dilepton_pairs){
        for (auto const& pj:dijet_pairs){
          if (std::abs((pl.first->p4()+pl.second->p4()+pj.first->p4()+pj.second->p4()).M()-125.)<0.05){
            has_llqq = true;
            break;
          }
        }
      }
      for (auto const& pl:dilepton_pairs){
        for (auto const& pn:dinu_pairs){
          if (std::abs((pl.first->p4()+pl.second->p4()+pn.first->p4()+pn.second->p4()).M()-125.)<0.05){
            has_llqq = true;
            break;
          }
        }
      }
      pass_LHEFinalState = (has_llqq || has_llnn);
    }

    if (!pass_LHEFinalState) continue;

    cms3_id_t motherVid = -9000;
    if (!isGG && mothers.size()==2) motherVid = PDGHelpers::getCoupledVertex(mothers.front().pdgId(), mothers.back().pdgId());

    bool pass_leptons_pt_eta = true;
    bool pass_leptons_min_mll = true;
    bool pass_leptons_extraSel = true;
    double min_mll = -1;
    unsigned char nleps_pt10 = 0;
    unsigned char nleps_pt20 = 0;
    double smallestZdiff = -1;
    double mll_smallestZdiff = -1;
    char idx_lepZ11 = -1;
    char idx_lepZ12 = -1;
    char idx_lepZ21 = -1;
    char idx_lepZ22 = -1;
    int nleps = 0;
    for (unsigned int il=0; il<leptons.size(); il++){
      auto const& lep1 = leptons.at(il);
      if (!use4lSel){
        if (lep1.pt()<7. || std::abs(lep1.eta())>=2.4){
          pass_leptons_pt_eta = false;
          break;
        }
      }
      else{
        bool const isMuon = (std::abs(lep1.pdgId())==13);
        if (lep1.pt()<(isMuon ? 5. : 7.) || std::abs(lep1.eta())>=(isMuon ? 2.4 : 2.5)){
          pass_leptons_pt_eta = false;
          break;
        }
      }
      if (lep1.pt()>=10.){
        nleps_pt10++;
        if (lep1.pt()>=20.){
          nleps_pt20++;
        }
      }

      nleps++;

      for (unsigned int jl=il+1; jl<leptons.size(); jl++){
        auto const& lep2 = leptons.at(jl);
        if (lep1.pdgId()==-lep2.pdgId() && std::abs(lep1.pdgId())%2==1){
          double mll = (lep1.p4() + lep2.p4()).M();
          if (min_mll<0. || mll<min_mll) min_mll = mll;
          double Zdiff = std::abs(mll - PDGHelpers::Zmass);
          if (smallestZdiff<0. || Zdiff<smallestZdiff){
            smallestZdiff = Zdiff;
            mll_smallestZdiff = mll;
            idx_lepZ11 = il;
            idx_lepZ12 = jl;
            if (lep1.pdgId()<lep2.pdgId()) std::swap(idx_lepZ11, idx_lepZ12);
          }
        }
      }
    }
    pass_leptons_min_mll = (min_mll>=4.);
    if (!pass_leptons_pt_eta || !pass_leptons_min_mll) continue;

    {
      std::vector<std::pair<unsigned int, unsigned int>> Z2pairs;
      for (unsigned int il=0; il<leptons.size(); il++){
        if (static_cast<char>(il)==idx_lepZ11 || static_cast<char>(il)==idx_lepZ12) continue;
        auto const& lep1 = leptons.at(il);
        for (unsigned int jl=il+1; jl<leptons.size(); jl++){
          if (static_cast<char>(jl)==idx_lepZ11 || static_cast<char>(jl)==idx_lepZ12) continue;
          auto const& lep2 = leptons.at(jl);
          if (lep1.pdgId()==-lep2.pdgId() && std::abs(lep1.pdgId())%2==1){
            double mll = (lep1.p4() + lep2.p4()).M();
            if (!use4lSel || mll>12.) Z2pairs.emplace_back(il, jl);
          }
        }
      }
      if (Z2pairs.empty()) continue;
      double highest_ptsum = -1;
      for (auto const& pp:Z2pairs){
        double ptsum = leptons.at(pp.first).pt() + leptons.at(pp.second).pt();
        if (highest_ptsum<ptsum){
          highest_ptsum = ptsum;
          idx_lepZ21 = pp.first;
          idx_lepZ22 = pp.second;
        }
      }
    }

    pass_leptons_extraSel = !use4lSel || (idx_lepZ11*idx_lepZ12>=0 && idx_lepZ21*idx_lepZ22>=0 && mll_smallestZdiff>=40. && nleps_pt10>=2 && nleps_pt20>=1);
    if (!pass_leptons_extraSel) continue;
    if (use4lSel && (nleps_req>0 && nleps!=nleps_req)) continue;

    double const invmass = (!use4lSel ? leptons_sump4.M() : (leptons.at(idx_lepZ11).p4()+leptons.at(idx_lepZ12).p4()+leptons.at(idx_lepZ21).p4()+leptons.at(idx_lepZ22).p4()).M());

    bool pass_jets_pt_eta = true;
    bool pass_jets_mjj = true;
    bool pass_jets_dRjj = true;
    bool pass_jets_extraSel = true;
    double min_mjj = -1;
    double min_dRjj = -1;
    if (!isGG){
      if (jets.size()==0){
        pass_jets_pt_eta = pass_jets_mjj = false;
      }
      if (useGenJetExtraCuts) pass_jets_extraSel = false;
      for (unsigned int ip=0; ip<jets.size(); ip++){
        auto const& jet1 = jets.at(ip);
        if (jet1.pt()<30. || std::abs(jet1.eta())>=4.7){
          pass_jets_pt_eta = false;
          break;
        }
        for (unsigned int jp=ip+1; jp<jets.size(); jp++){
          auto const& jet2 = jets.at(jp);
          double mjj = (jet1.p4() + jet2.p4()).M();
          if (min_mjj<0. || mjj<min_mjj) min_mjj = mjj;

          double dRjj = jet1.deltaR(jet2);
          if (min_dRjj<0. || dRjj<min_dRjj) min_dRjj = dRjj;

          if (useGenJetExtraCuts){
            cms3_id_t assocVid = PDGHelpers::getCoupledVertex(jet1.pdgId(), jet2.pdgId());
            if (
                mjj>130.
                ||
                (
                  motherVid==assocVid
                  &&
                  (
                    (PDGHelpers::isAWBoson(assocVid) && std::abs(mjj-PDGHelpers::Wmass)<10.)
                    ||
                    (PDGHelpers::isAZBoson(assocVid) && std::abs(mjj-PDGHelpers::Zmass)<10.)
                  )
                )
              ) pass_jets_extraSel = true;
          }
        }
      }
      if (min_mjj>=0.) pass_jets_mjj = (min_mjj>=50.);
      if (min_dRjj>=0.) pass_jets_dRjj = (min_dRjj>=0.4);
    }

    // If we are using 4l cuts, we are intested in the jet-inclusive yield, so we do not look at jets.
    if (!use4lSel && useJetCuts){
      if (!pass_jets_pt_eta || !pass_jets_mjj || !pass_jets_dRjj || !pass_jets_extraSel) continue;
    }

    htmp.Fill(invmass, wgt);

    if (invmass<minMass || invmass>=maxMass) continue;

    sum_wgts_accepted += wgt;
  }

  IVYout << "Sum of weights: " << sum_wgts_accepted << " / " << sum_wgts_passGenCompSelection << " / " << sum_wgts << endl;
  IVYout << "\t- Fraction of events passing gen. cuts: " << sum_wgts_passGenCompSelection / sum_wgts << endl;
  IVYout << "\t- Fraction of events accepted after gen. cuts: " << sum_wgts_accepted / sum_wgts_passGenCompSelection << endl;

  delete tin;

  foutput->WriteTObject(&htmp);
  foutput->Close();
}
