#include <cctype>
#include <algorithm>

#include <CommonTools/UtilAlgos/interface/TFileService.h>

#include <CMS3/NtupleMaker/interface/plugins/MCUtilities.h>
#include <CMS3/NtupleMaker/interface/plugins/CMS3Ntuplizer.h>
#include <CMS3/NtupleMaker/interface/FSRCandidateInfo.h>
#include "CMS3/NtupleMaker/interface/VertexSelectionHelpers.h"
#include "CMS3/NtupleMaker/interface/MuonSelectionHelpers.h"
#include "CMS3/NtupleMaker/interface/ElectronSelectionHelpers.h"
#include "CMS3/NtupleMaker/interface/PhotonSelectionHelpers.h"
#include "CMS3/NtupleMaker/interface/FSRSelectionHelpers.h"
#include "CMS3/NtupleMaker/interface/AK4JetSelectionHelpers.h"
#include "CMS3/NtupleMaker/interface/AK8JetSelectionHelpers.h"
#include "CMS3/NtupleMaker/interface/IsotrackSelectionHelpers.h"
#include <CMS3/NtupleMaker/interface/CMS3ObjectHelpers.h>

#include <CMSDataTools/AnalysisTree/interface/HelperFunctionsCore.h>

#include "MELAStreamHelpers.hh"


using namespace std;
using namespace edm;
using namespace MELAStreamHelpers;


CMS3Ntuplizer::CMS3Ntuplizer(const edm::ParameterSet& pset_) :
  pset(pset_),
  outtree(nullptr),
  firstEvent(true),

  year(pset.getParameter<int>("year")),
  treename(pset.getUntrackedParameter<std::string>("treename")),
  isMC(pset.getParameter<bool>("isMC")),
  is80X(pset.getParameter<bool>("is80X")),

  prefiringWeightsTag(pset.getUntrackedParameter<std::string>("prefiringWeightsTag")),
  applyPrefiringWeights(prefiringWeightsTag!=""),

  keepGenParticles(CMS3Ntuplizer::getParticleRecordLevel(pset.getUntrackedParameter<std::string>("keepGenParticles"))),
  keepGenJets(pset.getParameter<bool>("keepGenJets")),

  minNmuons(pset.getParameter<int>("minNmuons")),
  minNelectrons(pset.getParameter<int>("minNelectrons")),
  minNleptons(pset.getParameter<int>("minNleptons")),
  minNphotons(pset.getParameter<int>("minNphotons")),
  minNak4jets(pset.getParameter<int>("minNak4jets")),
  minNak8jets(pset.getParameter<int>("minNak8jets"))
{
  if (year!=2016 && year!=2017 && year!=2018) throw cms::Exception("CMS3Ntuplizer::CMS3Ntuplizer: Year is undefined!");

  muonsToken  = consumes< edm::View<pat::Muon> >(pset.getParameter<edm::InputTag>("muonSrc"));
  electronsToken  = consumes< edm::View<pat::Electron> >(pset.getParameter<edm::InputTag>("electronSrc"));
  photonsToken  = consumes< edm::View<pat::Photon> >(pset.getParameter<edm::InputTag>("photonSrc"));
  ak4jetsToken  = consumes< edm::View<pat::Jet> >(pset.getParameter<edm::InputTag>("ak4jetSrc"));
  ak8jetsToken  = consumes< edm::View<pat::Jet> >(pset.getParameter<edm::InputTag>("ak8jetSrc"));
  isotracksToken  = consumes< edm::View<IsotrackInfo> >(pset.getParameter<edm::InputTag>("isotrackSrc"));
  pfcandsToken  = consumes< edm::View<pat::PackedCandidate> >(pset.getParameter<edm::InputTag>("pfcandSrc"));

  pfmetToken = consumes< METInfo >(pset.getParameter<edm::InputTag>("pfmetSrc"));
  puppimetToken = consumes< METInfo >(pset.getParameter<edm::InputTag>("puppimetSrc"));

  vtxToken = consumes<reco::VertexCollection>(pset.getParameter<edm::InputTag>("vtxSrc"));

  rhoToken  = consumes< double >(pset.getParameter<edm::InputTag>("rhoSrc"));
  triggerInfoToken = consumes< edm::View<TriggerInfo> >(pset.getParameter<edm::InputTag>("triggerInfoSrc"));
  puInfoToken = consumes< std::vector<PileupSummaryInfo> >(pset.getParameter<edm::InputTag>("puInfoSrc"));
  metFilterInfoToken = consumes< METFilterInfo >(pset.getParameter<edm::InputTag>("metFilterInfoSrc"));

  if (applyPrefiringWeights){
    prefiringWeightToken = consumes< double >(edm::InputTag(prefiringWeightsTag, "nonPrefiringProb"));
    prefiringWeightToken_Dn = consumes< double >(edm::InputTag(prefiringWeightsTag, "nonPrefiringProbDown"));
    prefiringWeightToken_Up = consumes< double >(edm::InputTag(prefiringWeightsTag, "nonPrefiringProbUp"));
  }

  genInfoToken = consumes< GenInfo >(pset.getParameter<edm::InputTag>("genInfoSrc"));
  prunedGenParticlesToken = consumes< reco::GenParticleCollection >(pset.getParameter<edm::InputTag>("prunedGenParticlesSrc"));
  packedGenParticlesToken = consumes< pat::PackedGenParticleCollection >(pset.getParameter<edm::InputTag>("packedGenParticlesSrc"));
  genAK4JetsToken = consumes< edm::View<reco::GenJet> >(pset.getParameter<edm::InputTag>("genAK4JetsSrc"));
  genAK8JetsToken = consumes< edm::View<reco::GenJet> >(pset.getParameter<edm::InputTag>("genAK8JetsSrc"));

  this->usesResource("TFileService");
}
CMS3Ntuplizer::~CMS3Ntuplizer(){
  //delete pileUpReweight;
  //delete metCorrHandler;
}


void CMS3Ntuplizer::beginJob(){
  edm::Service<TFileService> fs;
  TTree* tout = fs->make<TTree>(treename, "Selected event summary");
  outtree = std::make_shared<BaseTree>(nullptr, tout, nullptr, nullptr, false);
  outtree->setAcquireTreePossession(false);
  outtree->setAutoSave(0);
}
void CMS3Ntuplizer::endJob(){}


// Convenience macros to easily make and push vector values
#define MAKE_VECTOR_WITH_RESERVE(type_, name_, size_) std::vector<type_> name_; name_.reserve(size_);
#define PUSH_USERINT_INTO_VECTOR(name_) name_.push_back(obj->userInt(#name_));
#define PUSH_USERFLOAT_INTO_VECTOR(name_) name_.push_back(obj->userFloat(#name_));
#define PUSH_VECTOR_WITH_NAME(name_, var_) commonEntry.setNamedVal(TString(name_)+"_"+#var_, var_);


void CMS3Ntuplizer::analyze(edm::Event const& iEvent, const edm::EventSetup& iSetup){
  bool isSelected = true;

  /********************************/
  /* Set the communicator entries */
  /********************************/
  /*
  When naeing variables, try to be conscious of the nanoAOD naming conventions, but do not make a big fuss about them either!
  The latest list of variables are documented at https://cms-nanoaod-integration.web.cern.ch/integration/master-102X/mc102X_doc.html
  */

  // Gen. variables
  std::vector<reco::GenParticle const*> filledPrunedGenParts;
  std::vector<pat::PackedGenParticle const*> filledPackedGenParts;
  std::vector<reco::GenJet const*> filledGenAK4Jets;
  std::vector<reco::GenJet const*> filledGenAK8Jets;
  if (this->isMC) isSelected &= this->fillGenVariables(
    iEvent,
    &filledPrunedGenParts, &filledPackedGenParts,
    &filledGenAK4Jets, &filledGenAK8Jets
  );

  // Vertices
  size_t n_vtxs = this->fillVertices(iEvent, nullptr);
  isSelected &= (n_vtxs>0);

  // Muons
  std::vector<pat::Muon const*> filledMuons;
  size_t n_muons = this->fillMuons(iEvent, &filledMuons);

  // Electrons
  std::vector<pat::Electron const*> filledElectrons;
  size_t n_electrons = this->fillElectrons(iEvent, &filledElectrons);

  // Photons
  std::vector<pat::Photon const*> filledPhotons;
  size_t n_photons = this->fillPhotons(iEvent, &filledPhotons);

  // PF candidates, including FSR information
  /*size_t n_pfcands = */this->fillPFCandidates(
    iEvent,
    &filledMuons, &filledElectrons, &filledPhotons,
    nullptr
  );

  // ak4 jets
  size_t n_ak4jets = this->fillAK4Jets(iEvent, nullptr);

  // ak8 jets
  size_t n_ak8jets = this->fillAK8Jets(iEvent, nullptr);

  // Isolated tracks
  /*size_t n_isotracks = */this->fillIsotracks(iEvent, nullptr);

  // The (data) event should have at least one electron, muon, or photon.
  // If all cuts are -1, passNobjects is true and no filtering on the number of objects is done.
  bool passNobjects = (minNmuons<0 && minNelectrons<0 && minNleptons<0 && minNphotons<0 && minNak4jets<0 && minNak8jets<0);
#define passNobjects_or_statement(n_objs, min_n_objs) if (min_n_objs>=0) passNobjects |= (n_objs>=static_cast<size_t const>(min_n_objs));
  passNobjects_or_statement(n_muons, minNmuons);
  passNobjects_or_statement(n_electrons, minNelectrons);
  passNobjects_or_statement((n_muons+n_electrons), minNleptons);
  passNobjects_or_statement(n_photons, minNphotons);
  passNobjects_or_statement(n_ak4jets, minNak4jets);
  passNobjects_or_statement(n_ak8jets, minNak8jets);
#undef passNobjects_or_statement
  isSelected &= passNobjects;

  // MET info
  isSelected &= this->fillMETVariables(iEvent);

  // Event info
  isSelected &= this->fillEventVariables(iEvent);

  // Trigger info
  isSelected &= this->fillTriggerInfo(iEvent);

  // MET filters
  isSelected &= this->fillMETFilterVariables(iEvent);


  /************************************************************/
  /************************************************************/
  /* NO MORE CALLS TO THE FILL SUBROUTINES BEYOND THIS POINT! */
  /************************************************************/
  /************************************************************/

  commonEntry.setNamedVal("passCommonSkim", isSelected); // Can use this flag to match data and MC selections

  /************************************************/
  /* Record the communicator values into the tree */
  /************************************************/

  // If this is the first event, create the tree branches based on what is available in the commonEntry.
  if (firstEvent){
#define SIMPLE_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=commonEntry.named##name_t##s.begin(); itb!=commonEntry.named##name_t##s.end(); itb++) outtree->putBranch(itb->first, itb->second);
#define VECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=commonEntry.namedV##name_t##s.begin(); itb!=commonEntry.namedV##name_t##s.end(); itb++) outtree->putBranch(itb->first, &(itb->second));
#define DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=commonEntry.namedVV##name_t##s.begin(); itb!=commonEntry.namedVV##name_t##s.end(); itb++) outtree->putBranch(itb->first, &(itb->second));
    SIMPLE_DATA_OUTPUT_DIRECTIVES
    VECTOR_DATA_OUTPUT_DIRECTIVES
    DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVES
#undef SIMPLE_DATA_OUTPUT_DIRECTIVE
#undef VECTOR_DATA_OUTPUT_DIRECTIVE
#undef DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE

    outtree->getSelectedTree()->SetBasketSize("triggers_*", 16384*23);
  }

  // Record whatever is in commonEntry into the tree.
#define SIMPLE_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=commonEntry.named##name_t##s.begin(); itb!=commonEntry.named##name_t##s.end(); itb++) outtree->setVal(itb->first, itb->second);
#define VECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=commonEntry.namedV##name_t##s.begin(); itb!=commonEntry.namedV##name_t##s.end(); itb++) outtree->setVal(itb->first, &(itb->second));
#define DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=commonEntry.namedVV##name_t##s.begin(); itb!=commonEntry.namedVV##name_t##s.end(); itb++) outtree->setVal(itb->first, &(itb->second));
  SIMPLE_DATA_OUTPUT_DIRECTIVES
  VECTOR_DATA_OUTPUT_DIRECTIVES
  DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVES
#undef SIMPLE_DATA_OUTPUT_DIRECTIVE
#undef VECTOR_DATA_OUTPUT_DIRECTIVE
#undef DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE

  // Fill the tree
  if (this->isMC || isSelected) outtree->fill();

  // No longer the first event...
  if (firstEvent) firstEvent = false;
}

void CMS3Ntuplizer::recordGenInfo(edm::Event const& iEvent){
  edm::Handle< GenInfo > genInfoHandle;
  iEvent.getByToken(genInfoToken, genInfoHandle);
  if (!genInfoHandle.isValid()) throw cms::Exception("CMS3Ntuplizer::recordGenInfo: Error getting the gen. info. from the event...");
  const GenInfo& genInfo = *genInfoHandle;

#define SET_GENINFO_VARIABLE(var) commonEntry.setNamedVal(#var, genInfo.var);

  SET_GENINFO_VARIABLE(xsec);
  SET_GENINFO_VARIABLE(xsecerr);

  SET_GENINFO_VARIABLE(xsec_lhe);

  SET_GENINFO_VARIABLE(qscale);
  SET_GENINFO_VARIABLE(alphaS);

  SET_GENINFO_VARIABLE(genmet_met);
  SET_GENINFO_VARIABLE(genmet_metPhi);

  SET_GENINFO_VARIABLE(sumEt);
  SET_GENINFO_VARIABLE(pThat);

  // LHE variations
  SET_GENINFO_VARIABLE(genHEPMCweight_default);
  SET_GENINFO_VARIABLE(genHEPMCweight_NNPDF30);

  SET_GENINFO_VARIABLE(LHEweight_QCDscale_muR1_muF1);
  SET_GENINFO_VARIABLE(LHEweight_QCDscale_muR1_muF2);
  SET_GENINFO_VARIABLE(LHEweight_QCDscale_muR1_muF0p5);
  SET_GENINFO_VARIABLE(LHEweight_QCDscale_muR2_muF1);
  SET_GENINFO_VARIABLE(LHEweight_QCDscale_muR2_muF2);
  SET_GENINFO_VARIABLE(LHEweight_QCDscale_muR2_muF0p5);
  SET_GENINFO_VARIABLE(LHEweight_QCDscale_muR0p5_muF1);
  SET_GENINFO_VARIABLE(LHEweight_QCDscale_muR0p5_muF2);
  SET_GENINFO_VARIABLE(LHEweight_QCDscale_muR0p5_muF0p5);

  SET_GENINFO_VARIABLE(LHEweight_PDFVariation_Up_default);
  SET_GENINFO_VARIABLE(LHEweight_PDFVariation_Dn_default);
  SET_GENINFO_VARIABLE(LHEweight_AsMZ_Up_default);
  SET_GENINFO_VARIABLE(LHEweight_AsMZ_Dn_default);

  SET_GENINFO_VARIABLE(LHEweight_PDFVariation_Up_NNPDF30);
  SET_GENINFO_VARIABLE(LHEweight_PDFVariation_Dn_NNPDF30);
  SET_GENINFO_VARIABLE(LHEweight_AsMZ_Up_NNPDF30);
  SET_GENINFO_VARIABLE(LHEweight_AsMZ_Dn_NNPDF30);

  // Pythis PS weights
  SET_GENINFO_VARIABLE(PythiaWeight_isr_muRoneoversqrt2);
  SET_GENINFO_VARIABLE(PythiaWeight_fsr_muRoneoversqrt2);
  SET_GENINFO_VARIABLE(PythiaWeight_isr_muRsqrt2);
  SET_GENINFO_VARIABLE(PythiaWeight_fsr_muRsqrt2);
  SET_GENINFO_VARIABLE(PythiaWeight_isr_muR0p5);
  SET_GENINFO_VARIABLE(PythiaWeight_fsr_muR0p5);
  SET_GENINFO_VARIABLE(PythiaWeight_isr_muR2);
  SET_GENINFO_VARIABLE(PythiaWeight_fsr_muR2);
  SET_GENINFO_VARIABLE(PythiaWeight_isr_muR0p25);
  SET_GENINFO_VARIABLE(PythiaWeight_fsr_muR0p25);
  SET_GENINFO_VARIABLE(PythiaWeight_isr_muR4);
  SET_GENINFO_VARIABLE(PythiaWeight_fsr_muR4);

  // LHE particles
  SET_GENINFO_VARIABLE(lheparticles_px);
  SET_GENINFO_VARIABLE(lheparticles_py);
  SET_GENINFO_VARIABLE(lheparticles_pz);
  SET_GENINFO_VARIABLE(lheparticles_E);
  SET_GENINFO_VARIABLE(lheparticles_id);
  SET_GENINFO_VARIABLE(lheparticles_status);
  SET_GENINFO_VARIABLE(lheparticles_mother0_index);
  SET_GENINFO_VARIABLE(lheparticles_mother1_index);

#undef SET_GENINFO_VARIABLE

  for (auto const& it:genInfo.LHE_ME_weights) commonEntry.setNamedVal(it.first, it.second);
  for (auto const& it:genInfo.Kfactors) commonEntry.setNamedVal(it.first, it.second);
}
void CMS3Ntuplizer::recordGenParticles(edm::Event const& iEvent, std::vector<reco::GenParticle const*>* filledGenParts, std::vector<pat::PackedGenParticle const*>* filledPackedGenParts){
  const char colName[] = "genparticles";

  edm::Handle<reco::GenParticleCollection> prunedGenParticlesHandle;
  iEvent.getByToken(prunedGenParticlesToken, prunedGenParticlesHandle);
  if (!prunedGenParticlesHandle.isValid()) throw cms::Exception("CMS3Ntuplizer::recordGenParticles: Error getting the pruned gen. particles from the event...");
  std::vector<reco::GenParticle> const* prunedGenParticles = prunedGenParticlesHandle.product();

  edm::Handle<pat::PackedGenParticleCollection> packedGenParticlesHandle;
  iEvent.getByToken(packedGenParticlesToken, packedGenParticlesHandle);
  if (!prunedGenParticlesHandle.isValid()) throw cms::Exception("CMS3Ntuplizer::recordGenParticles: Error getting the packed gen. particles from the event...");
  std::vector<pat::PackedGenParticle> const* packedGenParticles = packedGenParticlesHandle.product();

  // Make a collection of all unique reco::GenParticle pointers
  std::vector<reco::GenParticle const*> allGenParticles; allGenParticles.reserve(prunedGenParticles->size() + packedGenParticles->size());
  // Fill with the pruned collection first
  for (reco::GenParticle const& part:(*prunedGenParticles)){
    int st = part.status();
    int id = std::abs(part.pdgId());
    if (
      (this->keepGenParticles==kReducedFinalStatesAndHardProcesses && !part.isHardProcess() && st!=1)
      ||
      ((this->keepGenParticles==kAllFinalStates || this->keepGenParticles==kReducedFinalStates) && st!=1)
      ) continue;
    else if (
      (this->keepGenParticles==kReducedFinalStates || this->keepGenParticles==kReducedFinalStatesAndHardProcesses)
      &&
      st==1
      && !(
        (id>=11 && id<=16) // Neutral or charged leptons
        ||
        id==22 // Photons
        )
      ) continue;
    allGenParticles.push_back(&part);
  }

  // Get the packed gen. particles unique from the pruned collection
  // Adapted from GeneratorInterface/RivetInterface/plugins/MergedGenParticleProducer.cc
  std::vector<pat::PackedGenParticle const*> uniquePackedGenParticles; uniquePackedGenParticles.reserve(packedGenParticles->size());
  for (pat::PackedGenParticle const& packedGenParticle:(*packedGenParticlesHandle)){
    int st = packedGenParticle.status();
    int id_signed = packedGenParticle.pdgId();
    int id = std::abs(id_signed);

    double match_ref = -1;
    for (reco::GenParticle const& prunedGenParticle:(*prunedGenParticles)){
      if (prunedGenParticle.status() != 1 || packedGenParticle.pdgId() != id_signed) continue;

      double euc_dot_prod = packedGenParticle.px()*prunedGenParticle.px() + packedGenParticle.py()*prunedGenParticle.py() + packedGenParticle.pz()*prunedGenParticle.pz() + packedGenParticle.energy()*prunedGenParticle.energy();
      double comp_ref = prunedGenParticle.px()*prunedGenParticle.px() + prunedGenParticle.py()*prunedGenParticle.py() + prunedGenParticle.pz()*prunedGenParticle.pz() + prunedGenParticle.energy()*prunedGenParticle.energy();
      double match_ref_tmp = std::abs(euc_dot_prod/comp_ref - 1.);
      if (match_ref_tmp<1e-5 && (match_ref<0. || match_ref_tmp<match_ref)) match_ref = match_ref_tmp;
    }

    // Record if NOT matched to any pruned gen. particle.
    if (match_ref<0.){
      if (
        ((this->keepGenParticles==kAllFinalStates || this->keepGenParticles==kReducedFinalStates || this->keepGenParticles==kReducedFinalStatesAndHardProcesses) && st!=1)
        ) continue;
      else if (
        (this->keepGenParticles==kReducedFinalStates || this->keepGenParticles==kReducedFinalStatesAndHardProcesses)
        &&
        st==1
        && !(
          (id>=11 && id<=16) // Neutral or charged leptons
          /*
          // Disable photons in packed candidates for reduced final states because there are a lot of them.
          ||
          id==22 // Photons
          */
          )
        ) continue;
      uniquePackedGenParticles.push_back(&packedGenParticle);
    }
  }

  // Add the mothers of the packed gen. particles to the bigger collection
  if (this->keepGenParticles==kAll){
    for (pat::PackedGenParticle const* part:uniquePackedGenParticles) MCUtilities::getAllMothers(part, allGenParticles, false);
  }

  // Make the variables to record
  // Size of the variable collections are known at this point.
  size_t n_objects = allGenParticles.size() + uniquePackedGenParticles.size();
  if (filledGenParts) filledGenParts->reserve(allGenParticles.size());
  if (filledPackedGenParts) filledPackedGenParts->reserve(uniquePackedGenParticles.size());

  MAKE_VECTOR_WITH_RESERVE(float, pt, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, eta, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, phi, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, mass, n_objects);

  MAKE_VECTOR_WITH_RESERVE(bool, is_packed, n_objects);

  // a) isPromptFinalState(): is particle prompt (not from hadron, muon, or tau decay) and final state
  // b) isPromptDecayed(): is particle prompt (not from hadron, muon, or tau decay) and decayed
  // c) isDirectPromptTauDecayProductFinalState(): this particle is a direct decay product of a prompt tau and is final state
  //    (eg an electron or muon from a leptonic decay of a prompt tau)
  // d) isHardProcess(): this particle is part of the hard process
  // e) fromHardProcessFinalState(): this particle is the final state direct descendant of a hard process particle
  // f) fromHardProcessDecayed(): this particle is the decayed direct descendant of a hard process particle such as a tau from the hard process
  // g) isDirectHardProcessTauDecayProductFinalState(): this particle is a direct decay product of a hardprocess tau and is final state
  //    (eg an electron or muon from a leptonic decay of a tau from the hard process)
  // h) fromHardProcessBeforeFSR(): this particle is the direct descendant of a hard process particle of the same pdg id.
  //    For outgoing particles the kinematics are those before QCD or QED FSR
  //    This corresponds roughly to status code 3 in pythia 6
  //    This is the most complex and error prone of all the flags and you are strongly encouraged
  //    to consider using the others to fill your needs.
  // i) isLastCopy(): this particle is the last copy of the particle in the chain  with the same pdg id
  //    (and therefore is more likely, but not guaranteed, to carry the final physical momentum)
  // j) isLastCopyBeforeFSR(): this particle is the last copy of the particle in the chain with the same pdg id
  //    before QED or QCD FSR (and therefore is more likely, but not guaranteed, to carry the momentum after ISR)
  MAKE_VECTOR_WITH_RESERVE(bool, isPromptFinalState, n_objects); // (a)
  MAKE_VECTOR_WITH_RESERVE(bool, isPromptDecayed, n_objects); // (b)
  MAKE_VECTOR_WITH_RESERVE(bool, isDirectPromptTauDecayProductFinalState, n_objects); // (c)
  MAKE_VECTOR_WITH_RESERVE(bool, isHardProcess, n_objects); // (d)
  MAKE_VECTOR_WITH_RESERVE(bool, fromHardProcessFinalState, n_objects); // (e)
  MAKE_VECTOR_WITH_RESERVE(bool, fromHardProcessDecayed, n_objects); // (f)
  MAKE_VECTOR_WITH_RESERVE(bool, isDirectHardProcessTauDecayProductFinalState, n_objects); // (g)
  MAKE_VECTOR_WITH_RESERVE(bool, fromHardProcessBeforeFSR, n_objects); // (h)
  MAKE_VECTOR_WITH_RESERVE(bool, isLastCopy, n_objects); // (i)
  MAKE_VECTOR_WITH_RESERVE(bool, isLastCopyBeforeFSR, n_objects); // (j)

  MAKE_VECTOR_WITH_RESERVE(int, id, n_objects);
  MAKE_VECTOR_WITH_RESERVE(int, status, n_objects);
  MAKE_VECTOR_WITH_RESERVE(int, mom0_index, n_objects);
  MAKE_VECTOR_WITH_RESERVE(int, mom1_index, n_objects);

  // Record all reco::GenParticle objects
  for (reco::GenParticle const* obj:allGenParticles){
    if (filledGenParts && obj->status()==1) filledGenParts->push_back(obj);

    pt.push_back(obj->pt());
    eta.push_back(obj->eta());
    phi.push_back(obj->phi());
    mass.push_back(obj->mass());

    is_packed.push_back(false);

    isPromptFinalState.push_back(obj->isPromptFinalState()); // (a)
    isPromptDecayed.push_back(obj->isPromptDecayed()); // (b)
    isDirectPromptTauDecayProductFinalState.push_back(obj->isDirectPromptTauDecayProductFinalState()); // (c)
    isHardProcess.push_back(obj->isHardProcess()); // (d)
    fromHardProcessFinalState.push_back(obj->fromHardProcessFinalState()); // (e)
    fromHardProcessDecayed.push_back(obj->fromHardProcessDecayed()); // (f)
    isDirectHardProcessTauDecayProductFinalState.push_back(obj->isDirectHardProcessTauDecayProductFinalState()); // (g)
    fromHardProcessBeforeFSR.push_back(obj->fromHardProcessBeforeFSR()); // (h)
    isLastCopy.push_back(obj->isLastCopy()); // (i)
    isLastCopyBeforeFSR.push_back(obj->isLastCopyBeforeFSR()); // (j)

    id.push_back(obj->pdgId());
    status.push_back(obj->status());

    std::vector<const reco::GenParticle*> mothers;
    if (this->keepGenParticles==kAll) MCUtilities::getAllMothers(obj, mothers, false);
    if (mothers.size()>0){
      const reco::GenParticle* mom = mothers.at(0);
      int index=-1;
      for (reco::GenParticle const* tmpobj:allGenParticles){
        index++;
        if (tmpobj == obj) continue;
        if (mom == tmpobj) break;
      }
      mom0_index.push_back(index);
    }
    else mom0_index.push_back(-1);
    if (mothers.size()>1){
      const reco::GenParticle* mom = mothers.at(1);
      int index=-1;
      for (reco::GenParticle const* tmpobj:allGenParticles){
        index++;
        if (tmpobj == obj) continue;
        if (mom == tmpobj) break;
      }
      mom1_index.push_back(index);
    }
    else mom1_index.push_back(-1);
  }
  // Record the remaining unique pat::PackedGenParticle objects
  for (pat::PackedGenParticle const* obj:uniquePackedGenParticles){
    if (filledPackedGenParts && obj->status()==1) filledPackedGenParts->push_back(obj);

    pt.push_back(obj->pt());
    eta.push_back(obj->eta());
    phi.push_back(obj->phi());
    mass.push_back(obj->mass());

    is_packed.push_back(true);

    isPromptFinalState.push_back(false); // (a)
    isPromptDecayed.push_back(false); // (b)
    isDirectPromptTauDecayProductFinalState.push_back(false); // (c)
    isHardProcess.push_back(false); // (d)
    fromHardProcessFinalState.push_back(false); // (e)
    fromHardProcessDecayed.push_back(false); // (f)
    isDirectHardProcessTauDecayProductFinalState.push_back(false); // (g)
    fromHardProcessBeforeFSR.push_back(false); // (h)
    isLastCopy.push_back(false); // (i)
    isLastCopyBeforeFSR.push_back(false); // (j)

    id.push_back(obj->pdgId());
    status.push_back(obj->status());

    std::vector<const reco::GenParticle*> mothers;
    if (this->keepGenParticles==kAll) MCUtilities::getAllMothers(obj, mothers, false);
    if (mothers.size()>0){
      const reco::GenParticle* mom = mothers.at(0);
      int index=-1;
      for (reco::GenParticle const* tmpobj:allGenParticles){
        index++;
        if (mom == tmpobj) break;
      }
      mom0_index.push_back(index);
    }
    else mom0_index.push_back(-1);
    if (mothers.size()>1){
      const reco::GenParticle* mom = mothers.at(1);
      int index=-1;
      for (reco::GenParticle const* tmpobj:allGenParticles){
        index++;
        if (mom == tmpobj) break;
      }
      mom1_index.push_back(index);
    }
    else mom1_index.push_back(-1);
  }

  PUSH_VECTOR_WITH_NAME(colName, pt);
  PUSH_VECTOR_WITH_NAME(colName, eta);
  PUSH_VECTOR_WITH_NAME(colName, phi);
  PUSH_VECTOR_WITH_NAME(colName, mass);

  PUSH_VECTOR_WITH_NAME(colName, is_packed);

  PUSH_VECTOR_WITH_NAME(colName, isPromptFinalState); // (a)
  PUSH_VECTOR_WITH_NAME(colName, isPromptDecayed); // (b)
  PUSH_VECTOR_WITH_NAME(colName, isDirectPromptTauDecayProductFinalState); // (c)
  PUSH_VECTOR_WITH_NAME(colName, isHardProcess); // (d)
  PUSH_VECTOR_WITH_NAME(colName, fromHardProcessFinalState); // (e)
  PUSH_VECTOR_WITH_NAME(colName, fromHardProcessDecayed); // (f)
  PUSH_VECTOR_WITH_NAME(colName, isDirectHardProcessTauDecayProductFinalState); // (g)
  PUSH_VECTOR_WITH_NAME(colName, fromHardProcessBeforeFSR); // (h)
  PUSH_VECTOR_WITH_NAME(colName, isLastCopy); // (i)
  PUSH_VECTOR_WITH_NAME(colName, isLastCopyBeforeFSR); // (j)

  PUSH_VECTOR_WITH_NAME(colName, id);
  PUSH_VECTOR_WITH_NAME(colName, status);
  PUSH_VECTOR_WITH_NAME(colName, mom0_index);
  PUSH_VECTOR_WITH_NAME(colName, mom1_index);

}
void CMS3Ntuplizer::recordGenJets(edm::Event const& iEvent, bool const& isFatJet, std::vector<reco::GenJet const*>* filledObjects){
  std::string strColName = (isFatJet ? "genak4jets" : "genak8jets");
  const char* colName = strColName.data();
  edm::Handle< edm::View<reco::GenJet> > genJetsHandle;
  iEvent.getByToken((isFatJet ? genAK4JetsToken : genAK8JetsToken), genJetsHandle);
  if (!genJetsHandle.isValid()) throw cms::Exception((std::string("CMS3Ntuplizer::recordGenJets: Error getting the gen. ") + (isFatJet ? "ak4" : "ak8") + " jets from the event...").data());

  size_t n_objects = genJetsHandle->size();
  if (filledObjects) filledObjects->reserve(n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, pt, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, eta, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, phi, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, mass, n_objects);

  for (edm::View<reco::GenJet>::const_iterator obj = genJetsHandle->begin(); obj != genJetsHandle->end(); obj++){
    if (std::abs(obj->eta())>5. || obj->pt()<(isFatJet ? 100. : 20.)) continue;

    if (filledObjects) filledObjects->push_back(&(*obj));

    pt.push_back(obj->pt());
    eta.push_back(obj->eta());
    phi.push_back(obj->phi());
    mass.push_back(obj->mass());
  }

  // Pass collections to the communicator
  PUSH_VECTOR_WITH_NAME(colName, pt);
  PUSH_VECTOR_WITH_NAME(colName, eta);
  PUSH_VECTOR_WITH_NAME(colName, phi);
  PUSH_VECTOR_WITH_NAME(colName, mass);
}

size_t CMS3Ntuplizer::fillMuons(edm::Event const& iEvent, std::vector<pat::Muon const*>* filledObjects){
  const char colName[] = "muons";
  edm::Handle< edm::View<pat::Muon> > muonsHandle;
  iEvent.getByToken(muonsToken, muonsHandle);
  if (!muonsHandle.isValid()) throw cms::Exception("CMS3Ntuplizer::fillMuons: Error getting the muon collection from the event...");
  size_t n_objects = muonsHandle->size();

  if (filledObjects) filledObjects->reserve(n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, pt, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, eta, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, phi, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, mass, n_objects);

  MAKE_VECTOR_WITH_RESERVE(int, charge, n_objects);

  MAKE_VECTOR_WITH_RESERVE(unsigned int, POG_selector_bits, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, pfIso03_sum_charged_nofsr, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, pfIso03_sum_neutral_nofsr, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, pfIso03_comb_nofsr, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, pfIso04_sum_charged_nofsr, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, pfIso04_sum_neutral_nofsr, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, pfIso04_comb_nofsr, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, miniIso_sum_charged_nofsr, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, miniIso_sum_neutral_nofsr, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, miniIso_comb_nofsr, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, miniIso_comb_nofsr_uncorrected, n_objects);

  MAKE_VECTOR_WITH_RESERVE(int, time_comb_ndof, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, time_comb_IPInOut, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, time_comb_IPOutIn, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, time_comb_IPInOutError, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, time_comb_IPOutInError, n_objects);
  MAKE_VECTOR_WITH_RESERVE(int, time_rpc_ndof, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, time_rpc_IPInOut, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, time_rpc_IPOutIn, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, time_rpc_IPInOutError, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, time_rpc_IPOutInError, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, pull_dxdz_noArb_DT, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, pull_dxdz_noArb_CSC, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, scale_smear_pt_corr, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, scale_smear_pt_corr_scale_totalUp, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, scale_smear_pt_corr_scale_totalDn, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, scale_smear_pt_corr_smear_totalUp, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, scale_smear_pt_corr_smear_totalDn, n_objects);

  size_t n_skimmed_objects=0;
  for (edm::View<pat::Muon>::const_iterator obj = muonsHandle->begin(); obj != muonsHandle->end(); obj++){
    if (!MuonSelectionHelpers::testSkimMuon(*obj, this->year)) continue;

    // Core particle quantities
    pt.push_back(obj->pt());
    eta.push_back(obj->eta());
    phi.push_back(obj->phi());
    mass.push_back(obj->mass());

    PUSH_USERINT_INTO_VECTOR(charge);

    PUSH_USERINT_INTO_VECTOR(POG_selector_bits);

    PUSH_USERFLOAT_INTO_VECTOR(pfIso03_sum_charged_nofsr);
    PUSH_USERFLOAT_INTO_VECTOR(pfIso03_sum_neutral_nofsr);
    PUSH_USERFLOAT_INTO_VECTOR(pfIso03_comb_nofsr);

    PUSH_USERFLOAT_INTO_VECTOR(pfIso04_sum_charged_nofsr);
    PUSH_USERFLOAT_INTO_VECTOR(pfIso04_sum_neutral_nofsr);
    PUSH_USERFLOAT_INTO_VECTOR(pfIso04_comb_nofsr);

    PUSH_USERFLOAT_INTO_VECTOR(miniIso_sum_charged_nofsr);
    PUSH_USERFLOAT_INTO_VECTOR(miniIso_sum_neutral_nofsr);
    PUSH_USERFLOAT_INTO_VECTOR(miniIso_comb_nofsr);
    PUSH_USERFLOAT_INTO_VECTOR(miniIso_comb_nofsr_uncorrected);

    PUSH_USERINT_INTO_VECTOR(time_comb_ndof);
    PUSH_USERFLOAT_INTO_VECTOR(time_comb_IPInOut);
    PUSH_USERFLOAT_INTO_VECTOR(time_comb_IPOutIn);
    PUSH_USERFLOAT_INTO_VECTOR(time_comb_IPInOutError);
    PUSH_USERFLOAT_INTO_VECTOR(time_comb_IPOutInError);
    PUSH_USERINT_INTO_VECTOR(time_rpc_ndof);
    PUSH_USERFLOAT_INTO_VECTOR(time_rpc_IPInOut);
    PUSH_USERFLOAT_INTO_VECTOR(time_rpc_IPOutIn);
    PUSH_USERFLOAT_INTO_VECTOR(time_rpc_IPInOutError);
    PUSH_USERFLOAT_INTO_VECTOR(time_rpc_IPOutInError);

    PUSH_USERFLOAT_INTO_VECTOR(pull_dxdz_noArb_DT);
    PUSH_USERFLOAT_INTO_VECTOR(pull_dxdz_noArb_CSC);

    PUSH_USERFLOAT_INTO_VECTOR(scale_smear_pt_corr);
    PUSH_USERFLOAT_INTO_VECTOR(scale_smear_pt_corr_scale_totalUp);
    PUSH_USERFLOAT_INTO_VECTOR(scale_smear_pt_corr_scale_totalDn);
    PUSH_USERFLOAT_INTO_VECTOR(scale_smear_pt_corr_smear_totalUp);
    PUSH_USERFLOAT_INTO_VECTOR(scale_smear_pt_corr_smear_totalDn);

    if (filledObjects) filledObjects->push_back(&(*obj));
    n_skimmed_objects++;
  }

  // Pass collections to the communicator
  PUSH_VECTOR_WITH_NAME(colName, pt);
  PUSH_VECTOR_WITH_NAME(colName, eta);
  PUSH_VECTOR_WITH_NAME(colName, phi);
  PUSH_VECTOR_WITH_NAME(colName, mass);

  PUSH_VECTOR_WITH_NAME(colName, charge);

  PUSH_VECTOR_WITH_NAME(colName, POG_selector_bits);

  PUSH_VECTOR_WITH_NAME(colName, pfIso03_sum_charged_nofsr);
  PUSH_VECTOR_WITH_NAME(colName, pfIso03_sum_neutral_nofsr);
  PUSH_VECTOR_WITH_NAME(colName, pfIso03_comb_nofsr);

  PUSH_VECTOR_WITH_NAME(colName, pfIso04_sum_charged_nofsr);
  PUSH_VECTOR_WITH_NAME(colName, pfIso04_sum_neutral_nofsr);
  PUSH_VECTOR_WITH_NAME(colName, pfIso04_comb_nofsr);

  PUSH_VECTOR_WITH_NAME(colName, miniIso_sum_charged_nofsr);
  PUSH_VECTOR_WITH_NAME(colName, miniIso_sum_neutral_nofsr);
  PUSH_VECTOR_WITH_NAME(colName, miniIso_comb_nofsr);
  PUSH_VECTOR_WITH_NAME(colName, miniIso_comb_nofsr_uncorrected);

  PUSH_VECTOR_WITH_NAME(colName, time_comb_ndof);
  PUSH_VECTOR_WITH_NAME(colName, time_comb_IPInOut);
  PUSH_VECTOR_WITH_NAME(colName, time_comb_IPOutIn);
  PUSH_VECTOR_WITH_NAME(colName, time_comb_IPInOutError);
  PUSH_VECTOR_WITH_NAME(colName, time_comb_IPOutInError);
  PUSH_VECTOR_WITH_NAME(colName, time_rpc_ndof);
  PUSH_VECTOR_WITH_NAME(colName, time_rpc_IPInOut);
  PUSH_VECTOR_WITH_NAME(colName, time_rpc_IPOutIn);
  PUSH_VECTOR_WITH_NAME(colName, time_rpc_IPInOutError);
  PUSH_VECTOR_WITH_NAME(colName, time_rpc_IPOutInError);

  PUSH_VECTOR_WITH_NAME(colName, pull_dxdz_noArb_DT);
  PUSH_VECTOR_WITH_NAME(colName, pull_dxdz_noArb_CSC);

  PUSH_VECTOR_WITH_NAME(colName, scale_smear_pt_corr);
  PUSH_VECTOR_WITH_NAME(colName, scale_smear_pt_corr_scale_totalUp);
  PUSH_VECTOR_WITH_NAME(colName, scale_smear_pt_corr_scale_totalDn);
  PUSH_VECTOR_WITH_NAME(colName, scale_smear_pt_corr_smear_totalUp);
  PUSH_VECTOR_WITH_NAME(colName, scale_smear_pt_corr_smear_totalDn);

  return n_skimmed_objects;
}
size_t CMS3Ntuplizer::fillElectrons(edm::Event const& iEvent, std::vector<pat::Electron const*>* filledObjects){
  const char colName[] = "electrons";
  edm::Handle< edm::View<pat::Electron> > electronsHandle;
  iEvent.getByToken(electronsToken, electronsHandle);
  if (!electronsHandle.isValid()) throw cms::Exception("CMS3Ntuplizer::fillElectrons: Error getting the electron collection from the event...");
  size_t n_objects = electronsHandle->size();

  if (filledObjects) filledObjects->reserve(n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, pt, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, eta, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, phi, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, mass, n_objects);

  MAKE_VECTOR_WITH_RESERVE(int, charge, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, etaSC, n_objects);

  // Has no convention correspondence in nanoAOD
  MAKE_VECTOR_WITH_RESERVE(float, scale_smear_corr, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, scale_smear_corr_scale_totalUp, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, scale_smear_corr_scale_totalDn, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, scale_smear_corr_smear_totalUp, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, scale_smear_corr_smear_totalDn, n_objects);

  MAKE_VECTOR_WITH_RESERVE(bool, conv_vtx_flag, n_objects);
  MAKE_VECTOR_WITH_RESERVE(int, n_missing_inner_hits, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, id_MVA_Fall17V2_Iso_Val, n_objects);
  MAKE_VECTOR_WITH_RESERVE(unsigned int, id_MVA_Fall17V2_Iso_Cat, n_objects);
  MAKE_VECTOR_WITH_RESERVE(bool, id_MVA_Fall17V2_Iso_pass_wpLoose, n_objects);
  MAKE_VECTOR_WITH_RESERVE(bool, id_MVA_Fall17V2_Iso_pass_wp90, n_objects);
  MAKE_VECTOR_WITH_RESERVE(bool, id_MVA_Fall17V2_Iso_pass_wp80, n_objects);
  MAKE_VECTOR_WITH_RESERVE(bool, id_MVA_Fall17V2_Iso_pass_wpHZZ, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, id_MVA_Fall17V2_NoIso_Val, n_objects);
  MAKE_VECTOR_WITH_RESERVE(unsigned int, id_MVA_Fall17V2_NoIso_Cat, n_objects);
  MAKE_VECTOR_WITH_RESERVE(bool, id_MVA_Fall17V2_NoIso_pass_wpLoose, n_objects);
  MAKE_VECTOR_WITH_RESERVE(bool, id_MVA_Fall17V2_NoIso_pass_wp90, n_objects);
  MAKE_VECTOR_WITH_RESERVE(bool, id_MVA_Fall17V2_NoIso_pass_wp80, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, id_MVA_HZZRun2Legacy_Iso_Val, n_objects);
  MAKE_VECTOR_WITH_RESERVE(unsigned int, id_MVA_HZZRun2Legacy_Iso_Cat, n_objects);
  MAKE_VECTOR_WITH_RESERVE(bool, id_MVA_HZZRun2Legacy_Iso_pass_wpHZZ, n_objects);

  MAKE_VECTOR_WITH_RESERVE(unsigned int, id_cutBased_Fall17V2_Veto_Bits, n_objects);
  MAKE_VECTOR_WITH_RESERVE(unsigned int, id_cutBased_Fall17V2_Loose_Bits, n_objects);
  MAKE_VECTOR_WITH_RESERVE(unsigned int, id_cutBased_Fall17V2_Medium_Bits, n_objects);
  MAKE_VECTOR_WITH_RESERVE(unsigned int, id_cutBased_Fall17V2_Tight_Bits, n_objects);

  MAKE_VECTOR_WITH_RESERVE(unsigned int, id_cutBased_Fall17V1_Veto_Bits, n_objects);
  MAKE_VECTOR_WITH_RESERVE(unsigned int, id_cutBased_Fall17V1_Loose_Bits, n_objects);
  MAKE_VECTOR_WITH_RESERVE(unsigned int, id_cutBased_Fall17V1_Medium_Bits, n_objects);
  MAKE_VECTOR_WITH_RESERVE(unsigned int, id_cutBased_Fall17V1_Tight_Bits, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, pfIso03_sum_charged_nofsr, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, pfIso03_sum_neutral_nofsr, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, pfIso03_comb_nofsr, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, pfIso04_sum_charged_nofsr, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, pfIso04_sum_neutral_nofsr, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, pfIso04_comb_nofsr, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, miniIso_sum_charged_nofsr, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, miniIso_sum_neutral_nofsr, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, miniIso_comb_nofsr, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, miniIso_comb_nofsr_uncorrected, n_objects);

  MAKE_VECTOR_WITH_RESERVE(unsigned int, n_associated_pfcands, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, associated_pfcands_sum_sc_pt, n_objects);

  MAKE_VECTOR_WITH_RESERVE(unsigned int, fid_mask, n_objects);
  MAKE_VECTOR_WITH_RESERVE(unsigned int, type_mask, n_objects);

  size_t n_skimmed_objects=0;
  for (edm::View<pat::Electron>::const_iterator obj = electronsHandle->begin(); obj != electronsHandle->end(); obj++){
    if (
      !ElectronSelectionHelpers::testSkimElectron(
        *obj, this->year,
        {
          "id_cutBased_Fall17V2_Veto_Bits", "id_cutBased_Fall17V2_Loose_Bits", "id_cutBased_Fall17V2_Medium_Bits", "id_cutBased_Fall17V2_Tight_Bits",
          "id_cutBased_Fall17V1_Veto_Bits", "id_cutBased_Fall17V1_Loose_Bits", "id_cutBased_Fall17V1_Medium_Bits", "id_cutBased_Fall17V1_Tight_Bits"
        },
        {
          "id_MVA_Fall17V2_Iso_pass_wpLoose", "id_MVA_Fall17V2_Iso_pass_wp90", "id_MVA_Fall17V2_Iso_pass_wp80", "id_MVA_Fall17V2_Iso_pass_wpHZZ",
          "id_MVA_Fall17V2_NoIso_pass_wpLoose", "id_MVA_Fall17V2_NoIso_pass_wp90", "id_MVA_Fall17V2_NoIso_pass_wp80",
          "id_MVA_HZZRun2Legacy_Iso_pass_wpHZZ"
        }
      )
      ) continue;

    // Core particle quantities
    // Uncorrected p4
    pt.push_back(obj->pt());
    eta.push_back(obj->eta());
    phi.push_back(obj->phi());
    mass.push_back(obj->mass());

    // Charge: Can obtain pdgId from this, so no need to record pdgId again
    PUSH_USERINT_INTO_VECTOR(charge);
    PUSH_USERFLOAT_INTO_VECTOR(etaSC);

    // Scale and smear
    // Nominal value: Needs to multiply the uncorrected p4 at analysis level
    PUSH_USERFLOAT_INTO_VECTOR(scale_smear_corr);
    // Uncertainties: Only store total up/dn for the moment
    PUSH_USERFLOAT_INTO_VECTOR(scale_smear_corr_scale_totalUp);
    PUSH_USERFLOAT_INTO_VECTOR(scale_smear_corr_scale_totalDn);
    PUSH_USERFLOAT_INTO_VECTOR(scale_smear_corr_smear_totalUp);
    PUSH_USERFLOAT_INTO_VECTOR(scale_smear_corr_smear_totalDn);

    PUSH_USERINT_INTO_VECTOR(conv_vtx_flag);
    PUSH_USERINT_INTO_VECTOR(n_missing_inner_hits);

    // Id variables
    // Fall17V2_Iso MVA id
    PUSH_USERFLOAT_INTO_VECTOR(id_MVA_Fall17V2_Iso_Val);
    PUSH_USERINT_INTO_VECTOR(id_MVA_Fall17V2_Iso_Cat);
    PUSH_USERINT_INTO_VECTOR(id_MVA_Fall17V2_Iso_pass_wpLoose);
    PUSH_USERINT_INTO_VECTOR(id_MVA_Fall17V2_Iso_pass_wp90);
    PUSH_USERINT_INTO_VECTOR(id_MVA_Fall17V2_Iso_pass_wp80);
    PUSH_USERINT_INTO_VECTOR(id_MVA_Fall17V2_Iso_pass_wpHZZ);

    // Fall17V2_NoIso MVA id
    PUSH_USERFLOAT_INTO_VECTOR(id_MVA_Fall17V2_NoIso_Val);
    PUSH_USERINT_INTO_VECTOR(id_MVA_Fall17V2_NoIso_Cat);
    PUSH_USERINT_INTO_VECTOR(id_MVA_Fall17V2_NoIso_pass_wpLoose);
    PUSH_USERINT_INTO_VECTOR(id_MVA_Fall17V2_NoIso_pass_wp90);
    PUSH_USERINT_INTO_VECTOR(id_MVA_Fall17V2_NoIso_pass_wp80);

    // HZZ Run 2 legacy electron MVA id
    PUSH_USERFLOAT_INTO_VECTOR(id_MVA_HZZRun2Legacy_Iso_Val);
    PUSH_USERINT_INTO_VECTOR(id_MVA_HZZRun2Legacy_Iso_Cat);
    PUSH_USERINT_INTO_VECTOR(id_MVA_HZZRun2Legacy_Iso_pass_wpHZZ);

    // Fall17V2 cut-based ids
    PUSH_USERINT_INTO_VECTOR(id_cutBased_Fall17V2_Veto_Bits);
    PUSH_USERINT_INTO_VECTOR(id_cutBased_Fall17V2_Loose_Bits);
    PUSH_USERINT_INTO_VECTOR(id_cutBased_Fall17V2_Medium_Bits);
    PUSH_USERINT_INTO_VECTOR(id_cutBased_Fall17V2_Tight_Bits);

    // Fall17V1 cut-based ids
    PUSH_USERINT_INTO_VECTOR(id_cutBased_Fall17V1_Veto_Bits);
    PUSH_USERINT_INTO_VECTOR(id_cutBased_Fall17V1_Loose_Bits);
    PUSH_USERINT_INTO_VECTOR(id_cutBased_Fall17V1_Medium_Bits);
    PUSH_USERINT_INTO_VECTOR(id_cutBased_Fall17V1_Tight_Bits);

    // Isolation variables
    PUSH_USERFLOAT_INTO_VECTOR(pfIso03_sum_charged_nofsr);
    PUSH_USERFLOAT_INTO_VECTOR(pfIso03_sum_neutral_nofsr);
    PUSH_USERFLOAT_INTO_VECTOR(pfIso03_comb_nofsr);

    PUSH_USERFLOAT_INTO_VECTOR(pfIso04_sum_charged_nofsr);
    PUSH_USERFLOAT_INTO_VECTOR(pfIso04_sum_neutral_nofsr);
    PUSH_USERFLOAT_INTO_VECTOR(pfIso04_comb_nofsr);

    PUSH_USERFLOAT_INTO_VECTOR(miniIso_sum_charged_nofsr);
    PUSH_USERFLOAT_INTO_VECTOR(miniIso_sum_neutral_nofsr);
    PUSH_USERFLOAT_INTO_VECTOR(miniIso_comb_nofsr);
    PUSH_USERFLOAT_INTO_VECTOR(miniIso_comb_nofsr_uncorrected);

    PUSH_USERINT_INTO_VECTOR(n_associated_pfcands);
    PUSH_USERFLOAT_INTO_VECTOR(associated_pfcands_sum_sc_pt);

    // Masks
    PUSH_USERINT_INTO_VECTOR(fid_mask);
    PUSH_USERINT_INTO_VECTOR(type_mask);

    if (filledObjects) filledObjects->push_back(&(*obj));
    n_skimmed_objects++;
  }

  // Pass collections to the communicator
  PUSH_VECTOR_WITH_NAME(colName, pt);
  PUSH_VECTOR_WITH_NAME(colName, eta);
  PUSH_VECTOR_WITH_NAME(colName, phi);
  PUSH_VECTOR_WITH_NAME(colName, mass);

  PUSH_VECTOR_WITH_NAME(colName, charge);
  PUSH_VECTOR_WITH_NAME(colName, etaSC);

  PUSH_VECTOR_WITH_NAME(colName, conv_vtx_flag);
  PUSH_VECTOR_WITH_NAME(colName, n_missing_inner_hits);

  // Has no convention correspondence in nanoAOD
  PUSH_VECTOR_WITH_NAME(colName, scale_smear_corr);
  PUSH_VECTOR_WITH_NAME(colName, scale_smear_corr_scale_totalUp);
  PUSH_VECTOR_WITH_NAME(colName, scale_smear_corr_scale_totalDn);
  PUSH_VECTOR_WITH_NAME(colName, scale_smear_corr_smear_totalUp);
  PUSH_VECTOR_WITH_NAME(colName, scale_smear_corr_smear_totalDn);

  PUSH_VECTOR_WITH_NAME(colName, id_MVA_Fall17V2_Iso_Val);
  PUSH_VECTOR_WITH_NAME(colName, id_MVA_Fall17V2_Iso_Cat);
  PUSH_VECTOR_WITH_NAME(colName, id_MVA_Fall17V2_Iso_pass_wpLoose);
  PUSH_VECTOR_WITH_NAME(colName, id_MVA_Fall17V2_Iso_pass_wp90);
  PUSH_VECTOR_WITH_NAME(colName, id_MVA_Fall17V2_Iso_pass_wp80);
  PUSH_VECTOR_WITH_NAME(colName, id_MVA_Fall17V2_Iso_pass_wpHZZ);

  PUSH_VECTOR_WITH_NAME(colName, id_MVA_Fall17V2_NoIso_Val);
  PUSH_VECTOR_WITH_NAME(colName, id_MVA_Fall17V2_NoIso_Cat);
  PUSH_VECTOR_WITH_NAME(colName, id_MVA_Fall17V2_NoIso_pass_wpLoose);
  PUSH_VECTOR_WITH_NAME(colName, id_MVA_Fall17V2_NoIso_pass_wp90);
  PUSH_VECTOR_WITH_NAME(colName, id_MVA_Fall17V2_NoIso_pass_wp80);

  PUSH_VECTOR_WITH_NAME(colName, id_MVA_HZZRun2Legacy_Iso_Val);
  PUSH_VECTOR_WITH_NAME(colName, id_MVA_HZZRun2Legacy_Iso_Cat);
  PUSH_VECTOR_WITH_NAME(colName, id_MVA_HZZRun2Legacy_Iso_pass_wpHZZ);

  PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V2_Veto_Bits);
  PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V2_Loose_Bits);
  PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V2_Medium_Bits);
  PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V2_Tight_Bits);

  PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V1_Veto_Bits);
  PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V1_Loose_Bits);
  PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V1_Medium_Bits);
  PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V1_Tight_Bits);

  PUSH_VECTOR_WITH_NAME(colName, pfIso03_sum_charged_nofsr);
  PUSH_VECTOR_WITH_NAME(colName, pfIso03_sum_neutral_nofsr);
  PUSH_VECTOR_WITH_NAME(colName, pfIso03_comb_nofsr);

  PUSH_VECTOR_WITH_NAME(colName, pfIso04_sum_charged_nofsr);
  PUSH_VECTOR_WITH_NAME(colName, pfIso04_sum_neutral_nofsr);
  PUSH_VECTOR_WITH_NAME(colName, pfIso04_comb_nofsr);

  PUSH_VECTOR_WITH_NAME(colName, miniIso_sum_charged_nofsr);
  PUSH_VECTOR_WITH_NAME(colName, miniIso_sum_neutral_nofsr);
  PUSH_VECTOR_WITH_NAME(colName, miniIso_comb_nofsr);
  PUSH_VECTOR_WITH_NAME(colName, miniIso_comb_nofsr_uncorrected);

  PUSH_VECTOR_WITH_NAME(colName, n_associated_pfcands);
  PUSH_VECTOR_WITH_NAME(colName, associated_pfcands_sum_sc_pt);

  PUSH_VECTOR_WITH_NAME(colName, fid_mask);
  PUSH_VECTOR_WITH_NAME(colName, type_mask);

  return n_skimmed_objects;
}
size_t CMS3Ntuplizer::fillPhotons(edm::Event const& iEvent, std::vector<pat::Photon const*>* filledObjects){
  const char colName[] = "photons";
  edm::Handle< edm::View<pat::Photon> > photonsHandle;
  iEvent.getByToken(photonsToken, photonsHandle);
  if (!photonsHandle.isValid()) throw cms::Exception("CMS3Ntuplizer::fillPhotons: Error getting the photon collection from the event...");
  size_t n_objects = photonsHandle->size();

  if (filledObjects) filledObjects->reserve(n_objects);

  // Begin filling the objects
  MAKE_VECTOR_WITH_RESERVE(float, pt, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, eta, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, phi, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, mass, n_objects);

  // Has no convention correspondence in nanoAOD
  MAKE_VECTOR_WITH_RESERVE(float, scale_smear_corr, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, scale_smear_corr_scale_totalUp, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, scale_smear_corr_scale_totalDn, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, scale_smear_corr_smear_totalUp, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, scale_smear_corr_smear_totalDn, n_objects);

  // These two are needed for the MVA id
  MAKE_VECTOR_WITH_RESERVE(bool, hasPixelSeed, n_objects);
  MAKE_VECTOR_WITH_RESERVE(bool, passElectronVeto, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, id_MVA_Fall17V2_Val, n_objects);
  MAKE_VECTOR_WITH_RESERVE(unsigned int, id_MVA_Fall17V2_Cat, n_objects);
  MAKE_VECTOR_WITH_RESERVE(bool, id_MVA_Fall17V2_pass_wp90, n_objects);
  MAKE_VECTOR_WITH_RESERVE(bool, id_MVA_Fall17V2_pass_wp80, n_objects);

  MAKE_VECTOR_WITH_RESERVE(unsigned int, id_cutBased_Fall17V2_Loose_Bits, n_objects);
  MAKE_VECTOR_WITH_RESERVE(unsigned int, id_cutBased_Fall17V2_Medium_Bits, n_objects);
  MAKE_VECTOR_WITH_RESERVE(unsigned int, id_cutBased_Fall17V2_Tight_Bits, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, pfIso_comb, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, pfChargedHadronIso_EAcorr, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, pfNeutralHadronIso_EAcorr, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, pfEMIso_EAcorr, n_objects);

  MAKE_VECTOR_WITH_RESERVE(unsigned int, n_associated_pfcands, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, associated_pfcands_sum_sc_pt, n_objects);

  //MAKE_VECTOR_WITH_RESERVE(bool, pass_fsr_preselection, n_objects);
  //MAKE_VECTOR_WITH_RESERVE(float, fsrIso, n_objects);

  size_t n_skimmed_objects=0;
  for (edm::View<pat::Photon>::const_iterator obj = photonsHandle->begin(); obj != photonsHandle->end(); obj++){
    //bool passStandardSkim = PhotonSelectionHelpers::testSkimPhoton(*obj, this->year);
    //bool passFSRSkim = (HelperFunctions::checkListVariable(allFSRCandidates, &(*obj)) && PhotonSelectionHelpers::testSkimFSRPhoton(*obj, fsr_mindr_map[&(*obj)], this->year));
    //if (!passStandardSkim && !passFSRSkim) continue;
    if (
      !PhotonSelectionHelpers::testSkimPhoton(
        *obj,
        this->year,
        { "id_cutBased_Fall17V2_Loose_Bits", "id_cutBased_Fall17V2_Medium_Bits", "id_cutBased_Fall17V2_Tight_Bits" },
        { "id_MVA_Fall17V2_pass_wp90", "id_MVA_Fall17V2_pass_wp80" }
      )
      ) continue;

    // Core particle quantities
    // Uncorrected p4
    pt.push_back(obj->pt());
    eta.push_back(obj->eta());
    phi.push_back(obj->phi());
    mass.push_back(obj->mass());

    // Flag to identify FSR-preselected candidates
    //pass_fsr_preselection.push_back(passFSRSkim);

    // Scale and smear
    // Nominal value: Needs to multiply the uncorrected p4 at analysis level
    PUSH_USERFLOAT_INTO_VECTOR(scale_smear_corr);
    // Uncertainties: Only store total up/dn for the moment
    PUSH_USERFLOAT_INTO_VECTOR(scale_smear_corr_scale_totalUp);
    PUSH_USERFLOAT_INTO_VECTOR(scale_smear_corr_scale_totalDn);
    PUSH_USERFLOAT_INTO_VECTOR(scale_smear_corr_smear_totalUp);
    PUSH_USERFLOAT_INTO_VECTOR(scale_smear_corr_smear_totalDn);

    // Id variables
    PUSH_USERINT_INTO_VECTOR(hasPixelSeed);
    PUSH_USERINT_INTO_VECTOR(passElectronVeto);

    // Fall17V2 MVA id
    PUSH_USERFLOAT_INTO_VECTOR(id_MVA_Fall17V2_Val);
    PUSH_USERINT_INTO_VECTOR(id_MVA_Fall17V2_Cat);
    PUSH_USERINT_INTO_VECTOR(id_MVA_Fall17V2_pass_wp90);
    PUSH_USERINT_INTO_VECTOR(id_MVA_Fall17V2_pass_wp80);

    // Fall17V2 cut-based ids
    PUSH_USERINT_INTO_VECTOR(id_cutBased_Fall17V2_Loose_Bits);
    PUSH_USERINT_INTO_VECTOR(id_cutBased_Fall17V2_Medium_Bits);
    PUSH_USERINT_INTO_VECTOR(id_cutBased_Fall17V2_Tight_Bits);

    PUSH_USERFLOAT_INTO_VECTOR(pfIso_comb);
    PUSH_USERFLOAT_INTO_VECTOR(pfChargedHadronIso_EAcorr);
    PUSH_USERFLOAT_INTO_VECTOR(pfNeutralHadronIso_EAcorr);
    PUSH_USERFLOAT_INTO_VECTOR(pfEMIso_EAcorr);

    PUSH_USERINT_INTO_VECTOR(n_associated_pfcands);
    PUSH_USERFLOAT_INTO_VECTOR(associated_pfcands_sum_sc_pt);

    //PUSH_USERFLOAT_INTO_VECTOR(fsrIso);

    if (filledObjects) filledObjects->push_back(&(*obj));
    n_skimmed_objects++;
  }

  // Pass collections to the communicator
  PUSH_VECTOR_WITH_NAME(colName, pt);
  PUSH_VECTOR_WITH_NAME(colName, eta);
  PUSH_VECTOR_WITH_NAME(colName, phi);
  PUSH_VECTOR_WITH_NAME(colName, mass);

  // Has no convention correspondence in nanoAOD
  PUSH_VECTOR_WITH_NAME(colName, scale_smear_corr);
  PUSH_VECTOR_WITH_NAME(colName, scale_smear_corr_scale_totalUp);
  PUSH_VECTOR_WITH_NAME(colName, scale_smear_corr_scale_totalDn);
  PUSH_VECTOR_WITH_NAME(colName, scale_smear_corr_smear_totalUp);
  PUSH_VECTOR_WITH_NAME(colName, scale_smear_corr_smear_totalDn);

  PUSH_VECTOR_WITH_NAME(colName, hasPixelSeed);
  PUSH_VECTOR_WITH_NAME(colName, passElectronVeto);

  PUSH_VECTOR_WITH_NAME(colName, id_MVA_Fall17V2_Val);
  PUSH_VECTOR_WITH_NAME(colName, id_MVA_Fall17V2_Cat);
  PUSH_VECTOR_WITH_NAME(colName, id_MVA_Fall17V2_pass_wp90);
  PUSH_VECTOR_WITH_NAME(colName, id_MVA_Fall17V2_pass_wp80);

  PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V2_Loose_Bits);
  PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V2_Medium_Bits);
  PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V2_Tight_Bits);

  PUSH_VECTOR_WITH_NAME(colName, pfIso_comb);
  PUSH_VECTOR_WITH_NAME(colName, pfChargedHadronIso_EAcorr);
  PUSH_VECTOR_WITH_NAME(colName, pfNeutralHadronIso_EAcorr);
  PUSH_VECTOR_WITH_NAME(colName, pfEMIso_EAcorr);

  PUSH_VECTOR_WITH_NAME(colName, n_associated_pfcands);
  PUSH_VECTOR_WITH_NAME(colName, associated_pfcands_sum_sc_pt);

  //PUSH_VECTOR_WITH_NAME(colName, pass_fsr_preselection);
  //PUSH_VECTOR_WITH_NAME(colName, fsrIso);

  return n_skimmed_objects;
}
size_t CMS3Ntuplizer::fillAK4Jets(edm::Event const& iEvent, std::vector<pat::Jet const*>* filledObjects){
  constexpr AK4JetSelectionHelpers::AK4JetType jetType = AK4JetSelectionHelpers::AK4PFCHS;

  const char colName[] = "ak4jets";
  edm::Handle< edm::View<pat::Jet> > ak4jetsHandle;
  iEvent.getByToken(ak4jetsToken, ak4jetsHandle);
  if (!ak4jetsHandle.isValid()) throw cms::Exception("CMS3Ntuplizer::fillAK4Jets: Error getting the ak4 jet collection from the event...");
  size_t n_objects = ak4jetsHandle->size();

  if (filledObjects) filledObjects->reserve(n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, pt, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, eta, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, phi, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, mass, n_objects);

  MAKE_VECTOR_WITH_RESERVE(bool, pass_looseId, n_objects);
  MAKE_VECTOR_WITH_RESERVE(bool, pass_tightId, n_objects);
  MAKE_VECTOR_WITH_RESERVE(bool, pass_leptonVetoId, n_objects);
  MAKE_VECTOR_WITH_RESERVE(bool, pass_puId, n_objects);

  MAKE_VECTOR_WITH_RESERVE(size_t, n_pfcands, n_objects);
  MAKE_VECTOR_WITH_RESERVE(size_t, n_mucands, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, area, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, pt_resolution, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, deepFlavourprobb, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, deepFlavourprobbb, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, deepFlavourprobc, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, deepFlavourprobg, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, deepFlavourproblepb, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, deepFlavourprobuds, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, deepCSVprobb, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, deepCSVprobbb, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, deepCSVprobc, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, deepCSVprobudsg, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, ptDistribution, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, totalMultiplicity, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, axis1, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, axis2, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, JECNominal, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, JECUp, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, JECDn, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, JERNominal, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, JERUp, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, JERDn, n_objects);

  MAKE_VECTOR_WITH_RESERVE(int, partonFlavour, n_objects);
  MAKE_VECTOR_WITH_RESERVE(int, hadronFlavour, n_objects);

  size_t n_skimmed_objects=0;
  for (edm::View<pat::Jet>::const_iterator obj = ak4jetsHandle->begin(); obj != ak4jetsHandle->end(); obj++){
    if (!AK4JetSelectionHelpers::testSkimAK4Jet(*obj, this->year, jetType)) continue;

    // Core particle quantities
    // These are the uncorrected momentum components!
    pt.push_back(obj->pt());
    eta.push_back(obj->eta());
    phi.push_back(obj->phi());
    mass.push_back(obj->mass());

    pass_looseId.push_back(AK4JetSelectionHelpers::testLooseAK4Jet(*obj, this->year, jetType));
    pass_tightId.push_back(AK4JetSelectionHelpers::testTightAK4Jet(*obj, this->year, jetType));
    pass_leptonVetoId.push_back(AK4JetSelectionHelpers::testLeptonVetoAK4Jet(*obj, this->year, jetType));
    pass_puId.push_back(AK4JetSelectionHelpers::testPileUpAK4Jet(*obj, this->year, jetType));

    PUSH_USERINT_INTO_VECTOR(n_pfcands);
    PUSH_USERINT_INTO_VECTOR(n_mucands);

    PUSH_USERFLOAT_INTO_VECTOR(area);
    PUSH_USERFLOAT_INTO_VECTOR(pt_resolution);

    PUSH_USERFLOAT_INTO_VECTOR(deepFlavourprobb);
    PUSH_USERFLOAT_INTO_VECTOR(deepFlavourprobbb);
    PUSH_USERFLOAT_INTO_VECTOR(deepFlavourprobc);
    PUSH_USERFLOAT_INTO_VECTOR(deepFlavourprobg);
    PUSH_USERFLOAT_INTO_VECTOR(deepFlavourproblepb);
    PUSH_USERFLOAT_INTO_VECTOR(deepFlavourprobuds);

    PUSH_USERFLOAT_INTO_VECTOR(deepCSVprobb);
    PUSH_USERFLOAT_INTO_VECTOR(deepCSVprobbb);
    PUSH_USERFLOAT_INTO_VECTOR(deepCSVprobc);
    PUSH_USERFLOAT_INTO_VECTOR(deepCSVprobudsg);

    PUSH_USERFLOAT_INTO_VECTOR(ptDistribution);
    PUSH_USERFLOAT_INTO_VECTOR(totalMultiplicity);
    PUSH_USERFLOAT_INTO_VECTOR(axis1);
    PUSH_USERFLOAT_INTO_VECTOR(axis2);

    PUSH_USERFLOAT_INTO_VECTOR(JECNominal);
    PUSH_USERFLOAT_INTO_VECTOR(JECUp);
    PUSH_USERFLOAT_INTO_VECTOR(JECDn);

    PUSH_USERFLOAT_INTO_VECTOR(JERNominal);
    PUSH_USERFLOAT_INTO_VECTOR(JERUp);
    PUSH_USERFLOAT_INTO_VECTOR(JERDn);

    PUSH_USERINT_INTO_VECTOR(partonFlavour);
    PUSH_USERINT_INTO_VECTOR(hadronFlavour);

    if (filledObjects) filledObjects->push_back(&(*obj));
    n_skimmed_objects++;
  }

  // Pass collections to the communicator
  PUSH_VECTOR_WITH_NAME(colName, pt);
  PUSH_VECTOR_WITH_NAME(colName, eta);
  PUSH_VECTOR_WITH_NAME(colName, phi);
  PUSH_VECTOR_WITH_NAME(colName, mass);

  PUSH_VECTOR_WITH_NAME(colName, pass_looseId);
  PUSH_VECTOR_WITH_NAME(colName, pass_tightId);
  PUSH_VECTOR_WITH_NAME(colName, pass_leptonVetoId);
  PUSH_VECTOR_WITH_NAME(colName, pass_puId);

  PUSH_VECTOR_WITH_NAME(colName, n_pfcands);
  PUSH_VECTOR_WITH_NAME(colName, n_mucands);

  PUSH_VECTOR_WITH_NAME(colName, area);
  PUSH_VECTOR_WITH_NAME(colName, pt_resolution);

  PUSH_VECTOR_WITH_NAME(colName, deepFlavourprobb);
  PUSH_VECTOR_WITH_NAME(colName, deepFlavourprobbb);
  PUSH_VECTOR_WITH_NAME(colName, deepFlavourprobc);
  PUSH_VECTOR_WITH_NAME(colName, deepFlavourprobg);
  PUSH_VECTOR_WITH_NAME(colName, deepFlavourproblepb);
  PUSH_VECTOR_WITH_NAME(colName, deepFlavourprobuds);

  PUSH_VECTOR_WITH_NAME(colName, deepCSVprobb);
  PUSH_VECTOR_WITH_NAME(colName, deepCSVprobbb);
  PUSH_VECTOR_WITH_NAME(colName, deepCSVprobc);
  PUSH_VECTOR_WITH_NAME(colName, deepCSVprobudsg);

  PUSH_VECTOR_WITH_NAME(colName, ptDistribution);
  PUSH_VECTOR_WITH_NAME(colName, totalMultiplicity);
  PUSH_VECTOR_WITH_NAME(colName, axis1);
  PUSH_VECTOR_WITH_NAME(colName, axis2);

  PUSH_VECTOR_WITH_NAME(colName, JECNominal);
  PUSH_VECTOR_WITH_NAME(colName, JECUp);
  PUSH_VECTOR_WITH_NAME(colName, JECDn);

  PUSH_VECTOR_WITH_NAME(colName, JERNominal);
  PUSH_VECTOR_WITH_NAME(colName, JERUp);
  PUSH_VECTOR_WITH_NAME(colName, JERDn);

  PUSH_VECTOR_WITH_NAME(colName, partonFlavour);
  PUSH_VECTOR_WITH_NAME(colName, hadronFlavour);

  return n_skimmed_objects;
}
size_t CMS3Ntuplizer::fillAK8Jets(edm::Event const& iEvent, std::vector<pat::Jet const*>* filledObjects){
  const char colName[] = "ak8jets";
  edm::Handle< edm::View<pat::Jet> > ak8jetsHandle;
  iEvent.getByToken(ak8jetsToken, ak8jetsHandle);
  if (!ak8jetsHandle.isValid()) throw cms::Exception("CMS3Ntuplizer::fillAK8Jets: Error getting the ak8 jet collection from the event...");
  size_t n_objects = ak8jetsHandle->size();

  if (filledObjects) filledObjects->reserve(n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, pt, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, eta, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, phi, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, mass, n_objects);

  //MAKE_VECTOR_WITH_RESERVE(bool, pass_looseId, n_objects);
  //MAKE_VECTOR_WITH_RESERVE(bool, pass_tightId, n_objects);

  MAKE_VECTOR_WITH_RESERVE(size_t, n_pfcands, n_objects);
  MAKE_VECTOR_WITH_RESERVE(size_t, n_mucands, n_objects);
  MAKE_VECTOR_WITH_RESERVE(size_t, n_softdrop_subjets, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, area, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, pt_resolution, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, softdrop_pt, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, softdrop_eta, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, softdrop_phi, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, softdrop_mass, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, softdrop_subjet0_pt, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, softdrop_subjet0_eta, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, softdrop_subjet0_phi, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, softdrop_subjet0_mass, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, softdrop_subjet1_pt, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, softdrop_subjet1_eta, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, softdrop_subjet1_phi, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, softdrop_subjet1_mass, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, ptDistribution, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, totalMultiplicity, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, axis1, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, axis2, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, JECNominal, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, JECUp, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, JECDn, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, JERNominal, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, JERUp, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, JERDn, n_objects);

  MAKE_VECTOR_WITH_RESERVE(int, partonFlavour, n_objects);
  MAKE_VECTOR_WITH_RESERVE(int, hadronFlavour, n_objects);

  size_t n_skimmed_objects=0;
  for (edm::View<pat::Jet>::const_iterator obj = ak8jetsHandle->begin(); obj != ak8jetsHandle->end(); obj++){
    if (!AK8JetSelectionHelpers::testSkimAK8Jet(*obj, this->year)) continue;

    // Core particle quantities
    // These are the uncorrected momentum components!
    pt.push_back(obj->pt());
    eta.push_back(obj->eta());
    phi.push_back(obj->phi());
    mass.push_back(obj->mass());

    //pass_looseId.push_back(AK8JetSelectionHelpers::testLooseAK8Jet(*obj, this->year));
    //pass_tightId.push_back(AK8JetSelectionHelpers::testTightAK8Jet(*obj, this->year));

    PUSH_USERINT_INTO_VECTOR(n_pfcands);
    PUSH_USERINT_INTO_VECTOR(n_mucands);
    PUSH_USERINT_INTO_VECTOR(n_softdrop_subjets);

    PUSH_USERFLOAT_INTO_VECTOR(area);
    PUSH_USERFLOAT_INTO_VECTOR(pt_resolution);

    PUSH_USERFLOAT_INTO_VECTOR(softdrop_pt);
    PUSH_USERFLOAT_INTO_VECTOR(softdrop_eta);
    PUSH_USERFLOAT_INTO_VECTOR(softdrop_phi);
    PUSH_USERFLOAT_INTO_VECTOR(softdrop_mass);

    PUSH_USERFLOAT_INTO_VECTOR(softdrop_subjet0_pt);
    PUSH_USERFLOAT_INTO_VECTOR(softdrop_subjet0_eta);
    PUSH_USERFLOAT_INTO_VECTOR(softdrop_subjet0_phi);
    PUSH_USERFLOAT_INTO_VECTOR(softdrop_subjet0_mass);

    PUSH_USERFLOAT_INTO_VECTOR(softdrop_subjet1_pt);
    PUSH_USERFLOAT_INTO_VECTOR(softdrop_subjet1_eta);
    PUSH_USERFLOAT_INTO_VECTOR(softdrop_subjet1_phi);
    PUSH_USERFLOAT_INTO_VECTOR(softdrop_subjet1_mass);

    PUSH_USERFLOAT_INTO_VECTOR(ptDistribution);
    PUSH_USERFLOAT_INTO_VECTOR(totalMultiplicity);
    PUSH_USERFLOAT_INTO_VECTOR(axis1);
    PUSH_USERFLOAT_INTO_VECTOR(axis2);

    PUSH_USERFLOAT_INTO_VECTOR(JECNominal);
    PUSH_USERFLOAT_INTO_VECTOR(JECUp);
    PUSH_USERFLOAT_INTO_VECTOR(JECDn);

    PUSH_USERFLOAT_INTO_VECTOR(JERNominal);
    PUSH_USERFLOAT_INTO_VECTOR(JERUp);
    PUSH_USERFLOAT_INTO_VECTOR(JERDn);

    PUSH_USERINT_INTO_VECTOR(partonFlavour);
    PUSH_USERINT_INTO_VECTOR(hadronFlavour);

    if (filledObjects) filledObjects->push_back(&(*obj));
    n_skimmed_objects++;
  }

  // Pass collections to the communicator
  PUSH_VECTOR_WITH_NAME(colName, pt);
  PUSH_VECTOR_WITH_NAME(colName, eta);
  PUSH_VECTOR_WITH_NAME(colName, phi);
  PUSH_VECTOR_WITH_NAME(colName, mass);

  //PUSH_VECTOR_WITH_NAME(colName, pass_looseId);
  //PUSH_VECTOR_WITH_NAME(colName, pass_tightId);

  PUSH_VECTOR_WITH_NAME(colName, n_pfcands);
  PUSH_VECTOR_WITH_NAME(colName, n_mucands);
  PUSH_VECTOR_WITH_NAME(colName, n_softdrop_subjets);

  PUSH_VECTOR_WITH_NAME(colName, area);
  PUSH_VECTOR_WITH_NAME(colName, pt_resolution);

  // Disable softdrop for now
  /*
  PUSH_VECTOR_WITH_NAME(colName, softdrop_pt);
  PUSH_VECTOR_WITH_NAME(colName, softdrop_eta);
  PUSH_VECTOR_WITH_NAME(colName, softdrop_phi);
  PUSH_VECTOR_WITH_NAME(colName, softdrop_mass);

  PUSH_VECTOR_WITH_NAME(colName, softdrop_subjet0_pt);
  PUSH_VECTOR_WITH_NAME(colName, softdrop_subjet0_eta);
  PUSH_VECTOR_WITH_NAME(colName, softdrop_subjet0_phi);
  PUSH_VECTOR_WITH_NAME(colName, softdrop_subjet0_mass);

  PUSH_VECTOR_WITH_NAME(colName, softdrop_subjet1_pt);
  PUSH_VECTOR_WITH_NAME(colName, softdrop_subjet1_eta);
  PUSH_VECTOR_WITH_NAME(colName, softdrop_subjet1_phi);
  PUSH_VECTOR_WITH_NAME(colName, softdrop_subjet1_mass);
  */

  PUSH_VECTOR_WITH_NAME(colName, ptDistribution);
  PUSH_VECTOR_WITH_NAME(colName, totalMultiplicity);
  PUSH_VECTOR_WITH_NAME(colName, axis1);
  PUSH_VECTOR_WITH_NAME(colName, axis2);

  PUSH_VECTOR_WITH_NAME(colName, JECNominal);
  PUSH_VECTOR_WITH_NAME(colName, JECUp);
  PUSH_VECTOR_WITH_NAME(colName, JECDn);

  PUSH_VECTOR_WITH_NAME(colName, JERNominal);
  PUSH_VECTOR_WITH_NAME(colName, JERUp);
  PUSH_VECTOR_WITH_NAME(colName, JERDn);

  PUSH_VECTOR_WITH_NAME(colName, partonFlavour);
  PUSH_VECTOR_WITH_NAME(colName, hadronFlavour);

  return n_skimmed_objects;
}
size_t CMS3Ntuplizer::fillIsotracks(edm::Event const& iEvent, std::vector<IsotrackInfo const*>* filledObjects){
#define PUSH_ISOTRACK_VARIABLE(NAME) NAME.push_back(obj->NAME);

  if (this->is80X) return 0;

  const char colName[] = "isotracks";
  edm::Handle< edm::View<IsotrackInfo> > isotracksHandle;
  iEvent.getByToken(isotracksToken, isotracksHandle);
  if (!isotracksHandle.isValid()) throw cms::Exception("CMS3Ntuplizer::fillIsotracks: Error getting the isotrack collection from the event...");
  size_t n_objects = isotracksHandle->size();

  if (filledObjects) filledObjects->reserve(n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, pt, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, eta, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, phi, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, mass, n_objects);

  MAKE_VECTOR_WITH_RESERVE(int, id, n_objects);
  MAKE_VECTOR_WITH_RESERVE(int, charge, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, pfIso03_ch, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, pfIso03_nh, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, pfIso03_em, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, pfIso03_db, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, pfIso03_comb_nofsr, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, miniIso_ch, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, miniIso_nh, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, miniIso_em, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, miniIso_db, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, miniIso_comb_nofsr, n_objects);

  MAKE_VECTOR_WITH_RESERVE(bool, fromPV, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, dxy, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, dxyerr, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, dz, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, dzerr, n_objects);

  MAKE_VECTOR_WITH_RESERVE(bool, is_pfCand, n_objects);
  MAKE_VECTOR_WITH_RESERVE(bool, is_lostTrack, n_objects);
  MAKE_VECTOR_WITH_RESERVE(bool, lepOverlap, n_objects);
  MAKE_VECTOR_WITH_RESERVE(bool, is_highPurityTrack, n_objects);
  MAKE_VECTOR_WITH_RESERVE(bool, is_tightTrack, n_objects);

  MAKE_VECTOR_WITH_RESERVE(int, nearestPFcand_id, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, nearestPFcand_deltaR, n_objects);

  size_t n_skimmed_objects=0;
  for (edm::View<IsotrackInfo>::const_iterator obj = isotracksHandle->begin(); obj != isotracksHandle->end(); obj++){
    if (!IsotrackSelectionHelpers::testSkimIsotrack(*obj, this->year)) continue;

    // Core particle quantities
    // Uncorrected p4
    pt.push_back(obj->p4.Pt());
    eta.push_back(obj->p4.Eta());
    phi.push_back(obj->p4.Phi());
    mass.push_back(obj->p4.M());

    PUSH_ISOTRACK_VARIABLE(id);
    PUSH_ISOTRACK_VARIABLE(charge);

    PUSH_ISOTRACK_VARIABLE(pfIso03_ch);
    PUSH_ISOTRACK_VARIABLE(pfIso03_nh);
    PUSH_ISOTRACK_VARIABLE(pfIso03_em);
    PUSH_ISOTRACK_VARIABLE(pfIso03_db);
    PUSH_ISOTRACK_VARIABLE(pfIso03_comb_nofsr);
    PUSH_ISOTRACK_VARIABLE(miniIso_ch);
    PUSH_ISOTRACK_VARIABLE(miniIso_nh);
    PUSH_ISOTRACK_VARIABLE(miniIso_em);
    PUSH_ISOTRACK_VARIABLE(miniIso_db);
    PUSH_ISOTRACK_VARIABLE(miniIso_comb_nofsr);

    PUSH_ISOTRACK_VARIABLE(fromPV);
    PUSH_ISOTRACK_VARIABLE(dxy);
    PUSH_ISOTRACK_VARIABLE(dxyerr);
    PUSH_ISOTRACK_VARIABLE(dz);
    PUSH_ISOTRACK_VARIABLE(dzerr);

    PUSH_ISOTRACK_VARIABLE(is_pfCand);
    PUSH_ISOTRACK_VARIABLE(is_lostTrack);
    PUSH_ISOTRACK_VARIABLE(lepOverlap);
    PUSH_ISOTRACK_VARIABLE(is_highPurityTrack);
    PUSH_ISOTRACK_VARIABLE(is_tightTrack);

    PUSH_ISOTRACK_VARIABLE(nearestPFcand_id);
    PUSH_ISOTRACK_VARIABLE(nearestPFcand_deltaR);

    if (filledObjects) filledObjects->push_back(&(*obj));
    n_skimmed_objects++;
  }

  // Pass collections to the communicator
  PUSH_VECTOR_WITH_NAME(colName, pt);
  PUSH_VECTOR_WITH_NAME(colName, eta);
  PUSH_VECTOR_WITH_NAME(colName, phi);
  PUSH_VECTOR_WITH_NAME(colName, mass);

  PUSH_VECTOR_WITH_NAME(colName, id);
  PUSH_VECTOR_WITH_NAME(colName, charge);

  PUSH_VECTOR_WITH_NAME(colName, pfIso03_ch);
  PUSH_VECTOR_WITH_NAME(colName, pfIso03_nh);
  PUSH_VECTOR_WITH_NAME(colName, pfIso03_em);
  PUSH_VECTOR_WITH_NAME(colName, pfIso03_db);
  PUSH_VECTOR_WITH_NAME(colName, pfIso03_comb_nofsr);
  PUSH_VECTOR_WITH_NAME(colName, miniIso_ch);
  PUSH_VECTOR_WITH_NAME(colName, miniIso_nh);
  PUSH_VECTOR_WITH_NAME(colName, miniIso_em);
  PUSH_VECTOR_WITH_NAME(colName, miniIso_db);
  PUSH_VECTOR_WITH_NAME(colName, miniIso_comb_nofsr);

  PUSH_VECTOR_WITH_NAME(colName, fromPV);
  PUSH_VECTOR_WITH_NAME(colName, dxy);
  PUSH_VECTOR_WITH_NAME(colName, dxyerr);
  PUSH_VECTOR_WITH_NAME(colName, dz);
  PUSH_VECTOR_WITH_NAME(colName, dzerr);

  PUSH_VECTOR_WITH_NAME(colName, is_pfCand);
  PUSH_VECTOR_WITH_NAME(colName, is_lostTrack);
  PUSH_VECTOR_WITH_NAME(colName, lepOverlap);
  PUSH_VECTOR_WITH_NAME(colName, is_highPurityTrack);
  PUSH_VECTOR_WITH_NAME(colName, is_tightTrack);

  PUSH_VECTOR_WITH_NAME(colName, nearestPFcand_id);
  PUSH_VECTOR_WITH_NAME(colName, nearestPFcand_deltaR);

  return n_skimmed_objects;

#undef PUSH_ISOTRACK_VARIABLE
}
size_t CMS3Ntuplizer::fillPFCandidates(edm::Event const& iEvent, std::vector<pat::Muon const*> const* filledMuons, std::vector<pat::Electron const*> const* filledElectrons, std::vector<pat::Photon const*> const* filledPhotons, std::vector<pat::PackedCandidate const*>* filledObjects){
  //const char colName[] = "pfcands";
  const char colNameFSR[] = "fsrcands";
  edm::Handle< edm::View<pat::PackedCandidate> > pfcandsHandle;
  iEvent.getByToken(pfcandsToken, pfcandsHandle);
  if (!pfcandsHandle.isValid()) throw cms::Exception("CMS3Ntuplizer::fillPFCandidates: Error getting the PF candidate collection from the event...");
  size_t n_objects = pfcandsHandle->size();
  size_t n_skimmed_objects=0;

  if (filledObjects) filledObjects->reserve(n_objects);

  // FSR preselection
  std::vector<FSRCandidateInfo> preselectedFSRCandidates; preselectedFSRCandidates.reserve(n_objects);
  edm::View<pat::PackedCandidate>::const_iterator it_pfcands_begin = pfcandsHandle->begin();
  edm::View<pat::PackedCandidate>::const_iterator it_pfcands_end = pfcandsHandle->end();
  for (edm::View<pat::PackedCandidate>::const_iterator obj = it_pfcands_begin; obj != it_pfcands_end; obj++){
    if (obj->pdgId()!=22) continue; // Check only photons

    if (!FSRSelectionHelpers::testSkimFSR_PtEta(*obj, this->year)) continue;

    double fsrInfo = FSRSelectionHelpers::fsrIso(*obj, this->year, it_pfcands_begin, it_pfcands_end);
    if (!FSRSelectionHelpers::testSkimFSR_Iso(*obj, this->year, fsrInfo)) continue;

    FSRCandidateInfo fsrinfo;
    fsrinfo.obj = &(*obj);
    fsrinfo.fsrIso =  fsrInfo;
    if (filledElectrons){ for (auto const& electron:(*filledElectrons)){ if (FSRSelectionHelpers::testSCVeto(&(*obj), electron)){ fsrinfo.veto_electron_list.push_back(electron); } } }
    if (filledPhotons){ for (auto const& photon:(*filledPhotons)){ if (FSRSelectionHelpers::testSCVeto(&(*obj), photon)){ fsrinfo.veto_photon_list.push_back(photon); } } }

    preselectedFSRCandidates.emplace_back(fsrinfo);
  }
  // Match FSR candidates to leptons
  std::vector<reco::LeafCandidate const*> leptons;
  if (filledMuons){ for (auto const& muon:(*filledMuons)) leptons.push_back(muon); }
  if (filledElectrons){ for (auto const& electron:(*filledElectrons)) leptons.push_back(electron); }
  std::unordered_map< FSRCandidateInfo const*, std::vector<reco::LeafCandidate const*> > fsrcand_lepton_map;
  CMS3ObjectHelpers::matchParticles_OneToMany(
    CMS3ObjectHelpers::kMatchBy_DeltaR, FSRSelectionHelpers::selection_match_fsr_deltaR,
    preselectedFSRCandidates.begin(), preselectedFSRCandidates.end(),
    leptons.cbegin(), leptons.cend(),
    fsrcand_lepton_map
  );
  std::vector<FSRCandidateInfo*> writableFSRCandidates;
  for (auto& it:fsrcand_lepton_map){
    FSRCandidateInfo const* fsrcand_const = it.first;
    FSRCandidateInfo* fsrcand = nullptr;
    for (auto& cand:preselectedFSRCandidates){ if (&cand == fsrcand_const) fsrcand = &cand; }
    reco::LeafCandidate const* bestLepton = nullptr;
    for (auto const& lepton:it.second){
      pat::Muon const* muon = dynamic_cast<pat::Muon const*>(lepton);
      pat::Electron const* electron = dynamic_cast<pat::Electron const*>(lepton);
      if (muon){
        if (!bestLepton) bestLepton = muon;
        fsrcand->matched_muon_list.push_back(muon);
      }
      else if (electron){
        if (HelperFunctions::checkListVariable(fsrcand->veto_electron_list, electron)) continue;
        if (!bestLepton) bestLepton = electron;
        fsrcand->matched_electron_list.push_back(electron);
      }
    }
    if (!bestLepton) continue;

    double mindr = reco::deltaR(fsrcand->p4(), bestLepton->p4());
    if (!FSRSelectionHelpers::testSkimFSR_MinDeltaR(*(fsrcand->obj), this->year, mindr)) continue;

    writableFSRCandidates.push_back(fsrcand);
  }

  // Fill FSR candidates
  MAKE_VECTOR_WITH_RESERVE(float, pt, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, eta, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, phi, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, mass, n_objects);
  MAKE_VECTOR_WITH_RESERVE(std::vector<unsigned int>, fsrMatch_muon_index_list, n_objects);
  MAKE_VECTOR_WITH_RESERVE(std::vector<unsigned int>, fsrMatch_electron_index_list, n_objects);
  MAKE_VECTOR_WITH_RESERVE(std::vector<unsigned int>, photonVeto_index_list, n_objects);
  for (auto const& obj:writableFSRCandidates){
    pt.push_back(obj->pt());
    eta.push_back(obj->eta());
    phi.push_back(obj->phi());
    mass.push_back(obj->mass());

    // Get the dR-ordered matched lepton index lists
    std::vector<unsigned int> muon_indices;
    for (auto const& muon:obj->matched_muon_list){
      unsigned int imuon=0;
      for (auto const& obj:(*filledMuons)){
        if (obj==muon){
          muon_indices.push_back(imuon);
          break;
        }
        imuon++;
      }
    }
    fsrMatch_muon_index_list.push_back(muon_indices);

    std::vector<unsigned int> electron_indices;
    for (auto const& electron:obj->matched_electron_list){
      unsigned int ielectron=0;
      for (auto const& obj:(*filledElectrons)){
        if (obj==electron){
          electron_indices.push_back(ielectron);
          break;
        }
        ielectron++;
      }
    }
    fsrMatch_electron_index_list.push_back(electron_indices);

    std::vector<unsigned int> photonVeto_indices;
    for (auto const& photon:obj->veto_photon_list){
      unsigned int iphoton=0;
      for (auto const& obj:(*filledPhotons)){
        if (obj==photon){
          photonVeto_indices.push_back(iphoton);
          break;
        }
        iphoton++;
      }
    }
    photonVeto_index_list.push_back(photonVeto_indices);

    n_skimmed_objects++;
  }
  PUSH_VECTOR_WITH_NAME(colNameFSR, pt);
  PUSH_VECTOR_WITH_NAME(colNameFSR, eta);
  PUSH_VECTOR_WITH_NAME(colNameFSR, phi);
  PUSH_VECTOR_WITH_NAME(colNameFSR, mass);
  PUSH_VECTOR_WITH_NAME(colNameFSR, fsrMatch_muon_index_list);
  PUSH_VECTOR_WITH_NAME(colNameFSR, fsrMatch_electron_index_list);
  PUSH_VECTOR_WITH_NAME(colNameFSR, photonVeto_index_list);

  return n_skimmed_objects;
}

size_t CMS3Ntuplizer::fillVertices(edm::Event const& iEvent, std::vector<reco::Vertex const*>* filledObjects){
  const char colName[] = "vtxs";
  edm::Handle< reco::VertexCollection > vtxHandle;
  iEvent.getByToken(vtxToken, vtxHandle);
  if (!vtxHandle.isValid()) throw cms::Exception("CMS3Ntuplizer::fillVertices: Error getting the vertex collection from the event...");
  size_t n_objects = vtxHandle->size();

  if (filledObjects) filledObjects->reserve(n_objects);

  MAKE_VECTOR_WITH_RESERVE(bool, is_fake, n_objects);
  MAKE_VECTOR_WITH_RESERVE(bool, is_valid, n_objects);
  MAKE_VECTOR_WITH_RESERVE(bool, is_good, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, ndof, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, pos_x, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, pos_y, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, pos_z, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, pos_t, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, pos_dx, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, pos_dy, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, pos_dz, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, pos_dt, n_objects);

  bool didFirstVertex = false;
  bool didFirstGoodVertex = false;
  unsigned int nvtxs=0, nvtxs_good=0;
  for (reco::VertexCollection::const_iterator obj = vtxHandle->begin(); obj != vtxHandle->end(); obj++){
    bool isGoodVtx = VertexSelectionHelpers::testGoodVertex(*obj);

    // Avoid recording all vertices
    if (!didFirstVertex || (!didFirstGoodVertex && isGoodVtx)){
      auto const& pos = obj->position();

      is_fake.push_back(obj->isFake());
      is_valid.push_back(obj->isValid());
      is_good.push_back(isGoodVtx);

      ndof.push_back(obj->ndof());

      pos_x.push_back(pos.x());
      pos_y.push_back(pos.y());
      pos_z.push_back(pos.z());
      pos_t.push_back(obj->t());

      pos_dx.push_back(obj->xError());
      pos_dy.push_back(obj->yError());
      pos_dz.push_back(obj->zError());
      pos_dt.push_back(obj->tError());

      if (filledObjects) filledObjects->push_back(&(*obj));

      if (!didFirstVertex) didFirstVertex=true;
      if (isGoodVtx) didFirstGoodVertex=true;
    }

    if (isGoodVtx) nvtxs_good++;
    nvtxs++;
  }

  // Record the counts
  commonEntry.setNamedVal(TString(colName)+"_nvtxs", nvtxs);
  commonEntry.setNamedVal(TString(colName)+"_nvtxs_good", nvtxs_good);

  // Pass collections to the communicator
  PUSH_VECTOR_WITH_NAME(colName, is_fake);
  PUSH_VECTOR_WITH_NAME(colName, is_valid);
  PUSH_VECTOR_WITH_NAME(colName, is_good);

  PUSH_VECTOR_WITH_NAME(colName, ndof);

  PUSH_VECTOR_WITH_NAME(colName, pos_x);
  PUSH_VECTOR_WITH_NAME(colName, pos_y);
  PUSH_VECTOR_WITH_NAME(colName, pos_z);
  PUSH_VECTOR_WITH_NAME(colName, pos_t);

  PUSH_VECTOR_WITH_NAME(colName, pos_dx);
  PUSH_VECTOR_WITH_NAME(colName, pos_dy);
  PUSH_VECTOR_WITH_NAME(colName, pos_dz);
  PUSH_VECTOR_WITH_NAME(colName, pos_dt);

  return nvtxs_good;
}

bool CMS3Ntuplizer::fillEventVariables(edm::Event const& iEvent){
  edm::Handle< double > rhoHandle;
  iEvent.getByToken(rhoToken, rhoHandle);
  if (!rhoHandle.isValid()) throw cms::Exception("CMS3Ntuplizer::fillEventVariables: Error getting the rho collection from the event...");

  // Simple event-level variables
  commonEntry.setNamedVal("EventNumber", iEvent.id().event());
  commonEntry.setNamedVal("RunNumber", iEvent.id().run());
  commonEntry.setNamedVal("LuminosityBlock", iEvent.luminosityBlock());
  commonEntry.setNamedVal("event_rho", (float) (*rhoHandle));
  if (isMC){
    edm::Handle< std::vector<PileupSummaryInfo> > puInfoHandle;
    iEvent.getByToken(puInfoToken, puInfoHandle);
    if (!puInfoHandle.isValid()) throw cms::Exception("CMS3Ntuplizer::fillEventVariables: Error getting the PU info from the event...");

    commonEntry.setNamedVal("n_vtxs_PU", (int) ((*puInfoHandle)[0].getPU_NumInteractions()));
    commonEntry.setNamedVal("n_true_int", (float) ((*puInfoHandle)[0].getTrueNumInteractions()));
  }
  else{
    commonEntry.setNamedVal("n_vtxs_PU", (int) -1);
    commonEntry.setNamedVal("n_true_int", (float) -1);
  }

  if (applyPrefiringWeights){
    edm::Handle< double > prefiringweight;
    float prefiringweightval=0;

    iEvent.getByToken(prefiringWeightToken, prefiringweight);
    if (!prefiringweight.isValid()) throw cms::Exception("CMS3Ntuplizer::fillEventVariables: Error getting the nominal prefiring weight from the event...");
    prefiringweightval = (*prefiringweight);
    commonEntry.setNamedVal("prefiringWeight_Nominal", prefiringweightval);

    iEvent.getByToken(prefiringWeightToken_Dn, prefiringweight);
    if (!prefiringweight.isValid()) throw cms::Exception("CMS3Ntuplizer::fillEventVariables: Error getting the prefiring weight down variation from the event...");
    prefiringweightval = (*prefiringweight);
    commonEntry.setNamedVal("prefiringWeight_Dn", prefiringweightval);

    iEvent.getByToken(prefiringWeightToken_Up, prefiringweight);
    if (!prefiringweight.isValid()) throw cms::Exception("CMS3Ntuplizer::fillEventVariables: Error getting the prefiring weight up variation from the event...");
    prefiringweightval = (*prefiringweight);
    commonEntry.setNamedVal("prefiringWeight_Up", prefiringweightval);
  }

  return true;
}
bool CMS3Ntuplizer::fillTriggerInfo(edm::Event const& iEvent){
  const char colName[] = "triggers";
  edm::Handle< edm::View<TriggerInfo> > triggerInfoHandle;
  iEvent.getByToken(triggerInfoToken, triggerInfoHandle);
  if (!triggerInfoHandle.isValid()) throw cms::Exception("CMS3Ntuplizer::fillTriggerInfo: Error getting the trigger infos. from the event...");
  size_t n_triggers = triggerInfoHandle->size();

  MAKE_VECTOR_WITH_RESERVE(std::string, name, n_triggers);
  MAKE_VECTOR_WITH_RESERVE(bool, passTrigger, n_triggers);
  MAKE_VECTOR_WITH_RESERVE(int, L1prescale, n_triggers);
  MAKE_VECTOR_WITH_RESERVE(int, HLTprescale, n_triggers);

  bool passAtLeastOneTrigger = false;
  for (edm::View<TriggerInfo>::const_iterator obj = triggerInfoHandle->begin(); obj != triggerInfoHandle->end(); obj++){
    name.emplace_back(obj->name);
    passTrigger.emplace_back(obj->passTrigger);
    L1prescale.emplace_back(obj->L1prescale);
    HLTprescale.emplace_back(obj->HLTprescale);

    passAtLeastOneTrigger |= obj->passTrigger;
  }

  PUSH_VECTOR_WITH_NAME(colName, name);
  PUSH_VECTOR_WITH_NAME(colName, passTrigger);
  PUSH_VECTOR_WITH_NAME(colName, L1prescale);
  PUSH_VECTOR_WITH_NAME(colName, HLTprescale);

  // If the (data) event does not pass any triggers, do not record it.
  return passAtLeastOneTrigger;
}
bool CMS3Ntuplizer::fillMETFilterVariables(edm::Event const& iEvent){
  // See https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2 for recommendations
  // See also PhysicsTools/PatAlgos/python/slimming/metFilterPaths_cff.py for the collection names
  const char metFiltCollName[] = "metfilter";
  edm::Handle<METFilterInfo> metFilterInfoHandle;
  iEvent.getByToken(metFilterInfoToken, metFilterInfoHandle);
  if (!metFilterInfoHandle.isValid()) throw cms::Exception("CMS3Ntuplizer::fillMETFilterVariables: Error getting the MET filter handle from the event...");

  static const std::vector<std::string> metfilterflags{
    "CSCTightHaloFilter",
    "CSCTightHalo2015Filter",
    "globalTightHalo2016Filter",
    "globalSuperTightHalo2016Filter",
    "HBHENoiseFilter",
    "EcalDeadCellTriggerPrimitiveFilter",
    "hcalLaserEventFilter",
    "trackingFailureFilter",
    "chargedHadronTrackResolutionFilter",
    "eeBadScFilter",
    "ecalLaserCorrFilter",
    "METFilters",
    "goodVertices",
    "trkPOGFilters",
    "trkPOG_logErrorTooManyClusters",
    "trkPOG_manystripclus53X",
    "trkPOG_toomanystripclus53X",
    "HBHENoiseIsoFilter",
    "CSCTightHaloTrkMuUnvetoFilter",
    "HcalStripHaloFilter",
    "EcalDeadCellBoundaryEnergyFilter",
    "muonBadTrackFilter",
    "BadPFMuonFilter",
    "BadChargedCandidateFilter",
    "ecalBadCalibFilter",
    "ecalBadCalibFilterUpdated"
  };
  for (auto const& flagname:metfilterflags){
    bool flag = false;
    auto it_flag = metFilterInfoHandle->flag_accept_map.find(flagname);
    if (it_flag!=metFilterInfoHandle->flag_accept_map.cend()) flag = it_flag->second;
    commonEntry.setNamedVal((std::string(metFiltCollName) + "_" + flagname).data(), flag);
  }

  return true;
}
bool CMS3Ntuplizer::fillMETVariables(edm::Event const& iEvent){
#define SET_MET_VARIABLE(HANDLE, NAME, COLLNAME) commonEntry.setNamedVal((std::string(COLLNAME) + "_" + #NAME).data(), HANDLE->NAME);

  edm::Handle<METInfo> metHandle;

  // PF MET
  const char pfmetCollName[] = "pfmet";
  iEvent.getByToken(pfmetToken, metHandle);
  if (!metHandle.isValid()) throw cms::Exception("CMS3Ntuplizer::fillMETVariables: Error getting the PF MET handle from the event...");

  SET_MET_VARIABLE(metHandle, met_Nominal, pfmetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_Nominal, pfmetCollName);
  SET_MET_VARIABLE(metHandle, sumEt_Nominal, pfmetCollName);
  SET_MET_VARIABLE(metHandle, metSignificance, pfmetCollName);
  SET_MET_VARIABLE(metHandle, met_over_sqrtSumEt, pfmetCollName);

  SET_MET_VARIABLE(metHandle, met_Raw, pfmetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_Raw, pfmetCollName);
  SET_MET_VARIABLE(metHandle, sumEt_Raw, pfmetCollName);

  SET_MET_VARIABLE(metHandle, met_JERUp, pfmetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_JERUp, pfmetCollName);
  SET_MET_VARIABLE(metHandle, met_JERDn, pfmetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_JERDn, pfmetCollName);

  SET_MET_VARIABLE(metHandle, met_JECUp, pfmetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_JECUp, pfmetCollName);
  SET_MET_VARIABLE(metHandle, met_JECDn, pfmetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_JECDn, pfmetCollName);

  /*
  SET_MET_VARIABLE(metHandle, met_MuonEnUp, pfmetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_MuonEnUp, pfmetCollName);
  SET_MET_VARIABLE(metHandle, met_MuonEnDn, pfmetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_MuonEnDn, pfmetCollName);

  SET_MET_VARIABLE(metHandle, met_ElectronEnUp, pfmetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_ElectronEnUp, pfmetCollName);
  SET_MET_VARIABLE(metHandle, met_ElectronEnDn, pfmetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_ElectronEnDn, pfmetCollName);

  SET_MET_VARIABLE(metHandle, met_TauEnUp, pfmetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_TauEnUp, pfmetCollName);
  SET_MET_VARIABLE(metHandle, met_TauEnDn, pfmetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_TauEnDn, pfmetCollName);

  SET_MET_VARIABLE(metHandle, met_UnclusteredEnUp, pfmetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_UnclusteredEnUp, pfmetCollName);
  SET_MET_VARIABLE(metHandle, met_UnclusteredEnDn, pfmetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_UnclusteredEnDn, pfmetCollName);

  SET_MET_VARIABLE(metHandle, met_PhotonEnUp, pfmetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_PhotonEnUp, pfmetCollName);
  SET_MET_VARIABLE(metHandle, met_PhotonEnDn, pfmetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_PhotonEnDn, pfmetCollName);
  */

  SET_MET_VARIABLE(metHandle, calo_met, pfmetCollName);
  SET_MET_VARIABLE(metHandle, calo_metPhi, pfmetCollName);

  /*
  SET_MET_VARIABLE(metHandle, gen_met, pfmetCollName);
  SET_MET_VARIABLE(metHandle, gen_metPhi, pfmetCollName);
  */

  // PUPPI MET
  const char puppimetCollName[] = "puppimet";
  iEvent.getByToken(puppimetToken, metHandle);
  if (!metHandle.isValid()) throw cms::Exception("CMS3Ntuplizer::fillMETVariables: Error getting the PUPPI MET handle from the event...");

  SET_MET_VARIABLE(metHandle, met_Nominal, puppimetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_Nominal, puppimetCollName);
  SET_MET_VARIABLE(metHandle, sumEt_Nominal, puppimetCollName);
  SET_MET_VARIABLE(metHandle, metSignificance, puppimetCollName);
  SET_MET_VARIABLE(metHandle, met_over_sqrtSumEt, puppimetCollName);

  SET_MET_VARIABLE(metHandle, met_Raw, puppimetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_Raw, puppimetCollName);
  SET_MET_VARIABLE(metHandle, sumEt_Raw, puppimetCollName);

  SET_MET_VARIABLE(metHandle, met_JERUp, puppimetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_JERUp, puppimetCollName);
  SET_MET_VARIABLE(metHandle, met_JERDn, puppimetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_JERDn, puppimetCollName);

  SET_MET_VARIABLE(metHandle, met_JECUp, puppimetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_JECUp, puppimetCollName);
  SET_MET_VARIABLE(metHandle, met_JECDn, puppimetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_JECDn, puppimetCollName);

  /*
  SET_MET_VARIABLE(metHandle, met_MuonEnUp, puppimetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_MuonEnUp, puppimetCollName);
  SET_MET_VARIABLE(metHandle, met_MuonEnDn, puppimetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_MuonEnDn, puppimetCollName);

  SET_MET_VARIABLE(metHandle, met_ElectronEnUp, puppimetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_ElectronEnUp, puppimetCollName);
  SET_MET_VARIABLE(metHandle, met_ElectronEnDn, puppimetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_ElectronEnDn, puppimetCollName);

  SET_MET_VARIABLE(metHandle, met_TauEnUp, puppimetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_TauEnUp, puppimetCollName);
  SET_MET_VARIABLE(metHandle, met_TauEnDn, puppimetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_TauEnDn, puppimetCollName);

  SET_MET_VARIABLE(metHandle, met_UnclusteredEnUp, puppimetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_UnclusteredEnUp, puppimetCollName);
  SET_MET_VARIABLE(metHandle, met_UnclusteredEnDn, puppimetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_UnclusteredEnDn, puppimetCollName);

  SET_MET_VARIABLE(metHandle, met_PhotonEnUp, puppimetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_PhotonEnUp, puppimetCollName);
  SET_MET_VARIABLE(metHandle, met_PhotonEnDn, puppimetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_PhotonEnDn, puppimetCollName);
  */

  /*
  SET_MET_VARIABLE(metHandle, calo_met, puppimetCollName);
  SET_MET_VARIABLE(metHandle, calo_metPhi, puppimetCollName);

  SET_MET_VARIABLE(metHandle, gen_met, puppimetCollName);
  SET_MET_VARIABLE(metHandle, gen_metPhi, puppimetCollName);
  */

#undef SET_MET_VARIABLE

  return true;
}

bool CMS3Ntuplizer::fillGenVariables(
  edm::Event const& iEvent,
  std::vector<reco::GenParticle const*>* filledPrunedGenParts,
  std::vector<pat::PackedGenParticle const*>* filledPackedGenParts,
  std::vector<reco::GenJet const*>* filledGenAK4Jets,
  std::vector<reco::GenJet const*>* filledGenAK8Jets
){
  if (!this->isMC) return true;

  // Gen. info.
  recordGenInfo(iEvent);

  // Gen. particles
  if (this->keepGenParticles!=kNone) recordGenParticles(iEvent, filledPrunedGenParts, filledPackedGenParts);

  // Gen. jets
  if (this->keepGenJets){
    recordGenJets(iEvent, false, filledGenAK4Jets); // ak4jets
    recordGenJets(iEvent, true, filledGenAK8Jets); // ak8jets
  }
  return true;
}


// Undefine the convenience macros
#undef PUSH_VECTOR_WITH_NAME
#undef PUSH_USERFLOAT_INTO_VECTOR
#undef PUSH_USERINT_INTO_VECTOR
#undef MAKE_VECTOR_WITH_RESERVE


CMS3Ntuplizer::ParticleRecordLevel CMS3Ntuplizer::getParticleRecordLevel(std::string str){
  std::string word;
  HelperFunctions::lowercase(str, word);

  ParticleRecordLevel res = nParticleRecordLevels;
  for (unsigned char i=kNone; i!=nParticleRecordLevels; i++){
    if (
      (i==kNone && word=="none")
      ||
      (i==kReducedFinalStates && word=="reducedfinalstates")
      ||
      (i==kAllFinalStates && word=="allfinalstates")
      ||
      (i==kReducedFinalStatesAndHardProcesses && word=="reducedfinalstatesandhardprocesses")
      ||
      (i==kAll && word=="all")
      ) res = static_cast<ParticleRecordLevel>(i);
  }
  if (res==nParticleRecordLevels) throw cms::Exception((std::string("CMS3Ntuplizer::getParticleRecordLevel: Word ") + word + " is not recognized!").data());
  return res;
}


DEFINE_FWK_MODULE(CMS3Ntuplizer);
