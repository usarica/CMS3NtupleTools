#include <cctype>
#include <algorithm>

#include <CommonTools/UtilAlgos/interface/TFileService.h>
#include <DataFormats/EgammaReco/interface/SuperCluster.h>
#include <DataFormats/EcalDetId/interface/EBDetId.h>
#include <DataFormats/EcalDetId/interface/EEDetId.h>
#include <DataFormats/EcalDetId/interface/EcalSubdetector.h>
#include <DataFormats/ForwardDetId/interface/ForwardSubdetector.h>
#include <DataFormats/ForwardDetId/interface/HGCalDetId.h>

#include <RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h>
#include <RecoEcal/EgammaCoreTools/interface/EcalTools.h>

#include <CMS3/NtupleMaker/interface/plugins/CMS3Ntuplizer.h>
#include "CMS3/NtupleMaker/interface/VertexSelectionHelpers.h"
#include "CMS3/NtupleMaker/interface/MuonSelectionHelpers.h"
#include "CMS3/NtupleMaker/interface/ElectronSelectionHelpers.h"
#include "CMS3/NtupleMaker/interface/PhotonSelectionHelpers.h"
#include "CMS3/NtupleMaker/interface/FSRSelectionHelpers.h"
#include "CMS3/NtupleMaker/interface/AK4JetSelectionHelpers.h"
#include "CMS3/NtupleMaker/interface/AK8JetSelectionHelpers.h"
#include "CMS3/NtupleMaker/interface/IsotrackSelectionHelpers.h"
#include <CMS3/NtupleMaker/interface/CMS3ObjectHelpers.h>
#include <CMS3/NtupleMaker/interface/MCUtilities.h>

#include <CMS3/Dictionaries/interface/CommonTypedefs.h>
#include <CMS3/Dictionaries/interface/EgammaFiduciality.h>

#include <CMSDataTools/AnalysisTree/interface/HelperFunctionsCore.h>

#include "MELAStreamHelpers.hh"


using namespace std;
using namespace edm;
using namespace MELAStreamHelpers;


// Collection names
const std::string CMS3Ntuplizer::colName_muons = "muons";
const std::string CMS3Ntuplizer::colName_electrons = "electrons";
const std::string CMS3Ntuplizer::colName_photons = "photons";
const std::string CMS3Ntuplizer::colName_fsrcands = "fsrcands";
const std::string CMS3Ntuplizer::colName_superclusters = "superclusters";
const std::string CMS3Ntuplizer::colName_isotracks = "isotracks";
const std::string CMS3Ntuplizer::colName_ak4jets = "ak4jets";
const std::string CMS3Ntuplizer::colName_ak8jets = "ak8jets";
const std::string CMS3Ntuplizer::colName_pfcands = "pfcands";
const std::string CMS3Ntuplizer::colName_vtxs = "vtxs";
const std::string CMS3Ntuplizer::colName_triggerinfos = "triggers";
const std::string CMS3Ntuplizer::colName_triggerobjects = "triggerObjects";
const std::string CMS3Ntuplizer::colName_genparticles = "genparticles";

CMS3Ntuplizer::CMS3Ntuplizer(const edm::ParameterSet& pset_) :
  pset(pset_),
  outtree(nullptr),
  firstEvent(true),

  year(pset.getParameter<int>("year")),
  treename(pset.getUntrackedParameter<std::string>("treename")),
  isMC(pset.getParameter<bool>("isMC")),
  is80X(pset.getParameter<bool>("is80X")),

  processTriggerObjectInfos(pset.getParameter<bool>("processTriggerObjectInfos")),

  keepMuonTimingInfo(pset.getParameter<bool>("keepMuonTimingInfo")),
  keepMuonPullInfo(pset.getParameter<bool>("keepMuonPullInfo")),

  keepElectronMVAInfo(pset.getParameter<bool>("keepElectronMVAInfo")),

  keepExtraSuperclusters(pset.getParameter<bool>("keepExtraSuperclusters")),

  prefiringWeightsTag(pset.getUntrackedParameter<std::string>("prefiringWeightsTag")),
  applyPrefiringWeights(prefiringWeightsTag!=""),

  keepGenParticles(CMS3Ntuplizer::getParticleRecordLevel(pset.getUntrackedParameter<std::string>("keepGenParticles"))),
  keepGenJets(pset.getParameter<bool>("keepGenJets")),

  includeLJetsSelection(pset.getParameter<bool>("includeLJetsSelection")),
  minNmuons(pset.getParameter<int>("minNmuons")),
  minNelectrons(pset.getParameter<int>("minNelectrons")),
  minNleptons(pset.getParameter<int>("minNleptons")),
  minNphotons(pset.getParameter<int>("minNphotons")),
  minNak4jets(pset.getParameter<int>("minNak4jets")),
  minNak8jets(pset.getParameter<int>("minNak8jets"))
{
  if (year!=2016 && year!=2017 && year!=2018) throw cms::Exception("CMS3Ntuplizer::CMS3Ntuplizer: Year is undefined!");

  muonsToken = consumes< edm::View<pat::Muon> >(pset.getParameter<edm::InputTag>("muonSrc"));
  electronsToken = consumes< edm::View<pat::Electron> >(pset.getParameter<edm::InputTag>("electronSrc"));
  photonsToken = consumes< edm::View<pat::Photon> >(pset.getParameter<edm::InputTag>("photonSrc"));
  ak4jetsToken = consumes< edm::View<pat::Jet> >(pset.getParameter<edm::InputTag>("ak4jetSrc"));
  ak8jetsToken = consumes< edm::View<pat::Jet> >(pset.getParameter<edm::InputTag>("ak8jetSrc"));
  isotracksToken = consumes< edm::View<IsotrackInfo> >(pset.getParameter<edm::InputTag>("isotrackSrc"));
  pfcandsToken = consumes< edm::View<pat::PackedCandidate> >(pset.getParameter<edm::InputTag>("pfcandSrc"));

  if (keepExtraSuperclusters){
    reducedSuperclusterToken = consumes< reco::SuperClusterCollection >(pset.getParameter<edm::InputTag>("reducedSuperclusterSrc"));
  }

  vtxToken = consumes<reco::VertexCollection>(pset.getParameter<edm::InputTag>("vtxSrc"));

  rhoToken  = consumes< double >(pset.getParameter<edm::InputTag>("rhoSrc"));
  triggerInfoToken = consumes< edm::View<TriggerInfo> >(pset.getParameter<edm::InputTag>("triggerInfoSrc"));
  if (processTriggerObjectInfos) triggerObjectInfoToken = consumes< edm::View<TriggerObjectInfo> >(pset.getParameter<edm::InputTag>("triggerObjectInfoSrc"));
  puInfoToken = consumes< std::vector<PileupSummaryInfo> >(pset.getParameter<edm::InputTag>("puInfoSrc"));
  metFilterInfoToken = consumes< METFilterInfo >(pset.getParameter<edm::InputTag>("metFilterInfoSrc"));

  if (applyPrefiringWeights){
    prefiringWeightToken = consumes< double >(edm::InputTag(prefiringWeightsTag, "nonPrefiringProb"));
    prefiringWeightToken_Dn = consumes< double >(edm::InputTag(prefiringWeightsTag, "nonPrefiringProbDown"));
    prefiringWeightToken_Up = consumes< double >(edm::InputTag(prefiringWeightsTag, "nonPrefiringProbUp"));
  }

  pfmetToken = consumes< METInfo >(pset.getParameter<edm::InputTag>("pfmetSrc"));
  puppimetToken = consumes< METInfo >(pset.getParameter<edm::InputTag>("puppimetSrc"));
  if (isMC){
    pfmetshiftToken_JERNominal = consumes< reco::Particle::LorentzVector >(pset.getParameter<edm::InputTag>("pfmetShiftSrc_JERNominal"));
    pfmetshiftToken_JERUp = consumes< reco::Particle::LorentzVector >(pset.getParameter<edm::InputTag>("pfmetShiftSrc_JERUp"));
    pfmetshiftToken_JERDn = consumes< reco::Particle::LorentzVector >(pset.getParameter<edm::InputTag>("pfmetShiftSrc_JERDn"));

    genInfoToken = consumes< GenInfo >(pset.getParameter<edm::InputTag>("genInfoSrc"));
    prunedGenParticlesToken = consumes< reco::GenParticleCollection >(pset.getParameter<edm::InputTag>("prunedGenParticlesSrc"));
    packedGenParticlesToken = consumes< pat::PackedGenParticleCollection >(pset.getParameter<edm::InputTag>("packedGenParticlesSrc"));
    genAK4JetsToken = consumes< edm::View<reco::GenJet> >(pset.getParameter<edm::InputTag>("genAK4JetsSrc"));
    genAK8JetsToken = consumes< edm::View<reco::GenJet> >(pset.getParameter<edm::InputTag>("genAK8JetsSrc"));
  }

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
#define MAKE_VECTOR_WITHOUT_RESERVE(type_, name_) std::vector<type_> name_;
#define RESERVE_VECTOR(name_, size_) name_.reserve(size_);
#define MAKE_VECTOR_WITH_RESERVE(type_, name_, size_) std::vector<type_> name_; name_.reserve(size_);
#define MAKE_VECTOR_WITH_DEFAULT_ASSIGN(type_, name_, size_) std::vector<type_> name_(size_);
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

  // Vertices
  size_t n_vtxs = this->fillVertices(iEvent, nullptr);
  isSelected &= (n_vtxs>0);

  // Muons
  std::vector<pat::Muon const*> filledMuons;
  size_t n_muons = this->fillMuons(iEvent, &filledMuons);

  // Electrons
  std::vector<pat::Electron const*> filledElectrons;
  size_t n_electrons = this->fillElectrons(iEvent, &filledElectrons);

  // FSR candidates
  std::vector<FSRCandidateInfo> filledFSRInfos;
  /*size_t n_fsrcands = */this->fillFSRCandidates(
    iEvent,
    &filledMuons, &filledElectrons,
    &filledFSRInfos
  );

  // Photons
  std::vector<pat::Photon const*> filledPhotons;
  size_t n_photons = this->fillPhotons(iEvent, &filledFSRInfos, &filledPhotons);

  //std::vector<reco::SuperCluster const*>* filledReducedSuperclusters;
  /*size_t n_reducedSuperclusters = */this->fillReducedSuperclusters(iEvent, &filledElectrons, &filledPhotons, /*filledReducedSuperclusters*/nullptr);

  // ak4 jets
  //std::vector<pat::Jet const*> filledAK4Jets;
  size_t n_ak4jets = this->fillAK4Jets(iEvent, /*&filledAK4Jets*/nullptr);

  // ak8 jets
  //std::vector<pat::Jet const*> filledAK8Jets;
  size_t n_ak8jets = this->fillAK8Jets(iEvent, /*&filledAK8Jets*/nullptr);

  // Isolated tracks
  /*size_t n_isotracks = */this->fillIsotracks(iEvent, nullptr);

  // Gen. variables
  std::vector<reco::GenParticle const*> filledPrunedGenParts;
  std::vector<pat::PackedGenParticle const*> filledPackedGenParts;
  std::vector<reco::GenJet const*> filledGenAK4Jets;
  std::vector<reco::GenJet const*> filledGenAK8Jets;
  if (this->isMC) isSelected &= this->fillGenVariables(
    iEvent,
    &filledMuons, &filledElectrons, &filledPhotons,
    // No need to pass reco jets since gen.-matching info is already filled.
    &filledPrunedGenParts, &filledPackedGenParts,
    &filledGenAK4Jets, &filledGenAK8Jets
  );

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
  if (includeLJetsSelection) passNobjects |= ((n_muons+n_electrons)>=1 && (n_ak4jets+n_ak8jets)>=1);
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

    outtree->getSelectedTree()->SetBasketSize((CMS3Ntuplizer::colName_triggerinfos+"_*").data(), 16384*23);
    //outtree->getSelectedTree()->SetBasketSize((CMS3Ntuplizer::colName_triggerinfos+"_*").data(), 21846*32);
    outtree->getSelectedTree()->SetBasketSize((CMS3Ntuplizer::colName_triggerobjects+"_passedTriggers").data(), 64000);
    outtree->getSelectedTree()->SetBasketSize((CMS3Ntuplizer::colName_triggerobjects+"_associatedTriggers").data(), 64000);
  }

  // Record whatever is in commonEntry into the tree.
#define SIMPLE_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=commonEntry.named##name_t##s.begin(); itb!=commonEntry.named##name_t##s.end(); itb++) outtree->setVal(itb->first, itb->second);
#define VECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=commonEntry.namedV##name_t##s.begin(); itb!=commonEntry.namedV##name_t##s.end(); itb++){ cleanUnusedCollection(isSelected, itb->first, itb->second); outtree->setVal(itb->first, &(itb->second)); }
#define DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=commonEntry.namedVV##name_t##s.begin(); itb!=commonEntry.namedVV##name_t##s.end(); itb++){ cleanUnusedCollection(isSelected, itb->first, itb->second); outtree->setVal(itb->first, &(itb->second)); }
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

  SET_GENINFO_VARIABLE(genjets_HT);
  SET_GENINFO_VARIABLE(genjets_MHT);

  SET_GENINFO_VARIABLE(genmet_met);
  SET_GENINFO_VARIABLE(genmet_metPhi);

  SET_GENINFO_VARIABLE(sumEt);
  SET_GENINFO_VARIABLE(pThat);

  // Number of shower gluons that decay to charms and bottoms
  // Useful to assign a scaling or systematic for g->cc/bb
  SET_GENINFO_VARIABLE(n_shower_gluons_to_bottom);
  SET_GENINFO_VARIABLE(n_shower_gluons_to_charm);

  // LHE variations
  SET_GENINFO_VARIABLE(genHEPMCweight_default);
  SET_GENINFO_VARIABLE(genHEPMCweight_NNPDF30);

  SET_GENINFO_VARIABLE(LHEweight_scaledOriginalWeight_default);
  SET_GENINFO_VARIABLE(LHEweight_scaledOriginalWeight_NNPDF30);

  SET_GENINFO_VARIABLE(LHEweight_defaultMemberZero);

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
void CMS3Ntuplizer::recordGenParticles(
  edm::Event const& iEvent,
  std::vector<pat::Muon const*>* filledMuons,
  std::vector<pat::Electron const*>* filledElectrons,
  std::vector<pat::Photon const*>* filledPhotons,
  std::vector<reco::GenParticle const*>* filledGenParts,
  std::vector<pat::PackedGenParticle const*>* filledPackedGenParts
){
  std::string const& colName = CMS3Ntuplizer::colName_genparticles;

  edm::Handle<reco::GenParticleCollection> prunedGenParticlesHandle;
  iEvent.getByToken(prunedGenParticlesToken, prunedGenParticlesHandle);
  if (!prunedGenParticlesHandle.isValid()) throw cms::Exception("CMS3Ntuplizer::recordGenParticles: Error getting the pruned gen. particles from the event...");
  std::vector<reco::GenParticle> const* prunedGenParticles = prunedGenParticlesHandle.product();

  {
    std::vector<reco::GenParticle const*> genFSMuons, genFSElectrons, genFSPhotons;
    genFSMuons.reserve(prunedGenParticles->size());
    genFSElectrons.reserve(prunedGenParticles->size());
    genFSPhotons.reserve(prunedGenParticles->size());
    for (auto const& part:(*prunedGenParticles)){
      if (part.status()!=1) continue;
      unsigned int id = std::abs(part.pdgId());
      if (id==13) genFSMuons.push_back(&part);
      else if (id==11) genFSElectrons.push_back(&part);
      else if (id==22) genFSPhotons.push_back(&part);
    }
    if (filledMuons){
      std::unordered_map<pat::Muon const*, reco::GenParticle const*> reco_gen_map;
      CMS3ObjectHelpers::matchParticles(
        CMS3ObjectHelpers::kMatchBy_DeltaR, 0.2,
        filledMuons->begin(), filledMuons->end(),
        genFSMuons.begin(), genFSMuons.end(),
        reco_gen_map
      );
      size_t n_objects_reco = filledMuons->size();
      MAKE_VECTOR_WITH_RESERVE(bool, is_genMatched, n_objects_reco);
      MAKE_VECTOR_WITH_RESERVE(bool, is_genMatched_prompt, n_objects_reco);
      auto it_end = reco_gen_map.end();
      for (auto const& part:(*filledMuons)){
        bool val_is_genMatched = false;
        bool val_is_genMatched_prompt = false;
        auto it_tmp = reco_gen_map.find(part);
        if (it_tmp!=it_end && it_tmp->second){
          val_is_genMatched = true;
          val_is_genMatched_prompt = it_tmp->second->isPromptFinalState();
        }
        is_genMatched.push_back(val_is_genMatched);
        is_genMatched_prompt.push_back(val_is_genMatched_prompt);
      }
      PUSH_VECTOR_WITH_NAME(CMS3Ntuplizer::colName_muons, is_genMatched);
      PUSH_VECTOR_WITH_NAME(CMS3Ntuplizer::colName_muons, is_genMatched_prompt);
    }
    if (filledElectrons){
      std::unordered_map<pat::Electron const*, reco::GenParticle const*> reco_gen_map;
      CMS3ObjectHelpers::matchParticles(
        CMS3ObjectHelpers::kMatchBy_DeltaR, 0.2,
        filledElectrons->begin(), filledElectrons->end(),
        genFSElectrons.begin(), genFSElectrons.end(),
        reco_gen_map
      );
      size_t n_objects_reco = filledElectrons->size();
      MAKE_VECTOR_WITH_RESERVE(bool, is_genMatched, n_objects_reco);
      MAKE_VECTOR_WITH_RESERVE(bool, is_genMatched_prompt, n_objects_reco);
      auto it_end = reco_gen_map.end();
      for (auto const& part:(*filledElectrons)){
        bool val_is_genMatched = false;
        bool val_is_genMatched_prompt = false;
        auto it_tmp = reco_gen_map.find(part);
        if (it_tmp!=it_end && it_tmp->second){
          val_is_genMatched = true;
          val_is_genMatched_prompt = it_tmp->second->isPromptFinalState();
        }
        is_genMatched.push_back(val_is_genMatched);
        is_genMatched_prompt.push_back(val_is_genMatched_prompt);
      }
      PUSH_VECTOR_WITH_NAME(CMS3Ntuplizer::colName_electrons, is_genMatched);
      PUSH_VECTOR_WITH_NAME(CMS3Ntuplizer::colName_electrons, is_genMatched_prompt);
    }
    if (filledPhotons){
      std::unordered_map<pat::Photon const*, reco::GenParticle const*> reco_gen_map;
      CMS3ObjectHelpers::matchParticles(
        CMS3ObjectHelpers::kMatchBy_DeltaR, 0.2,
        filledPhotons->begin(), filledPhotons->end(),
        genFSPhotons.begin(), genFSPhotons.end(),
        reco_gen_map
      );
      size_t n_objects_reco = filledPhotons->size();
      MAKE_VECTOR_WITH_RESERVE(bool, is_genMatched, n_objects_reco);
      MAKE_VECTOR_WITH_RESERVE(bool, is_genMatched_prompt, n_objects_reco);
      auto it_end = reco_gen_map.end();
      for (auto const& part:(*filledPhotons)){
        bool val_is_genMatched = false;
        bool val_is_genMatched_prompt = false;
        auto it_tmp = reco_gen_map.find(part);
        if (it_tmp!=it_end && it_tmp->second){
          val_is_genMatched = true;
          val_is_genMatched_prompt = it_tmp->second->isPromptFinalState();
        }
        is_genMatched.push_back(val_is_genMatched);
        is_genMatched_prompt.push_back(val_is_genMatched_prompt);
      }
      PUSH_VECTOR_WITH_NAME(CMS3Ntuplizer::colName_photons, is_genMatched);
      PUSH_VECTOR_WITH_NAME(CMS3Ntuplizer::colName_photons, is_genMatched_prompt);
    }
  }

  edm::Handle<pat::PackedGenParticleCollection> packedGenParticlesHandle;
  iEvent.getByToken(packedGenParticlesToken, packedGenParticlesHandle);
  if (!prunedGenParticlesHandle.isValid()) throw cms::Exception("CMS3Ntuplizer::recordGenParticles: Error getting the packed gen. particles from the event...");
  std::vector<pat::PackedGenParticle> const* packedGenParticles = packedGenParticlesHandle.product();

  // Make a collection of all unique reco::GenParticle pointers
  std::vector<reco::GenParticle const*> allGenParticles; allGenParticles.reserve(prunedGenParticles->size() + packedGenParticles->size());
  // Fill with the pruned collection first
  for (reco::GenParticle const& part:(*prunedGenParticles)){
    int st = part.status();
    unsigned int id = std::abs(part.pdgId());
    if (
      (this->keepGenParticles==kReducedFinalStatesAndHardProcesses && !part.isHardProcess() && st!=1)
      ||
      ((this->keepGenParticles==kAllFinalStates || this->keepGenParticles==kReducedFinalStates || this->keepGenParticles==kPromptFinalStatePhotons) && st!=1)
      ) continue;
    else if (
      this->keepGenParticles==kReducedFinalStatesAndHardProcessFinalStates
      && (
        (st!=1 && !(part.isHardProcess() && st!=21 && st!=22))
        ||
        (st==1 && !(part.isHardProcess() || (id>=11 && id<=16) || id==22))
        )
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
    else if (
      this->keepGenParticles==kPromptFinalStatePhotons
      &&
      // Only photons for this flag
      !(id==22 && part.isPromptFinalState())
      ) continue;
    else if (
      (this->keepGenParticles==kHardProcesses && !part.isHardProcess())
      ||
      (this->keepGenParticles==kHardProcessFinalStates && (!part.isHardProcess() || st==21 || st==22))
      ) continue;
    allGenParticles.push_back(&part);
  }

  // Get the packed gen. particles unique from the pruned collection
  // Adapted from GeneratorInterface/RivetInterface/plugins/MergedGenParticleProducer.cc
  std::vector<pat::PackedGenParticle const*> uniquePackedGenParticles; uniquePackedGenParticles.reserve(packedGenParticles->size());
  for (pat::PackedGenParticle const& packedGenParticle:(*packedGenParticlesHandle)){
    int st = packedGenParticle.status();
    int id_signed = packedGenParticle.pdgId();
    unsigned int id = std::abs(id_signed);

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
        (
          this->keepGenParticles==kAllFinalStates
          ||
          this->keepGenParticles==kReducedFinalStates
          ||
          this->keepGenParticles==kReducedFinalStatesAndHardProcesses || this->keepGenParticles==kReducedFinalStatesAndHardProcessFinalStates
          )
        &&
        st!=1
        ) continue;
      else if (this->keepGenParticles==kPromptFinalStatePhotons) continue; // Disable photons from packed candidates since they are not prompt photons
      else if (this->keepGenParticles==kHardProcesses || this->keepGenParticles==kHardProcessFinalStates) continue; // Disable packed candidates since they are not hard process particles
      else if (
        (this->keepGenParticles==kReducedFinalStates || this->keepGenParticles==kReducedFinalStatesAndHardProcesses || this->keepGenParticles==kReducedFinalStatesAndHardProcessFinalStates)
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

  MAKE_VECTOR_WITH_RESERVE(cms3_id_t, id, n_objects);
  MAKE_VECTOR_WITH_RESERVE(cms3_genstatus_t, status, n_objects);

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

  MAKE_VECTOR_WITH_RESERVE(int, mom0_index, n_objects);
  MAKE_VECTOR_WITH_RESERVE(int, mom1_index, n_objects);

  // Record all reco::GenParticle objects
  for (reco::GenParticle const* obj:allGenParticles){
    if (filledGenParts && obj->status()==1) filledGenParts->push_back(obj);

    pt.push_back(obj->pt());
    eta.push_back(obj->eta());
    phi.push_back(obj->phi());
    mass.push_back(obj->mass());

    if (obj->pdgId() < std::numeric_limits<cms3_id_t>::min() || obj->pdgId() > std::numeric_limits<cms3_id_t>::max()){
      cms::Exception numexcpt("NumericLimits");
      numexcpt << "Particle id " << obj->pdgId() << " exceeds numerical bounds.";
      throw numexcpt;
    }
    if (obj->status() < std::numeric_limits<cms3_genstatus_t>::min() || obj->status() > std::numeric_limits<cms3_genstatus_t>::max()){
      cms::Exception numexcpt("NumericLimits");
      numexcpt << "Particle status " << obj->status() << " exceeds numerical bounds.";
      throw numexcpt;
    }

    id.push_back(obj->pdgId());
    status.push_back(obj->status());

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

    id.push_back(obj->pdgId());
    status.push_back(obj->status());

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

  PUSH_VECTOR_WITH_NAME(colName, id);
  PUSH_VECTOR_WITH_NAME(colName, status);

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
  std::string const& colName = CMS3Ntuplizer::colName_muons;

  edm::Handle< edm::View<pat::Muon> > muonsHandle;
  iEvent.getByToken(muonsToken, muonsHandle);
  if (!muonsHandle.isValid()) throw cms::Exception("CMS3Ntuplizer::fillMuons: Error getting the muon collection from the event...");
  size_t n_objects = muonsHandle->size();

  if (filledObjects) filledObjects->reserve(n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, pt, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, eta, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, phi, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, mass, n_objects);

  MAKE_VECTOR_WITH_RESERVE(cms3_charge_t, charge, n_objects);

  MAKE_VECTOR_WITH_RESERVE(cms3_muon_pogselectorbits_t, POG_selector_bits, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, pfIso03_comb_nofsr, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, pfIso03_sum_charged_nofsr, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, pfIso03_sum_neutral_nofsr, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, pfIso03_sum_neutral_EAcorr_nofsr, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, pfIso04_comb_nofsr, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, pfIso04_sum_charged_nofsr, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, pfIso04_sum_neutral_nofsr, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, pfIso04_sum_neutral_EAcorr_nofsr, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, miniIso_sum_charged_nofsr, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, miniIso_sum_neutral_nofsr, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, miniIso_comb_nofsr, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, miniIso_comb_nofsr_uncorrected, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, trkIso03_trackerSumPt, n_objects);

  MAKE_VECTOR_WITH_RESERVE(bool, pass_tightCharge, n_objects);
  MAKE_VECTOR_WITH_RESERVE(bool, is_probeForTnP, n_objects);
  MAKE_VECTOR_WITH_RESERVE(bool, is_probeForTnP_STA, n_objects);

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
  MAKE_VECTOR_WITH_RESERVE(bool, pass_muon_timing, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, pull_dxdz_noArb_DT, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, pull_dxdz_noArb_CSC, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, dxy_bestTrack_firstPV, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, dz_bestTrack_firstPV, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, SIP3D, n_objects);

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

    PUSH_USERFLOAT_INTO_VECTOR(pfIso03_comb_nofsr);
    PUSH_USERFLOAT_INTO_VECTOR(pfIso03_sum_charged_nofsr);
    PUSH_USERFLOAT_INTO_VECTOR(pfIso03_sum_neutral_nofsr);
    PUSH_USERFLOAT_INTO_VECTOR(pfIso03_sum_neutral_EAcorr_nofsr);

    PUSH_USERFLOAT_INTO_VECTOR(pfIso04_comb_nofsr);
    PUSH_USERFLOAT_INTO_VECTOR(pfIso04_sum_charged_nofsr);
    PUSH_USERFLOAT_INTO_VECTOR(pfIso04_sum_neutral_nofsr);
    PUSH_USERFLOAT_INTO_VECTOR(pfIso04_sum_neutral_EAcorr_nofsr);

    PUSH_USERFLOAT_INTO_VECTOR(miniIso_sum_charged_nofsr);
    PUSH_USERFLOAT_INTO_VECTOR(miniIso_sum_neutral_nofsr);
    PUSH_USERFLOAT_INTO_VECTOR(miniIso_comb_nofsr);
    PUSH_USERFLOAT_INTO_VECTOR(miniIso_comb_nofsr_uncorrected);

    PUSH_USERFLOAT_INTO_VECTOR(trkIso03_trackerSumPt);

    PUSH_USERINT_INTO_VECTOR(pass_tightCharge);
    PUSH_USERINT_INTO_VECTOR(is_probeForTnP);
    PUSH_USERINT_INTO_VECTOR(is_probeForTnP_STA);

    if (keepMuonTimingInfo){
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
    }
    else{
      PUSH_USERINT_INTO_VECTOR(pass_muon_timing);
    }

    if (keepMuonPullInfo){
      PUSH_USERFLOAT_INTO_VECTOR(pull_dxdz_noArb_DT);
      PUSH_USERFLOAT_INTO_VECTOR(pull_dxdz_noArb_CSC);
    }

    PUSH_USERFLOAT_INTO_VECTOR(dxy_bestTrack_firstPV);
    PUSH_USERFLOAT_INTO_VECTOR(dz_bestTrack_firstPV);
    PUSH_USERFLOAT_INTO_VECTOR(SIP3D);

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

  PUSH_VECTOR_WITH_NAME(colName, pfIso03_comb_nofsr);
  PUSH_VECTOR_WITH_NAME(colName, pfIso03_sum_charged_nofsr);
  PUSH_VECTOR_WITH_NAME(colName, pfIso03_sum_neutral_nofsr);
  PUSH_VECTOR_WITH_NAME(colName, pfIso03_sum_neutral_EAcorr_nofsr);

  PUSH_VECTOR_WITH_NAME(colName, pfIso04_comb_nofsr);
  PUSH_VECTOR_WITH_NAME(colName, pfIso04_sum_charged_nofsr);
  PUSH_VECTOR_WITH_NAME(colName, pfIso04_sum_neutral_nofsr);
  PUSH_VECTOR_WITH_NAME(colName, pfIso04_sum_neutral_EAcorr_nofsr);

  PUSH_VECTOR_WITH_NAME(colName, miniIso_sum_charged_nofsr);
  PUSH_VECTOR_WITH_NAME(colName, miniIso_sum_neutral_nofsr);
  PUSH_VECTOR_WITH_NAME(colName, miniIso_comb_nofsr);
  PUSH_VECTOR_WITH_NAME(colName, miniIso_comb_nofsr_uncorrected);

  PUSH_VECTOR_WITH_NAME(colName, trkIso03_trackerSumPt);

  PUSH_VECTOR_WITH_NAME(colName, pass_tightCharge);
  PUSH_VECTOR_WITH_NAME(colName, is_probeForTnP);
  PUSH_VECTOR_WITH_NAME(colName, is_probeForTnP_STA);

  if (keepMuonTimingInfo){
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
  }
  else{
    PUSH_VECTOR_WITH_NAME(colName, pass_muon_timing);
  }

  if (keepMuonPullInfo){
    PUSH_VECTOR_WITH_NAME(colName, pull_dxdz_noArb_DT);
    PUSH_VECTOR_WITH_NAME(colName, pull_dxdz_noArb_CSC);
  }

  PUSH_VECTOR_WITH_NAME(colName, dxy_bestTrack_firstPV);
  PUSH_VECTOR_WITH_NAME(colName, dz_bestTrack_firstPV);
  PUSH_VECTOR_WITH_NAME(colName, SIP3D);

  PUSH_VECTOR_WITH_NAME(colName, scale_smear_pt_corr);
  PUSH_VECTOR_WITH_NAME(colName, scale_smear_pt_corr_scale_totalUp);
  PUSH_VECTOR_WITH_NAME(colName, scale_smear_pt_corr_scale_totalDn);
  PUSH_VECTOR_WITH_NAME(colName, scale_smear_pt_corr_smear_totalUp);
  PUSH_VECTOR_WITH_NAME(colName, scale_smear_pt_corr_smear_totalDn);

  return n_skimmed_objects;
}
size_t CMS3Ntuplizer::fillElectrons(edm::Event const& iEvent, std::vector<pat::Electron const*>* filledObjects){
  std::string const& colName = CMS3Ntuplizer::colName_electrons;

  edm::Handle< edm::View<pat::Electron> > electronsHandle;
  iEvent.getByToken(electronsToken, electronsHandle);
  if (!electronsHandle.isValid()) throw cms::Exception("CMS3Ntuplizer::fillElectrons: Error getting the electron collection from the event...");
  size_t n_objects = electronsHandle->size();

  if (filledObjects) filledObjects->reserve(n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, pt, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, eta, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, phi, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, mass, n_objects);

  MAKE_VECTOR_WITH_RESERVE(cms3_charge_t, charge, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, etaSC, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, ecalEnergy, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, sinTheta_SC_pos, n_objects);

  // Has no convention correspondence in nanoAOD
  MAKE_VECTOR_WITH_RESERVE(float, scale_smear_corr, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, scale_smear_corr_scale_totalUp, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, scale_smear_corr_scale_totalDn, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, scale_smear_corr_smear_totalUp, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, scale_smear_corr_smear_totalDn, n_objects);

  MAKE_VECTOR_WITH_RESERVE(cms3_electron_charge_consistency_bits_t, charge_consistency_bits, n_objects);

  MAKE_VECTOR_WITH_RESERVE(bool, conv_vtx_flag, n_objects);
  MAKE_VECTOR_WITH_RESERVE(cms3_electron_missinghits_t, n_missing_inner_hits, n_objects);
  MAKE_VECTOR_WITH_RESERVE(cms3_electron_missinghits_t, n_all_missing_inner_hits, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, id_MVA_Fall17V2_Iso_Val, n_objects);
  MAKE_VECTOR_WITH_RESERVE(cms3_electron_mvacat_t, id_MVA_Fall17V2_Iso_Cat, n_objects);
  MAKE_VECTOR_WITH_RESERVE(bool, id_MVA_Fall17V2_Iso_pass_wpLoose, n_objects);
  MAKE_VECTOR_WITH_RESERVE(bool, id_MVA_Fall17V2_Iso_pass_wp90, n_objects);
  MAKE_VECTOR_WITH_RESERVE(bool, id_MVA_Fall17V2_Iso_pass_wp80, n_objects);
  MAKE_VECTOR_WITH_RESERVE(bool, id_MVA_Fall17V2_Iso_pass_wpHZZ, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, id_MVA_Fall17V2_NoIso_Val, n_objects);
  MAKE_VECTOR_WITH_RESERVE(cms3_electron_mvacat_t, id_MVA_Fall17V2_NoIso_Cat, n_objects);
  MAKE_VECTOR_WITH_RESERVE(bool, id_MVA_Fall17V2_NoIso_pass_wpLoose, n_objects);
  MAKE_VECTOR_WITH_RESERVE(bool, id_MVA_Fall17V2_NoIso_pass_wp90, n_objects);
  MAKE_VECTOR_WITH_RESERVE(bool, id_MVA_Fall17V2_NoIso_pass_wp80, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, id_MVA_HZZRun2Legacy_Iso_Val, n_objects);
  MAKE_VECTOR_WITH_RESERVE(cms3_electron_mvacat_t, id_MVA_HZZRun2Legacy_Iso_Cat, n_objects);
  MAKE_VECTOR_WITH_RESERVE(bool, id_MVA_HZZRun2Legacy_Iso_pass_wpHZZ, n_objects);

  MAKE_VECTOR_WITH_RESERVE(cms3_electron_cutbasedbits_t, id_cutBased_Fall17V2_Veto_Bits, n_objects);
  MAKE_VECTOR_WITH_RESERVE(cms3_electron_cutbasedbits_t, id_cutBased_Fall17V2_Loose_Bits, n_objects);
  MAKE_VECTOR_WITH_RESERVE(cms3_electron_cutbasedbits_t, id_cutBased_Fall17V2_Medium_Bits, n_objects);
  MAKE_VECTOR_WITH_RESERVE(cms3_electron_cutbasedbits_t, id_cutBased_Fall17V2_Tight_Bits, n_objects);
  MAKE_VECTOR_WITH_RESERVE(cms3_electron_cutbasedbits_triggeremulation_t, id_cutBased_triggerEmulationV1_Bits, n_objects);
  MAKE_VECTOR_WITH_RESERVE(cms3_electron_cutbasedbits_triggeremulation_t, id_cutBased_triggerEmulationV2_Bits, n_objects);

  /*
  MAKE_VECTOR_WITH_RESERVE(cms3_electron_cutbasedbits_t, id_cutBased_Fall17V1_Veto_Bits, n_objects);
  MAKE_VECTOR_WITH_RESERVE(cms3_electron_cutbasedbits_t, id_cutBased_Fall17V1_Loose_Bits, n_objects);
  MAKE_VECTOR_WITH_RESERVE(cms3_electron_cutbasedbits_t, id_cutBased_Fall17V1_Medium_Bits, n_objects);
  MAKE_VECTOR_WITH_RESERVE(cms3_electron_cutbasedbits_t, id_cutBased_Fall17V1_Tight_Bits, n_objects);
  */

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

  MAKE_VECTOR_WITH_RESERVE(float, dxy_firstPV, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, dz_firstPV, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, SIP3D, n_objects);

  /*
  MAKE_VECTOR_WITH_RESERVE(unsigned int, n_associated_pfcands, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, associated_pfcands_sum_sc_pt, n_objects);
  */

  // Fiduciality and type masks use bits from interface/EgammaFiduciality.h.
  // These are needed to define gap regions (fid_mask ISGAP bit) or being ECAL-driven (type mask ISECALDRIVEN and ISCUTPRESELECTED bits)
  // Do not comment these out!
  MAKE_VECTOR_WITH_RESERVE(cms3_egamma_fid_type_mask_t, fid_mask, n_objects);
  MAKE_VECTOR_WITH_RESERVE(cms3_egamma_fid_type_mask_t, type_mask, n_objects);

  size_t n_skimmed_objects=0;
  for (edm::View<pat::Electron>::const_iterator obj = electronsHandle->begin(); obj != electronsHandle->end(); obj++){
    if (!ElectronSelectionHelpers::testSkimElectron(*obj, this->year)) continue;

    // Core particle quantities
    // Uncorrected p4
    pt.push_back(obj->pt());
    eta.push_back(obj->eta());
    phi.push_back(obj->phi());
    mass.push_back(obj->mass());

    // Charge: Can obtain pdgId from this, so no need to record pdgId again
    PUSH_USERINT_INTO_VECTOR(charge);
    PUSH_USERFLOAT_INTO_VECTOR(etaSC);

    PUSH_USERFLOAT_INTO_VECTOR(ecalEnergy);
    sinTheta_SC_pos.push_back(!obj->superCluster().isNull() ? std::sin(obj->superCluster()->position().theta()) : -99.);

    // Scale and smear
    // Nominal value: Needs to multiply the uncorrected p4 at analysis level
    PUSH_USERFLOAT_INTO_VECTOR(scale_smear_corr);
    // Uncertainties: Only store total up/dn for the moment
    PUSH_USERFLOAT_INTO_VECTOR(scale_smear_corr_scale_totalUp);
    PUSH_USERFLOAT_INTO_VECTOR(scale_smear_corr_scale_totalDn);
    PUSH_USERFLOAT_INTO_VECTOR(scale_smear_corr_smear_totalUp);
    PUSH_USERFLOAT_INTO_VECTOR(scale_smear_corr_smear_totalDn);

    PUSH_USERINT_INTO_VECTOR(charge_consistency_bits);

    PUSH_USERINT_INTO_VECTOR(conv_vtx_flag);
    PUSH_USERINT_INTO_VECTOR(n_missing_inner_hits);
    PUSH_USERINT_INTO_VECTOR(n_all_missing_inner_hits);

    // Id variables
    // Fall17V2_Iso MVA id
    if (keepElectronMVAInfo){
      PUSH_USERFLOAT_INTO_VECTOR(id_MVA_Fall17V2_Iso_Val);
      PUSH_USERINT_INTO_VECTOR(id_MVA_Fall17V2_Iso_Cat);
    }
    PUSH_USERINT_INTO_VECTOR(id_MVA_Fall17V2_Iso_pass_wpLoose);
    PUSH_USERINT_INTO_VECTOR(id_MVA_Fall17V2_Iso_pass_wp90);
    PUSH_USERINT_INTO_VECTOR(id_MVA_Fall17V2_Iso_pass_wp80);
    PUSH_USERINT_INTO_VECTOR(id_MVA_Fall17V2_Iso_pass_wpHZZ);

    // Fall17V2_NoIso MVA id
    if (keepElectronMVAInfo){
      PUSH_USERFLOAT_INTO_VECTOR(id_MVA_Fall17V2_NoIso_Val);
      PUSH_USERINT_INTO_VECTOR(id_MVA_Fall17V2_NoIso_Cat);
    }
    PUSH_USERINT_INTO_VECTOR(id_MVA_Fall17V2_NoIso_pass_wpLoose);
    PUSH_USERINT_INTO_VECTOR(id_MVA_Fall17V2_NoIso_pass_wp90);
    PUSH_USERINT_INTO_VECTOR(id_MVA_Fall17V2_NoIso_pass_wp80);

    // HZZ Run 2 legacy electron MVA id
    if (keepElectronMVAInfo){
      PUSH_USERFLOAT_INTO_VECTOR(id_MVA_HZZRun2Legacy_Iso_Val);
      PUSH_USERINT_INTO_VECTOR(id_MVA_HZZRun2Legacy_Iso_Cat);
    }
    PUSH_USERINT_INTO_VECTOR(id_MVA_HZZRun2Legacy_Iso_pass_wpHZZ);

    // Fall17V2 cut-based ids
    PUSH_USERINT_INTO_VECTOR(id_cutBased_Fall17V2_Veto_Bits);
    PUSH_USERINT_INTO_VECTOR(id_cutBased_Fall17V2_Loose_Bits);
    PUSH_USERINT_INTO_VECTOR(id_cutBased_Fall17V2_Medium_Bits);
    PUSH_USERINT_INTO_VECTOR(id_cutBased_Fall17V2_Tight_Bits);

    PUSH_USERINT_INTO_VECTOR(id_cutBased_triggerEmulationV1_Bits);
    PUSH_USERINT_INTO_VECTOR(id_cutBased_triggerEmulationV2_Bits);

    /*
    // Fall17V1 cut-based ids
    PUSH_USERINT_INTO_VECTOR(id_cutBased_Fall17V1_Veto_Bits);
    PUSH_USERINT_INTO_VECTOR(id_cutBased_Fall17V1_Loose_Bits);
    PUSH_USERINT_INTO_VECTOR(id_cutBased_Fall17V1_Medium_Bits);
    PUSH_USERINT_INTO_VECTOR(id_cutBased_Fall17V1_Tight_Bits);
    */

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

    PUSH_USERFLOAT_INTO_VECTOR(dxy_firstPV);
    PUSH_USERFLOAT_INTO_VECTOR(dz_firstPV);
    PUSH_USERFLOAT_INTO_VECTOR(SIP3D);

    /*
    // PF candidate properties
    PUSH_USERINT_INTO_VECTOR(n_associated_pfcands);
    PUSH_USERFLOAT_INTO_VECTOR(associated_pfcands_sum_sc_pt);
    */

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

  PUSH_VECTOR_WITH_NAME(colName, ecalEnergy);
  PUSH_VECTOR_WITH_NAME(colName, sinTheta_SC_pos);

  PUSH_VECTOR_WITH_NAME(colName, charge_consistency_bits);

  PUSH_VECTOR_WITH_NAME(colName, conv_vtx_flag);
  PUSH_VECTOR_WITH_NAME(colName, n_missing_inner_hits);
  PUSH_VECTOR_WITH_NAME(colName, n_all_missing_inner_hits);

  // Has no convention correspondence in nanoAOD
  PUSH_VECTOR_WITH_NAME(colName, scale_smear_corr);
  PUSH_VECTOR_WITH_NAME(colName, scale_smear_corr_scale_totalUp);
  PUSH_VECTOR_WITH_NAME(colName, scale_smear_corr_scale_totalDn);
  PUSH_VECTOR_WITH_NAME(colName, scale_smear_corr_smear_totalUp);
  PUSH_VECTOR_WITH_NAME(colName, scale_smear_corr_smear_totalDn);

  if (keepElectronMVAInfo){
    PUSH_VECTOR_WITH_NAME(colName, id_MVA_Fall17V2_Iso_Val);
    PUSH_VECTOR_WITH_NAME(colName, id_MVA_Fall17V2_Iso_Cat);
  }
  PUSH_VECTOR_WITH_NAME(colName, id_MVA_Fall17V2_Iso_pass_wpLoose);
  PUSH_VECTOR_WITH_NAME(colName, id_MVA_Fall17V2_Iso_pass_wp90);
  PUSH_VECTOR_WITH_NAME(colName, id_MVA_Fall17V2_Iso_pass_wp80);
  PUSH_VECTOR_WITH_NAME(colName, id_MVA_Fall17V2_Iso_pass_wpHZZ);

  if (keepElectronMVAInfo){
    PUSH_VECTOR_WITH_NAME(colName, id_MVA_Fall17V2_NoIso_Val);
    PUSH_VECTOR_WITH_NAME(colName, id_MVA_Fall17V2_NoIso_Cat);
  }
  PUSH_VECTOR_WITH_NAME(colName, id_MVA_Fall17V2_NoIso_pass_wpLoose);
  PUSH_VECTOR_WITH_NAME(colName, id_MVA_Fall17V2_NoIso_pass_wp90);
  PUSH_VECTOR_WITH_NAME(colName, id_MVA_Fall17V2_NoIso_pass_wp80);

  if (keepElectronMVAInfo){
    PUSH_VECTOR_WITH_NAME(colName, id_MVA_HZZRun2Legacy_Iso_Val);
    PUSH_VECTOR_WITH_NAME(colName, id_MVA_HZZRun2Legacy_Iso_Cat);
  }
  PUSH_VECTOR_WITH_NAME(colName, id_MVA_HZZRun2Legacy_Iso_pass_wpHZZ);

  PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V2_Veto_Bits);
  PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V2_Loose_Bits);
  PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V2_Medium_Bits);
  PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V2_Tight_Bits);

  PUSH_VECTOR_WITH_NAME(colName, id_cutBased_triggerEmulationV1_Bits);
  PUSH_VECTOR_WITH_NAME(colName, id_cutBased_triggerEmulationV2_Bits);

  /*
  PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V1_Veto_Bits);
  PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V1_Loose_Bits);
  PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V1_Medium_Bits);
  PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V1_Tight_Bits);
  */

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

  PUSH_VECTOR_WITH_NAME(colName, dxy_firstPV);
  PUSH_VECTOR_WITH_NAME(colName, dz_firstPV);
  PUSH_VECTOR_WITH_NAME(colName, SIP3D);

  /*
  PUSH_VECTOR_WITH_NAME(colName, n_associated_pfcands);
  PUSH_VECTOR_WITH_NAME(colName, associated_pfcands_sum_sc_pt);
  */

  PUSH_VECTOR_WITH_NAME(colName, fid_mask);
  PUSH_VECTOR_WITH_NAME(colName, type_mask);

  return n_skimmed_objects;
}
size_t CMS3Ntuplizer::fillFSRCandidates(
  edm::Event const& iEvent,
  std::vector<pat::Muon const*> const* filledMuons, std::vector<pat::Electron const*> const* filledElectrons,
  std::vector<FSRCandidateInfo>* filledObjects
){
  std::string const& colNameFSR = CMS3Ntuplizer::colName_fsrcands;

  edm::Handle< edm::View<pat::PackedCandidate> > pfcandsHandle;
  iEvent.getByToken(pfcandsToken, pfcandsHandle);
  if (!pfcandsHandle.isValid()) throw cms::Exception("CMS3Ntuplizer::fillFSRCandidates: Error getting the PF candidate collection from the event...");
  size_t n_objects = pfcandsHandle->size();
  size_t n_skimmed_objects=0;

  // FSR preselection
  std::vector<FSRCandidateInfo> preselectedFSRCandidates; preselectedFSRCandidates.reserve(n_objects);
  edm::View<pat::PackedCandidate>::const_iterator it_pfcands_begin = pfcandsHandle->begin();
  edm::View<pat::PackedCandidate>::const_iterator it_pfcands_end = pfcandsHandle->end();
  for (edm::View<pat::PackedCandidate>::const_iterator obj = it_pfcands_begin; obj != it_pfcands_end; obj++){
    if (obj->pdgId()!=22) continue; // Check only photons

    if (!FSRSelectionHelpers::testSkimFSR_PtEta(*obj, this->year)) continue;

    double fsrIso = FSRSelectionHelpers::fsrIso(*obj, this->year, it_pfcands_begin, it_pfcands_end);
    if (!FSRSelectionHelpers::testSkimFSR_Iso(*obj, this->year, fsrIso)) continue;

    FSRCandidateInfo fsrInfo;
    fsrInfo.obj = &(*obj);
    fsrInfo.fsrIso = fsrIso;
    if (filledElectrons){ for (auto const& electron:(*filledElectrons)){ if (FSRSelectionHelpers::testSCVeto(&(*obj), electron)){ fsrInfo.veto_electron_list.push_back(electron); } } }

    preselectedFSRCandidates.emplace_back(fsrInfo);
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
  // Do the rel. min. dR check between leptons and FSR candidates and make the final, writable collection
  std::vector<FSRCandidateInfo*> writableFSRCandidates;
  for (auto& it:fsrcand_lepton_map){
    FSRCandidateInfo const* fsrcand_const = it.first;
    FSRCandidateInfo* fsrcand = nullptr;
    for (auto& cand:preselectedFSRCandidates){ if (&cand == fsrcand_const) fsrcand = &cand; }
    for (auto const& lepton:it.second){
      // Skip leptons that do not satisfy min. rel. dR requirements (dR/(pTfsr^2)<0.012)
      double const dR_fsr_lepton = reco::deltaR(fsrcand->p4(), lepton->p4());
      if (!FSRSelectionHelpers::testSkimFSR_MinDeltaR(*(fsrcand->obj), this->year, dR_fsr_lepton)) continue;

      pat::Muon const* muon = dynamic_cast<pat::Muon const*>(lepton);
      pat::Electron const* electron = dynamic_cast<pat::Electron const*>(lepton);
      if (muon){
        fsrcand->matched_muon_list.push_back(muon);
      }
      else if (electron){
        if (HelperFunctions::checkListVariable(fsrcand->veto_electron_list, electron)) continue;
        fsrcand->matched_electron_list.push_back(electron);
      }
    }
    if (!fsrcand->matched_muon_list.empty() || !fsrcand->matched_electron_list.empty()) writableFSRCandidates.push_back(fsrcand);
  }

  if (filledObjects) filledObjects->reserve(writableFSRCandidates.size());

  // Fill FSR candidates
  MAKE_VECTOR_WITH_RESERVE(float, pt, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, eta, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, phi, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, mass, n_objects);
  MAKE_VECTOR_WITH_RESERVE(std::vector<unsigned int>, fsrMatch_muon_index_list, n_objects);
  MAKE_VECTOR_WITH_RESERVE(std::vector<unsigned int>, fsrMatch_electron_index_list, n_objects);
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

    // Do not fill photonVeto_index_list here. This list is supposed to be filled in CMS3Ntuplizer::fillPhotons because those photons also need to be recorded even if they don't pass the photon skim selection.

    if (filledObjects) filledObjects->emplace_back(*obj);

    n_skimmed_objects++;
  }
  PUSH_VECTOR_WITH_NAME(colNameFSR, pt);
  PUSH_VECTOR_WITH_NAME(colNameFSR, eta);
  PUSH_VECTOR_WITH_NAME(colNameFSR, phi);
  PUSH_VECTOR_WITH_NAME(colNameFSR, mass);
  PUSH_VECTOR_WITH_NAME(colNameFSR, fsrMatch_muon_index_list);
  PUSH_VECTOR_WITH_NAME(colNameFSR, fsrMatch_electron_index_list);

  return n_skimmed_objects;
}
size_t CMS3Ntuplizer::fillPhotons(edm::Event const& iEvent, std::vector<FSRCandidateInfo>* filledFSRCandidates, std::vector<pat::Photon const*>* filledObjects){
  std::string const& colName = CMS3Ntuplizer::colName_photons;

  edm::Handle< edm::View<pat::Photon> > photonsHandle;
  iEvent.getByToken(photonsToken, photonsHandle);
  if (!photonsHandle.isValid()) throw cms::Exception("CMS3Ntuplizer::fillPhotons: Error getting the photon collection from the event...");
  size_t n_objects = photonsHandle->size();

  std::vector<pat::Photon const*> filledPhotons;
  filledPhotons.reserve(n_objects);

  // Begin filling the objects
  MAKE_VECTOR_WITH_RESERVE(float, pt, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, eta, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, phi, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, mass, n_objects);

  // etaSC
  MAKE_VECTOR_WITH_RESERVE(float, etaSC, n_objects);

  // Has no convention correspondence in nanoAOD
  MAKE_VECTOR_WITH_RESERVE(float, scale_smear_corr, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, scale_smear_corr_scale_totalUp, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, scale_smear_corr_scale_totalDn, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, scale_smear_corr_smear_totalUp, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, scale_smear_corr_smear_totalDn, n_objects);

  // These two are needed for the MVA id
  MAKE_VECTOR_WITH_RESERVE(bool, hasPixelSeed, n_objects);
  MAKE_VECTOR_WITH_RESERVE(bool, passElectronVeto, n_objects);

  /*
  MAKE_VECTOR_WITH_RESERVE(float, id_MVA_Fall17V2_Val, n_objects);
  MAKE_VECTOR_WITH_RESERVE(cms3_photon_mvacat_t, id_MVA_Fall17V2_Cat, n_objects);
  */
  MAKE_VECTOR_WITH_RESERVE(bool, id_MVA_Fall17V2_pass_wp90, n_objects);
  MAKE_VECTOR_WITH_RESERVE(bool, id_MVA_Fall17V2_pass_wp80, n_objects);

  MAKE_VECTOR_WITH_RESERVE(cms3_photon_cutbasedbits_t, id_cutBased_Fall17V2_Loose_Bits, n_objects);
  MAKE_VECTOR_WITH_RESERVE(cms3_photon_cutbasedbits_t, id_cutBased_Fall17V2_Medium_Bits, n_objects);
  MAKE_VECTOR_WITH_RESERVE(cms3_photon_cutbasedbits_t, id_cutBased_Fall17V2_Tight_Bits, n_objects);
  MAKE_VECTOR_WITH_RESERVE(cms3_photon_cutbasedbits_hgg_t, id_cutBased_HGG_Bits, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, full5x5_r9, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, full5x5_sigmaIEtaIEta, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, full5x5_sigmaIPhiIPhi, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, MIPTotalEnergy, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, E4overE1, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, seedTime, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, pfIso_comb, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, pfChargedHadronIso_EAcorr, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, pfNeutralHadronIso_EAcorr, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, pfEMIso_EAcorr, n_objects);

  /*
  MAKE_VECTOR_WITH_RESERVE(unsigned int, n_associated_pfcands, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, associated_pfcands_sum_sc_pt, n_objects);
  */

  MAKE_VECTOR_WITH_RESERVE(cms3_egamma_fid_type_mask_t, fid_mask, n_objects);

  size_t n_skimmed_objects=0;
  for (edm::View<pat::Photon>::const_iterator obj = photonsHandle->begin(); obj != photonsHandle->end(); obj++){
    bool const passStandardSkim = PhotonSelectionHelpers::testSkimPhoton(*obj, this->year);
    bool isFSRSCVetoed = false;
    if (filledFSRCandidates){
      for (auto& fsrInfo:(*filledFSRCandidates)){
        if (FSRSelectionHelpers::testSCVeto(fsrInfo.obj, &(*obj))){
          isFSRSCVetoed = true;
          fsrInfo.veto_photon_list.push_back(&(*obj)); // We can insert the photon here because this photon is going to be recorded now, per the if ... continue line below
        }
      }
    }
    if (!passStandardSkim && !isFSRSCVetoed) continue;

    // Core particle quantities
    // Uncorrected p4
    pt.push_back(obj->pt());
    eta.push_back(obj->eta());
    phi.push_back(obj->phi());
    mass.push_back(obj->mass());

    PUSH_USERFLOAT_INTO_VECTOR(etaSC);

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
    /*
    PUSH_USERFLOAT_INTO_VECTOR(id_MVA_Fall17V2_Val);
    PUSH_USERINT_INTO_VECTOR(id_MVA_Fall17V2_Cat);
    */
    PUSH_USERINT_INTO_VECTOR(id_MVA_Fall17V2_pass_wp90);
    PUSH_USERINT_INTO_VECTOR(id_MVA_Fall17V2_pass_wp80);

    // Fall17V2 cut-based ids
    PUSH_USERINT_INTO_VECTOR(id_cutBased_Fall17V2_Loose_Bits);
    PUSH_USERINT_INTO_VECTOR(id_cutBased_Fall17V2_Medium_Bits);
    PUSH_USERINT_INTO_VECTOR(id_cutBased_Fall17V2_Tight_Bits);
    PUSH_USERINT_INTO_VECTOR(id_cutBased_HGG_Bits);

    PUSH_USERFLOAT_INTO_VECTOR(full5x5_r9);
    PUSH_USERFLOAT_INTO_VECTOR(full5x5_sigmaIEtaIEta);
    PUSH_USERFLOAT_INTO_VECTOR(full5x5_sigmaIPhiIPhi);

    PUSH_USERFLOAT_INTO_VECTOR(MIPTotalEnergy);
    PUSH_USERFLOAT_INTO_VECTOR(E4overE1);
    PUSH_USERFLOAT_INTO_VECTOR(seedTime);

    PUSH_USERFLOAT_INTO_VECTOR(pfIso_comb);
    PUSH_USERFLOAT_INTO_VECTOR(pfChargedHadronIso_EAcorr);
    PUSH_USERFLOAT_INTO_VECTOR(pfNeutralHadronIso_EAcorr);
    PUSH_USERFLOAT_INTO_VECTOR(pfEMIso_EAcorr);

    /*
    PUSH_USERINT_INTO_VECTOR(n_associated_pfcands);
    PUSH_USERFLOAT_INTO_VECTOR(associated_pfcands_sum_sc_pt);
    */

    // Masks
    PUSH_USERINT_INTO_VECTOR(fid_mask);

    filledPhotons.push_back(&(*obj));
    n_skimmed_objects++;
  }

  // Pass collections to the communicator
  PUSH_VECTOR_WITH_NAME(colName, pt);
  PUSH_VECTOR_WITH_NAME(colName, eta);
  PUSH_VECTOR_WITH_NAME(colName, phi);
  PUSH_VECTOR_WITH_NAME(colName, mass);

  PUSH_VECTOR_WITH_NAME(colName, etaSC);

  // Has no convention correspondence in nanoAOD
  PUSH_VECTOR_WITH_NAME(colName, scale_smear_corr);
  PUSH_VECTOR_WITH_NAME(colName, scale_smear_corr_scale_totalUp);
  PUSH_VECTOR_WITH_NAME(colName, scale_smear_corr_scale_totalDn);
  PUSH_VECTOR_WITH_NAME(colName, scale_smear_corr_smear_totalUp);
  PUSH_VECTOR_WITH_NAME(colName, scale_smear_corr_smear_totalDn);

  PUSH_VECTOR_WITH_NAME(colName, hasPixelSeed);
  PUSH_VECTOR_WITH_NAME(colName, passElectronVeto);

  /*
  PUSH_VECTOR_WITH_NAME(colName, id_MVA_Fall17V2_Val);
  PUSH_VECTOR_WITH_NAME(colName, id_MVA_Fall17V2_Cat);
  */
  PUSH_VECTOR_WITH_NAME(colName, id_MVA_Fall17V2_pass_wp90);
  PUSH_VECTOR_WITH_NAME(colName, id_MVA_Fall17V2_pass_wp80);

  PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V2_Loose_Bits);
  PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V2_Medium_Bits);
  PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V2_Tight_Bits);
  PUSH_VECTOR_WITH_NAME(colName, id_cutBased_HGG_Bits);

  PUSH_VECTOR_WITH_NAME(colName, full5x5_r9);
  PUSH_VECTOR_WITH_NAME(colName, full5x5_sigmaIEtaIEta);
  PUSH_VECTOR_WITH_NAME(colName, full5x5_sigmaIPhiIPhi);

  PUSH_VECTOR_WITH_NAME(colName, MIPTotalEnergy);
  PUSH_VECTOR_WITH_NAME(colName, E4overE1);
  PUSH_VECTOR_WITH_NAME(colName, seedTime);

  PUSH_VECTOR_WITH_NAME(colName, pfIso_comb);
  PUSH_VECTOR_WITH_NAME(colName, pfChargedHadronIso_EAcorr);
  PUSH_VECTOR_WITH_NAME(colName, pfNeutralHadronIso_EAcorr);
  PUSH_VECTOR_WITH_NAME(colName, pfEMIso_EAcorr);

  /*
  PUSH_VECTOR_WITH_NAME(colName, n_associated_pfcands);
  PUSH_VECTOR_WITH_NAME(colName, associated_pfcands_sum_sc_pt);
  */

  PUSH_VECTOR_WITH_NAME(colName, fid_mask);


  /********************************************/
  /* Fill the FSR photon SC veto indices here */
  /********************************************/
  {
    std::string const& colNameFSR = CMS3Ntuplizer::colName_fsrcands;

    MAKE_VECTOR_WITH_DEFAULT_ASSIGN(std::vector<unsigned int>, photonVeto_index_list, (filledFSRCandidates ? filledFSRCandidates->size() : 0));

    if (filledFSRCandidates){
      size_t i_fsrInfo=0;
      for (auto const& fsrInfo:(*filledFSRCandidates)){
        std::vector<unsigned int>& photonVeto_indices = photonVeto_index_list.at(i_fsrInfo);
        photonVeto_indices.reserve(fsrInfo.veto_photon_list.size());
        for (auto const& photon:fsrInfo.veto_photon_list){
          unsigned int iphoton=0;
          for (auto const& obj:filledPhotons){
            if (obj==photon){
              photonVeto_indices.push_back(iphoton);
              break;
            }
            iphoton++;
          }
        }
        // No need to push photonVeto_indices
        i_fsrInfo++;
      }
    }

    PUSH_VECTOR_WITH_NAME(colNameFSR, photonVeto_index_list);
  }
  /*****************/
  /* End FSR block */
  /*****************/


  if (filledObjects) std::swap(filledPhotons, *filledObjects);
  return n_skimmed_objects;
}

size_t CMS3Ntuplizer::fillReducedSuperclusters(
  edm::Event const& iEvent,
  std::vector<pat::Electron const*> const* filledElectrons, std::vector<pat::Photon const*> const* filledPhotons,
  std::vector<reco::SuperCluster const*>* filledObjects
){
  if (!keepExtraSuperclusters) return 0;

  std::string const& colName = CMS3Ntuplizer::colName_superclusters;

  edm::Handle< reco::SuperClusterCollection > reducedSuperclusterHandle;
  iEvent.getByToken(reducedSuperclusterToken, reducedSuperclusterHandle);
  if (!reducedSuperclusterHandle.isValid()) throw cms::Exception("CMS3Ntuplizer::fillReducedSuperclusters: Error getting the reduced supercluster collection from the event...");
  size_t n_objects = reducedSuperclusterHandle->size();

  if (filledObjects) filledObjects->reserve(n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, eta, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, phi, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, energy, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, correctedEnergy, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, sinTheta_SC_pos, n_objects);

  MAKE_VECTOR_WITH_RESERVE(cms3_egamma_fid_type_mask_t, fid_mask, n_objects);

  size_t n_skimmed_objects=0;
  for (reco::SuperClusterCollection::const_iterator obj = reducedSuperclusterHandle->begin(); obj != reducedSuperclusterHandle->end(); obj++){
    if (
      std::max(obj->energy(), obj->correctedEnergy())/std::cosh(obj->eta())<ElectronSelectionHelpers::selection_skim_pt
      ||
      std::abs(obj->eta())>=ElectronSelectionHelpers::selection_skim_eta
      ) continue;

    bool doSkip = false;
    if (filledElectrons && !doSkip){
      for (auto const& part:(*filledElectrons)){
        double deltaR;
        HelperFunctions::deltaR<double>(obj->eta(), obj->phi(), part->eta(), part->phi(), deltaR);
        if (deltaR<0.4 || &(*(part->superCluster()))==&(*obj)){
          doSkip = true;
          break;
        }
      }
    }
    if (filledPhotons && !doSkip){
      for (auto const& part:(*filledPhotons)){
        double deltaR;
        HelperFunctions::deltaR<double>(obj->eta(), obj->phi(), part->eta(), part->phi(), deltaR);
        if (deltaR<0.4 || &(*(part->superCluster()))==&(*obj)){
          doSkip = true;
          break;
        }
      }
    }
    if (doSkip) continue;

    // Core particle quantities
    eta.push_back(obj->eta());
    phi.push_back(obj->phi());
    energy.push_back(obj->energy());
    correctedEnergy.push_back(obj->correctedEnergy());
    sinTheta_SC_pos.push_back(std::sin(obj->position().theta()));

    // Masks
    cms3_egamma_fid_type_mask_t fiducialityMask = 0;  // The enums are in interface/EgammaFiduciality.h
    // This part is from RecoEgamma/EgammaElectronAlgos/src/GsfElectronAlgo.cc
    auto const& seedHitsAndFractions = obj->seed()->hitsAndFractions();
    if (!seedHitsAndFractions.empty()){
      double abs_pos_eta = std::abs(obj->position().eta());
      DetId seedXtalId = seedHitsAndFractions.at(0).first;
      int detector = seedXtalId.subdetId();
      if (detector == EcalBarrel){
        fiducialityMask |= 1 << ISEB;
        EBDetId ebdetid(seedXtalId);
        if (EBDetId::isNextToEtaBoundary(ebdetid)){
          if (ebdetid.ietaAbs() == 85) fiducialityMask |= 1 << ISEBEEGAP;
          else fiducialityMask |= 1 << ISEBETAGAP;
        }
        if (EBDetId::isNextToPhiBoundary(ebdetid)) fiducialityMask |= 1 << ISEBPHIGAP;
      }
      else if (detector == EcalEndcap){
        fiducialityMask |= 1 << ISEE;
        EEDetId eedetid(seedXtalId);
        if (EEDetId::isNextToRingBoundary(eedetid)){
          if (abs_pos_eta < 2.) fiducialityMask |= 1 << ISEBEEGAP;
          else fiducialityMask |= 1 << ISEERINGGAP;
        }
        if (EEDetId::isNextToDBoundary(eedetid)) fiducialityMask |= 1 << ISEEDEEGAP;
      }
#if CMSSW_VERSION_MAJOR>=10 && CMSSW_VERSION_MINOR>=3 
      else if (EcalTools::isHGCalDet((DetId::Detector) seedXtalId.det())) fiducialityMask |= 1 << ISEE;
#endif
      else throw cms::Exception("UnknownXtalRegion") << "Seed cluster region unknown.";

      if (HelperFunctions::test_bit(fiducialityMask, ISEERINGGAP) || HelperFunctions::test_bit(fiducialityMask, ISEEDEEGAP)){
        HelperFunctions::set_bit(fiducialityMask, ISEEGAP);
        HelperFunctions::set_bit(fiducialityMask, ISGAP);
      }
      else if (HelperFunctions::test_bit(fiducialityMask, ISEBETAGAP) || HelperFunctions::test_bit(fiducialityMask, ISEBPHIGAP)){
        HelperFunctions::set_bit(fiducialityMask, ISEBGAP);
        HelperFunctions::set_bit(fiducialityMask, ISGAP);
      }
      else if (HelperFunctions::test_bit(fiducialityMask, ISEBEEGAP)) HelperFunctions::set_bit(fiducialityMask, ISGAP);
    }
    fid_mask.push_back(fiducialityMask);

    if (filledObjects) filledObjects->push_back(&(*obj));
    n_skimmed_objects++;
  }

  // Pass collections to the communicator
  PUSH_VECTOR_WITH_NAME(colName, eta);
  PUSH_VECTOR_WITH_NAME(colName, phi);
  PUSH_VECTOR_WITH_NAME(colName, energy);
  PUSH_VECTOR_WITH_NAME(colName, correctedEnergy);
  PUSH_VECTOR_WITH_NAME(colName, sinTheta_SC_pos);

  PUSH_VECTOR_WITH_NAME(colName, fid_mask);

  return n_skimmed_objects;
}

size_t CMS3Ntuplizer::fillAK4Jets(edm::Event const& iEvent, std::vector<pat::Jet const*>* filledObjects){
  constexpr AK4JetSelectionHelpers::AK4JetType jetType = AK4JetSelectionHelpers::AK4PFCHS;
  std::string const& colName = CMS3Ntuplizer::colName_ak4jets;

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

  // PU jet id should really be for pT<50 GeV, |eta|<5. Note that the order of indices is tight (1), medium (2), loose (4).
  MAKE_VECTOR_WITH_RESERVE(cms3_jet_pujetid_t, pileupJetId, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, pileupJetIdScore, n_objects);
  MAKE_VECTOR_WITH_RESERVE(cms3_jet_pujetid_t, pileupJetId_default, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, pileupJetIdScore_default, n_objects);

  MAKE_VECTOR_WITH_RESERVE(bool, is_genMatched, n_objects);
  MAKE_VECTOR_WITH_RESERVE(bool, is_genMatched_fullCone, n_objects);

  /*
  MAKE_VECTOR_WITH_RESERVE(size_t, n_pfcands, n_objects);
  MAKE_VECTOR_WITH_RESERVE(size_t, n_mucands, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, area, n_objects);
  */

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

  MAKE_VECTOR_WITH_RESERVE(float, NEMF, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, CEMF, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, JECNominal, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, JECUp, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, JECDn, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, JECL1Nominal, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, mucands_sump4_px, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, mucands_sump4_py, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, JERNominal, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, JERUp, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, JERDn, n_objects);

  MAKE_VECTOR_WITH_RESERVE(cms3_jet_genflavor_t, partonFlavour, n_objects);
  MAKE_VECTOR_WITH_RESERVE(cms3_jet_genflavor_t, hadronFlavour, n_objects);

  size_t n_skimmed_objects=0;
  for (edm::View<pat::Jet>::const_iterator obj = ak4jetsHandle->begin(); obj != ak4jetsHandle->end(); obj++){
    if (!AK4JetSelectionHelpers::testSkimAK4Jet(*obj, this->year, jetType)) continue;

    const double uncorrected_energy = AK4JetSelectionHelpers::getUncorrectedJetEnergy(*obj);

    // Core particle quantities
    // These are the uncorrected momentum components!
    pt.push_back(obj->pt());
    eta.push_back(obj->eta());
    phi.push_back(obj->phi());
    mass.push_back(obj->mass());

    pass_looseId.push_back(AK4JetSelectionHelpers::testLooseAK4Jet(*obj, this->year, jetType));
    pass_tightId.push_back(AK4JetSelectionHelpers::testTightAK4Jet(*obj, this->year, jetType));
    pass_leptonVetoId.push_back(AK4JetSelectionHelpers::testLeptonVetoAK4Jet(*obj, this->year, jetType));

    PUSH_USERINT_INTO_VECTOR(pileupJetId);
    PUSH_USERFLOAT_INTO_VECTOR(pileupJetIdScore);
    PUSH_USERINT_INTO_VECTOR(pileupJetId_default);
    PUSH_USERFLOAT_INTO_VECTOR(pileupJetIdScore_default);

    /*
    PUSH_USERINT_INTO_VECTOR(n_pfcands);
    PUSH_USERINT_INTO_VECTOR(n_mucands);
    PUSH_USERFLOAT_INTO_VECTOR(area);
    */

    PUSH_USERINT_INTO_VECTOR(is_genMatched);
    PUSH_USERINT_INTO_VECTOR(is_genMatched_fullCone);

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

    NEMF.push_back(obj->neutralEmEnergy() / uncorrected_energy);
    CEMF.push_back(obj->chargedEmEnergy() / uncorrected_energy);

    PUSH_USERFLOAT_INTO_VECTOR(JECNominal);
    PUSH_USERFLOAT_INTO_VECTOR(JECUp);
    PUSH_USERFLOAT_INTO_VECTOR(JECDn);

    PUSH_USERFLOAT_INTO_VECTOR(JECL1Nominal);
    PUSH_USERFLOAT_INTO_VECTOR(mucands_sump4_px);
    PUSH_USERFLOAT_INTO_VECTOR(mucands_sump4_py);

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

  PUSH_VECTOR_WITH_NAME(colName, is_genMatched);
  PUSH_VECTOR_WITH_NAME(colName, is_genMatched_fullCone);

  PUSH_VECTOR_WITH_NAME(colName, pass_looseId);
  PUSH_VECTOR_WITH_NAME(colName, pass_tightId);
  PUSH_VECTOR_WITH_NAME(colName, pass_leptonVetoId);

  PUSH_VECTOR_WITH_NAME(colName, pileupJetId);
  PUSH_VECTOR_WITH_NAME(colName, pileupJetIdScore);
  PUSH_VECTOR_WITH_NAME(colName, pileupJetId_default);
  PUSH_VECTOR_WITH_NAME(colName, pileupJetIdScore_default);

  /*
  PUSH_VECTOR_WITH_NAME(colName, n_pfcands);
  PUSH_VECTOR_WITH_NAME(colName, n_mucands);
  PUSH_VECTOR_WITH_NAME(colName, area);
  */

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

  PUSH_VECTOR_WITH_NAME(colName, NEMF);
  PUSH_VECTOR_WITH_NAME(colName, CEMF);

  PUSH_VECTOR_WITH_NAME(colName, JECNominal);
  PUSH_VECTOR_WITH_NAME(colName, JECUp);
  PUSH_VECTOR_WITH_NAME(colName, JECDn);

  PUSH_VECTOR_WITH_NAME(colName, JECL1Nominal);
  PUSH_VECTOR_WITH_NAME(colName, mucands_sump4_px);
  PUSH_VECTOR_WITH_NAME(colName, mucands_sump4_py);

  PUSH_VECTOR_WITH_NAME(colName, JERNominal);
  PUSH_VECTOR_WITH_NAME(colName, JERUp);
  PUSH_VECTOR_WITH_NAME(colName, JERDn);

  PUSH_VECTOR_WITH_NAME(colName, partonFlavour);
  PUSH_VECTOR_WITH_NAME(colName, hadronFlavour);

  return n_skimmed_objects;
}
size_t CMS3Ntuplizer::fillAK8Jets(edm::Event const& iEvent, std::vector<pat::Jet const*>* filledObjects){
  const AK8JetSelectionHelpers::AK8JetType jetType = (this->is80X ? AK8JetSelectionHelpers::AK8PFCHS : AK8JetSelectionHelpers::AK8PFPUPPI);
  std::string const& colName = CMS3Ntuplizer::colName_ak8jets;

  edm::Handle< edm::View<pat::Jet> > ak8jetsHandle;
  iEvent.getByToken(ak8jetsToken, ak8jetsHandle);
  if (!ak8jetsHandle.isValid()) throw cms::Exception("CMS3Ntuplizer::fillAK8Jets: Error getting the ak8 jet collection from the event...");
  size_t n_objects = ak8jetsHandle->size();

  if (filledObjects) filledObjects->reserve(n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, pt, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, eta, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, phi, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, mass, n_objects);

  MAKE_VECTOR_WITH_RESERVE(bool, pass_looseId, n_objects);
  MAKE_VECTOR_WITH_RESERVE(bool, pass_tightId, n_objects);
  MAKE_VECTOR_WITH_RESERVE(bool, pass_leptonVetoId, n_objects);

  MAKE_VECTOR_WITH_RESERVE(bool, is_genMatched, n_objects);
  MAKE_VECTOR_WITH_RESERVE(bool, is_genMatched_fullCone, n_objects);

  /*
  MAKE_VECTOR_WITH_RESERVE(size_t, n_pfcands, n_objects);
  MAKE_VECTOR_WITH_RESERVE(size_t, n_mucands, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, area, n_objects);
  */

  MAKE_VECTOR_WITH_RESERVE(size_t, n_softdrop_subjets, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, pt_resolution, n_objects);

  /*
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
  */

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

  MAKE_VECTOR_WITH_RESERVE(cms3_jet_genflavor_t, partonFlavour, n_objects);
  MAKE_VECTOR_WITH_RESERVE(cms3_jet_genflavor_t, hadronFlavour, n_objects);

  size_t n_skimmed_objects=0;
  for (edm::View<pat::Jet>::const_iterator obj = ak8jetsHandle->begin(); obj != ak8jetsHandle->end(); obj++){
    if (!AK8JetSelectionHelpers::testSkimAK8Jet(*obj, this->year)) continue;

    // Core particle quantities
    // These are the uncorrected momentum components!
    pt.push_back(obj->pt());
    eta.push_back(obj->eta());
    phi.push_back(obj->phi());
    mass.push_back(obj->mass());

    pass_looseId.push_back(AK8JetSelectionHelpers::testLooseAK8Jet(*obj, this->year, jetType));
    pass_tightId.push_back(AK8JetSelectionHelpers::testTightAK8Jet(*obj, this->year, jetType));
    pass_leptonVetoId.push_back(AK8JetSelectionHelpers::testLeptonVetoAK8Jet(*obj, this->year, jetType));

    PUSH_USERINT_INTO_VECTOR(is_genMatched);
    PUSH_USERINT_INTO_VECTOR(is_genMatched_fullCone);

    /*
    PUSH_USERINT_INTO_VECTOR(n_pfcands);
    PUSH_USERINT_INTO_VECTOR(n_mucands);
    PUSH_USERFLOAT_INTO_VECTOR(area);
    */

    PUSH_USERINT_INTO_VECTOR(n_softdrop_subjets);

    PUSH_USERFLOAT_INTO_VECTOR(pt_resolution);

    /*
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
    */

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

  PUSH_VECTOR_WITH_NAME(colName, is_genMatched);
  PUSH_VECTOR_WITH_NAME(colName, is_genMatched_fullCone);

  /*
  PUSH_VECTOR_WITH_NAME(colName, n_pfcands);
  PUSH_VECTOR_WITH_NAME(colName, n_mucands);
  PUSH_VECTOR_WITH_NAME(colName, area);
  */

  PUSH_VECTOR_WITH_NAME(colName, n_softdrop_subjets);

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

  std::string const& colName = CMS3Ntuplizer::colName_isotracks;

  edm::Handle< edm::View<IsotrackInfo> > isotracksHandle;
  iEvent.getByToken(isotracksToken, isotracksHandle);
  if (!isotracksHandle.isValid()) throw cms::Exception("CMS3Ntuplizer::fillIsotracks: Error getting the isotrack collection from the event...");
  size_t n_objects = isotracksHandle->size();

  if (filledObjects) filledObjects->reserve(n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, pt, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, eta, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, phi, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, mass, n_objects);

  MAKE_VECTOR_WITH_RESERVE(cms3_id_t, id, n_objects);
  MAKE_VECTOR_WITH_RESERVE(cms3_charge_t, charge, n_objects);

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

  /*
  MAKE_VECTOR_WITH_RESERVE(int, nearestPFcand_id, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, nearestPFcand_deltaR, n_objects);
  */

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

    /*
    PUSH_ISOTRACK_VARIABLE(nearestPFcand_id);
    PUSH_ISOTRACK_VARIABLE(nearestPFcand_deltaR);
    */

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

  /*
  PUSH_VECTOR_WITH_NAME(colName, nearestPFcand_id);
  PUSH_VECTOR_WITH_NAME(colName, nearestPFcand_deltaR);
  */

  return n_skimmed_objects;

#undef PUSH_ISOTRACK_VARIABLE
}

size_t CMS3Ntuplizer::fillVertices(edm::Event const& iEvent, std::vector<reco::Vertex const*>* filledObjects){
  std::string const& colName = CMS3Ntuplizer::colName_vtxs;

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
  std::string const& colName = CMS3Ntuplizer::colName_triggerinfos;

  edm::Handle< edm::View<TriggerInfo> > triggerInfoHandle;
  iEvent.getByToken(triggerInfoToken, triggerInfoHandle);
  if (!triggerInfoHandle.isValid()) throw cms::Exception("CMS3Ntuplizer::fillTriggerInfo: Error getting the trigger infos. from the event...");
  size_t n_triggers = triggerInfoHandle->size();

  MAKE_VECTOR_WITH_RESERVE(std::string, name, n_triggers);
#if TRIGGEROBJECTINFO_INDEX_BY_ORIGINAL == 0
#else
  MAKE_VECTOR_WITH_RESERVE(unsigned int, index, n_triggers);
#endif
  MAKE_VECTOR_WITH_RESERVE(bool, passTrigger, n_triggers);
  MAKE_VECTOR_WITH_RESERVE(int, L1prescale, n_triggers);
  MAKE_VECTOR_WITH_RESERVE(int, HLTprescale, n_triggers);

  bool passAtLeastOneTrigger = false;
  for (edm::View<TriggerInfo>::const_iterator trigInfo = triggerInfoHandle->begin(); trigInfo != triggerInfoHandle->end(); trigInfo++){
    name.emplace_back(trigInfo->name);
#if TRIGGEROBJECTINFO_INDEX_BY_ORIGINAL == 0
#else
    index.emplace_back(trigInfo->index);
#endif
    passTrigger.emplace_back(trigInfo->passTrigger);
    L1prescale.emplace_back(trigInfo->L1prescale);
    HLTprescale.emplace_back(trigInfo->HLTprescale);

    passAtLeastOneTrigger |= trigInfo->passTrigger;
  }

  PUSH_VECTOR_WITH_NAME(colName, name);
#if TRIGGEROBJECTINFO_INDEX_BY_ORIGINAL == 0
#else
  // No need to record indices since matching to position is done below
  /*PUSH_VECTOR_WITH_NAME(colName, index);*/
#endif
  PUSH_VECTOR_WITH_NAME(colName, passTrigger);
  PUSH_VECTOR_WITH_NAME(colName, L1prescale);
  PUSH_VECTOR_WITH_NAME(colName, HLTprescale);

  // Trigger objects
  MAKE_VECTOR_WITHOUT_RESERVE(cms3_id_t, type);
  MAKE_VECTOR_WITHOUT_RESERVE(float, pt);
  MAKE_VECTOR_WITHOUT_RESERVE(float, eta);
  MAKE_VECTOR_WITHOUT_RESERVE(float, phi);
  MAKE_VECTOR_WITHOUT_RESERVE(float, mass);
  MAKE_VECTOR_WITHOUT_RESERVE(std::vector<unsigned int>, associatedTriggers);
  MAKE_VECTOR_WITHOUT_RESERVE(std::vector<unsigned int>, passedTriggers);
  if (processTriggerObjectInfos){
    edm::Handle< edm::View<TriggerObjectInfo> > triggerObjectInfoHandle;
    iEvent.getByToken(triggerObjectInfoToken, triggerObjectInfoHandle);
    if (!triggerObjectInfoHandle.isValid()) throw cms::Exception("CMS3Ntuplizer::fillTriggerInfo: Error getting the trigger object infos. from the event...");
    size_t n_triggerobjects = triggerObjectInfoHandle->size();

    RESERVE_VECTOR(type, n_triggerobjects);
    RESERVE_VECTOR(pt, n_triggerobjects);
    RESERVE_VECTOR(eta, n_triggerobjects);
    RESERVE_VECTOR(phi, n_triggerobjects);
    RESERVE_VECTOR(mass, n_triggerobjects);
    RESERVE_VECTOR(associatedTriggers, n_triggerobjects);
    RESERVE_VECTOR(passedTriggers, n_triggerobjects);

    for (edm::View<TriggerObjectInfo>::const_iterator trigObj = triggerObjectInfoHandle->begin(); trigObj != triggerObjectInfoHandle->end(); trigObj++){
      type.push_back(trigObj->bestType());
      pt.push_back(trigObj->p4.Pt());
      eta.push_back(trigObj->p4.Eta());
      phi.push_back(trigObj->p4.Phi());
      mass.push_back(trigObj->p4.M());

      size_t const nAssociatedTriggers = trigObj->triggerCollectionIndices.size();
      assert(nAssociatedTriggers == trigObj->passAllTriggerFiltersList.size());

      associatedTriggers.emplace_back(std::vector<unsigned int>());
      std::vector<unsigned int>& trigObj_associatedTriggers = associatedTriggers.back();
      trigObj_associatedTriggers.reserve(nAssociatedTriggers);

      passedTriggers.emplace_back(std::vector<unsigned int>());
      std::vector<unsigned int>& trigObj_passedTriggers = passedTriggers.back();
      trigObj_passedTriggers.reserve(nAssociatedTriggers);

      std::vector<unsigned int>::const_iterator triggerCollectionIndices_end = trigObj->triggerCollectionIndices.cend();
      //std::vector<bool>::const_iterator passAllTriggerFiltersList_end = trigObj->passAllTriggerFiltersList.cend();
      std::vector<unsigned int>::const_iterator it_triggerCollectionIndices = trigObj->triggerCollectionIndices.cbegin();
      std::vector<bool>::const_iterator it_passAllTriggerFiltersList = trigObj->passAllTriggerFiltersList.cbegin();
      while (it_triggerCollectionIndices!=triggerCollectionIndices_end){
#if TRIGGEROBJECTINFO_INDEX_BY_ORIGINAL == 0
        unsigned int const& pos = *it_triggerCollectionIndices;
#else
        unsigned int pos=0;
        for (auto const& trigIndex:index){
          if (*it_triggerCollectionIndices == trigIndex) break;
          pos++;
        }
#endif
        if (pos>=name.size()) throw cms::Exception("CMS3Ntuplizer::fillTriggerInfo: Trigger object position index reached trigger list size!");

        trigObj_associatedTriggers.emplace_back(pos);
        if (*it_passAllTriggerFiltersList) trigObj_passedTriggers.emplace_back(pos);

        it_triggerCollectionIndices++;
        it_passAllTriggerFiltersList++;
      }
    }
  }
  PUSH_VECTOR_WITH_NAME(colName_triggerobjects, type);
  PUSH_VECTOR_WITH_NAME(colName_triggerobjects, pt);
  PUSH_VECTOR_WITH_NAME(colName_triggerobjects, eta);
  PUSH_VECTOR_WITH_NAME(colName_triggerobjects, phi);
  PUSH_VECTOR_WITH_NAME(colName_triggerobjects, mass);
  PUSH_VECTOR_WITH_NAME(colName_triggerobjects, associatedTriggers);
  PUSH_VECTOR_WITH_NAME(colName_triggerobjects, passedTriggers);

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
#define SET_MET_SHIFT(NAME, COLLNAME, VAL) commonEntry.setNamedVal((std::string(COLLNAME) + "_" + #NAME).data(), VAL);

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

  SET_MET_VARIABLE(metHandle, met_JECUp, pfmetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_JECUp, pfmetCollName);
  SET_MET_VARIABLE(metHandle, met_JECDn, pfmetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_JECDn, pfmetCollName);

  SET_MET_VARIABLE(metHandle, met_UnclusteredEnUp, pfmetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_UnclusteredEnUp, pfmetCollName);
  SET_MET_VARIABLE(metHandle, met_UnclusteredEnDn, pfmetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_UnclusteredEnDn, pfmetCollName);

  if (isMC){
    edm::Handle< reco::Particle::LorentzVector > pfmetshiftHandle_JERNominal;
    iEvent.getByToken(pfmetshiftToken_JERNominal, pfmetshiftHandle_JERNominal);
    if (!pfmetshiftHandle_JERNominal.isValid()) throw cms::Exception("CMS3Ntuplizer::fillMETVariables: Error getting the PF MET JER nominal shift handle from the event...");
    SET_MET_SHIFT(metShift_px_JERNominal, pfmetCollName, float(pfmetshiftHandle_JERNominal->Px()));
    SET_MET_SHIFT(metShift_py_JERNominal, pfmetCollName, float(pfmetshiftHandle_JERNominal->Py()));

    edm::Handle< reco::Particle::LorentzVector > pfmetshiftHandle_JERUp;
    iEvent.getByToken(pfmetshiftToken_JERUp, pfmetshiftHandle_JERUp);
    if (!pfmetshiftHandle_JERUp.isValid()) throw cms::Exception("CMS3Ntuplizer::fillMETVariables: Error getting the PF MET JER up shift handle from the event...");
    SET_MET_SHIFT(metShift_px_JERUp, pfmetCollName, float(pfmetshiftHandle_JERUp->Px()));
    SET_MET_SHIFT(metShift_py_JERUp, pfmetCollName, float(pfmetshiftHandle_JERUp->Py()));

    edm::Handle< reco::Particle::LorentzVector > pfmetshiftHandle_JERDn;
    iEvent.getByToken(pfmetshiftToken_JERDn, pfmetshiftHandle_JERDn);
    if (!pfmetshiftHandle_JERDn.isValid()) throw cms::Exception("CMS3Ntuplizer::fillMETVariables: Error getting the PF MET JER down shift handle from the event...");
    SET_MET_SHIFT(metShift_px_JERDn, pfmetCollName, float(pfmetshiftHandle_JERDn->Px()));
    SET_MET_SHIFT(metShift_py_JERDn, pfmetCollName, float(pfmetshiftHandle_JERDn->Py()));
  }
  else{
    SET_MET_SHIFT(metShift_px_JERNominal, pfmetCollName, float(0));
    SET_MET_SHIFT(metShift_py_JERNominal, pfmetCollName, float(0));
    SET_MET_SHIFT(metShift_px_JERUp, pfmetCollName, float(0));
    SET_MET_SHIFT(metShift_py_JERUp, pfmetCollName, float(0));
    SET_MET_SHIFT(metShift_px_JERDn, pfmetCollName, float(0));
    SET_MET_SHIFT(metShift_py_JERDn, pfmetCollName, float(0));
  }

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

  SET_MET_VARIABLE(metHandle, met_JECUp, puppimetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_JECUp, puppimetCollName);
  SET_MET_VARIABLE(metHandle, met_JECDn, puppimetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_JECDn, puppimetCollName);

  SET_MET_VARIABLE(metHandle, met_UnclusteredEnUp, puppimetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_UnclusteredEnUp, puppimetCollName);
  SET_MET_VARIABLE(metHandle, met_UnclusteredEnDn, puppimetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_UnclusteredEnDn, puppimetCollName);

  /*
  if (isMC){
    edm::Handle< reco::Particle::LorentzVector > puppimetshiftHandle_JERNominal;
    iEvent.getByToken(puppimetshiftToken_JERNominal, puppimetshiftHandle_JERNominal);
    if (!puppimetshiftHandle_JERNominal.isValid()) throw cms::Exception("CMS3Ntuplizer::fillMETVariables: Error getting the PF MET JER nominal shift handle from the event...");
    SET_MET_SHIFT(metShift_px_JERNominal, puppimetCollName, float(puppimetshiftHandle_JERNominal->Px()));
    SET_MET_SHIFT(metShift_py_JERNominal, puppimetCollName, float(puppimetshiftHandle_JERNominal->Py()));

    edm::Handle< reco::Particle::LorentzVector > puppimetshiftHandle_JERUp;
    iEvent.getByToken(puppimetshiftToken_JERUp, puppimetshiftHandle_JERUp);
    if (!puppimetshiftHandle_JERUp.isValid()) throw cms::Exception("CMS3Ntuplizer::fillMETVariables: Error getting the PF MET JER up shift handle from the event...");
    SET_MET_SHIFT(metShift_px_JERUp, puppimetCollName, float(puppimetshiftHandle_JERUp->Px()));
    SET_MET_SHIFT(metShift_py_JERUp, puppimetCollName, float(puppimetshiftHandle_JERUp->Py()));

    edm::Handle< reco::Particle::LorentzVector > puppimetshiftHandle_JERDn;
    iEvent.getByToken(puppimetshiftToken_JERDn, puppimetshiftHandle_JERDn);
    if (!puppimetshiftHandle_JERDn.isValid()) throw cms::Exception("CMS3Ntuplizer::fillMETVariables: Error getting the PF MET JER down shift handle from the event...");
    SET_MET_SHIFT(metShift_px_JERDn, puppimetCollName, float(puppimetshiftHandle_JERDn->Px()));
    SET_MET_SHIFT(metShift_py_JERDn, puppimetCollName, float(puppimetshiftHandle_JERDn->Py()));
  }
  else{
    SET_MET_SHIFT(metShift_px_JERNominal, puppimetCollName, float(0));
    SET_MET_SHIFT(metShift_py_JERNominal, puppimetCollName, float(0));
    SET_MET_SHIFT(metShift_px_JERUp, puppimetCollName, float(0));
    SET_MET_SHIFT(metShift_py_JERUp, puppimetCollName, float(0));
    SET_MET_SHIFT(metShift_px_JERDn, puppimetCollName, float(0));
    SET_MET_SHIFT(metShift_py_JERDn, puppimetCollName, float(0));
  }
  */

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

#undef SET_MET_SHIFT
#undef SET_MET_VARIABLE

  return true;
}

bool CMS3Ntuplizer::fillGenVariables(
  edm::Event const& iEvent,
  std::vector<pat::Muon const*>* filledMuons,
  std::vector<pat::Electron const*>* filledElectrons,
  std::vector<pat::Photon const*>* filledPhotons,
  std::vector<reco::GenParticle const*>* filledPrunedGenParts,
  std::vector<pat::PackedGenParticle const*>* filledPackedGenParts,
  std::vector<reco::GenJet const*>* filledGenAK4Jets,
  std::vector<reco::GenJet const*>* filledGenAK8Jets
){
  if (!this->isMC) return true;

  // Gen. info.
  recordGenInfo(iEvent);

  // Gen. particles
  if (this->keepGenParticles!=kNone) recordGenParticles(
    iEvent,
    filledMuons, filledElectrons, filledPhotons,
    filledPrunedGenParts, filledPackedGenParts
  );

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
      (i==kPromptFinalStatePhotons && word=="promptfinalstatephotons")
      ||
      (i==kReducedFinalStates && word=="reducedfinalstates")
      ||
      (i==kAllFinalStates && word=="allfinalstates")
      ||
      (i==kReducedFinalStatesAndHardProcessFinalStates && word=="reducedfinalstatesandhardprocessfinalstates")
      ||
      (i==kReducedFinalStatesAndHardProcesses && word=="reducedfinalstatesandhardprocesses")
      ||
      (i==kHardProcesses && word=="hardprocesses")
      ||
      (i==kHardProcessFinalStates && word=="hardprocessfinalstates")
      ||
      (i==kAll && word=="all")
      ) res = static_cast<ParticleRecordLevel>(i);
  }
  if (res==nParticleRecordLevels) throw cms::Exception((std::string("CMS3Ntuplizer::getParticleRecordLevel: Word ") + word + " is not recognized!").data());
  return res;
}


DEFINE_FWK_MODULE(CMS3Ntuplizer);
