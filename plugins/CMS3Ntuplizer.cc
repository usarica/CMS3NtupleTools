#include <CommonTools/UtilAlgos/interface/TFileService.h>
#include <CMS3/NtupleMaker/interface/plugins/CMS3Ntuplizer.h>
#include "CMS3/NtupleMaker/interface/VertexSelectionHelpers.h"
#include "CMS3/NtupleMaker/interface/AK4JetSelectionHelpers.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;
using namespace edm;


CMS3Ntuplizer::CMS3Ntuplizer(const edm::ParameterSet& pset_) :
  pset(pset_),
  outtree(nullptr),
  firstEvent(true),

  year(pset.getParameter<int>("year")),
  treename(pset.getUntrackedParameter<std::string>("treename")),
  isMC(pset.getParameter<bool>("isMC"))
{
  if (year!=2016 && year!=2017 && year!=2018) throw cms::Exception("CMS3Ntuplizer::CMS3Ntuplizer: Year is undefined!");

  electronsToken  = consumes< edm::View<pat::Electron> >(pset.getParameter<edm::InputTag>("electronSrc"));
  photonsToken  = consumes< edm::View<pat::Photon> >(pset.getParameter<edm::InputTag>("photonSrc"));
  muonsToken  = consumes< edm::View<pat::Muon> >(pset.getParameter<edm::InputTag>("muonSrc"));
  ak4jetsToken  = consumes< edm::View<pat::Jet> >(pset.getParameter<edm::InputTag>("ak4jetSrc"));

  pfmetToken = consumes< METInfo >(pset.getParameter<edm::InputTag>("pfmetSrc"));
  puppimetToken = consumes< METInfo >(pset.getParameter<edm::InputTag>("puppimetSrc"));

  vtxToken = consumes<reco::VertexCollection>(pset.getParameter<edm::InputTag>("vtxSrc"));

  rhoToken  = consumes< double >(pset.getParameter<edm::InputTag>("rhoSrc"));
  triggerInfoToken = consumes< edm::View<TriggerInfo> >(pset.getParameter<edm::InputTag>("triggerInfoSrc"));
  puInfoToken = consumes< std::vector<PileupSummaryInfo> >(pset.getParameter<edm::InputTag>("puInfoSrc"));
  metFilterInfoToken = consumes< METFilterInfo >(pset.getParameter<edm::InputTag>("metFilterInfoSrc"));

  genInfoToken = consumes<GenInfo>(pset.getParameter<edm::InputTag>("genInfoSrc"));
  prunedGenParticlesToken = consumes<reco::GenParticleCollection>(pset.getParameter<edm::InputTag>("prunedGenParticlesSrc"));
  genJetsToken = consumes< edm::View<reco::GenJet> >(pset.getParameter<edm::InputTag>("genJetsSrc"));

}
CMS3Ntuplizer::~CMS3Ntuplizer(){
  //delete pileUpReweight;
  //delete metCorrHandler;
}


void CMS3Ntuplizer::beginJob(){
  edm::Service<TFileService> fs;
  TTree* tout = fs->make<TTree>(treename, "Selected event summary");
  outtree = new BaseTree(nullptr, tout, nullptr, nullptr, false);
  outtree->setAcquireTreePossession(false);
}
void CMS3Ntuplizer::endJob(){
  delete outtree;
}

void CMS3Ntuplizer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){}
void CMS3Ntuplizer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){}

void CMS3Ntuplizer::beginRun(edm::Run const&, edm::EventSetup const&){}
void CMS3Ntuplizer::endRun(edm::Run const&, edm::EventSetup const&){}


// Convenience macros to easily make and push vector values
#define MAKE_VECTOR_WITH_RESERVE(type_, name_, size_) std::vector<type_> name_; name_.reserve(size_);
#define PUSH_USERINT_INTO_VECTOR(name_) name_.push_back(obj->userInt(#name_));
#define PUSH_USERFLOAT_INTO_VECTOR(name_) name_.push_back(obj->userFloat(#name_));
#define PUSH_VECTOR_WITH_NAME(name_, var_) commonEntry.setNamedVal(TString(name_)+"_"+#var_, var_);


void CMS3Ntuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
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

  // Electrons
  size_t n_electrons = this->fillElectrons(iEvent, nullptr);

  // Photons
  size_t n_photons = this->fillPhotons(iEvent, nullptr);

  // Muons
  size_t n_muons = this->fillMuons(iEvent, nullptr);

  // ak4 jets
  /*size_t n_ak4jets = */this->fillAK4Jets(iEvent, nullptr);

  // The (data) event should have at least one electron, muon, or photon.
  isSelected &= ((n_muons + n_electrons + n_photons)>0);

  // MET info
  isSelected &= this->fillMETVariables(iEvent);

  // Event info
  isSelected &= this->fillEventVariables(iEvent);

  // Trigger info
  isSelected &= this->fillTriggerInfo(iEvent);

  // MET filters
  isSelected &= this->fillMETFilterVariables(iEvent);

  // GenInfo
  isSelected &= this->fillGenVariables(iEvent);


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
  if (isMC || isSelected) outtree->fill();

  // No longer the first event...
  if (firstEvent) firstEvent = false;
}

void CMS3Ntuplizer::recordGenInfo(GenInfo const& genInfo){
#define SET_GENINFO_VARIABLE(var) commonEntry.setNamedVal(#var, genInfo.var);

  SET_GENINFO_VARIABLE(xsec)
  SET_GENINFO_VARIABLE(xsecerr)

  SET_GENINFO_VARIABLE(qscale)
  SET_GENINFO_VARIABLE(alphaS)

  SET_GENINFO_VARIABLE(genMET)
  SET_GENINFO_VARIABLE(genMETPhi)

  SET_GENINFO_VARIABLE(sumEt)
  SET_GENINFO_VARIABLE(pThat)

  // LHE variations
  SET_GENINFO_VARIABLE(genHEPMCweight_default)
  SET_GENINFO_VARIABLE(genHEPMCweight_NNPDF30)

  SET_GENINFO_VARIABLE(LHEweight_QCDscale_muR1_muF1)
  SET_GENINFO_VARIABLE(LHEweight_QCDscale_muR1_muF2)
  SET_GENINFO_VARIABLE(LHEweight_QCDscale_muR1_muF0p5)
  SET_GENINFO_VARIABLE(LHEweight_QCDscale_muR2_muF1)
  SET_GENINFO_VARIABLE(LHEweight_QCDscale_muR2_muF2)
  SET_GENINFO_VARIABLE(LHEweight_QCDscale_muR2_muF0p5)
  SET_GENINFO_VARIABLE(LHEweight_QCDscale_muR0p5_muF1)
  SET_GENINFO_VARIABLE(LHEweight_QCDscale_muR0p5_muF2)
  SET_GENINFO_VARIABLE(LHEweight_QCDscale_muR0p5_muF0p5)

  SET_GENINFO_VARIABLE(LHEweight_PDFVariation_Up_default)
  SET_GENINFO_VARIABLE(LHEweight_PDFVariation_Dn_default)
  SET_GENINFO_VARIABLE(LHEweight_AsMZ_Up_default)
  SET_GENINFO_VARIABLE(LHEweight_AsMZ_Dn_default)

  SET_GENINFO_VARIABLE(LHEweight_PDFVariation_Up_NNPDF30)
  SET_GENINFO_VARIABLE(LHEweight_PDFVariation_Dn_NNPDF30)
  SET_GENINFO_VARIABLE(LHEweight_AsMZ_Up_NNPDF30)
  SET_GENINFO_VARIABLE(LHEweight_AsMZ_Dn_NNPDF30)

  // Pythis PS weights
  SET_GENINFO_VARIABLE(PythiaWeight_isr_muRoneoversqrt2)
  SET_GENINFO_VARIABLE(PythiaWeight_fsr_muRoneoversqrt2)
  SET_GENINFO_VARIABLE(PythiaWeight_isr_muRsqrt2)
  SET_GENINFO_VARIABLE(PythiaWeight_fsr_muRsqrt2)
  SET_GENINFO_VARIABLE(PythiaWeight_isr_muR0p5)
  SET_GENINFO_VARIABLE(PythiaWeight_fsr_muR0p5)
  SET_GENINFO_VARIABLE(PythiaWeight_isr_muR2)
  SET_GENINFO_VARIABLE(PythiaWeight_fsr_muR2)
  SET_GENINFO_VARIABLE(PythiaWeight_isr_muR0p25)
  SET_GENINFO_VARIABLE(PythiaWeight_fsr_muR0p25)
  SET_GENINFO_VARIABLE(PythiaWeight_isr_muR4)
  SET_GENINFO_VARIABLE(PythiaWeight_fsr_muR4)

#undef SET_GENINFO_VARIABLE

  for (auto const it:genInfo.LHE_ME_weights) commonEntry.setNamedVal(it.first, it.second);
}

size_t CMS3Ntuplizer::fillElectrons(const edm::Event& iEvent, std::vector<pat::Electron const*>* filledObjects){
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

  // Has no convention correspondence in nanoAOD
  MAKE_VECTOR_WITH_RESERVE(float, scale_smear_corr, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, scale_smear_corr_scale_totalUp, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, scale_smear_corr_scale_totalDn, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, scale_smear_corr_smear_totalUp, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, scale_smear_corr_smear_totalDn, n_objects);

  MAKE_VECTOR_WITH_RESERVE(float, id_MVA_Fall17V2_Iso_Val, n_objects);
  MAKE_VECTOR_WITH_RESERVE(unsigned int, id_MVA_Fall17V2_Iso_Cat, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, id_MVA_Fall17V2_NoIso_Val, n_objects);
  MAKE_VECTOR_WITH_RESERVE(unsigned int, id_MVA_Fall17V2_NoIso_Cat, n_objects);

  MAKE_VECTOR_WITH_RESERVE(unsigned int, id_cutBased_Fall17V2_Veto_Bits, n_objects);
  MAKE_VECTOR_WITH_RESERVE(unsigned int, id_cutBased_Fall17V2_Loose_Bits, n_objects);
  MAKE_VECTOR_WITH_RESERVE(unsigned int, id_cutBased_Fall17V2_Medium_Bits, n_objects);
  MAKE_VECTOR_WITH_RESERVE(unsigned int, id_cutBased_Fall17V2_Tight_Bits, n_objects);

  MAKE_VECTOR_WITH_RESERVE(unsigned int, id_cutBased_Fall17V1_Veto_Bits, n_objects);
  MAKE_VECTOR_WITH_RESERVE(unsigned int, id_cutBased_Fall17V1_Loose_Bits, n_objects);
  MAKE_VECTOR_WITH_RESERVE(unsigned int, id_cutBased_Fall17V1_Medium_Bits, n_objects);
  MAKE_VECTOR_WITH_RESERVE(unsigned int, id_cutBased_Fall17V1_Tight_Bits, n_objects);

  MAKE_VECTOR_WITH_RESERVE(unsigned int, fid_mask, n_objects);
  MAKE_VECTOR_WITH_RESERVE(unsigned int, type_mask, n_objects);

  for (View<pat::Electron>::const_iterator obj = electronsHandle->begin(); obj != electronsHandle->end(); obj++){
    // Core particle quantities
    // Uncorrected p4
    pt.push_back(obj->pt());
    eta.push_back(obj->eta());
    phi.push_back(obj->phi());
    mass.push_back(obj->mass());

    // Charge: Can obtain pdgId from this, so no need to record pdgId again
    PUSH_USERINT_INTO_VECTOR(charge);

    // Scale and smear
    // Nominal value: Needs to multiply the uncorrected p4 at analysis level
    PUSH_USERFLOAT_INTO_VECTOR(scale_smear_corr);
    // Uncertainties: Only store total up/dn for the moment
    PUSH_USERFLOAT_INTO_VECTOR(scale_smear_corr_scale_totalUp);
    PUSH_USERFLOAT_INTO_VECTOR(scale_smear_corr_scale_totalDn);
    PUSH_USERFLOAT_INTO_VECTOR(scale_smear_corr_smear_totalUp);
    PUSH_USERFLOAT_INTO_VECTOR(scale_smear_corr_smear_totalDn);

    // Id variables
    // Fall17V2_Iso MVA id
    PUSH_USERFLOAT_INTO_VECTOR(id_MVA_Fall17V2_Iso_Val);
    PUSH_USERINT_INTO_VECTOR(id_MVA_Fall17V2_Iso_Cat);

    // Fall17V2_NoIso MVA id
    PUSH_USERFLOAT_INTO_VECTOR(id_MVA_Fall17V2_NoIso_Val);
    PUSH_USERINT_INTO_VECTOR(id_MVA_Fall17V2_NoIso_Cat);

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

    // Masks
    PUSH_USERINT_INTO_VECTOR(fid_mask);
    PUSH_USERINT_INTO_VECTOR(type_mask);

    if (filledObjects) filledObjects->push_back(&(*obj));
  }

  // Pass collections to the communicator
  PUSH_VECTOR_WITH_NAME(colName, pt);
  PUSH_VECTOR_WITH_NAME(colName, eta);
  PUSH_VECTOR_WITH_NAME(colName, phi);
  PUSH_VECTOR_WITH_NAME(colName, mass);

  PUSH_VECTOR_WITH_NAME(colName, charge);

  // Has no convention correspondence in nanoAOD
  PUSH_VECTOR_WITH_NAME(colName, scale_smear_corr);
  PUSH_VECTOR_WITH_NAME(colName, scale_smear_corr_scale_totalUp);
  PUSH_VECTOR_WITH_NAME(colName, scale_smear_corr_scale_totalDn);
  PUSH_VECTOR_WITH_NAME(colName, scale_smear_corr_smear_totalUp);
  PUSH_VECTOR_WITH_NAME(colName, scale_smear_corr_smear_totalDn);

  PUSH_VECTOR_WITH_NAME(colName, id_MVA_Fall17V2_Iso_Val);
  PUSH_VECTOR_WITH_NAME(colName, id_MVA_Fall17V2_Iso_Cat);
  PUSH_VECTOR_WITH_NAME(colName, id_MVA_Fall17V2_NoIso_Val);
  PUSH_VECTOR_WITH_NAME(colName, id_MVA_Fall17V2_NoIso_Cat);

  PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V2_Veto_Bits);
  PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V2_Loose_Bits);
  PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V2_Medium_Bits);
  PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V2_Tight_Bits);

  PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V1_Veto_Bits);
  PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V1_Loose_Bits);
  PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V1_Medium_Bits);
  PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V1_Tight_Bits);

  PUSH_VECTOR_WITH_NAME(colName, fid_mask);
  PUSH_VECTOR_WITH_NAME(colName, type_mask);

  return n_objects;
}
size_t CMS3Ntuplizer::fillPhotons(const edm::Event& iEvent, std::vector<pat::Photon const*>* filledObjects){
  const char colName[] = "photons";
  edm::Handle< edm::View<pat::Photon> > photonsHandle;
  iEvent.getByToken(photonsToken, photonsHandle);
  if (!photonsHandle.isValid()) throw cms::Exception("CMS3Ntuplizer::fillPhotons: Error getting the photon collection from the event...");
  size_t n_objects = photonsHandle->size();

  if (filledObjects) filledObjects->reserve(n_objects);

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

  MAKE_VECTOR_WITH_RESERVE(unsigned int, id_cutBased_Fall17V2_Loose_Bits, n_objects);
  MAKE_VECTOR_WITH_RESERVE(unsigned int, id_cutBased_Fall17V2_Medium_Bits, n_objects);
  MAKE_VECTOR_WITH_RESERVE(unsigned int, id_cutBased_Fall17V2_Tight_Bits, n_objects);

  for (View<pat::Photon>::const_iterator obj = photonsHandle->begin(); obj != photonsHandle->end(); obj++){
    // Core particle quantities
    // Uncorrected p4
    pt.push_back(obj->pt());
    eta.push_back(obj->eta());
    phi.push_back(obj->phi());
    mass.push_back(obj->mass());

    // Scale and smear
    // Nominal value: Needs to multiply the uncorrected p4 at analysis level
    PUSH_USERFLOAT_INTO_VECTOR(scale_smear_corr);
    // Uncertainties: Only store total up/dn for the moment
    PUSH_USERFLOAT_INTO_VECTOR(scale_smear_corr_scale_totalUp);
    PUSH_USERFLOAT_INTO_VECTOR(scale_smear_corr_scale_totalDn);
    PUSH_USERFLOAT_INTO_VECTOR(scale_smear_corr_smear_totalUp);
    PUSH_USERFLOAT_INTO_VECTOR(scale_smear_corr_smear_totalDn);

    // Id variables
    // Fall17V2 cut-based ids
    PUSH_USERINT_INTO_VECTOR(id_cutBased_Fall17V2_Loose_Bits);
    PUSH_USERINT_INTO_VECTOR(id_cutBased_Fall17V2_Medium_Bits);
    PUSH_USERINT_INTO_VECTOR(id_cutBased_Fall17V2_Tight_Bits);

    if (filledObjects) filledObjects->push_back(&(*obj));
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

  PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V2_Loose_Bits);
  PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V2_Medium_Bits);
  PUSH_VECTOR_WITH_NAME(colName, id_cutBased_Fall17V2_Tight_Bits);

  return n_objects;
}
size_t CMS3Ntuplizer::fillMuons(const edm::Event& iEvent, std::vector<pat::Muon const*>* filledObjects){
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

  MAKE_VECTOR_WITH_RESERVE(float, scale_smear_pt_corr, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, scale_smear_pt_corr_scale_totalUp, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, scale_smear_pt_corr_scale_totalDn, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, scale_smear_pt_corr_smear_totalUp, n_objects);
  MAKE_VECTOR_WITH_RESERVE(float, scale_smear_pt_corr_smear_totalDn, n_objects);

  for (View<pat::Muon>::const_iterator obj = muonsHandle->begin(); obj != muonsHandle->end(); obj++){
    // Core particle quantities
    pt.push_back(obj->pt());
    eta.push_back(obj->eta());
    phi.push_back(obj->phi());
    mass.push_back(obj->mass());

    PUSH_USERINT_INTO_VECTOR(charge);

    PUSH_USERFLOAT_INTO_VECTOR(scale_smear_pt_corr);
    PUSH_USERFLOAT_INTO_VECTOR(scale_smear_pt_corr_scale_totalUp);
    PUSH_USERFLOAT_INTO_VECTOR(scale_smear_pt_corr_scale_totalDn);
    PUSH_USERFLOAT_INTO_VECTOR(scale_smear_pt_corr_smear_totalUp);
    PUSH_USERFLOAT_INTO_VECTOR(scale_smear_pt_corr_smear_totalDn);

    if (filledObjects) filledObjects->push_back(&(*obj));
  }

  // Pass collections to the communicator
  PUSH_VECTOR_WITH_NAME(colName, pt);
  PUSH_VECTOR_WITH_NAME(colName, eta);
  PUSH_VECTOR_WITH_NAME(colName, phi);
  PUSH_VECTOR_WITH_NAME(colName, mass);

  PUSH_VECTOR_WITH_NAME(colName, charge);

  PUSH_VECTOR_WITH_NAME(colName, scale_smear_pt_corr);
  PUSH_VECTOR_WITH_NAME(colName, scale_smear_pt_corr_scale_totalUp);
  PUSH_VECTOR_WITH_NAME(colName, scale_smear_pt_corr_scale_totalDn);
  PUSH_VECTOR_WITH_NAME(colName, scale_smear_pt_corr_smear_totalUp);
  PUSH_VECTOR_WITH_NAME(colName, scale_smear_pt_corr_smear_totalDn);

  return n_objects;
}
size_t CMS3Ntuplizer::fillAK4Jets(const edm::Event& iEvent, std::vector<pat::Jet const*>* filledObjects){
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
  MAKE_VECTOR_WITH_RESERVE(bool, pass_puId, n_objects);

  for (View<pat::Jet>::const_iterator obj = ak4jetsHandle->begin(); obj != ak4jetsHandle->end(); obj++){
    // Core particle quantities
    // These are the uncorrected momentum components!
    pt.push_back(obj->pt());
    eta.push_back(obj->eta());
    phi.push_back(obj->phi());
    mass.push_back(obj->mass());

    pass_looseId.push_back(AK4JetSelectionHelpers::testLooseAK4Jet(*obj, this->year));
    pass_tightId.push_back(AK4JetSelectionHelpers::testTightAK4Jet(*obj, this->year));
    pass_puId.push_back(AK4JetSelectionHelpers::testPileUpAK4Jet(*obj, this->year));

    if (filledObjects) filledObjects->push_back(&(*obj));
  }

  // Pass collections to the communicator
  PUSH_VECTOR_WITH_NAME(colName, pt);
  PUSH_VECTOR_WITH_NAME(colName, eta);
  PUSH_VECTOR_WITH_NAME(colName, phi);
  PUSH_VECTOR_WITH_NAME(colName, mass);

  PUSH_VECTOR_WITH_NAME(colName, pass_looseId);
  PUSH_VECTOR_WITH_NAME(colName, pass_tightId);
  PUSH_VECTOR_WITH_NAME(colName, pass_puId);

  return n_objects;
}
size_t CMS3Ntuplizer::fillVertices(const edm::Event& iEvent, std::vector<reco::Vertex const*>* filledObjects){
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

  for (reco::VertexCollection::const_iterator obj = vtxHandle->begin(); obj != vtxHandle->end(); obj++){
    auto const& pos = obj->position();

    is_fake.push_back(obj->isFake());
    is_valid.push_back(obj->isValid());
    is_good.push_back(VertexSelectionHelpers::testGoodVertex(*obj));

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
  }

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

  return n_objects;
}

bool CMS3Ntuplizer::fillEventVariables(const edm::Event& iEvent){
  edm::Handle< double > rhoHandle;
  iEvent.getByToken(rhoToken, rhoHandle);
  if (!rhoHandle.isValid()) throw cms::Exception("CMS3Ntuplizer::fillEventVariables: Error getting the rho collection from the event...");

  edm::Handle< std::vector<PileupSummaryInfo> > puInfoHandle;
  iEvent.getByToken(puInfoToken, puInfoHandle);
  if (!puInfoHandle.isValid()) throw cms::Exception("CMS3Ntuplizer::fillEventVariables: Error getting the PU info from the event...");

  // Simple event-level variables
  commonEntry.setNamedVal("EventNumber", iEvent.id().event());
  commonEntry.setNamedVal("RunNumber", iEvent.id().run());
  commonEntry.setNamedVal("LuminosityBlock", iEvent.luminosityBlock());
  commonEntry.setNamedVal("event_rho", (float) (*rhoHandle));
  if (isMC){
    commonEntry.setNamedVal("n_vtxs_PU", (int) ((*puInfoHandle)[0].getPU_NumInteractions()));
    commonEntry.setNamedVal("n_true_int", (float) ((*puInfoHandle)[0].getTrueNumInteractions()));
  }
  else{
    commonEntry.setNamedVal("n_vtxs_PU", (int) 0);
    commonEntry.setNamedVal("n_true_int", (float) 0);
  }

  return true;
}
bool CMS3Ntuplizer::fillTriggerInfo(const edm::Event& iEvent){
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
  for (View<TriggerInfo>::const_iterator obj = triggerInfoHandle->begin(); obj != triggerInfoHandle->end(); obj++){
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
bool CMS3Ntuplizer::fillMETFilterVariables(const edm::Event& iEvent){
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
bool CMS3Ntuplizer::fillMETVariables(const edm::Event& iEvent){
#define SET_MET_VARIABLE(HANDLE, NAME, COLLNAME) commonEntry.setNamedVal((std::string(COLLNAME) + "_" + #NAME).data(), HANDLE->NAME);

  edm::Handle<METInfo> metHandle;

  // PF MET
  const char pfmetCollName[] = "pfmet";
  iEvent.getByToken(pfmetToken, metHandle);
  if (!metHandle.isValid()) throw cms::Exception("CMS3Ntuplizer::fillMETVariables: Error getting the PF MET handle from the event...");

  SET_MET_VARIABLE(metHandle, met, pfmetCollName);
  SET_MET_VARIABLE(metHandle, metPhi, pfmetCollName);
  SET_MET_VARIABLE(metHandle, sumEt, pfmetCollName);
  SET_MET_VARIABLE(metHandle, metSignificance, pfmetCollName);
  SET_MET_VARIABLE(metHandle, met_over_sqrtSumEt, pfmetCollName);

  SET_MET_VARIABLE(metHandle, met_raw, pfmetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_raw, pfmetCollName);
  SET_MET_VARIABLE(metHandle, sumEt_raw, pfmetCollName);

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

  SET_MET_VARIABLE(metHandle, gen_met, pfmetCollName);
  SET_MET_VARIABLE(metHandle, gen_metPhi, pfmetCollName);

  // PUPPI MET
  const char puppimetCollName[] = "puppimet";
  iEvent.getByToken(puppimetToken, metHandle);
  if (!metHandle.isValid()) throw cms::Exception("CMS3Ntuplizer::fillMETVariables: Error getting the PUPPI MET handle from the event...");

  SET_MET_VARIABLE(metHandle, met, puppimetCollName);
  SET_MET_VARIABLE(metHandle, metPhi, puppimetCollName);
  SET_MET_VARIABLE(metHandle, sumEt, puppimetCollName);
  SET_MET_VARIABLE(metHandle, metSignificance, puppimetCollName);
  SET_MET_VARIABLE(metHandle, met_over_sqrtSumEt, puppimetCollName);

  SET_MET_VARIABLE(metHandle, met_raw, puppimetCollName);
  SET_MET_VARIABLE(metHandle, metPhi_raw, puppimetCollName);
  SET_MET_VARIABLE(metHandle, sumEt_raw, puppimetCollName);

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

  SET_MET_VARIABLE(metHandle, calo_met, puppimetCollName);
  SET_MET_VARIABLE(metHandle, calo_metPhi, puppimetCollName);

  SET_MET_VARIABLE(metHandle, gen_met, puppimetCollName);
  SET_MET_VARIABLE(metHandle, gen_metPhi, puppimetCollName);


#undef SET_MET_VARIABLE

  return true;
}

bool CMS3Ntuplizer::fillGenVariables(const edm::Event& iEvent){
  if (!this->isMC) return true;

  // Gen info.
  edm::Handle< GenInfo > genInfo;
  iEvent.getByToken(genInfoToken, genInfo);
  if (!genInfo.isValid()) throw cms::Exception("CMS3Ntuplizer::fillGenVariables: Error getting the gen. info. from the event...");
  recordGenInfo(*genInfo);

  return true;
}


// Undefine the convenience macros
#undef PUSH_VECTOR_WITH_NAME
#undef PUSH_USERFLOAT_INTO_VECTOR
#undef PUSH_USERINT_INTO_VECTOR
#undef MAKE_VECTOR_WITH_RESERVE


DEFINE_FWK_MODULE(CMS3Ntuplizer);
