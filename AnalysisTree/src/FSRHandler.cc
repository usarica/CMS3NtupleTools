#include <cassert>
#include <algorithm>
#include <utility>

#include <CMS3/Dictionaries/interface/GlobalCollectionNames.h>

#include "ParticleObjectHelpers.h"
#include "FSRHandler.h"
#include "MuonSelectionHelpers.h"
#include "ElectronSelectionHelpers.h"
#include "ParticleSelectionHelpers.h"
#include <CMS3/Dictionaries/interface/CMS3StreamHelpers.h>


using namespace std;
using namespace IvyStreamHelpers;


#define VECTOR_ITERATOR_HANDLER_DIRECTIVES \
FSR_VARIABLE(float, pt, 0) \
FSR_VARIABLE(float, eta, 0) \
FSR_VARIABLE(float, phi, 0) \
FSR_VARIABLE(float, mass, 0) \
FSR_VARIABLES


const std::string FSRHandler::colName = GlobalCollectionNames::colName_fsrcands;

FSRHandler::FSRHandler() : IvyBase()
{
#define FSR_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(FSRHandler::colName + "_" + #NAME);
#define FSR_VECTOR_VARIABLE(TYPE, NAME) this->addConsumed<std::vector<TYPE>*>(FSRHandler::colName + "_" + #NAME);
  VECTOR_ITERATOR_HANDLER_DIRECTIVES;
#undef FSR_VECTOR_VARIABLE
#undef FSR_VARIABLE
}

void FSRHandler::clear(){
  this->resetCache();

  muons_postFSR.clear();
  electrons_postFSR.clear();
  photons_postFSR.clear();

  for (auto*& part:electrons_owned) delete part;
  electrons_owned.clear();
  for (auto*& part:muons_owned) delete part;
  muons_owned.clear();
  for (auto*& part:fsrCandidates) delete part;
  fsrCandidates.clear();
}

bool FSRHandler::constructPostFSRParticles(std::vector<MuonObject*> const* muons, std::vector<ElectronObject*> const* electrons, std::vector<PhotonObject*> const* photons, std::vector<PFCandidateObject*> const* pfcandidates){
  if (this->isAlreadyCached()) return true;

  clear();
  if (!currentTree) return false;

  if (!muons){
    if (this->verbosity>=MiscUtils::ERROR) IVYerr << "FSRHandler::constructPostFSRParticles: The input muon collection cannot be empty!" << endl;
    assert(0);
  }
  if (!electrons){
    if (this->verbosity>=MiscUtils::ERROR) IVYerr << "FSRHandler::constructPostFSRParticles: The input electron collection cannot be empty!" << endl;
    assert(0);
  }
  if (!photons){
    if (this->verbosity>=MiscUtils::ERROR) IVYerr << "FSRHandler::constructPostFSRParticles: The input photon collection cannot be empty!" << endl;
    assert(0);
  }
  bool res = constructFSRObjects() && associatePFCandidates(pfcandidates) && reconstructPostFSRObjects(muons, electrons, photons);

  if (res) this->cacheEvent();
  return res;
}

bool FSRHandler::constructFSRObjects(){
#define FSR_VARIABLE(TYPE, NAME, DEFVAL) std::vector<TYPE>::const_iterator itBegin_##NAME, itEnd_##NAME;
#define FSR_VECTOR_VARIABLE(TYPE, NAME) std::vector<TYPE>::const_iterator itBegin_##NAME, itEnd_##NAME;
  VECTOR_ITERATOR_HANDLER_DIRECTIVES;
#undef FSR_VECTOR_VARIABLE
#undef FSR_VARIABLE

  // Beyond this point starts checks and selection
  bool allVariablesPresent = true;
#define FSR_VARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumedCIterators<std::vector<TYPE>>(FSRHandler::colName + "_" + #NAME, &itBegin_##NAME, &itEnd_##NAME);
#define FSR_VECTOR_VARIABLE(TYPE, NAME) allVariablesPresent &= this->getConsumedCIterators<std::vector<TYPE>>(FSRHandler::colName + "_" + #NAME, &itBegin_##NAME, &itEnd_##NAME);
  VECTOR_ITERATOR_HANDLER_DIRECTIVES;
#undef FSR_VECTOR_VARIABLE
#undef FSR_VARIABLE

  if (!allVariablesPresent){
    if (this->verbosity>=MiscUtils::ERROR) IVYerr << "FSRHandler::constructFSRObjects: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=MiscUtils::DEBUG) IVYout << "FSRHandler::constructFSRObjects: All variables are set up!" << endl;

  if (itBegin_pt == itEnd_pt) return true; // Construction is successful, it is just that no photons exist.

  size_t nProducts = (itEnd_pt - itBegin_pt);
  fsrCandidates.reserve(nProducts);
#define FSR_VARIABLE(TYPE, NAME, DEFVAL) auto it_##NAME = itBegin_##NAME;
#define FSR_VECTOR_VARIABLE(TYPE, NAME) auto it_##NAME = itBegin_##NAME;
  VECTOR_ITERATOR_HANDLER_DIRECTIVES;
#undef FSR_VECTOR_VARIABLE
#undef FSR_VARIABLE
  {
    size_t ip=0;
    while (it_pt != itEnd_pt){
      if (this->verbosity>=MiscUtils::DEBUG) IVYout << "FSRHandler::constructFSRObjects: Attempting photon " << ip << "..." << endl;

      ParticleObject::LorentzVector_t momentum;
      momentum = ParticleObject::PolarLorentzVector_t(*it_pt, *it_eta, *it_phi, *it_mass); // Yes you have to do this on a separate line because CMSSW...
      fsrCandidates.push_back(new FSRObject(momentum));
      FSRObject*& obj = fsrCandidates.back();

      // Set extras
#define FSR_VARIABLE(TYPE, NAME, DEFVAL) obj->extras.NAME = *it_##NAME;
#define FSR_VECTOR_VARIABLE(TYPE, NAME) std::copy(it_##NAME->cbegin(), it_##NAME->cend(), obj->extras.NAME.begin());
      FSR_VARIABLES;
#undef FSR_VECTOR_VARIABLE
#undef FSR_VARIABLE

      // Set particle index as its unique identifier
      obj->setUniqueIdentifier(ip);

      if (this->verbosity>=MiscUtils::DEBUG) IVYout << "\t- Success!" << endl;

      ip++;
#define FSR_VARIABLE(TYPE, NAME, DEFVAL) it_##NAME++;
#define FSR_VECTOR_VARIABLE(TYPE, NAME) it_##NAME++;
      VECTOR_ITERATOR_HANDLER_DIRECTIVES;
#undef FSR_VECTOR_VARIABLE
#undef FSR_VARIABLE
    }
  }

  return true;
}

bool FSRHandler::associatePFCandidates(std::vector<PFCandidateObject*> const* pfcandidates) const{
  if (!pfcandidates) return true;

  for (auto const& pfcand:(*pfcandidates)){
    for (auto const& part:fsrCandidates){
      if (pfcand->extras.matched_FSRCandidate_index>=0 && (ParticleObject::UniqueId_t) pfcand->extras.matched_FSRCandidate_index == part->getUniqueIdentifier()) part->addDaughter(pfcand);
    }
  }

  return true;
}

bool FSRHandler::reconstructPostFSRObjects(std::vector<MuonObject*> const* muons, std::vector<ElectronObject*> const* electrons, std::vector<PhotonObject*> const* photons){
  if (fsrCandidates.empty()) return true; // No need to do anything if there are no FSR candidates to merge

  // Temporary struct for intermediate purposes
  struct fsr_lepton_matchSpecs{
    FSRObject* fsrObj;
    ParticleObject* lepton;
    ParticleObject::LorentzVector_t::Scalar dR;

    fsr_lepton_matchSpecs() :
      fsrObj(nullptr),
      lepton(nullptr),
      dR(9999)
    {}
    fsr_lepton_matchSpecs(FSRObject* fsrObj_, ParticleObject* lepton_) :
      fsrObj(fsrObj_),
      lepton(lepton_),
      dR((fsrObj && lepton ? ParticleObject::LorentzVector_t::Scalar(fsrObj->deltaR(lepton)) : ParticleObject::LorentzVector_t::Scalar(9999)))
    {}
    fsr_lepton_matchSpecs(fsr_lepton_matchSpecs const& other) :
      fsrObj(other.fsrObj),
      lepton(other.lepton),
      dR(other.dR)
    {}

    bool operator == (fsr_lepton_matchSpecs const& other) const{ return (dR == other.dR); }
    bool operator > (fsr_lepton_matchSpecs const& other) const{ return (dR > other.dR); }
    bool operator < (fsr_lepton_matchSpecs const& other) const{ return (dR < other.dR); }
    bool operator >= (fsr_lepton_matchSpecs const& other) const{ return (*this > other || *this == other); }
    bool operator <= (fsr_lepton_matchSpecs const& other) const{ return (*this < other || *this == other); }
  };

  muons_owned.reserve(fsrCandidates.size());
  electrons_owned.reserve(fsrCandidates.size());

  if (muons) muons_postFSR.reserve(muons->size());
  if (electrons) electrons_postFSR.reserve(electrons->size());
  if (photons) photons_postFSR.reserve(photons->size());

  // Final map to be used
  std::vector<fsr_lepton_matchSpecs> fsr_lepton_final_map;
  {
    // Intermediate map to collect and sort
    std::vector<fsr_lepton_matchSpecs> fsr_lepton_match_list; fsr_lepton_match_list.reserve(fsrCandidates.size()*((muons ? muons->size() : 0)+(electrons ? electrons->size() : 0)));
    if (muons){
      for (auto const& part:(*muons)){
        if (ParticleSelectionHelpers::isFSRSuitable(part)){
          auto const& uid = part->getUniqueIdentifier();
          for (auto const& fsrCand:fsrCandidates){
            if (HelperFunctions::checkListVariable(fsrCand->extras.fsrMatch_muon_index_list, uid)) HelperFunctions::addByLowest(fsr_lepton_match_list, fsr_lepton_matchSpecs(fsrCand, part), false);
          }
        }
      }
    }
    if (electrons){
      for (auto const& part:(*electrons)){
        if (ParticleSelectionHelpers::isFSRSuitable(part)){
          auto const& uid = part->getUniqueIdentifier();
          for (auto const& fsrCand:fsrCandidates){
            if (HelperFunctions::checkListVariable(fsrCand->extras.fsrMatch_electron_index_list, uid)) HelperFunctions::addByLowest(fsr_lepton_match_list, fsr_lepton_matchSpecs(fsrCand, part), false);
          }
        }
      }
    }

    // Ensure that there is a 1-1 matching by smallest dR
    fsr_lepton_final_map.reserve(fsr_lepton_match_list.size());
    for (auto const& fsr_lepton_match:fsr_lepton_match_list){
      bool doSkip = false;
      for (auto const& fsr_lepton_final_match:fsr_lepton_final_map){
        if (fsr_lepton_match.fsrObj == fsr_lepton_final_match.fsrObj || fsr_lepton_match.lepton == fsr_lepton_final_match.lepton){ doSkip=true; break; }
      }
      if (!doSkip) fsr_lepton_final_map.emplace_back(fsr_lepton_match);
    }
  }

  // Add unpaired muons and electrons to the post-FSR collections without cloning them
  if (muons){
    for (auto const& part:(*muons)){
      bool isFound=false;
      for (auto const& fsr_lepton_match:fsr_lepton_final_map){ if (fsr_lepton_match.lepton == part){ isFound=true; break; } }
      if (!isFound) muons_postFSR.push_back(part);
    }
  }
  if (electrons){
    for (auto const& part:(*electrons)){
      bool isFound=false;
      for (auto const& fsr_lepton_match:fsr_lepton_final_map){ if (fsr_lepton_match.lepton == part){ isFound=true; break; } }
      if (!isFound) electrons_postFSR.push_back(part);
    }
  }
  // Clean photons from FSR candidates
  if (photons){
    for (auto const& part:(*photons)){
      auto const& uid = part->getUniqueIdentifier();
      bool isRejected=false;
      for (auto const& fsr_lepton_match:fsr_lepton_final_map){
        auto const& fsrObj = fsr_lepton_match.fsrObj;
        if (HelperFunctions::checkListVariable(fsrObj->extras.photonVeto_index_list, uid)){
          float dR = fsrObj->deltaR(part);
          isRejected = (
            part->extras.n_associated_pfphotons>=1
            &&
            (
              dR<0.1f
              ||
              (fsrObj->px() == part->extras.closestPFPhoton_associated_px && fsrObj->py() == part->extras.closestPFPhoton_associated_py)
              )
            );
          if (isRejected) fsrObj->setAssociatedPhoton(part);
        }
      }
      if (!isRejected) photons_postFSR.push_back(part);
    }
  }

  // Now reconstruct dressed leptons by cloning undressed leptons and applying corrections and selection
  for (auto const& fsr_lepton_match:fsr_lepton_final_map){
    auto const& fsrObj = fsr_lepton_match.fsrObj;
    auto const& lepton = fsr_lepton_match.lepton;
    auto const& dR_fsr_lepton = fsr_lepton_match.dR;

    bool isMuon = std::abs(lepton->pdgId())==13; // Either a muon or an electron
    if (isMuon){
      MuonObject* lepton_postFSR = new MuonObject(*(dynamic_cast<MuonObject const*>(lepton)));
      lepton_postFSR->resetMothers();
      lepton_postFSR->resetDaughters();

      // Apply object corrections
      lepton_postFSR->applyFSRIsoCorr(dR_fsr_lepton, fsrObj->pt());

      // Add p4 of the fsrObj
      lepton_postFSR->p4() += fsrObj->p4();

      // Add the original lepton and the FSR cand. as the daughters of the dressed lepton
      lepton_postFSR->addDaughter(lepton); lepton->addMother(lepton_postFSR);
      lepton_postFSR->addDaughter(fsrObj); fsrObj->addMother(lepton_postFSR);

      // Set selection flags
      MuonSelectionHelpers::setSelectionBits(*lepton_postFSR);

      muons_owned.emplace_back(lepton_postFSR);
      muons_postFSR.push_back(lepton_postFSR);
    }
    else{
      ElectronObject* lepton_postFSR = new ElectronObject(*(dynamic_cast<ElectronObject const*>(lepton)));
      lepton_postFSR->resetMothers();
      lepton_postFSR->resetDaughters();

      // Apply object corrections
      lepton_postFSR->applyFSRIsoCorr(dR_fsr_lepton, fsrObj->pt());

      // Add p4 of the fsrObj
      lepton_postFSR->p4() += fsrObj->p4();

      // Add the original lepton and the FSR cand. as the daughters of the dressed lepton
      lepton_postFSR->addDaughter(lepton); lepton->addMother(lepton_postFSR);
      lepton_postFSR->addDaughter(fsrObj); fsrObj->addMother(lepton_postFSR);

      // Set selection flags
      ElectronSelectionHelpers::setSelectionBits(*lepton_postFSR);

      electrons_owned.emplace_back(lepton_postFSR);
      electrons_postFSR.push_back(lepton_postFSR);
    }
  }

  // Sort particles
  ParticleObjectHelpers::sortByGreaterPt(muons_owned);
  ParticleObjectHelpers::sortByGreaterPt(muons_postFSR);
  ParticleObjectHelpers::sortByGreaterPt(electrons_owned);
  ParticleObjectHelpers::sortByGreaterPt(electrons_postFSR);
  // No need to sort photons: Rejecting some does not change the order.

  // Check collection sizes
  if (muons) assert(muons->size() == muons_postFSR.size());
  if (electrons) assert(electrons->size() == electrons_postFSR.size());
  if (photons) assert(photons->size() >= photons_postFSR.size());

  return true;
}


void FSRHandler::bookBranches(BaseTree* tree){
  if (!tree) return;

#define FSR_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<std::vector<TYPE>*>(FSRHandler::colName + "_" + #NAME, nullptr);
#define FSR_VECTOR_VARIABLE(TYPE, NAME) tree->bookBranch<std::vector<TYPE>*>(FSRHandler::colName + "_" + #NAME, nullptr);
  VECTOR_ITERATOR_HANDLER_DIRECTIVES
#undef FSR_VECTOR_VARIABLE
#undef FSR_VARIABLE
}


#undef VECTOR_ITERATOR_HANDLER_DIRECTIVES
