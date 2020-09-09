#include <cassert>

#include <CMS3/Dictionaries/interface/GlobalCollectionNames.h>

#include "ParticleObjectHelpers.h"
#include "PhotonHandler.h"
#include "PhotonSelectionHelpers.h"
#include "ParticleSelectionHelpers.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


#define PHOTON_MOMENTUM_VARIABLES \
PHOTON_VARIABLE(float, pt, 0) \
PHOTON_VARIABLE(float, eta, 0) \
PHOTON_VARIABLE(float, phi, 0) \
PHOTON_VARIABLE(float, mass, 0)


const std::string PhotonHandler::colName = GlobalCollectionNames::colName_photons;

PhotonHandler::PhotonHandler() :
  IvyBase(),
  has_genmatching(false),
  hasOverlapMaps(false),
  overlapMap_photons_ak4jets(nullptr),
  overlapMap_photons_ak8jets(nullptr)
{
#define PHOTON_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(PhotonHandler::colName + "_" + #NAME);
  PHOTON_MOMENTUM_VARIABLES;
  PHOTON_RECO_VARIABLES;
#undef PHOTON_VARIABLE
}


bool PhotonHandler::constructPhotons(SystematicsHelpers::SystematicVariationTypes const& syst, std::vector<PFCandidateObject*> const* pfcandidates){
  if (this->isAlreadyCached()) return true;

  clear();
  if (!currentTree) return false;

  bool res = (constructPhotonObjects(syst) && associatePFCandidates(pfcandidates) && linkOverlapElements());

  if (res) this->cacheEvent();
  return res;
}

bool PhotonHandler::associatePFCandidates(std::vector<PFCandidateObject*> const* pfcandidates) const{
  if (!pfcandidates) return true;

  for (auto const& pfcand:(*pfcandidates)){
    auto const& associated_particle_indices = pfcand->extras.matched_photon_index_list;
    for (auto const& part:productList){
      if (HelperFunctions::checkListVariable(associated_particle_indices, part->getUniqueIdentifier())){
        part->addDaughter(pfcand);
        pfcand->addMother(part);
      }
    }
  }

  return true;
}

bool PhotonHandler::linkOverlapElements() const{
  if (!hasOverlapMaps) return true;

  overlapMap_photons_ak4jets->constructOverlapMaps();
  for (auto const& ome:overlapMap_photons_ak4jets->getProducts()) ome->linkFirstElement(productList);

  overlapMap_photons_ak8jets->constructOverlapMaps();
  for (auto const& ome:overlapMap_photons_ak8jets->getProducts()) ome->linkFirstElement(productList);

  return true;
}

bool PhotonHandler::constructPhotonObjects(SystematicsHelpers::SystematicVariationTypes const& syst){
#define PHOTON_VARIABLE(TYPE, NAME, DEFVAL) std::vector<TYPE>::const_iterator itBegin_##NAME, itEnd_##NAME;
  PHOTON_MOMENTUM_VARIABLES;
  PHOTON_VARIABLES;
#undef PHOTON_VARIABLE

    // Beyond this point starts checks and selection
  bool allVariablesPresent = true;
#define PHOTON_VARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumedCIterators<std::vector<TYPE>>(PhotonHandler::colName + "_" + #NAME, &itBegin_##NAME, &itEnd_##NAME);
  PHOTON_MOMENTUM_VARIABLES;
  PHOTON_RECO_VARIABLES;
  if (this->has_genmatching){
    PHOTON_GENINFO_VARIABLES;
  }
#undef PHOTON_VARIABLE

  if (!allVariablesPresent){
    if (this->verbosity>=TVar::ERROR) MELAerr << "PhotonHandler::constructPhotonObjects: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=TVar::DEBUG) MELAout << "PhotonHandler::constructPhotonObjects: All variables are set up!" << endl;

  if (itBegin_pt == itEnd_pt) return true; // Construction is successful, it is just that no photons exist.

  size_t nProducts = (itEnd_pt - itBegin_pt);
  productList.reserve(nProducts);
#define PHOTON_VARIABLE(TYPE, NAME, DEFVAL) auto it_##NAME = itBegin_##NAME;
  PHOTON_MOMENTUM_VARIABLES;
  PHOTON_VARIABLES;
#undef PHOTON_VARIABLE
  {
    size_t ip=0;
    while (it_pt != itEnd_pt){
      if (this->verbosity>=TVar::DEBUG) MELAout << "PhotonHandler::constructPhotonObjects: Attempting photon " << ip << "..." << endl;

      ParticleObject::LorentzVector_t momentum;
      momentum = ParticleObject::PolarLorentzVector_t(*it_pt, *it_eta, *it_phi, *it_mass); // Yes you have to do this on a separate line because CMSSW...
      productList.push_back(new ProductType_t(momentum));
      PhotonObject*& obj = productList.back();

      // Set extras
#define PHOTON_VARIABLE(TYPE, NAME, DEFVAL) obj->extras.NAME = *it_##NAME;
      PHOTON_RECO_VARIABLES;
      if (this->has_genmatching){
        PHOTON_GENINFO_VARIABLES;
      }
#undef PHOTON_VARIABLE

      // Set particle index as its unique identifier
      obj->setUniqueIdentifier(ip);

      // Replace momentum
      obj->makeFinalMomentum(syst);

      // Set the selection bits
      PhotonSelectionHelpers::setSelectionBits(*obj);

      if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Success!" << endl;

      ip++;
#define PHOTON_VARIABLE(TYPE, NAME, DEFVAL) it_##NAME++;
      PHOTON_MOMENTUM_VARIABLES;
      PHOTON_RECO_VARIABLES;
      if (this->has_genmatching){
        PHOTON_GENINFO_VARIABLES;
      }
#undef PHOTON_VARIABLE
    }
  }
    // Sort particles
  ParticleObjectHelpers::sortByGreaterPt(productList);

  return true;
}

void PhotonHandler::checkOptionalInfo(BaseTree* tree, bool& flag_genmatching){
  flag_genmatching = true;

  std::vector<TString> bnames;
  tree->getValidBranchNamesWithoutAlias(bnames, false);

#define PHOTON_VARIABLE(TYPE, NAME, DEFVAL) flag_genmatching &= (std::find(bnames.cbegin(), bnames.cend(), PhotonHandler::colName + "_" + #NAME)!=bnames.cend());
  PHOTON_GENINFO_VARIABLES;
#undef PHOTON_VARIABLE
}

bool PhotonHandler::wrapTree(BaseTree* tree){
  if (!tree) return false;

  PhotonHandler::checkOptionalInfo(tree, this->has_genmatching);

  return IvyBase::wrapTree(tree);
}

void PhotonHandler::bookBranches(BaseTree* tree){
  if (!tree) return;

  PhotonHandler::checkOptionalInfo(tree, this->has_genmatching);
#define PHOTON_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(PhotonHandler::colName + "_" + #NAME); this->defineConsumedSloppy(#NAME);
  if (this->has_genmatching){
    PHOTON_GENINFO_VARIABLES;
  }
#undef PHOTON_VARIABLE

#define PHOTON_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<std::vector<TYPE>*>(PhotonHandler::colName + "_" + #NAME, nullptr);
  PHOTON_MOMENTUM_VARIABLES;
  PHOTON_RECO_VARIABLES;
  if (this->has_genmatching){
    PHOTON_GENINFO_VARIABLES;
  }
#undef PHOTON_VARIABLE
}

void PhotonHandler::registerOverlapMaps(
  OverlapMapHandler<PhotonObject, AK4JetObject>& overlapMap_photons_ak4jets_,
  OverlapMapHandler<PhotonObject, AK8JetObject>& overlapMap_photons_ak8jets_
){
  overlapMap_photons_ak4jets = &overlapMap_photons_ak4jets_;
  overlapMap_photons_ak8jets = &overlapMap_photons_ak8jets_;
  hasOverlapMaps = true;
}


#undef VECTOR_ITERATOR_HANDLER_DIRECTIVES
