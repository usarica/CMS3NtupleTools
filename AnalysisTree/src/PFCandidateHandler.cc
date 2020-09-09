#include <cassert>

#include <CMS3/Dictionaries/interface/GlobalCollectionNames.h>

#include "ParticleObjectHelpers.h"
#include "PFCandidateHandler.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


#define PFCANDIDATE_MOMENTUM_VARIABLES \
PFCANDIDATE_VARIABLE(float, pt, 0) \
PFCANDIDATE_VARIABLE(float, eta, 0) \
PFCANDIDATE_VARIABLE(float, phi, 0) \
PFCANDIDATE_VARIABLE(float, mass, 0) \
PFCANDIDATE_VARIABLE(cms3_id_t, id, 0)
#define PFCANDIDATE_DIRECTIVES \
PFCANDIDATE_MOMENTUM_VARIABLES \
PFCANDIDATE_VARIABLES


const std::string PFCandidateHandler::colName = GlobalCollectionNames::colName_pfcands;

PFCandidateHandler::PFCandidateHandler() :
  IvyBase()
{
#define PFCANDIDATE_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(PFCandidateHandler::colName + "_" + #NAME);
  PFCANDIDATE_DIRECTIVES;
#undef PFCANDIDATE_VARIABLE
}


bool PFCandidateHandler::constructPFCandidates(SystematicsHelpers::SystematicVariationTypes const& syst){
  if (this->isAlreadyCached()) return true;

  clear();
  if (!currentTree) return false;

#define PFCANDIDATE_VARIABLE(TYPE, NAME, DEFVAL) std::vector<TYPE>::const_iterator itBegin_##NAME, itEnd_##NAME;
  PFCANDIDATE_DIRECTIVES;
#undef PFCANDIDATE_VARIABLE

    // Beyond this point starts checks and selection
  bool allVariablesPresent = true;
#define PFCANDIDATE_VARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumedCIterators<std::vector<TYPE>>(PFCandidateHandler::colName + "_" + #NAME, &itBegin_##NAME, &itEnd_##NAME);
  PFCANDIDATE_DIRECTIVES;
#undef PFCANDIDATE_VARIABLE

  if (!allVariablesPresent){
    if (this->verbosity>=TVar::ERROR) MELAerr << "PFCandidateHandler::constructPFCandidates: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=TVar::DEBUG) MELAout << "PFCandidateHandler::constructPFCandidates: All variables are set up!" << endl;

  size_t n_products = (itEnd_id - itBegin_id);
  productList.reserve(n_products);
#define PFCANDIDATE_VARIABLE(TYPE, NAME, DEFVAL) auto it_##NAME = itBegin_##NAME;
  PFCANDIDATE_DIRECTIVES;
#undef PFCANDIDATE_VARIABLE
  {
    size_t ip = 0;
    while (ip != n_products){
      if (this->verbosity>=TVar::DEBUG) MELAout << "PFCandidateHandler::constructPFCandidates: Attempting muon " << ip << "..." << endl;

      ParticleObject::LorentzVector_t momentum;
      momentum = ParticleObject::PolarLorentzVector_t(*it_pt, *it_eta, *it_phi, *it_mass); // Yes you have to do this on a separate line because CMSSW...
      productList.push_back(new ProductType_t(*it_id, momentum));
      ProductType_t*& obj = productList.back();

      // Set extras
#define PFCANDIDATE_VARIABLE(TYPE, NAME, DEFVAL) obj->extras.NAME = *it_##NAME;
      PFCANDIDATE_VARIABLES;
#undef PFCANDIDATE_VARIABLE

      // Set particle index as its unique identifier
      obj->setUniqueIdentifier(ip);

      // Replace momentum
      obj->makeFinalMomentum(syst);

      if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Success!" << endl;

      ip++;
#define PFCANDIDATE_VARIABLE(TYPE, NAME, DEFVAL) it_##NAME++;
      PFCANDIDATE_DIRECTIVES;
#undef PFCANDIDATE_VARIABLE
    }
  }
  // Sort particles
  ParticleObjectHelpers::sortByGreaterPt(productList);

  this->cacheEvent();
  return true;
}

void PFCandidateHandler::bookBranches(BaseTree* tree){
  if (!tree) return;

#define PFCANDIDATE_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<std::vector<TYPE>*>(PFCandidateHandler::colName + "_" + #NAME, nullptr);
  PFCANDIDATE_DIRECTIVES;
#undef PFCANDIDATE_VARIABLE
}


#undef PFCANDIDATE_DIRECTIVES
