#include <cassert>
#include "ParticleObjectHelpers.h"
#include "PhotonHandler.h"
#include "PhotonSelectionHelpers.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


#define VECTOR_ITERATOR_HANDLER_DIRECTIVES \
PHOTON_VARIABLE(float, pt, 0) \
PHOTON_VARIABLE(float, eta, 0) \
PHOTON_VARIABLE(float, phi, 0) \
PHOTON_VARIABLE(float, mass, 0) \
PHOTON_VARIABLES


const std::string PhotonHandler::colName = "photons";

PhotonHandler::PhotonHandler() : IvyBase()
{
#define PHOTON_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(PhotonHandler::colName + "_" + #NAME);
  VECTOR_ITERATOR_HANDLER_DIRECTIVES;
#undef PHOTON_VARIABLE
}


bool PhotonHandler::constructPhotons(SystematicsHelpers::SystematicVariationTypes const& syst){
  clear();
  if (!currentTree) return false;

#define PHOTON_VARIABLE(TYPE, NAME, DEFVAL) std::vector<TYPE>::const_iterator itBegin_##NAME, itEnd_##NAME;
  VECTOR_ITERATOR_HANDLER_DIRECTIVES;
#undef PHOTON_VARIABLE

    // Beyond this point starts checks and selection
  bool allVariablesPresent = true;
#define PHOTON_VARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumedCIterators<std::vector<TYPE>>(PhotonHandler::colName + "_" + #NAME, &itBegin_##NAME, &itEnd_##NAME);
  VECTOR_ITERATOR_HANDLER_DIRECTIVES;
#undef PHOTON_VARIABLE

  if (!allVariablesPresent){
    if (this->verbosity>=TVar::ERROR) MELAerr << "PhotonHandler::constructPhotons: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=TVar::DEBUG) MELAout << "PhotonHandler::constructPhotons: All variables are set up!" << endl;

  if (itBegin_pt == itEnd_pt) return true; // Construction is successful, it is just that no photons exist.

  size_t nProducts = (itEnd_pt - itBegin_pt);
  productList.reserve(nProducts);
#define PHOTON_VARIABLE(TYPE, NAME, DEFVAL) auto it_##NAME = itBegin_##NAME;
  VECTOR_ITERATOR_HANDLER_DIRECTIVES;
#undef PHOTON_VARIABLE
  {
    size_t ip=0;
    while (it_pt != itEnd_pt){
      if (this->verbosity>=TVar::DEBUG) MELAout << "PhotonHandler::constructPhotons: Attempting photon " << ip << "..." << endl;

      ParticleObject::LorentzVector_t momentum;
      momentum = ParticleObject::PolarLorentzVector_t(*it_pt, *it_eta, *it_phi, *it_mass); // Yes you have to do this on a separate line because CMSSW...
      productList.push_back(new PhotonObject(momentum));
      PhotonObject*& obj = productList.back();

      // Set extras
#define PHOTON_VARIABLE(TYPE, NAME, DEFVAL) obj->extras.NAME = *it_##NAME;
      PHOTON_VARIABLES;
#undef PHOTON_VARIABLE

      // Replace momentum
      obj->makeFinalMomentum(syst);

      // Set the selection bits
      PhotonSelectionHelpers::setSelectionBits(*obj);

      if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Success!" << endl;

      ip++;
#define PHOTON_VARIABLE(TYPE, NAME, DEFVAL) it_##NAME++;
      VECTOR_ITERATOR_HANDLER_DIRECTIVES;
#undef PHOTON_VARIABLE
    }
  }
    // Sort particles
  ParticleObjectHelpers::sortByGreaterPt(productList);

  return true;
}

void PhotonHandler::bookBranches(BaseTree* tree){
  if (!tree) return;

#define PHOTON_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<std::vector<TYPE>*>(PhotonHandler::colName + "_" + #NAME, nullptr);
  VECTOR_ITERATOR_HANDLER_DIRECTIVES
#undef PHOTON_VARIABLE
}


#undef VECTOR_ITERATOR_HANDLER_DIRECTIVES
