#include <cassert>

#include <CMS3/Dictionaries/interface/GlobalCollectionNames.h>

#include "ParticleObjectHelpers.h"
#include "SuperclusterHandler.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


#define VECTOR_ITERATOR_HANDLER_DIRECTIVES \
SUPERCLUSTER_VARIABLE(float, correctedEnergy, 0) \
SUPERCLUSTER_VARIABLE(float, eta, 0) \
SUPERCLUSTER_VARIABLE(float, phi, 0) \
SUPERCLUSTER_VARIABLES


const std::string SuperclusterHandler::colName = GlobalCollectionNames::colName_superclusters;

SuperclusterHandler::SuperclusterHandler() : IvyBase()
{
#define SUPERCLUSTER_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(SuperclusterHandler::colName + "_" + #NAME);
  VECTOR_ITERATOR_HANDLER_DIRECTIVES;
#undef SUPERCLUSTER_VARIABLE
}


bool SuperclusterHandler::constructSuperclusters(SystematicsHelpers::SystematicVariationTypes const& syst){
  if (this->isAlreadyCached()) return true;

  clear();
  if (!currentTree) return false;

#define SUPERCLUSTER_VARIABLE(TYPE, NAME, DEFVAL) std::vector<TYPE>::const_iterator itBegin_##NAME, itEnd_##NAME;
  VECTOR_ITERATOR_HANDLER_DIRECTIVES;
#undef SUPERCLUSTER_VARIABLE

    // Beyond this point starts checks and selection
  bool allVariablesPresent = true;
#define SUPERCLUSTER_VARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumedCIterators<std::vector<TYPE>>(SuperclusterHandler::colName + "_" + #NAME, &itBegin_##NAME, &itEnd_##NAME);
  VECTOR_ITERATOR_HANDLER_DIRECTIVES;
#undef SUPERCLUSTER_VARIABLE

  if (!allVariablesPresent){
    if (this->verbosity>=TVar::ERROR) MELAerr << "SuperclusterHandler::constructSuperclusters: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=TVar::DEBUG) MELAout << "SuperclusterHandler::constructSuperclusters: All variables are set up!" << endl;

  if (itBegin_correctedEnergy == itEnd_correctedEnergy) return true; // Construction is successful, it is just that no photons exist.

  size_t nProducts = (itEnd_correctedEnergy - itBegin_correctedEnergy);
  productList.reserve(nProducts);
#define SUPERCLUSTER_VARIABLE(TYPE, NAME, DEFVAL) auto it_##NAME = itBegin_##NAME;
  VECTOR_ITERATOR_HANDLER_DIRECTIVES;
#undef SUPERCLUSTER_VARIABLE
  {
    size_t ip=0;
    while (it_correctedEnergy != itEnd_correctedEnergy){
      if (this->verbosity>=TVar::DEBUG) MELAout << "SuperclusterHandler::constructSuperclusters: Attempting supercluster " << ip << "..." << endl;

      ParticleObject::LorentzVector_t momentum;
      momentum = ParticleObject::PolarLorentzVector_t(*it_correctedEnergy/std::cosh(*it_eta), *it_eta, *it_phi, 0.); // Yes you have to do this on a separate line because CMSSW...
      productList.push_back(new ProductType_t(momentum));
      SuperclusterObject*& obj = productList.back();

      // Set extras
#define SUPERCLUSTER_VARIABLE(TYPE, NAME, DEFVAL) obj->extras.NAME = *it_##NAME;
      SUPERCLUSTER_VARIABLES;
#undef SUPERCLUSTER_VARIABLE

      // Set particle index as its unique identifier
      obj->setUniqueIdentifier(ip);

      // Replace momentum
      obj->makeFinalMomentum(syst);

      if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Success!" << endl;

      ip++;
#define SUPERCLUSTER_VARIABLE(TYPE, NAME, DEFVAL) it_##NAME++;
      VECTOR_ITERATOR_HANDLER_DIRECTIVES;
#undef SUPERCLUSTER_VARIABLE
    }
  }
    // Sort particles
  ParticleObjectHelpers::sortByGreaterPt(productList);

  this->cacheEvent();
  return true;
}

void SuperclusterHandler::bookBranches(BaseTree* tree){
  if (!tree) return;

#define SUPERCLUSTER_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<std::vector<TYPE>*>(SuperclusterHandler::colName + "_" + #NAME, nullptr);
  VECTOR_ITERATOR_HANDLER_DIRECTIVES
#undef SUPERCLUSTER_VARIABLE
}


#undef VECTOR_ITERATOR_HANDLER_DIRECTIVES
