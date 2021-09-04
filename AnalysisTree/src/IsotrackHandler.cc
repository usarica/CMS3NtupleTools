#include <cassert>

#include <CMS3/Dictionaries/interface/GlobalCollectionNames.h>

#include "ParticleObjectHelpers.h"
#include "IsotrackHandler.h"
#include "IsotrackSelectionHelpers.h"
#include "MuonSelectionHelpers.h"
#include "ElectronSelectionHelpers.h"
#include "ParticleSelectionHelpers.h"
#include <CMS3/Dictionaries/interface/CMS3StreamHelpers.h>


using namespace std;
using namespace IvyStreamHelpers;


#define VECTOR_ITERATOR_HANDLER_DIRECTIVES \
ISOTRACK_VARIABLE(float, pt, 0) \
ISOTRACK_VARIABLE(float, eta, 0) \
ISOTRACK_VARIABLE(float, phi, 0) \
ISOTRACK_VARIABLE(float, mass, 0) \
ISOTRACK_VARIABLE(cms3_id_t, id, 0) \
ISOTRACK_VARIABLES


const std::string IsotrackHandler::colName = GlobalCollectionNames::colName_isotracks;

IsotrackHandler::IsotrackHandler() : IvyBase()
{
#define ISOTRACK_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(IsotrackHandler::colName + "_" + #NAME);
  VECTOR_ITERATOR_HANDLER_DIRECTIVES;
#undef ISOTRACK_VARIABLE
}

bool IsotrackHandler::constructIsotracks(std::vector<MuonObject*> const* muons, std::vector<ElectronObject*> const* electrons){
  if (this->isAlreadyCached()) return true;

  clear();
  if (!currentTree) return false;

#define ISOTRACK_VARIABLE(TYPE, NAME, DEFVAL) std::vector<TYPE>::const_iterator itBegin_##NAME, itEnd_##NAME;
  VECTOR_ITERATOR_HANDLER_DIRECTIVES;
#undef ISOTRACK_VARIABLE

    // Beyond this point starts checks and selection
  bool allVariablesPresent = true;
#define ISOTRACK_VARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumedCIterators<std::vector<TYPE>>(IsotrackHandler::colName + "_" + #NAME, &itBegin_##NAME, &itEnd_##NAME);
  VECTOR_ITERATOR_HANDLER_DIRECTIVES;
#undef ISOTRACK_VARIABLE

  if (!allVariablesPresent){
    if (this->verbosity>=MiscUtils::ERROR) IVYerr << "IsotrackHandler::constructIsotracks: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=MiscUtils::DEBUG) IVYout << "IsotrackHandler::constructIsotracks: All variables are set up!" << endl;

  if (itBegin_id == itEnd_id) return true; // Construction is successful, it is just that no isotracks exist.

  size_t nProducts = (itEnd_id - itBegin_id);
  productList.reserve(nProducts);
#define ISOTRACK_VARIABLE(TYPE, NAME, DEFVAL) auto it_##NAME = itBegin_##NAME;
  VECTOR_ITERATOR_HANDLER_DIRECTIVES;
#undef ISOTRACK_VARIABLE
  {
    size_t ip=0;
    while (it_id != itEnd_id){
      if (this->verbosity>=MiscUtils::DEBUG) IVYout << "IsotrackHandler::constructIsotracks: Attempting isotrack " << ip << "..." << endl;

      ParticleObject::LorentzVector_t momentum;
      momentum = ParticleObject::PolarLorentzVector_t(*it_pt, *it_eta, *it_phi, *it_mass); // Yes you have to do this on a separate line because CMSSW...
      productList.push_back(new ProductType_t(*it_id, momentum));
      ProductType_t*& obj = productList.back();

      // Set extras
#define ISOTRACK_VARIABLE(TYPE, NAME, DEFVAL) obj->extras.NAME = *it_##NAME;
      ISOTRACK_VARIABLES;
#undef ISOTRACK_VARIABLE

      // Set the selection bits
      IsotrackSelectionHelpers::setSelectionBits(*obj);

      if (this->verbosity>=MiscUtils::DEBUG) IVYout << "\t- Success!" << endl;

      ip++;
#define ISOTRACK_VARIABLE(TYPE, NAME, DEFVAL) it_##NAME++;
      VECTOR_ITERATOR_HANDLER_DIRECTIVES;
#undef ISOTRACK_VARIABLE
    }
  }
    // Sort particles
  ParticleObjectHelpers::sortByGreaterPt(productList);

  bool res = applyCleaning(muons, electrons);

  if (res) this->cacheEvent();
  return res;
}

bool IsotrackHandler::applyCleaning(std::vector<MuonObject*> const* muons, std::vector<ElectronObject*> const* electrons){
  std::vector<ProductType_t*> productList_new; productList_new.reserve(productList.size());
  for (auto*& product:productList){
    bool doSkip=false;
    float const dR_isotrack = IsotrackSelectionHelpers::getIsolationDRmax(*product);
    if (muons){
      for (auto const* part:*(muons)){
        if (!ParticleSelectionHelpers::isParticleForIsotrackCleaning(part)) continue;
        float const separation_deltaR = std::max(dR_isotrack, MuonSelectionHelpers::getIsolationDRmax(*part));
        if (product->deltaR(part)<separation_deltaR){ doSkip=true; break; }
      }
    }
    if (electrons){
      for (auto const* part:*(electrons)){
        if (!ParticleSelectionHelpers::isParticleForIsotrackCleaning(part)) continue;
        float const separation_deltaR = std::max(dR_isotrack, ElectronSelectionHelpers::getIsolationDRmax(*part));
        if (product->deltaR(part)<separation_deltaR){ doSkip=true; break; }
      }
    }
    if (!doSkip) productList_new.push_back(product);
    else delete product;
  }
  productList = productList_new;

  return true;
}

void IsotrackHandler::bookBranches(BaseTree* tree){
  if (!tree) return;

#define ISOTRACK_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<std::vector<TYPE>*>(IsotrackHandler::colName + "_" + #NAME, nullptr);
  VECTOR_ITERATOR_HANDLER_DIRECTIVES
#undef ISOTRACK_VARIABLE
}


#undef VECTOR_ITERATOR_HANDLER_DIRECTIVES
