#include <cassert>
#include "ParticleObjectHelpers.h"
#include "ElectronHandler.h"
#include "SamplesCore.h"
#include "ElectronSelectionHelpers.h"
#include "ParticleSelectionHelpers.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


#define ELECTRON_MOMENTUM_VARIABLES \
ELECTRON_VARIABLE(float, pt, 0) \
ELECTRON_VARIABLE(float, eta, 0) \
ELECTRON_VARIABLE(float, phi, 0) \
ELECTRON_VARIABLE(float, mass, 0) \
ELECTRON_VARIABLE(cms3_charge_t, charge, 0)


const std::string ElectronHandler::colName = "electrons";

ElectronHandler::ElectronHandler() :
  IvyBase(),
  has_mvaid_extras(false)
{
#define ELECTRON_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(ElectronHandler::colName + "_" + #NAME);
  ELECTRON_MOMENTUM_VARIABLES;
  ELECTRON_COMMON_VARIABLES;
#undef ELECTRON_VARIABLE
}


bool ElectronHandler::constructElectrons(SystematicsHelpers::SystematicVariationTypes const& syst){
  clear();
  if (!currentTree) return false;

#define ELECTRON_VARIABLE(TYPE, NAME, DEFVAL) std::vector<TYPE>::const_iterator itBegin_##NAME, itEnd_##NAME;
  ELECTRON_MOMENTUM_VARIABLES;
  ELECTRON_VARIABLES;
#undef ELECTRON_VARIABLE

    // Beyond this point starts checks and selection
  bool allVariablesPresent = true;
#define ELECTRON_VARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumedCIterators<std::vector<TYPE>>(ElectronHandler::colName + "_" + #NAME, &itBegin_##NAME, &itEnd_##NAME);
  ELECTRON_MOMENTUM_VARIABLES;
  ELECTRON_COMMON_VARIABLES;
  if (this->has_mvaid_extras){
    ELECTRON_MVAID_EXTRA_VARIABLES;
  }
#undef ELECTRON_VARIABLE
  if (!allVariablesPresent){
    if (this->verbosity>=TVar::ERROR) MELAerr << "ElectronHandler::constructElectrons: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=TVar::DEBUG) MELAout << "ElectronHandler::constructElectrons: All variables are set up!" << endl;

  if (itBegin_charge == itEnd_charge) return true; // Construction is successful, it is just that no electrons exist.

  size_t nProducts = (itEnd_charge - itBegin_charge);
  productList.reserve(nProducts);
#define ELECTRON_VARIABLE(TYPE, NAME, DEFVAL) auto it_##NAME = itBegin_##NAME;
  ELECTRON_MOMENTUM_VARIABLES;
  ELECTRON_VARIABLES;
#undef ELECTRON_VARIABLE
  {
    size_t ip=0;
    while (it_charge != itEnd_charge){
      if (this->verbosity>=TVar::DEBUG) MELAout << "ElectronHandler::constructElectrons: Attempting electron " << ip << "..." << endl;

      ParticleObject::LorentzVector_t momentum;
      momentum = ParticleObject::PolarLorentzVector_t(*it_pt, *it_eta, *it_phi, *it_mass); // Yes you have to do this on a separate line because CMSSW...
      productList.push_back(new ProductType_t(-11*(*it_charge>0 ? 1 : -1), momentum));
      ElectronObject*& obj = productList.back();

      // Set extras
#define ELECTRON_VARIABLE(TYPE, NAME, DEFVAL) obj->extras.NAME = *it_##NAME;
      ELECTRON_COMMON_VARIABLES;
      if (this->has_mvaid_extras){
        ELECTRON_MVAID_EXTRA_VARIABLES;
      }
#undef ELECTRON_VARIABLE

      // Set particle index as its unique identifier
      obj->setUniqueIdentifier(ip);

      // Replace momentum
      obj->makeFinalMomentum(syst);

      // Set the selection bits
      ElectronSelectionHelpers::setSelectionBits(*obj);

      if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Success!" << endl;

      ip++;
#define ELECTRON_VARIABLE(TYPE, NAME, DEFVAL) it_##NAME++;
      ELECTRON_MOMENTUM_VARIABLES;
      ELECTRON_COMMON_VARIABLES;
      if (this->has_mvaid_extras){
        ELECTRON_MVAID_EXTRA_VARIABLES;
      }
#undef ELECTRON_VARIABLE
    }
  }
  // Sort particles
  ParticleObjectHelpers::sortByGreaterPt(productList);

  return true;
}

bool ElectronHandler::wrapTree(BaseTree* tree){
  if (!tree) return false;

  std::vector<TString> bnames;
  tree->getValidBranchNamesWithoutAlias(bnames, false);
  this->has_mvaid_extras = (std::find(bnames.cbegin(), bnames.cend(), ElectronHandler::colName + "_id_MVA_Fall17V2_NoIso_Val")!=bnames.cend());

  return IvyBase::wrapTree(tree);
}

void ElectronHandler::bookBranches(BaseTree* tree){
  if (!tree) return;

  std::vector<TString> bnames;
  tree->getValidBranchNamesWithoutAlias(bnames, false);
  this->has_mvaid_extras = (std::find(bnames.cbegin(), bnames.cend(), ElectronHandler::colName + "_id_MVA_Fall17V2_NoIso_Val")!=bnames.cend());
#define ELECTRON_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(ElectronHandler::colName + "_" + #NAME); this->defineConsumedSloppy(#NAME);
  if (this->has_mvaid_extras){
    ELECTRON_MVAID_EXTRA_VARIABLES;
  }
#undef ELECTRON_VARIABLE

#define ELECTRON_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<std::vector<TYPE>*>(ElectronHandler::colName + "_" + #NAME, nullptr);
  ELECTRON_MOMENTUM_VARIABLES;
  ELECTRON_COMMON_VARIABLES;
  if (this->has_mvaid_extras){
    ELECTRON_MVAID_EXTRA_VARIABLES;
  }
#undef ELECTRON_VARIABLE
}


#undef VECTOR_ITERATOR_HANDLER_DIRECTIVES
