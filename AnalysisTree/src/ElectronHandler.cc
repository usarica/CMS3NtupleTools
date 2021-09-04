#include <cassert>

#include <CMS3/Dictionaries/interface/GlobalCollectionNames.h>

#include "ParticleObjectHelpers.h"
#include "ElectronHandler.h"
#include "SamplesCore.h"
#include "ElectronSelectionHelpers.h"
#include "ParticleSelectionHelpers.h"
#include <CMS3/Dictionaries/interface/CMS3StreamHelpers.h>


using namespace std;
using namespace IvyStreamHelpers;


#define ELECTRON_MOMENTUM_VARIABLES \
ELECTRON_VARIABLE(float, pt, 0) \
ELECTRON_VARIABLE(float, eta, 0) \
ELECTRON_VARIABLE(float, phi, 0) \
ELECTRON_VARIABLE(float, mass, 0) \
ELECTRON_VARIABLE(cms3_charge_t, charge, 0)


const std::string ElectronHandler::colName = GlobalCollectionNames::colName_electrons;

ElectronHandler::ElectronHandler() :
  IvyBase(),
  has_mvaid_extras(false),
  has_genmatching(false),
  hasOverlapMaps(false),
  overlapMap_electrons_ak4jets(nullptr),
  overlapMap_electrons_ak8jets(nullptr)
{
#define ELECTRON_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(ElectronHandler::colName + "_" + #NAME);
  ELECTRON_MOMENTUM_VARIABLES;
  ELECTRON_COMMON_VARIABLES;
#undef ELECTRON_VARIABLE
}


bool ElectronHandler::constructElectrons(SystematicsHelpers::SystematicVariationTypes const& syst, std::vector<PFCandidateObject*> const* pfcandidates){
  if (this->isAlreadyCached()) return true;

  clear();
  if (!currentTree) return false;

  bool res = (constructElectronObjects(syst) && associatePFCandidates(pfcandidates) && linkOverlapElements());

  if (res) this->cacheEvent();
  return res;
}

bool ElectronHandler::associatePFCandidates(std::vector<PFCandidateObject*> const* pfcandidates) const{
  if (!pfcandidates) return true;

  for (auto const& pfcand:(*pfcandidates)){
    auto const& associated_particle_indices = pfcand->extras.matched_electron_index_list;
    for (auto const& part:productList){
      if (HelperFunctions::checkListVariable(associated_particle_indices, part->getUniqueIdentifier())){
        part->addDaughter(pfcand);
        pfcand->addMother(part);
      }
    }
  }

  return true;
}

bool ElectronHandler::linkOverlapElements() const{
  if (!hasOverlapMaps) return true;

  overlapMap_electrons_ak4jets->constructOverlapMaps();
  for (auto const& ome:overlapMap_electrons_ak4jets->getProducts()) ome->linkFirstElement(productList);

  overlapMap_electrons_ak8jets->constructOverlapMaps();
  for (auto const& ome:overlapMap_electrons_ak8jets->getProducts()) ome->linkFirstElement(productList);

  return true;
}

bool ElectronHandler::constructElectronObjects(SystematicsHelpers::SystematicVariationTypes const& syst){
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
  if (this->has_genmatching){
    ELECTRON_GENINFO_VARIABLES;
  }
#undef ELECTRON_VARIABLE
  if (!allVariablesPresent){
    if (this->verbosity>=MiscUtils::ERROR) IVYerr << "ElectronHandler::constructElectronObjects: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=MiscUtils::DEBUG) IVYout << "ElectronHandler::constructElectronObjects: All variables are set up!" << endl;

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
      if (this->verbosity>=MiscUtils::DEBUG) IVYout << "ElectronHandler::constructElectronObjects: Attempting electron " << ip << "..." << endl;

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
      if (this->has_genmatching){
        ELECTRON_GENINFO_VARIABLES;
      }
#undef ELECTRON_VARIABLE

      // Set particle index as its unique identifier
      obj->setUniqueIdentifier(ip);

      // Replace momentum
      obj->makeFinalMomentum(syst);

      // Set the selection bits
      ElectronSelectionHelpers::setSelectionBits(*obj);

      if (this->verbosity>=MiscUtils::DEBUG) IVYout << "\t- Success!" << endl;

      ip++;
#define ELECTRON_VARIABLE(TYPE, NAME, DEFVAL) it_##NAME++;
      ELECTRON_MOMENTUM_VARIABLES;
      ELECTRON_COMMON_VARIABLES;
      if (this->has_mvaid_extras){
        ELECTRON_MVAID_EXTRA_VARIABLES;
      }
      if (this->has_genmatching){
        ELECTRON_GENINFO_VARIABLES;
      }
#undef ELECTRON_VARIABLE
    }
  }
  // Sort particles
  ParticleObjectHelpers::sortByGreaterPt(productList);

  return true;
}

void ElectronHandler::checkOptionalInfo(BaseTree* tree, bool& flag_mvaid_extras, bool& flag_genmatching){
  flag_mvaid_extras = flag_genmatching = true;

  std::vector<TString> bnames;
  tree->getValidBranchNamesWithoutAlias(bnames, false);

#define ELECTRON_VARIABLE(TYPE, NAME, DEFVAL) flag_mvaid_extras &= (std::find(bnames.cbegin(), bnames.cend(), ElectronHandler::colName + "_" + #NAME)!=bnames.cend());
  ELECTRON_MVAID_EXTRA_VARIABLES;
#undef ELECTRON_VARIABLE
#define ELECTRON_VARIABLE(TYPE, NAME, DEFVAL) flag_genmatching &= (std::find(bnames.cbegin(), bnames.cend(), ElectronHandler::colName + "_" + #NAME)!=bnames.cend());
  ELECTRON_GENINFO_VARIABLES;
#undef ELECTRON_VARIABLE
}

bool ElectronHandler::wrapTree(BaseTree* tree){
  if (!tree) return false;

  ElectronHandler::checkOptionalInfo(tree, this->has_mvaid_extras, this->has_genmatching);

  return IvyBase::wrapTree(tree);
}

void ElectronHandler::bookBranches(BaseTree* tree){
  if (!tree) return;

  ElectronHandler::checkOptionalInfo(tree, this->has_mvaid_extras, this->has_genmatching);
#define ELECTRON_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(ElectronHandler::colName + "_" + #NAME); this->defineConsumedSloppy(ElectronHandler::colName + "_" + #NAME);
  if (this->has_mvaid_extras){
    ELECTRON_MVAID_EXTRA_VARIABLES;
  }
  if (this->has_genmatching){
    ELECTRON_GENINFO_VARIABLES;
  }
#undef ELECTRON_VARIABLE

#define ELECTRON_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<std::vector<TYPE>*>(ElectronHandler::colName + "_" + #NAME, nullptr);
  ELECTRON_MOMENTUM_VARIABLES;
  ELECTRON_COMMON_VARIABLES;
  if (this->has_mvaid_extras){
    ELECTRON_MVAID_EXTRA_VARIABLES;
  }
  if (this->has_genmatching){
    ELECTRON_GENINFO_VARIABLES;
  }
#undef ELECTRON_VARIABLE
}

void ElectronHandler::registerOverlapMaps(
  OverlapMapHandler<ElectronObject, AK4JetObject>& overlapMap_electrons_ak4jets_,
  OverlapMapHandler<ElectronObject, AK8JetObject>& overlapMap_electrons_ak8jets_
){
  overlapMap_electrons_ak4jets = &overlapMap_electrons_ak4jets_;
  overlapMap_electrons_ak8jets = &overlapMap_electrons_ak8jets_;
  hasOverlapMaps = true;
}


#undef VECTOR_ITERATOR_HANDLER_DIRECTIVES
