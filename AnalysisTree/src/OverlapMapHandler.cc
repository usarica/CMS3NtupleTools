#include <cassert>

#include <CMS3/Dictionaries/interface/GlobalCollectionNames.h>

#include "OverlapMapHandler.h"
#include <CMS3/Dictionaries/interface/CMS3StreamHelpers.h>


using namespace std;
using namespace IvyStreamHelpers;


template<> const std::string OverlapMapHandler<MuonObject, AK4JetObject>::colName = GlobalCollectionNames::colName_overlapMap + "_" + GlobalCollectionNames::colName_muons + "_" + GlobalCollectionNames::colName_ak4jets;
template<> OverlapMapHandler<MuonObject, AK4JetObject>::OverlapMapHandler() :
  IvyBase()
{
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(this->colName + "_" + #NAME);
  OVERLAPMAP_MUONS_JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE
}
template<> void OverlapMapHandler<MuonObject, AK4JetObject>::bookBranches(BaseTree* tree){
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<std::vector<TYPE>*>(this->colName + "_" + #NAME, nullptr);
  OVERLAPMAP_MUONS_JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE
}
template<> bool OverlapMapHandler<MuonObject, AK4JetObject>::constructOverlapMaps(){
  if (this->isAlreadyCached()) return true;

  clear();
  if (!currentTree) return false;

  bool allVariablesPresent = true;

#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) \
  std::vector<TYPE>::const_iterator itBegin_##NAME, itEnd_##NAME; \
  allVariablesPresent &= this->getConsumedCIterators<std::vector<TYPE>>(this->colName + "_" + #NAME, &itBegin_##NAME, &itEnd_##NAME);
  OVERLAPMAP_MUONS_JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE

  if (!allVariablesPresent){
    if (this->verbosity>=MiscUtils::ERROR) IVYerr << "OverlapMapHandler<MuonObject, AK4JetObject>::constructOverlapMaps: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=MiscUtils::DEBUG) IVYout << "OverlapMapHandler<MuonObject, AK4JetObject>::constructOverlapMaps: All variables are set up!" << endl;

  size_t n_products = (itEnd_particle_match_index - itBegin_particle_match_index);
  productList.reserve(n_products);
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) auto it_##NAME = itBegin_##NAME;
  OVERLAPMAP_MUONS_JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE

  {
    size_t ip=0;
    while (ip != n_products){
      if (this->verbosity>=MiscUtils::DEBUG) IVYout << "OverlapMapHandler<MuonObject, AK4JetObject>::constructOverlapMaps: Attempting overlap map " << ip << "..." << endl;

      productList.push_back(new ProductType_t);
      ProductType_t*& obj = productList.back();

      obj->setFirstParticleIndex(*it_particle_match_index);
      obj->setSecondParticleIndex(*it_jet_match_index);

      // Set extras
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) obj->NAME = *it_##NAME;
      OVERLAPMAP_MUONS_JETS_LINKED_VARIABLES;
#undef OVERLAPMAP_VARIABLE

      if (this->verbosity>=MiscUtils::DEBUG) IVYout << "\t- Success!" << endl;

      ip++;
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) it_##NAME++;
      OVERLAPMAP_MUONS_JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE
    }
  }

  this->cacheEvent();
  return true;
}

template<> const std::string OverlapMapHandler<MuonObject, AK8JetObject>::colName = GlobalCollectionNames::colName_overlapMap + "_" + GlobalCollectionNames::colName_muons + "_" + GlobalCollectionNames::colName_ak8jets;
template<> OverlapMapHandler<MuonObject, AK8JetObject>::OverlapMapHandler() :
  IvyBase()
{
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(this->colName + "_" + #NAME);
  OVERLAPMAP_MUONS_JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE
}
template<> void OverlapMapHandler<MuonObject, AK8JetObject>::bookBranches(BaseTree* tree){
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<std::vector<TYPE>*>(this->colName + "_" + #NAME, nullptr);
  OVERLAPMAP_MUONS_JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE
}
template<> bool OverlapMapHandler<MuonObject, AK8JetObject>::constructOverlapMaps(){
  if (this->isAlreadyCached()) return true;

  clear();
  if (!currentTree) return false;

  bool allVariablesPresent = true;

#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) \
  std::vector<TYPE>::const_iterator itBegin_##NAME, itEnd_##NAME; \
  allVariablesPresent &= this->getConsumedCIterators<std::vector<TYPE>>(this->colName + "_" + #NAME, &itBegin_##NAME, &itEnd_##NAME);
  OVERLAPMAP_MUONS_JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE

  if (!allVariablesPresent){
    if (this->verbosity>=MiscUtils::ERROR) IVYerr << "OverlapMapHandler<MuonObject, AK8JetObject>::constructOverlapMaps: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=MiscUtils::DEBUG) IVYout << "OverlapMapHandler<MuonObject, AK8JetObject>::constructOverlapMaps: All variables are set up!" << endl;

  size_t n_products = (itEnd_particle_match_index - itBegin_particle_match_index);
  productList.reserve(n_products);
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) auto it_##NAME = itBegin_##NAME;
  OVERLAPMAP_MUONS_JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE

  {
    size_t ip=0;
    while (ip != n_products){
      if (this->verbosity>=MiscUtils::DEBUG) IVYout << "OverlapMapHandler<MuonObject, AK8JetObject>::constructOverlapMaps: Attempting overlap map " << ip << "..." << endl;

      productList.push_back(new ProductType_t);
      ProductType_t*& obj = productList.back();

      obj->setFirstParticleIndex(*it_particle_match_index);
      obj->setSecondParticleIndex(*it_jet_match_index);

      // Set extras
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) obj->NAME = *it_##NAME;
      OVERLAPMAP_MUONS_JETS_LINKED_VARIABLES;
#undef OVERLAPMAP_VARIABLE

      if (this->verbosity>=MiscUtils::DEBUG) IVYout << "\t- Success!" << endl;

      ip++;
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) it_##NAME++;
      OVERLAPMAP_MUONS_JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE
    }
  }

  this->cacheEvent();
  return true;
}

template<> const std::string OverlapMapHandler<ElectronObject, AK4JetObject>::colName = GlobalCollectionNames::colName_overlapMap + "_" + GlobalCollectionNames::colName_electrons + "_" + GlobalCollectionNames::colName_ak4jets;
template<> OverlapMapHandler<ElectronObject, AK4JetObject>::OverlapMapHandler() :
  IvyBase()
{
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(this->colName + "_" + #NAME);
  OVERLAPMAP_ELECTRONS_AK4JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE
}
template<> void OverlapMapHandler<ElectronObject, AK4JetObject>::bookBranches(BaseTree* tree){
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<std::vector<TYPE>*>(this->colName + "_" + #NAME, nullptr);
  OVERLAPMAP_ELECTRONS_AK4JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE
}
template<> bool OverlapMapHandler<ElectronObject, AK4JetObject>::constructOverlapMaps(){
  if (this->isAlreadyCached()) return true;

  clear();
  if (!currentTree) return false;

  bool allVariablesPresent = true;

#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) \
  std::vector<TYPE>::const_iterator itBegin_##NAME, itEnd_##NAME; \
  allVariablesPresent &= this->getConsumedCIterators<std::vector<TYPE>>(this->colName + "_" + #NAME, &itBegin_##NAME, &itEnd_##NAME);
  OVERLAPMAP_ELECTRONS_AK4JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE

  if (!allVariablesPresent){
    if (this->verbosity>=MiscUtils::ERROR) IVYerr << "OverlapMapHandler<ElectronObject, AK4JetObject>::constructOverlapMaps: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=MiscUtils::DEBUG) IVYout << "OverlapMapHandler<ElectronObject, AK4JetObject>::constructOverlapMaps: All variables are set up!" << endl;

  size_t n_products = (itEnd_particle_match_index - itBegin_particle_match_index);
  productList.reserve(n_products);
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) auto it_##NAME = itBegin_##NAME;
  OVERLAPMAP_ELECTRONS_AK4JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE

  {
    size_t ip=0;
    while (ip != n_products){
      if (this->verbosity>=MiscUtils::DEBUG) IVYout << "OverlapMapHandler<ElectronObject, AK4JetObject>::constructOverlapMaps: Attempting overlap map " << ip << "..." << endl;

      productList.push_back(new ProductType_t);
      ProductType_t*& obj = productList.back();

      obj->setFirstParticleIndex(*it_particle_match_index);
      obj->setSecondParticleIndex(*it_jet_match_index);

      // Set extras
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) obj->NAME = *it_##NAME;
      OVERLAPMAP_ELECTRONS_AK4JETS_LINKED_VARIABLES;
#undef OVERLAPMAP_VARIABLE

      if (this->verbosity>=MiscUtils::DEBUG) IVYout << "\t- Success!" << endl;

      ip++;
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) it_##NAME++;
      OVERLAPMAP_ELECTRONS_AK4JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE
    }
  }

  this->cacheEvent();
  return true;
}

template<> const std::string OverlapMapHandler<ElectronObject, AK8JetObject>::colName = GlobalCollectionNames::colName_overlapMap + "_" + GlobalCollectionNames::colName_electrons + "_" + GlobalCollectionNames::colName_ak8jets;
template<> OverlapMapHandler<ElectronObject, AK8JetObject>::OverlapMapHandler() :
  IvyBase()
{
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(this->colName + "_" + #NAME);
  OVERLAPMAP_ELECTRONS_AK8JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE
}
template<> void OverlapMapHandler<ElectronObject, AK8JetObject>::bookBranches(BaseTree* tree){
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<std::vector<TYPE>*>(this->colName + "_" + #NAME, nullptr);
  OVERLAPMAP_ELECTRONS_AK8JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE
}
template<> bool OverlapMapHandler<ElectronObject, AK8JetObject>::constructOverlapMaps(){
  if (this->isAlreadyCached()) return true;

  clear();
  if (!currentTree) return false;

  bool allVariablesPresent = true;

#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) \
  std::vector<TYPE>::const_iterator itBegin_##NAME, itEnd_##NAME; \
  allVariablesPresent &= this->getConsumedCIterators<std::vector<TYPE>>(this->colName + "_" + #NAME, &itBegin_##NAME, &itEnd_##NAME);
  OVERLAPMAP_ELECTRONS_AK8JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE

  if (!allVariablesPresent){
    if (this->verbosity>=MiscUtils::ERROR) IVYerr << "OverlapMapHandler<ElectronObject, AK8JetObject>::constructOverlapMaps: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=MiscUtils::DEBUG) IVYout << "OverlapMapHandler<ElectronObject, AK8JetObject>::constructOverlapMaps: All variables are set up!" << endl;

  size_t n_products = (itEnd_particle_match_index - itBegin_particle_match_index);
  productList.reserve(n_products);
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) auto it_##NAME = itBegin_##NAME;
  OVERLAPMAP_ELECTRONS_AK8JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE

  {
    size_t ip=0;
    while (ip != n_products){
      if (this->verbosity>=MiscUtils::DEBUG) IVYout << "OverlapMapHandler<ElectronObject, AK8JetObject>::constructOverlapMaps: Attempting overlap map " << ip << "..." << endl;

      productList.push_back(new ProductType_t);
      ProductType_t*& obj = productList.back();

      obj->setFirstParticleIndex(*it_particle_match_index);
      obj->setSecondParticleIndex(*it_jet_match_index);

      // Set extras
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) obj->NAME = *it_##NAME;
      OVERLAPMAP_ELECTRONS_AK8JETS_LINKED_VARIABLES;
#undef OVERLAPMAP_VARIABLE

      if (this->verbosity>=MiscUtils::DEBUG) IVYout << "\t- Success!" << endl;

      ip++;
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) it_##NAME++;
      OVERLAPMAP_ELECTRONS_AK8JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE
    }
  }

  this->cacheEvent();
  return true;
}

template<> const std::string OverlapMapHandler<PhotonObject, AK4JetObject>::colName = GlobalCollectionNames::colName_overlapMap + "_" + GlobalCollectionNames::colName_photons + "_" + GlobalCollectionNames::colName_ak4jets;
template<> OverlapMapHandler<PhotonObject, AK4JetObject>::OverlapMapHandler() :
  IvyBase()
{
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(this->colName + "_" + #NAME);
  OVERLAPMAP_PHOTONS_AK4JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE
}
template<> void OverlapMapHandler<PhotonObject, AK4JetObject>::bookBranches(BaseTree* tree){
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<std::vector<TYPE>*>(this->colName + "_" + #NAME, nullptr);
  OVERLAPMAP_PHOTONS_AK4JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE
}
template<> bool OverlapMapHandler<PhotonObject, AK4JetObject>::constructOverlapMaps(){
  if (this->isAlreadyCached()) return true;

  clear();
  if (!currentTree) return false;

  bool allVariablesPresent = true;

#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) \
  std::vector<TYPE>::const_iterator itBegin_##NAME, itEnd_##NAME; \
  allVariablesPresent &= this->getConsumedCIterators<std::vector<TYPE>>(this->colName + "_" + #NAME, &itBegin_##NAME, &itEnd_##NAME);
  OVERLAPMAP_PHOTONS_AK4JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE

  if (!allVariablesPresent){
    if (this->verbosity>=MiscUtils::ERROR) IVYerr << "OverlapMapHandler<PhotonObject, AK4JetObject>::constructOverlapMaps: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=MiscUtils::DEBUG) IVYout << "OverlapMapHandler<PhotonObject, AK4JetObject>::constructOverlapMaps: All variables are set up!" << endl;

  size_t n_products = (itEnd_particle_match_index - itBegin_particle_match_index);
  productList.reserve(n_products);
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) auto it_##NAME = itBegin_##NAME;
  OVERLAPMAP_PHOTONS_AK4JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE

  {
    size_t ip=0;
    while (ip != n_products){
      if (this->verbosity>=MiscUtils::DEBUG) IVYout << "OverlapMapHandler<PhotonObject, AK4JetObject>::constructOverlapMaps: Attempting overlap map " << ip << "..." << endl;

      productList.push_back(new ProductType_t);
      ProductType_t*& obj = productList.back();

      obj->setFirstParticleIndex(*it_particle_match_index);
      obj->setSecondParticleIndex(*it_jet_match_index);

      // Set extras
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) obj->NAME = *it_##NAME;
      OVERLAPMAP_PHOTONS_AK4JETS_LINKED_VARIABLES;
#undef OVERLAPMAP_VARIABLE

      if (this->verbosity>=MiscUtils::DEBUG) IVYout << "\t- Success!" << endl;

      ip++;
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) it_##NAME++;
      OVERLAPMAP_PHOTONS_AK4JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE
    }
  }

  this->cacheEvent();
  return true;
}

template<> const std::string OverlapMapHandler<PhotonObject, AK8JetObject>::colName = GlobalCollectionNames::colName_overlapMap + "_" + GlobalCollectionNames::colName_photons + "_" + GlobalCollectionNames::colName_ak8jets;
template<> OverlapMapHandler<PhotonObject, AK8JetObject>::OverlapMapHandler() :
  IvyBase()
{
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(this->colName + "_" + #NAME);
  OVERLAPMAP_PHOTONS_AK8JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE
}
template<> void OverlapMapHandler<PhotonObject, AK8JetObject>::bookBranches(BaseTree* tree){
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<std::vector<TYPE>*>(this->colName + "_" + #NAME, nullptr);
  OVERLAPMAP_PHOTONS_AK8JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE
}
template<> bool OverlapMapHandler<PhotonObject, AK8JetObject>::constructOverlapMaps(){
  if (this->isAlreadyCached()) return true;

  clear();
  if (!currentTree) return false;

  bool allVariablesPresent = true;

#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) \
  std::vector<TYPE>::const_iterator itBegin_##NAME, itEnd_##NAME; \
  allVariablesPresent &= this->getConsumedCIterators<std::vector<TYPE>>(this->colName + "_" + #NAME, &itBegin_##NAME, &itEnd_##NAME);
  OVERLAPMAP_PHOTONS_AK8JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE

  if (!allVariablesPresent){
    if (this->verbosity>=MiscUtils::ERROR) IVYerr << "OverlapMapHandler<PhotonObject, AK8JetObject>::constructOverlapMaps: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=MiscUtils::DEBUG) IVYout << "OverlapMapHandler<PhotonObject, AK8JetObject>::constructOverlapMaps: All variables are set up!" << endl;

  size_t n_products = (itEnd_particle_match_index - itBegin_particle_match_index);
  productList.reserve(n_products);
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) auto it_##NAME = itBegin_##NAME;
  OVERLAPMAP_PHOTONS_AK8JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE

  {
    size_t ip=0;
    while (ip != n_products){
      if (this->verbosity>=MiscUtils::DEBUG) IVYout << "OverlapMapHandler<PhotonObject, AK8JetObject>::constructOverlapMaps: Attempting overlap map " << ip << "..." << endl;

      productList.push_back(new ProductType_t);
      ProductType_t*& obj = productList.back();

      obj->setFirstParticleIndex(*it_particle_match_index);
      obj->setSecondParticleIndex(*it_jet_match_index);

      // Set extras
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) obj->NAME = *it_##NAME;
      OVERLAPMAP_PHOTONS_AK8JETS_LINKED_VARIABLES;
#undef OVERLAPMAP_VARIABLE

      if (this->verbosity>=MiscUtils::DEBUG) IVYout << "\t- Success!" << endl;

      ip++;
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) it_##NAME++;
      OVERLAPMAP_PHOTONS_AK8JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE
    }
  }

  this->cacheEvent();
  return true;
}
