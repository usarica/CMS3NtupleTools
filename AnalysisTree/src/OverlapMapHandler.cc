#include <cassert>

#include "OverlapMapHandler.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


template<> const std::string OverlapMapHandler<MuonObject, AK4JetObject>::colName = "overlapMap_" + MuonObject::colName + "_" + AK4JetObject::colName;
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
  clear();
  if (!currentTree) return false;

  bool allVariablesPresent = true;

#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) \
  std::vector<TYPE>::const_iterator itBegin_##NAME, itEnd_##NAME; \
  allVariablesPresent &= this->getConsumedCIterators<std::vector<TYPE>>(this->colName + "_" + #NAME, &itBegin_##NAME, &itEnd_##NAME);
  OVERLAPMAP_MUONS_JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE

  if (!allVariablesPresent){
    if (this->verbosity>=TVar::ERROR) MELAerr << "OverlapMapHandler<MuonObject, AK4JetObject>::constructOverlapMaps: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=TVar::DEBUG) MELAout << "OverlapMapHandler<MuonObject, AK4JetObject>::constructOverlapMaps: All variables are set up!" << endl;

  size_t n_products = (itEnd_particle_match_index - itBegin_particle_match_index);
  productList.reserve(n_products);
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) auto it_##NAME = itBegin_##NAME;
  OVERLAPMAP_MUONS_JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE

  {
    size_t ip=0;
    while (ip != n_products){
      if (this->verbosity>=TVar::DEBUG) MELAout << "OverlapMapHandler<MuonObject, AK4JetObject>::constructOverlapMaps: Attempting overlap map " << ip << "..." << endl;

      productList.push_back(new ProductType_t);
      ProductType_t*& obj = productList.back();

      obj->setFirstParticleIndex(*it_particle_match_index);
      obj->setSecondParticleIndex(*it_jet_match_index);

      // Set extras
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) obj->NAME = *it_##NAME;
      OVERLAPMAP_MUONS_JETS_LINKED_VARIABLES;
#undef OVERLAPMAP_VARIABLE

      if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Success!" << endl;

      ip++;
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) it_##NAME++;
      OVERLAPMAP_MUONS_JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE
    }
  }

  return true;
}

template<> const std::string OverlapMapHandler<MuonObject, AK8JetObject>::colName = "overlapMap_" + MuonObject::colName + "_" + AK8JetObject::colName;
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
  clear();
  if (!currentTree) return false;

  bool allVariablesPresent = true;

#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) \
  std::vector<TYPE>::const_iterator itBegin_##NAME, itEnd_##NAME; \
  allVariablesPresent &= this->getConsumedCIterators<std::vector<TYPE>>(this->colName + "_" + #NAME, &itBegin_##NAME, &itEnd_##NAME);
  OVERLAPMAP_MUONS_JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE

  if (!allVariablesPresent){
    if (this->verbosity>=TVar::ERROR) MELAerr << "OverlapMapHandler<MuonObject, AK8JetObject>::constructOverlapMaps: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=TVar::DEBUG) MELAout << "OverlapMapHandler<MuonObject, AK8JetObject>::constructOverlapMaps: All variables are set up!" << endl;

  size_t n_products = (itEnd_particle_match_index - itBegin_particle_match_index);
  productList.reserve(n_products);
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) auto it_##NAME = itBegin_##NAME;
  OVERLAPMAP_MUONS_JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE

  {
    size_t ip=0;
    while (ip != n_products){
      if (this->verbosity>=TVar::DEBUG) MELAout << "OverlapMapHandler<MuonObject, AK8JetObject>::constructOverlapMaps: Attempting overlap map " << ip << "..." << endl;

      productList.push_back(new ProductType_t);
      ProductType_t*& obj = productList.back();

      obj->setFirstParticleIndex(*it_particle_match_index);
      obj->setSecondParticleIndex(*it_jet_match_index);

      // Set extras
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) obj->NAME = *it_##NAME;
      OVERLAPMAP_MUONS_JETS_LINKED_VARIABLES;
#undef OVERLAPMAP_VARIABLE

      if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Success!" << endl;

      ip++;
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) it_##NAME++;
      OVERLAPMAP_MUONS_JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE
    }
  }

  return true;
}

template<> const std::string OverlapMapHandler<ElectronObject, AK4JetObject>::colName = "overlapMap_" + ElectronObject::colName + "_" + AK4JetObject::colName;
template<> OverlapMapHandler<ElectronObject, AK4JetObject>::OverlapMapHandler() :
  IvyBase()
{
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(this->colName + "_" + #NAME);
  OVERLAPMAP_ELECTRONS_JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE
}
template<> void OverlapMapHandler<ElectronObject, AK4JetObject>::bookBranches(BaseTree* tree){
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<std::vector<TYPE>*>(this->colName + "_" + #NAME, nullptr);
  OVERLAPMAP_ELECTRONS_JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE
}
template<> bool OverlapMapHandler<ElectronObject, AK4JetObject>::constructOverlapMaps(){
  clear();
  if (!currentTree) return false;

  bool allVariablesPresent = true;

#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) \
  std::vector<TYPE>::const_iterator itBegin_##NAME, itEnd_##NAME; \
  allVariablesPresent &= this->getConsumedCIterators<std::vector<TYPE>>(this->colName + "_" + #NAME, &itBegin_##NAME, &itEnd_##NAME);
  OVERLAPMAP_ELECTRONS_JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE

  if (!allVariablesPresent){
    if (this->verbosity>=TVar::ERROR) MELAerr << "OverlapMapHandler<ElectronObject, AK4JetObject>::constructOverlapMaps: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=TVar::DEBUG) MELAout << "OverlapMapHandler<ElectronObject, AK4JetObject>::constructOverlapMaps: All variables are set up!" << endl;

  size_t n_products = (itEnd_particle_match_index - itBegin_particle_match_index);
  productList.reserve(n_products);
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) auto it_##NAME = itBegin_##NAME;
  OVERLAPMAP_ELECTRONS_JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE

  {
    size_t ip=0;
    while (ip != n_products){
      if (this->verbosity>=TVar::DEBUG) MELAout << "OverlapMapHandler<ElectronObject, AK4JetObject>::constructOverlapMaps: Attempting overlap map " << ip << "..." << endl;

      productList.push_back(new ProductType_t);
      ProductType_t*& obj = productList.back();

      obj->setFirstParticleIndex(*it_particle_match_index);
      obj->setSecondParticleIndex(*it_jet_match_index);

      // Set extras
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) obj->NAME = *it_##NAME;
      OVERLAPMAP_ELECTRONS_JETS_LINKED_VARIABLES;
#undef OVERLAPMAP_VARIABLE

      if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Success!" << endl;

      ip++;
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) it_##NAME++;
      OVERLAPMAP_ELECTRONS_JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE
    }
  }

  return true;
}

template<> const std::string OverlapMapHandler<ElectronObject, AK8JetObject>::colName = "overlapMap_" + ElectronObject::colName + "_" + AK8JetObject::colName;
template<> OverlapMapHandler<ElectronObject, AK8JetObject>::OverlapMapHandler() :
  IvyBase()
{
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(this->colName + "_" + #NAME);
  OVERLAPMAP_ELECTRONS_JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE
}
template<> void OverlapMapHandler<ElectronObject, AK8JetObject>::bookBranches(BaseTree* tree){
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<std::vector<TYPE>*>(this->colName + "_" + #NAME, nullptr);
  OVERLAPMAP_ELECTRONS_JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE
}
template<> bool OverlapMapHandler<ElectronObject, AK8JetObject>::constructOverlapMaps(){
  clear();
  if (!currentTree) return false;

  bool allVariablesPresent = true;

#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) \
  std::vector<TYPE>::const_iterator itBegin_##NAME, itEnd_##NAME; \
  allVariablesPresent &= this->getConsumedCIterators<std::vector<TYPE>>(this->colName + "_" + #NAME, &itBegin_##NAME, &itEnd_##NAME);
  OVERLAPMAP_ELECTRONS_JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE

  if (!allVariablesPresent){
    if (this->verbosity>=TVar::ERROR) MELAerr << "OverlapMapHandler<ElectronObject, AK8JetObject>::constructOverlapMaps: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=TVar::DEBUG) MELAout << "OverlapMapHandler<ElectronObject, AK8JetObject>::constructOverlapMaps: All variables are set up!" << endl;

  size_t n_products = (itEnd_particle_match_index - itBegin_particle_match_index);
  productList.reserve(n_products);
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) auto it_##NAME = itBegin_##NAME;
  OVERLAPMAP_ELECTRONS_JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE

  {
    size_t ip=0;
    while (ip != n_products){
      if (this->verbosity>=TVar::DEBUG) MELAout << "OverlapMapHandler<ElectronObject, AK8JetObject>::constructOverlapMaps: Attempting overlap map " << ip << "..." << endl;

      productList.push_back(new ProductType_t);
      ProductType_t*& obj = productList.back();

      obj->setFirstParticleIndex(*it_particle_match_index);
      obj->setSecondParticleIndex(*it_jet_match_index);

      // Set extras
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) obj->NAME = *it_##NAME;
      OVERLAPMAP_ELECTRONS_JETS_LINKED_VARIABLES;
#undef OVERLAPMAP_VARIABLE

      if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Success!" << endl;

      ip++;
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) it_##NAME++;
      OVERLAPMAP_ELECTRONS_JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE
    }
  }

  return true;
}

template<> const std::string OverlapMapHandler<PhotonObject, AK4JetObject>::colName = "overlapMap_" + PhotonObject::colName + "_" + AK4JetObject::colName;
template<> OverlapMapHandler<PhotonObject, AK4JetObject>::OverlapMapHandler() :
  IvyBase()
{
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(this->colName + "_" + #NAME);
  OVERLAPMAP_PHOTONS_JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE
}
template<> void OverlapMapHandler<PhotonObject, AK4JetObject>::bookBranches(BaseTree* tree){
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<std::vector<TYPE>*>(this->colName + "_" + #NAME, nullptr);
  OVERLAPMAP_PHOTONS_JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE
}
template<> bool OverlapMapHandler<PhotonObject, AK4JetObject>::constructOverlapMaps(){
  clear();
  if (!currentTree) return false;

  bool allVariablesPresent = true;

#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) \
  std::vector<TYPE>::const_iterator itBegin_##NAME, itEnd_##NAME; \
  allVariablesPresent &= this->getConsumedCIterators<std::vector<TYPE>>(this->colName + "_" + #NAME, &itBegin_##NAME, &itEnd_##NAME);
  OVERLAPMAP_PHOTONS_JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE

  if (!allVariablesPresent){
    if (this->verbosity>=TVar::ERROR) MELAerr << "OverlapMapHandler<PhotonObject, AK4JetObject>::constructOverlapMaps: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=TVar::DEBUG) MELAout << "OverlapMapHandler<PhotonObject, AK4JetObject>::constructOverlapMaps: All variables are set up!" << endl;

  size_t n_products = (itEnd_particle_match_index - itBegin_particle_match_index);
  productList.reserve(n_products);
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) auto it_##NAME = itBegin_##NAME;
  OVERLAPMAP_PHOTONS_JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE

  {
    size_t ip=0;
    while (ip != n_products){
      if (this->verbosity>=TVar::DEBUG) MELAout << "OverlapMapHandler<PhotonObject, AK4JetObject>::constructOverlapMaps: Attempting overlap map " << ip << "..." << endl;

      productList.push_back(new ProductType_t);
      ProductType_t*& obj = productList.back();

      obj->setFirstParticleIndex(*it_particle_match_index);
      obj->setSecondParticleIndex(*it_jet_match_index);

      // Set extras
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) obj->NAME = *it_##NAME;
      OVERLAPMAP_PHOTONS_JETS_LINKED_VARIABLES;
#undef OVERLAPMAP_VARIABLE

      if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Success!" << endl;

      ip++;
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) it_##NAME++;
      OVERLAPMAP_PHOTONS_JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE
    }
  }

  return true;
}

template<> const std::string OverlapMapHandler<PhotonObject, AK8JetObject>::colName = "overlapMap_" + PhotonObject::colName + "_" + AK8JetObject::colName;
template<> OverlapMapHandler<PhotonObject, AK8JetObject>::OverlapMapHandler() :
  IvyBase()
{
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(this->colName + "_" + #NAME);
  OVERLAPMAP_PHOTONS_JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE
}
template<> void OverlapMapHandler<PhotonObject, AK8JetObject>::bookBranches(BaseTree* tree){
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<std::vector<TYPE>*>(this->colName + "_" + #NAME, nullptr);
  OVERLAPMAP_PHOTONS_JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE
}
template<> bool OverlapMapHandler<PhotonObject, AK8JetObject>::constructOverlapMaps(){
  clear();
  if (!currentTree) return false;

  bool allVariablesPresent = true;

#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) \
  std::vector<TYPE>::const_iterator itBegin_##NAME, itEnd_##NAME; \
  allVariablesPresent &= this->getConsumedCIterators<std::vector<TYPE>>(this->colName + "_" + #NAME, &itBegin_##NAME, &itEnd_##NAME);
  OVERLAPMAP_PHOTONS_JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE

  if (!allVariablesPresent){
    if (this->verbosity>=TVar::ERROR) MELAerr << "OverlapMapHandler<PhotonObject, AK8JetObject>::constructOverlapMaps: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=TVar::DEBUG) MELAout << "OverlapMapHandler<PhotonObject, AK8JetObject>::constructOverlapMaps: All variables are set up!" << endl;

  size_t n_products = (itEnd_particle_match_index - itBegin_particle_match_index);
  productList.reserve(n_products);
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) auto it_##NAME = itBegin_##NAME;
  OVERLAPMAP_PHOTONS_JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE

  {
    size_t ip=0;
    while (ip != n_products){
      if (this->verbosity>=TVar::DEBUG) MELAout << "OverlapMapHandler<PhotonObject, AK8JetObject>::constructOverlapMaps: Attempting overlap map " << ip << "..." << endl;

      productList.push_back(new ProductType_t);
      ProductType_t*& obj = productList.back();

      obj->setFirstParticleIndex(*it_particle_match_index);
      obj->setSecondParticleIndex(*it_jet_match_index);

      // Set extras
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) obj->NAME = *it_##NAME;
      OVERLAPMAP_PHOTONS_JETS_LINKED_VARIABLES;
#undef OVERLAPMAP_VARIABLE

      if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Success!" << endl;

      ip++;
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) it_##NAME++;
      OVERLAPMAP_PHOTONS_JETS_VARIABLES;
#undef OVERLAPMAP_VARIABLE
    }
  }

  return true;
}
