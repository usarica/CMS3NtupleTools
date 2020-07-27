#include <algorithm>
#include <utility>
#include <CMSDataTools/AnalysisTree/interface/HelperFunctions.h>
#include <CMS3/Dictionaries/interface/EgammaFiduciality.h>
#include "SuperclusterObject.h"


SuperclusterVariables::SuperclusterVariables(){
#define SUPERCLUSTER_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=DEFVAL;
  SUPERCLUSTER_VARIABLES;
#undef SUPERCLUSTER_VARIABLE
}
SuperclusterVariables::SuperclusterVariables(SuperclusterVariables const& other){
#define SUPERCLUSTER_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=other.NAME;
  SUPERCLUSTER_VARIABLES;
#undef SUPERCLUSTER_VARIABLE
}
void SuperclusterVariables::swap(SuperclusterVariables& other){
#define SUPERCLUSTER_VARIABLE(TYPE, NAME, DEFVAL) std::swap(this->NAME, other.NAME);
  SUPERCLUSTER_VARIABLES;
#undef SUPERCLUSTER_VARIABLE
}
SuperclusterVariables& SuperclusterVariables::operator=(const SuperclusterVariables& other){
  SuperclusterVariables tmp(other);
  swap(tmp);
  return *this;
}

SuperclusterObject::SuperclusterObject() :
  ParticleObject(),
  extras()
{}
SuperclusterObject::SuperclusterObject(LorentzVector_t const& momentum_) :
  ParticleObject(0, momentum_),
  extras()
{}
SuperclusterObject::SuperclusterObject(const SuperclusterObject& other) :
  ParticleObject(other),
  extras(other.extras)
{}
void SuperclusterObject::swap(SuperclusterObject& other){
  ParticleObject::swap(other);
  extras.swap(other.extras);
}
SuperclusterObject& SuperclusterObject::operator=(const SuperclusterObject& other){
  SuperclusterObject tmp(other);
  swap(tmp);
  return *this;
}
SuperclusterObject::~SuperclusterObject(){}

bool SuperclusterObject::isEBEEGap() const{
  return HelperFunctions::test_bit(extras.fid_mask, EgammaFiduciality::ISEBEEGAP);
}
bool SuperclusterObject::isEB() const{
  return HelperFunctions::test_bit(extras.fid_mask, EgammaFiduciality::ISEB);
}
bool SuperclusterObject::isEE() const{
  return HelperFunctions::test_bit(extras.fid_mask, EgammaFiduciality::ISEE);
}
bool SuperclusterObject::isAnyGap() const{
  return HelperFunctions::test_bit(extras.fid_mask, EgammaFiduciality::ISGAP);
}
