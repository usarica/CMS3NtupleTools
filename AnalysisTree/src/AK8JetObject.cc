#include <algorithm>
#include <utility>
#include "AK8JetObject.h"


AK8JetVariables::AK8JetVariables(){
#define AK8JET_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=DEFVAL;
  AK8JET_VARIABLES;
#undef AK8JET_VARIABLE
}
AK8JetVariables::AK8JetVariables(AK8JetVariables const& other){
#define AK8JET_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=other.NAME;
  AK8JET_VARIABLES;
#undef AK8JET_VARIABLE
}
void AK8JetVariables::swap(AK8JetVariables& other){
#define AK8JET_VARIABLE(TYPE, NAME, DEFVAL) std::swap(this->NAME, other.NAME);
  AK8JET_VARIABLES;
#undef AK8JET_VARIABLE
}
AK8JetVariables& AK8JetVariables::operator=(const AK8JetVariables& other){
  AK8JetVariables tmp(other);
  swap(tmp);
  return *this;
}

AK8JetObject::AK8JetObject() :
  ParticleObject(),
  extras(),
  currentSystScale(1)
{}
AK8JetObject::AK8JetObject(LorentzVector_t const& momentum_) :
  ParticleObject(0, momentum_),
  extras(),
  currentSystScale(1)
{}
AK8JetObject::AK8JetObject(const AK8JetObject& other) :
  ParticleObject(other),
  extras(other.extras),
  currentSystScale(other.currentSystScale)
{}
void AK8JetObject::swap(AK8JetObject& other){
  std::swap(id, other.id);
  std::swap(selectionBits, other.selectionBits);
  std::swap(momentum, other.momentum);
  extras.swap(other.extras);
  std::swap(currentSystScale, other.currentSystScale);
}
AK8JetObject& AK8JetObject::operator=(const AK8JetObject& other){
  AK8JetObject tmp(other);
  swap(tmp);
  return *this;
}
AK8JetObject::~AK8JetObject(){}

void AK8JetObject::makeFinalMomentum(SystematicsHelpers::SystematicVariationTypes const& syst){
  using namespace SystematicsHelpers;

  float scale=1;
  switch (syst){
  case eJECDn:
    scale = extras.JECDn * extras.JERNominal;
    break;
  case eJECUp:
    scale = extras.JECUp * extras.JERNominal;
    break;
  case eJERDn:
    scale = extras.JECNominal * extras.JERDn;
    break;
  case eJERUp:
    scale = extras.JECNominal * extras.JERUp;
    break;
  case sUncorrected:
    break;
  default:
    scale = extras.JECNominal * extras.JERNominal;
    break;
  }
  momentum = momentum * (scale/currentSystScale);
  currentSystScale = scale;
}
