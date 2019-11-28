#include <algorithm>
#include <utility>
#include "AK4JetObject.h"


AK4JetVariables::AK4JetVariables(){
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=DEFVAL;
  AK4JET_VARIABLES;
#undef AK4JET_VARIABLE
}
AK4JetVariables::AK4JetVariables(AK4JetVariables const& other){
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=other.NAME;
  AK4JET_VARIABLES;
#undef AK4JET_VARIABLE
}
void AK4JetVariables::swap(AK4JetVariables& other){
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) std::swap(this->NAME, other.NAME);
  AK4JET_VARIABLES;
#undef AK4JET_VARIABLE
}
AK4JetVariables& AK4JetVariables::operator=(const AK4JetVariables& other){
  AK4JetVariables tmp(other);
  swap(tmp);
  return *this;
}

AK4JetObject::AK4JetObject() :
  ParticleObject(),
  extras(),
  currentSystScale(1)
{}
AK4JetObject::AK4JetObject(LorentzVector_t const& momentum_) :
  ParticleObject(0, momentum_),
  extras(),
  currentSystScale(1)
{}
AK4JetObject::AK4JetObject(const AK4JetObject& other) :
  ParticleObject(other),
  extras(other.extras),
  currentSystScale(other.currentSystScale)
{}
void AK4JetObject::swap(AK4JetObject& other){
  std::swap(id, other.id);
  std::swap(selectionBits, other.selectionBits);
  std::swap(momentum, other.momentum);
  extras.swap(other.extras);
  std::swap(currentSystScale, other.currentSystScale);
}
AK4JetObject& AK4JetObject::operator=(const AK4JetObject& other){
  AK4JetObject tmp(other);
  swap(tmp);
  return *this;
}
AK4JetObject::~AK4JetObject(){}

float AK4JetObject::getBtagValue() const{
  return extras.deepFlavourprobb + extras.deepFlavourprobbb;
}

void AK4JetObject::makeFinalMomentum(SystematicsHelpers::SystematicVariationTypes const& syst){
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
