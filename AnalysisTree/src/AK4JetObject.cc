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
  extras()
{}
AK4JetObject::AK4JetObject(int id_) :
  ParticleObject(id_),
  extras()
{}
AK4JetObject::AK4JetObject(int id_, LorentzVector_t const& momentum_) :
  ParticleObject(id_, momentum_),
  extras()
{}
AK4JetObject::AK4JetObject(const AK4JetObject& other) :
  ParticleObject(other),
  extras(other.extras)
{}
void AK4JetObject::swap(AK4JetObject& other){
  std::swap(id, other.id);
  std::swap(selectionBits, other.selectionBits);
  std::swap(momentum, other.momentum);
  extras.swap(other.extras);
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

  LorentzVector_t res=momentum;
  switch (syst){
  case eJECDn:
    res = res * extras.JECDn * extras.JERNominal;
    break;
  case eJECUp:
    res = res * extras.JECUp * extras.JERNominal;
    break;
  case eJERDn:
    res = res * extras.JECNominal * extras.JERDn;
    break;
  case eJERUp:
    res = res * extras.JECNominal * extras.JERUp;
    break;
  case sUncorrected:
    break;
  default:
    res = res * extras.JECNominal * extras.JERNominal;
    break;
  }

  momentum=res;
}
