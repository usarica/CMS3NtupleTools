#include <algorithm>
#include <utility>
#include "ElectronObject.h"


ElectronVariables::ElectronVariables(){
#define ELECTRON_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=DEFVAL;
  ELECTRON_VARIABLES;
#undef ELECTRON_VARIABLE
}
ElectronVariables::ElectronVariables(ElectronVariables const& other){
#define ELECTRON_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=other.NAME;
  ELECTRON_VARIABLES;
#undef ELECTRON_VARIABLE
}
void ElectronVariables::swap(ElectronVariables& other){
#define ELECTRON_VARIABLE(TYPE, NAME, DEFVAL) std::swap(this->NAME, other.NAME);
  ELECTRON_VARIABLES;
#undef ELECTRON_VARIABLE
}
ElectronVariables& ElectronVariables::operator=(const ElectronVariables& other){
  ElectronVariables tmp(other);
  swap(tmp);
  return *this;
}

ElectronObject::ElectronObject() :
  ParticleObject(),
  extras()
{}
ElectronObject::ElectronObject(int id_) :
  ParticleObject(id_),
  extras()
{}
ElectronObject::ElectronObject(int id_, LorentzVector_t const& momentum_) :
  ParticleObject(id_, momentum_),
  extras()
{}
ElectronObject::ElectronObject(const ElectronObject& other) :
  ParticleObject(other),
  extras(other.extras)
{}
void ElectronObject::swap(ElectronObject& other){
  std::swap(id, other.id);
  std::swap(selectionBits, other.selectionBits);
  std::swap(momentum, other.momentum);
  extras.swap(other.extras);
}
ElectronObject& ElectronObject::operator=(const ElectronObject& other){
  ElectronObject tmp(other);
  swap(tmp);
  return *this;
}
ElectronObject::~ElectronObject(){}

void ElectronObject::makeFinalMomentum(SystematicsHelpers::SystematicVariationTypes const& syst){
  using namespace SystematicsHelpers;

  LorentzVector_t res=momentum;
  switch (syst){
  case eEleScaleDn:
    res = res * extras.scale_smear_corr_scale_totalDn;
    break;
  case eEleScaleUp:
    res = res * extras.scale_smear_corr_scale_totalUp;
    break;
  case eEleResDn:
    res = res * extras.scale_smear_corr_smear_totalDn;
    break;
  case eEleResUp:
    res = res * extras.scale_smear_corr_smear_totalUp;
    break;
  case sUncorrected:
    break;
  default:
    res = res * extras.scale_smear_corr;
    break;
  }

  momentum=res;
}
