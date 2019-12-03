#include <algorithm>
#include <utility>
#include <cmath>
#include "METObject.h"


METVariables::METVariables(){
#define MET_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=DEFVAL;
  MET_VARIABLES;
#undef MET_VARIABLE
}
METVariables::METVariables(METVariables const& other){
#define MET_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=other.NAME;
  MET_VARIABLES;
#undef MET_VARIABLE
}
void METVariables::swap(METVariables& other){
#define MET_VARIABLE(TYPE, NAME, DEFVAL) std::swap(this->NAME, other.NAME);
  MET_VARIABLES;
#undef MET_VARIABLE
}
METVariables& METVariables::operator=(const METVariables& other){
  METVariables tmp(other);
  swap(tmp);
  return *this;
}

METObject::METObject() :
  extras(),
  currentSyst(SystematicsHelpers::sNominal)
{}
METObject::METObject(const METObject& other) :
  extras(other.extras),
  currentSyst(other.currentSyst)
{}
void METObject::swap(METObject& other){
  extras.swap(other.extras);
  std::swap(currentSyst, other.currentSyst);
}
METObject& METObject::operator=(const METObject& other){
  METObject tmp(other);
  swap(tmp);
  return *this;
}
METObject::~METObject(){}

void METObject::setSystematic(SystematicsHelpers::SystematicVariationTypes const& syst){ currentSyst = syst; }

void METObject::getPtPhi(float const*& pt, float const*& phi) const{
  using namespace SystematicsHelpers;

  switch (currentSyst){
  case eJECDn:
    pt = &(extras.met_JECDn);
    phi = &(extras.metPhi_JECDn);
    break;
  case eJECUp:
    pt = &(extras.met_JECUp);
    phi = &(extras.metPhi_JECUp);
    break;
  case eMETDn:
    pt = &(extras.met_METDn);
    phi = &(extras.metPhi_METDn);
    break;
  case eMETUp:
    pt = &(extras.met_METUp);
    phi = &(extras.metPhi_METUp);
    break;
  case sUncorrected:
    pt = &(extras.met_original);
    phi = &(extras.metPhi_original);
    break;
  default:
    pt = &(extras.met_Nominal);
    phi = &(extras.metPhi_Nominal);
    break;
  }
}
float const& METObject::met() const{
  float const* pt=nullptr;
  float const* phi=nullptr;
  getPtPhi(pt, phi);
  return *pt;
}
float const& METObject::phi() const{
  float const* pt=nullptr;
  float const* phi=nullptr;
  getPtPhi(pt, phi);
  return *phi;
}
float METObject::px(float phi_rot) const{
  float const* pt=nullptr;
  float const* phi=nullptr;
  getPtPhi(pt, phi);
  return (*pt) * std::cos((*phi) + phi_rot);
}
float METObject::py(float phi_rot) const{
  float const* pt=nullptr;
  float const* phi=nullptr;
  getPtPhi(pt, phi);
  return (*pt) * std::sin((*phi) + phi_rot);
}
ParticleObject::LorentzVector_t METObject::p4(float phi_rot) const{
  float const* pt=nullptr;
  float const* phi=nullptr;
  getPtPhi(pt, phi);
  return ParticleObject::LorentzVector_t((*pt) * std::cos((*phi) + phi_rot), (*pt) * std::sin((*phi) + phi_rot), 0., *pt);
}
