#include <algorithm>
#include <utility>
#include "PhotonObject.h"


PhotonVariables::PhotonVariables(){
#define PHOTON_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=DEFVAL;
  PHOTON_VARIABLES;
#undef PHOTON_VARIABLE
}
PhotonVariables::PhotonVariables(PhotonVariables const& other){
#define PHOTON_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=other.NAME;
  PHOTON_VARIABLES;
#undef PHOTON_VARIABLE
}
void PhotonVariables::swap(PhotonVariables& other){
#define PHOTON_VARIABLE(TYPE, NAME, DEFVAL) std::swap(this->NAME, other.NAME);
  PHOTON_VARIABLES;
#undef PHOTON_VARIABLE
}
PhotonVariables& PhotonVariables::operator=(const PhotonVariables& other){
  PhotonVariables tmp(other);
  swap(tmp);
  return *this;
}

PhotonObject::PhotonObject() :
  ParticleObject(),
  extras(),
  currentSystScale(1)
{}
PhotonObject::PhotonObject(LorentzVector_t const& momentum_) :
  ParticleObject(22, momentum_),
  extras(),
  currentSystScale(1)
{}
PhotonObject::PhotonObject(const PhotonObject& other) :
  ParticleObject(other),
  extras(other.extras),
  currentSystScale(other.currentSystScale)
{}
void PhotonObject::swap(PhotonObject& other){
  ParticleObject::swap(other);
  extras.swap(other.extras);
  std::swap(currentSystScale, other.currentSystScale);
}
PhotonObject& PhotonObject::operator=(const PhotonObject& other){
  PhotonObject tmp(other);
  swap(tmp);
  return *this;
}
PhotonObject::~PhotonObject(){}

void PhotonObject::makeFinalMomentum(SystematicsHelpers::SystematicVariationTypes const& syst){
  using namespace SystematicsHelpers;

  float scale=1;
  switch (syst){
  case ePhoScaleDn:
    scale = extras.scale_smear_corr_scale_totalDn;
    break;
  case ePhoScaleUp:
    scale = extras.scale_smear_corr_scale_totalUp;
    break;
  case ePhoResDn:
    scale = extras.scale_smear_corr_smear_totalDn;
    break;
  case ePhoResUp:
    scale = extras.scale_smear_corr_smear_totalUp;
    break;
  case sUncorrected:
    break;
  default:
    scale = extras.scale_smear_corr;
    break;
  }
  momentum = momentum * (scale/currentSystScale);
  currentSystScale = scale;
}

ParticleObject::LorentzVector_t::Scalar PhotonObject::uncorrected_pt() const{ return this->pt()/currentSystScale; }
