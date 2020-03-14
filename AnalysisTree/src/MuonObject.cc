#include <algorithm>
#include <utility>
#include "MuonObject.h"


MuonVariables::MuonVariables(){
#define MUON_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=DEFVAL;
  MUON_VARIABLES;
#undef MUON_VARIABLE
}
MuonVariables::MuonVariables(MuonVariables const& other){
#define MUON_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=other.NAME;
  MUON_VARIABLES;
#undef MUON_VARIABLE
}
void MuonVariables::swap(MuonVariables& other){
#define MUON_VARIABLE(TYPE, NAME, DEFVAL) std::swap(this->NAME, other.NAME);
  MUON_VARIABLES;
#undef MUON_VARIABLE
}
MuonVariables& MuonVariables::operator=(const MuonVariables& other){
  MuonVariables tmp(other);
  swap(tmp);
  return *this;
}

MuonObject::MuonObject() :
  ParticleObject(),
  extras(),
  currentSystScale(1)
{}
MuonObject::MuonObject(cms3_id_t id_) :
  ParticleObject(id_),
  extras(),
  currentSystScale(1)
{}
MuonObject::MuonObject(cms3_id_t id_, LorentzVector_t const& momentum_) :
  ParticleObject(id_, momentum_),
  extras(),
  currentSystScale(1)
{}
MuonObject::MuonObject(const MuonObject& other) :
  ParticleObject(other),
  extras(other.extras),
  currentSystScale(other.currentSystScale)
{}
void MuonObject::swap(MuonObject& other){
  ParticleObject::swap(other);
  extras.swap(other.extras);
  std::swap(currentSystScale, other.currentSystScale);
}
MuonObject& MuonObject::operator=(const MuonObject& other){
  MuonObject tmp(other);
  swap(tmp);
  return *this;
}
MuonObject::~MuonObject(){}

void MuonObject::makeFinalMomentum(SystematicsHelpers::SystematicVariationTypes const& syst){
  using namespace SystematicsHelpers;

  float scale=1;
  switch (syst){
  case eEleScaleDn:
    scale = extras.scale_smear_pt_corr_scale_totalDn;
    break;
  case eEleScaleUp:
    scale = extras.scale_smear_pt_corr_scale_totalUp;
    break;
  case eEleResDn:
    scale = extras.scale_smear_pt_corr_smear_totalDn;
    break;
  case eEleResUp:
    scale = extras.scale_smear_pt_corr_smear_totalUp;
    break;
  case sUncorrected:
    break;
  default:
    scale = extras.scale_smear_pt_corr;
    break;
  }
  momentum = PolarLorentzVector_t(momentum.Pt() * scale/currentSystScale, momentum.Eta(), momentum.Phi(), momentum.M());
  currentSystScale = scale;
}
