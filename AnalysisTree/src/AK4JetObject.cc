#include <cassert>
#include <algorithm>
#include <utility>
#include "AK4JetObject.h"
#include "BtagHelpers.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


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
  currentJEC_full(1),
  currentJEC_L1only(1),
  currentJER(1),
  currentSystScale(1)
{}
AK4JetObject::AK4JetObject(LorentzVector_t const& momentum_) :
  ParticleObject(0, momentum_),
  extras(),
  currentJEC_full(1),
  currentJEC_L1only(1),
  currentJER(1),
  currentSystScale(1)
{}
AK4JetObject::AK4JetObject(const AK4JetObject& other) :
  ParticleObject(other),
  extras(other.extras),
  currentJEC_full(other.currentJEC_full),
  currentJEC_L1only(other.currentJEC_L1only),
  currentJER(other.currentJER),
  currentSystScale(other.currentSystScale)
{}
void AK4JetObject::swap(AK4JetObject& other){
  ParticleObject::swap(other);
  extras.swap(other.extras);
  std::swap(currentJEC_full, other.currentJEC_full);
  std::swap(currentJEC_L1only, other.currentJEC_L1only);
  std::swap(currentJER, other.currentJER);
  std::swap(currentSystScale, other.currentSystScale);
}
AK4JetObject& AK4JetObject::operator=(const AK4JetObject& other){
  AK4JetObject tmp(other);
  swap(tmp);
  return *this;
}
AK4JetObject::~AK4JetObject(){}

float AK4JetObject::getBtagValue() const{
  switch (BtagHelpers::btagWPType){
  case BtagHelpers::kDeepFlav_Loose:
  case BtagHelpers::kDeepFlav_Medium:
  case BtagHelpers::kDeepFlav_Tight:
    return extras.deepFlavourprobb + extras.deepFlavourprobbb + extras.deepFlavourproblepb;
  case BtagHelpers::kDeepCSV_Loose:
  case BtagHelpers::kDeepCSV_Medium:
  case BtagHelpers::kDeepCSV_Tight:
    return extras.deepCSVprobb + extras.deepCSVprobbb;
  default:
    MELAerr << "AK4JetObject::getBtagValue: b-tagging WP type " << BtagHelpers::btagWPType << " is not implemented." << endl;
    assert(0);
    return -1;
  }
}

void AK4JetObject::makeFinalMomentum(SystematicsHelpers::SystematicVariationTypes const& syst){
  using namespace SystematicsHelpers;

  currentJEC_full = 1;
  currentJER = 1;
  currentJEC_L1only = 1;
  switch (syst){
  case eJECDn:
    currentJEC_full = extras.JECDn;
    currentJER = extras.JERNominal;
    currentJEC_L1only = currentJEC_full/extras.JECNominal*extras.JECL1Nominal;
    break;
  case eJECUp:
    currentJEC_full = extras.JECUp;
    currentJER = extras.JERNominal;
    currentJEC_L1only = currentJEC_full/extras.JECNominal*extras.JECL1Nominal;
    break;
  case eJERDn:
    currentJEC_full = extras.JECNominal;
    currentJER = extras.JERDn;
    currentJEC_L1only = extras.JECL1Nominal;
    break;
  case eJERUp:
    currentJEC_full = extras.JECNominal;
    currentJER = extras.JERUp;
    currentJEC_L1only = extras.JECL1Nominal;
    break;
  case sUncorrected:
    break;
  default:
    currentJEC_full = extras.JECNominal;
    currentJER = extras.JERNominal;
    currentJEC_L1only = extras.JECL1Nominal;
    break;
  }
  float scale = currentJEC_full * currentJER;
  // Test new pt
  float newpt = momentum.Pt() * (scale/currentSystScale);
  if (newpt<1e-5 && momentum.Pt()>0.f) scale = 1e-5 / momentum.Pt() * currentSystScale;
  momentum = momentum * (scale/currentSystScale);
  currentSystScale = scale;
}

ParticleObject::LorentzVector_t AK4JetObject::p4_nomu() const{
  if (momentum.Pt()<=1e-5) return momentum;
  else return momentum - LorentzVector_t(extras.mucands_sump4_px, extras.mucands_sump4_py, 0, 0)*currentJEC_full;
}
