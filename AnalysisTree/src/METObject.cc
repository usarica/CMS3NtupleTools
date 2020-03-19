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
  currentSyst(SystematicsHelpers::sNominal),
  currentXYshift(0, 0, 0, 0),
  currentJERShift(0, 0, 0, 0),
  particleMomentumCorrections(0, 0, 0, 0)
{}
METObject::METObject(const METObject& other) :
  extras(other.extras),
  currentSyst(other.currentSyst),
  currentXYshift(other.currentXYshift),
  currentJERShift(other.currentJERShift),
  particleMomentumCorrections(other.particleMomentumCorrections)
{}
void METObject::swap(METObject& other){
  extras.swap(other.extras);
  std::swap(currentSyst, other.currentSyst);
  std::swap(currentXYshift, other.currentXYshift);
  std::swap(currentJERShift, other.currentJERShift);
  std::swap(particleMomentumCorrections, other.particleMomentumCorrections);
}
METObject& METObject::operator=(const METObject& other){
  METObject tmp(other);
  swap(tmp);
  return *this;
}
METObject::~METObject(){}

void METObject::setSystematic(SystematicsHelpers::SystematicVariationTypes const& syst){
  currentSyst = syst;
  this->setJERShifts();
}

void METObject::setJERShifts(){
  using namespace SystematicsHelpers;

  float const* jerShift_px;
  float const* jerShift_py;
  switch (currentSyst){
  case eJERDn:
    jerShift_px = &(extras.metShift_px_JERDn);
    jerShift_py = &(extras.metShift_py_JERDn);
    break;
  case eJERUp:
    jerShift_px = &(extras.metShift_px_JERUp);
    jerShift_py = &(extras.metShift_py_JERUp);
    break;
  default:
    jerShift_px = &(extras.metShift_px_JERNominal);
    jerShift_py = &(extras.metShift_py_JERNominal);
    break;
  }
  currentJERShift = ParticleObject::LorentzVector_t(*jerShift_px, *jerShift_py, 0., 0.);
}

void METObject::getPtPhi(float& pt, float& phi, bool addXYShifts, bool addJERShifts, bool addParticleShifts) const{
  using namespace SystematicsHelpers;

  float const* pt_ref;
  float const* phi_ref;
  switch (currentSyst){
  case eJECDn:
    pt_ref = &(extras.met_JECDn);
    phi_ref = &(extras.metPhi_JECDn);
    break;
  case eJECUp:
    pt_ref = &(extras.met_JECUp);
    phi_ref = &(extras.metPhi_JECUp);
    break;
  case ePUDn:
    pt_ref = &(extras.met_PUDn);
    phi_ref = &(extras.metPhi_PUDn);
    break;
  case ePUUp:
    pt_ref = &(extras.met_PUUp);
    phi_ref = &(extras.metPhi_PUUp);
    break;
  case eMETDn:
    pt_ref = &(extras.met_METDn);
    phi_ref = &(extras.metPhi_METDn);
    break;
  case eMETUp:
    pt_ref = &(extras.met_METUp);
    phi_ref = &(extras.metPhi_METUp);
    break;
  case sUncorrected:
    pt_ref = &(extras.met_original);
    phi_ref = &(extras.metPhi_original);
    break;
  default:
    pt_ref = &(extras.met_Nominal);
    phi_ref = &(extras.metPhi_Nominal);
    break;
  }

  ParticleObject::LorentzVector_t tmp_p4((*pt_ref) * std::cos(*phi_ref), (*pt_ref) * std::sin(*phi_ref), 0., 0.);
  if (addXYShifts) tmp_p4 += currentXYshift;
  if (addJERShifts) tmp_p4 += currentJERShift;
  if (addParticleShifts) tmp_p4 += particleMomentumCorrections;

  pt = tmp_p4.Pt();
  phi = tmp_p4.Phi();
}
ParticleObject::LorentzVector_t::Scalar METObject::met(bool addXYShifts, bool addJERShifts, bool addParticleShifts) const{
  float pt, phi;
  getPtPhi(pt, phi, addXYShifts, addJERShifts, addParticleShifts);
  return pt;
}
ParticleObject::LorentzVector_t::Scalar METObject::phi(bool addXYShifts, bool addJERShifts, bool addParticleShifts) const{
  float pt, phi;
  getPtPhi(pt, phi, addXYShifts, addJERShifts, addParticleShifts);
  return phi;
}
ParticleObject::LorentzVector_t::Scalar METObject::px(bool addXYShifts, bool addJERShifts, bool addParticleShifts, float phi_rot) const{
  float pt, phi;
  getPtPhi(pt, phi, addXYShifts, addJERShifts, addParticleShifts);
  return pt * std::cos(phi + phi_rot);
}
ParticleObject::LorentzVector_t::Scalar METObject::py(bool addXYShifts, bool addJERShifts, bool addParticleShifts, float phi_rot) const{
  float pt, phi;
  getPtPhi(pt, phi, addXYShifts, addJERShifts, addParticleShifts);
  return pt * std::sin(phi + phi_rot);
}
ParticleObject::LorentzVector_t METObject::p4(bool addXYShifts, bool addJERShifts, bool addParticleShifts, float phi_rot) const{
  float pt, phi;
  getPtPhi(pt, phi, addXYShifts, addJERShifts, addParticleShifts);
  return ParticleObject::LorentzVector_t(pt * std::cos(phi + phi_rot), pt * std::sin(phi + phi_rot), 0., pt);
}
