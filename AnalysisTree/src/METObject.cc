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
  currentMETShift_noJER(0, 0, 0, 0),
  currentMETShift(0, 0, 0, 0),
  currentMETShift_p4Preserved_noJER(0, 0, 0, 0),
  currentMETShift_p4Preserved(0, 0, 0, 0),
  currentXYshift(0, 0, 0, 0),
  particleMomentumCorrections(0, 0, 0, 0)
{}
METObject::METObject(const METObject& other) :
  extras(other.extras),
  currentSyst(other.currentSyst),
  currentMETShift_noJER(other.currentMETShift_noJER),
  currentMETShift(other.currentMETShift),
  currentMETShift_p4Preserved_noJER(other.currentMETShift_p4Preserved_noJER),
  currentMETShift_p4Preserved(other.currentMETShift_p4Preserved),
  currentXYshift(other.currentXYshift),
  particleMomentumCorrections(other.particleMomentumCorrections),
  currentMETCorrections(other.currentMETCorrections),
  currentJetOverlapCorrections(other.currentJetOverlapCorrections)
{}
void METObject::swap(METObject& other){
  extras.swap(other.extras);
  std::swap(currentSyst, other.currentSyst);
  std::swap(currentMETShift_noJER, other.currentMETShift_noJER);
  std::swap(currentMETShift, other.currentMETShift);
  std::swap(currentMETShift_p4Preserved_noJER, other.currentMETShift_p4Preserved_noJER);
  std::swap(currentMETShift_p4Preserved, other.currentMETShift_p4Preserved);
  std::swap(currentXYshift, other.currentXYshift);
  std::swap(particleMomentumCorrections, other.particleMomentumCorrections);
  std::swap(currentMETCorrections, other.currentMETCorrections);
  std::swap(currentJetOverlapCorrections, other.currentJetOverlapCorrections);
}
METObject& METObject::operator=(const METObject& other){
  METObject tmp(other);
  swap(tmp);
  return *this;
}
METObject::~METObject(){}

void METObject::setSystematic(SystematicsHelpers::SystematicVariationTypes const& syst){
  currentSyst = syst;
  this->setMETShifts();
}

void METObject::setMETShifts(){
  using namespace SystematicsHelpers;

  switch (currentSyst){
  case eJECDn:
    currentMETShift_noJER = ParticleObject::LorentzVector_t(extras.metShift_px_JECDn, extras.metShift_py_JECDn, 0., 0.);
    currentMETShift = ParticleObject::LorentzVector_t(extras.metShift_px_JECDn_JERNominal, extras.metShift_py_JECDn_JERNominal, 0., 0.);
    currentMETShift_p4Preserved_noJER = ParticleObject::LorentzVector_t(extras.metShift_p4Preserved_px_JECDn, extras.metShift_p4Preserved_py_JECDn, 0., 0.);
    currentMETShift_p4Preserved = ParticleObject::LorentzVector_t(extras.metShift_p4Preserved_px_JECDn_JERNominal, extras.metShift_p4Preserved_py_JECDn_JERNominal, 0., 0.);
    break;
  case eJECUp:
    currentMETShift_noJER = ParticleObject::LorentzVector_t(extras.metShift_px_JECUp, extras.metShift_py_JECUp, 0., 0.);
    currentMETShift = ParticleObject::LorentzVector_t(extras.metShift_px_JECUp_JERNominal, extras.metShift_py_JECUp_JERNominal, 0., 0.);
    currentMETShift_p4Preserved_noJER = ParticleObject::LorentzVector_t(extras.metShift_p4Preserved_px_JECUp, extras.metShift_p4Preserved_py_JECUp, 0., 0.);
    currentMETShift_p4Preserved = ParticleObject::LorentzVector_t(extras.metShift_p4Preserved_px_JECUp_JERNominal, extras.metShift_p4Preserved_py_JECUp_JERNominal, 0., 0.);
    break;
  case eJERDn:
    currentMETShift_noJER = ParticleObject::LorentzVector_t(0., 0., 0., 0.);
    currentMETShift = ParticleObject::LorentzVector_t(extras.metShift_px_JECNominal_JERDn, extras.metShift_py_JECNominal_JERDn, 0., 0.);
    currentMETShift_p4Preserved_noJER = ParticleObject::LorentzVector_t(extras.metShift_p4Preserved_px_JECNominal, extras.metShift_p4Preserved_py_JECNominal, 0., 0.);
    currentMETShift_p4Preserved = ParticleObject::LorentzVector_t(extras.metShift_p4Preserved_px_JECNominal_JERDn, extras.metShift_p4Preserved_py_JECNominal_JERDn, 0., 0.);
    break;
  case eJERUp:
    currentMETShift_noJER = ParticleObject::LorentzVector_t(0., 0., 0., 0.);
    currentMETShift = ParticleObject::LorentzVector_t(extras.metShift_px_JECNominal_JERUp, extras.metShift_py_JECNominal_JERUp, 0., 0.);
    currentMETShift_p4Preserved_noJER = ParticleObject::LorentzVector_t(extras.metShift_p4Preserved_px_JECNominal, extras.metShift_p4Preserved_py_JECNominal, 0., 0.);
    currentMETShift_p4Preserved = ParticleObject::LorentzVector_t(extras.metShift_p4Preserved_px_JECNominal_JERUp, extras.metShift_p4Preserved_py_JECNominal_JERUp, 0., 0.);
    break;
  default:
    currentMETShift_noJER = ParticleObject::LorentzVector_t(0., 0., 0., 0.);
    currentMETShift = ParticleObject::LorentzVector_t(extras.metShift_px_JECNominal_JERNominal, extras.metShift_py_JECNominal_JERNominal, 0., 0.);
    currentMETShift_p4Preserved_noJER = ParticleObject::LorentzVector_t(extras.metShift_p4Preserved_px_JECNominal, extras.metShift_p4Preserved_py_JECNominal, 0., 0.);
    currentMETShift_p4Preserved = ParticleObject::LorentzVector_t(extras.metShift_p4Preserved_px_JECNominal_JERNominal, extras.metShift_p4Preserved_py_JECNominal_JERNominal, 0., 0.);
    break;
  }

  // p4-preserved momenta are recorded in the ntuples relative to default shifts, so add the default shifts back here.
  currentMETShift_p4Preserved_noJER += currentMETShift_noJER;
  currentMETShift_p4Preserved += currentMETShift;
}

void METObject::setMETCorrection(ParticleObject::LorentzVector_t const& corr, bool hasXYShifts, bool hasJERShifts, bool hasParticleShifts, bool preserveP4){
  if (currentMETCorrections.empty()) currentMETCorrections.assign(16, ParticleObject::LorentzVector_t(0, 0, 0, 0));
  currentMETCorrections.at(8*preserveP4 + 4*hasParticleShifts + 2*hasJERShifts + 1*hasXYShifts) = corr;
}

void METObject::setJetOverlapCorrection(ParticleObject::LorentzVector_t const& corr, bool hasJERShifts, bool preserveP4){
  if (currentJetOverlapCorrections.empty()) currentJetOverlapCorrections.assign(4, ParticleObject::LorentzVector_t(0, 0, 0, 0));
  currentJetOverlapCorrections.at(2*preserveP4 + 1*hasJERShifts) = corr;
}

void METObject::getPtPhi(float& pt, float& phi, bool addXYShifts, bool addJERShifts, bool addParticleShifts, bool preserveP4) const{
  using namespace SystematicsHelpers;

  ParticleObject::LorentzVector_t tmp_p4;
  if (currentSyst != sUncorrected){
    float const& pt_ref = extras.met_Nominal;
    float const& phi_ref = extras.metPhi_Nominal;
    tmp_p4 = ParticleObject::LorentzVector_t(pt_ref*std::cos(phi_ref), pt_ref*std::sin(phi_ref), 0., 0.);

    if (addJERShifts && preserveP4) tmp_p4 += currentMETShift_p4Preserved;
    else if (!addJERShifts && preserveP4) tmp_p4 += currentMETShift_p4Preserved_noJER;
    else if (addJERShifts && !preserveP4) tmp_p4 += currentMETShift;
    else tmp_p4 += currentMETShift_noJER;

    if (addXYShifts) tmp_p4 += currentXYshift;
    if (addParticleShifts) tmp_p4 += particleMomentumCorrections;

    if (!currentMETCorrections.empty()) tmp_p4 += currentMETCorrections.at(8*preserveP4 + 4*addParticleShifts + 2*addJERShifts + 1*addXYShifts);
    if (!currentJetOverlapCorrections.empty()) tmp_p4 += currentMETCorrections.at(2*preserveP4 + 1*addJERShifts);
  }
  else{
    float const& pt_ref = extras.met_Raw;
    float const& phi_ref = extras.metPhi_Raw;
    tmp_p4 = ParticleObject::LorentzVector_t(pt_ref*std::cos(phi_ref), pt_ref*std::sin(phi_ref), 0., 0.);
  }

  pt = tmp_p4.Pt();
  phi = tmp_p4.Phi();
}
ParticleObject::LorentzVector_t::Scalar METObject::met(bool addXYShifts, bool addJERShifts, bool addParticleShifts, bool preserveP4) const{
  float pt, phi;
  getPtPhi(pt, phi, addXYShifts, addJERShifts, addParticleShifts, preserveP4);
  return pt;
}
ParticleObject::LorentzVector_t::Scalar METObject::phi(bool addXYShifts, bool addJERShifts, bool addParticleShifts, bool preserveP4) const{
  float pt, phi;
  getPtPhi(pt, phi, addXYShifts, addJERShifts, addParticleShifts, preserveP4);
  return phi;
}
ParticleObject::LorentzVector_t::Scalar METObject::px(bool addXYShifts, bool addJERShifts, bool addParticleShifts, bool preserveP4, float phi_rot) const{
  float pt, phi;
  getPtPhi(pt, phi, addXYShifts, addJERShifts, addParticleShifts, preserveP4);
  return pt * std::cos(phi + phi_rot);
}
ParticleObject::LorentzVector_t::Scalar METObject::py(bool addXYShifts, bool addJERShifts, bool addParticleShifts, bool preserveP4, float phi_rot) const{
  float pt, phi;
  getPtPhi(pt, phi, addXYShifts, addJERShifts, addParticleShifts, preserveP4);
  return pt * std::sin(phi + phi_rot);
}
ParticleObject::LorentzVector_t METObject::p4(bool addXYShifts, bool addJERShifts, bool addParticleShifts, bool preserveP4, float phi_rot) const{
  float pt, phi;
  getPtPhi(pt, phi, addXYShifts, addJERShifts, addParticleShifts, preserveP4);
  return ParticleObject::LorentzVector_t(pt * std::cos(phi + phi_rot), pt * std::sin(phi + phi_rot), 0., pt);
}
