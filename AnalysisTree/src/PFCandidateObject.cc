#include <algorithm>
#include <utility>
#include "PFCandidateObject.h"


PFCandidateVariables::PFCandidateVariables(){
#define PFCANDIDATE_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=DEFVAL;
  PFCANDIDATE_VARIABLES;
#undef PFCANDIDATE_VARIABLE
}
PFCandidateVariables::PFCandidateVariables(PFCandidateVariables const& other){
#define PFCANDIDATE_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=other.NAME;
  PFCANDIDATE_VARIABLES;
#undef PFCANDIDATE_VARIABLE
}
void PFCandidateVariables::swap(PFCandidateVariables& other){
#define PFCANDIDATE_VARIABLE(TYPE, NAME, DEFVAL) std::swap(this->NAME, other.NAME);
  PFCANDIDATE_VARIABLES;
#undef PFCANDIDATE_VARIABLE
}
PFCandidateVariables& PFCandidateVariables::operator=(const PFCandidateVariables& other){
  PFCandidateVariables tmp(other);
  swap(tmp);
  return *this;
}


PFCandidateObject::PFCandidateObject() :
  ParticleObject(),
  extras()
{}
PFCandidateObject::PFCandidateObject(cms3_id_t id_) :
  ParticleObject(id_),
  extras()
{}
PFCandidateObject::PFCandidateObject(cms3_id_t id_, LorentzVector_t const& momentum_) :
  ParticleObject(id_, momentum_),
  extras()
{}
PFCandidateObject::PFCandidateObject(const PFCandidateObject& other) :
  ParticleObject(other),
  extras(other.extras)
{}
void PFCandidateObject::swap(PFCandidateObject& other){
  ParticleObject::swap(other);
  extras.swap(other.extras);
}
PFCandidateObject& PFCandidateObject::operator=(const PFCandidateObject& other){
  PFCandidateObject tmp(other);
  swap(tmp);
  return *this;
}
PFCandidateObject::~PFCandidateObject(){}

bool PFCandidateObject::isGlobalMuon() const{
  auto const& qflags = extras.qualityFlag;
  return ((qflags & muonFlagsMask) >> muonFlagsShift) & 2;
}
bool PFCandidateObject::isStandaloneMuon() const{
  auto const& qflags = extras.qualityFlag;
  return ((qflags & muonFlagsMask) >> muonFlagsShift) & 1;
}
