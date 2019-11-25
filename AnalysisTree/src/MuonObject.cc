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
  extras()
{}
MuonObject::MuonObject(int id_) :
  ParticleObject(id_),
  extras()
{}
MuonObject::MuonObject(int id_, LorentzVector_t const& momentum_) :
  ParticleObject(id_, momentum_),
  extras()
{}
MuonObject::MuonObject(const MuonObject& other) :
  ParticleObject(other),
  extras(other.extras)
{}
void MuonObject::swap(MuonObject& other){
  std::swap(id, other.id);
  std::swap(selectionBits, other.selectionBits);
  std::swap(momentum, other.momentum);
  extras.swap(other.extras);
}
MuonObject& MuonObject::operator=(const MuonObject& other){
  MuonObject tmp(other);
  swap(tmp);
  return *this;
}
MuonObject::~MuonObject(){}
