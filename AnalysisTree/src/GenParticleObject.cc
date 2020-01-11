#include <algorithm>
#include <utility>
#include "GenParticleObject.h"


GenParticleVariables::GenParticleVariables(){
#define GENPARTICLE_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=DEFVAL;
  GENPARTICLE_EXTRA_VARIABLES;
#undef GENPARTICLE_VARIABLE
}
GenParticleVariables::GenParticleVariables(GenParticleVariables const& other){
#define GENPARTICLE_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=other.NAME;
  GENPARTICLE_EXTRA_VARIABLES;
#undef GENPARTICLE_VARIABLE
}
void GenParticleVariables::swap(GenParticleVariables& other){
#define GENPARTICLE_VARIABLE(TYPE, NAME, DEFVAL) std::swap(this->NAME, other.NAME);
  GENPARTICLE_EXTRA_VARIABLES;
#undef GENPARTICLE_VARIABLE
}
GenParticleVariables& GenParticleVariables::operator=(const GenParticleVariables& other){
  GenParticleVariables tmp(other);
  swap(tmp);
  return *this;
}

GenParticleObject::GenParticleObject() :
  ParticleObject(),
  st(0),
  extras()
{}
GenParticleObject::GenParticleObject(int id_, int st_) :
  ParticleObject(id_),
  st(st_),
  extras()
{}
GenParticleObject::GenParticleObject(int id_, int st_, LorentzVector_t const& momentum_) :
  ParticleObject(id_, momentum_),
  st(st_),
  extras()
{}
GenParticleObject::GenParticleObject(const GenParticleObject& other) :
  ParticleObject(other),
  st(other.st),
  extras(other.extras)
{}
void GenParticleObject::swap(GenParticleObject& other){
  ParticleObject::swap(other);
  std::swap(st, other.st);
  extras.swap(other.extras);
}
GenParticleObject& GenParticleObject::operator=(const GenParticleObject& other){
  GenParticleObject tmp(other);
  swap(tmp);
  return *this;
}
GenParticleObject::~GenParticleObject(){}
