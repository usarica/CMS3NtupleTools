#include <algorithm>
#include <utility>
#include "LHEParticleObject.h"


LHEParticleObject::LHEParticleObject() :
  ParticleObject(),
  st(0)
{}
LHEParticleObject::LHEParticleObject(int id_, int st_) :
  ParticleObject(id_),
  st(st_)
{}
LHEParticleObject::LHEParticleObject(int id_, int st_, LorentzVector_t const& momentum_) :
  ParticleObject(id_, momentum_),
  st(st_)
{}
LHEParticleObject::LHEParticleObject(const LHEParticleObject& other) :
  ParticleObject(other),
  st(other.st)
{}
void LHEParticleObject::swap(LHEParticleObject& other){
  ParticleObject::swap(other);
  std::swap(st, other.st);
}
LHEParticleObject& LHEParticleObject::operator=(const LHEParticleObject& other){
  LHEParticleObject tmp(other);
  swap(tmp);
  return *this;
}
LHEParticleObject::~LHEParticleObject(){}
