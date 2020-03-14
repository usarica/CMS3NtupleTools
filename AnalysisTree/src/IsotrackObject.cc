#include <algorithm>
#include <utility>
#include "IsotrackObject.h"


IsotrackVariables::IsotrackVariables(){
#define ISOTRACK_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=DEFVAL;
  ISOTRACK_VARIABLES;
#undef ISOTRACK_VARIABLE
}
IsotrackVariables::IsotrackVariables(IsotrackVariables const& other){
#define ISOTRACK_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=other.NAME;
  ISOTRACK_VARIABLES;
#undef ISOTRACK_VARIABLE
}
void IsotrackVariables::swap(IsotrackVariables& other){
#define ISOTRACK_VARIABLE(TYPE, NAME, DEFVAL) std::swap(this->NAME, other.NAME);
  ISOTRACK_VARIABLES;
#undef ISOTRACK_VARIABLE
}
IsotrackVariables& IsotrackVariables::operator=(const IsotrackVariables& other){
  IsotrackVariables tmp(other);
  swap(tmp);
  return *this;
}

IsotrackObject::IsotrackObject() :
  ParticleObject(),
  extras()
{}
IsotrackObject::IsotrackObject(cms3_id_t id_) :
  ParticleObject(id_),
  extras()
{}
IsotrackObject::IsotrackObject(cms3_id_t id_, LorentzVector_t const& momentum_) :
  ParticleObject(id_, momentum_),
  extras()
{}
IsotrackObject::IsotrackObject(const IsotrackObject& other) :
  ParticleObject(other),
  extras(other.extras)
{}
void IsotrackObject::swap(IsotrackObject& other){
  ParticleObject::swap(other);
  extras.swap(other.extras);
}
IsotrackObject& IsotrackObject::operator=(const IsotrackObject& other){
  IsotrackObject tmp(other);
  swap(tmp);
  return *this;
}
IsotrackObject::~IsotrackObject(){}
