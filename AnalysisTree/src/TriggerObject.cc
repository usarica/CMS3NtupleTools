#include <algorithm>
#include <utility>
#include "TriggerObject.h"


TriggerObjectVariables::TriggerObjectVariables(){
#define TRIGGEROBJECT_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=DEFVAL;
  TRIGGEROBJECT_VARIABLES;
#undef TRIGGEROBJECT_VARIABLE
}
TriggerObjectVariables::TriggerObjectVariables(TriggerObjectVariables const& other){
#define TRIGGEROBJECT_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=other.NAME;
  TRIGGEROBJECT_VARIABLES;
#undef TRIGGEROBJECT_VARIABLE
}
void TriggerObjectVariables::swap(TriggerObjectVariables& other){
#define TRIGGEROBJECT_VARIABLE(TYPE, NAME, DEFVAL) std::swap(this->NAME, other.NAME);
  TRIGGEROBJECT_VARIABLES;
#undef TRIGGEROBJECT_VARIABLE
}
TriggerObjectVariables& TriggerObjectVariables::operator=(const TriggerObjectVariables& other){
  TriggerObjectVariables tmp(other);
  swap(tmp);
  return *this;
}

TriggerObject::TriggerObject() :
  ParticleObject(),
  extras()
{}
TriggerObject::TriggerObject(cms3_id_t id_) :
  ParticleObject(id_),
  extras()
{}
TriggerObject::TriggerObject(cms3_id_t id_, LorentzVector_t const& momentum_) :
  ParticleObject(id_, momentum_),
  extras()
{}
TriggerObject::TriggerObject(const TriggerObject& other) :
  ParticleObject(other),
  extras(other.extras)
{}
void TriggerObject::swap(TriggerObject& other){
  ParticleObject::swap(other);
  extras.swap(other.extras);
}
TriggerObject& TriggerObject::operator=(const TriggerObject& other){
  TriggerObject tmp(other);
  swap(tmp);
  return *this;
}
