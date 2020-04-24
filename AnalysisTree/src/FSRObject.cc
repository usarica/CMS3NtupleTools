#include <algorithm>
#include <utility>
#include "FSRObject.h"


FSRVariables::FSRVariables(){
#define FSR_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=DEFVAL;
#define FSR_VECTOR_VARIABLE(TYPE, NAME) this->NAME=TYPE();
  FSR_VARIABLES;
#undef FSR_VECTOR_VARIABLE
#undef FSR_VARIABLE
}
FSRVariables::FSRVariables(FSRVariables const& other){
#define FSR_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=other.NAME;
#define FSR_VECTOR_VARIABLE(TYPE, NAME) std::copy(other.NAME.cbegin(), other.NAME.cend(), this->NAME.begin());
  FSR_VARIABLES;
#undef FSR_VECTOR_VARIABLE
#undef FSR_VARIABLE
}
void FSRVariables::swap(FSRVariables& other){
#define FSR_VARIABLE(TYPE, NAME, DEFVAL) std::swap(this->NAME, other.NAME);
#define FSR_VECTOR_VARIABLE(TYPE, NAME) std::swap(this->NAME, other.NAME);
  FSR_VARIABLES;
#undef FSR_VECTOR_VARIABLE
#undef FSR_VARIABLE
}
FSRVariables& FSRVariables::operator=(const FSRVariables& other){
  FSRVariables tmp(other);
  swap(tmp);
  return *this;
}

FSRObject::FSRObject() :
  ParticleObject(),
  extras(),
  associatedPhoton(nullptr)
{}
FSRObject::FSRObject(LorentzVector_t const& momentum_) :
  ParticleObject(22, momentum_),
  extras(),
  associatedPhoton(nullptr)
{}
FSRObject::FSRObject(const FSRObject& other) :
  ParticleObject(other),
  extras(other.extras),
  associatedPhoton(other.associatedPhoton)
{}
void FSRObject::swap(FSRObject& other){
  ParticleObject::swap(other);
  extras.swap(other.extras);
  std::swap(associatedPhoton, other.associatedPhoton);
}
FSRObject& FSRObject::operator=(const FSRObject& other){
  FSRObject tmp(other);
  swap(tmp);
  return *this;
}
FSRObject::~FSRObject(){}
