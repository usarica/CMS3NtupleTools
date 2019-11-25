#include <algorithm>
#include <utility>
#include "PhotonObject.h"


PhotonVariables::PhotonVariables(){
#define PHOTON_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=DEFVAL;
  PHOTON_VARIABLES;
#undef PHOTON_VARIABLE
}
PhotonVariables::PhotonVariables(PhotonVariables const& other){
#define PHOTON_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=other.NAME;
  PHOTON_VARIABLES;
#undef PHOTON_VARIABLE
}
void PhotonVariables::swap(PhotonVariables& other){
#define PHOTON_VARIABLE(TYPE, NAME, DEFVAL) std::swap(this->NAME, other.NAME);
  PHOTON_VARIABLES;
#undef PHOTON_VARIABLE
}
PhotonVariables& PhotonVariables::operator=(const PhotonVariables& other){
  PhotonVariables tmp(other);
  swap(tmp);
  return *this;
}

PhotonObject::PhotonObject() :
  ParticleObject(),
  extras()
{}
PhotonObject::PhotonObject(int id_) :
  ParticleObject(id_),
  extras()
{}
PhotonObject::PhotonObject(int id_, LorentzVector_t const& momentum_) :
  ParticleObject(id_, momentum_),
  extras()
{}
PhotonObject::PhotonObject(const PhotonObject& other) :
  ParticleObject(other),
  extras(other.extras)
{}
void PhotonObject::swap(PhotonObject& other){
  std::swap(id, other.id);
  std::swap(selectionBits, other.selectionBits);
  std::swap(momentum, other.momentum);
  extras.swap(other.extras);
}
PhotonObject& PhotonObject::operator=(const PhotonObject& other){
  PhotonObject tmp(other);
  swap(tmp);
  return *this;
}
PhotonObject::~PhotonObject(){}
