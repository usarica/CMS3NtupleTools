#ifndef OVERLAPMAPELEMENT_H
#define OVERLAPMAPELEMENT_H

#include <utility>
#include <algorithm>

#include "MuonObject.h"
#include "ElectronObject.h"
#include "PhotonObject.h"
#include "AK4JetObject.h"
#include "AK8JetObject.h"


class OverlapMapElementBase{
protected:
  std::pair<cms3_listIndex_signed_long_t, cms3_listIndex_signed_long_t> index_pair;

public:
  OverlapMapElementBase() : index_pair(-1,-1){}
  OverlapMapElementBase(OverlapMapElementBase const& other) : index_pair(other.index_pair){}
  virtual ~OverlapMapElementBase(){}

  bool isValid() const{ return index_pair.first>=0 && index_pair.second>=0; }

  void setFirstParticleIndex(cms3_listIndex_short_t const& idx){ index_pair.first = idx; }
  void setSecondParticleIndex(cms3_listIndex_short_t const& idx){ index_pair.second = idx; }

  std::pair<cms3_listIndex_signed_long_t, cms3_listIndex_signed_long_t>& getIndices(){ return index_pair; }
  std::pair<cms3_listIndex_signed_long_t, cms3_listIndex_signed_long_t> const& getIndices() const{ return index_pair; }

};

// Overlap map element extras
template<typename T, typename U> struct OverlapMapElementExtras{
  OverlapMapElementExtras(){}
  OverlapMapElementExtras(OverlapMapElementExtras<T, U> const& other){}
  virtual ~OverlapMapElementExtras(){}

  void swap(OverlapMapElementExtras& other){}
  OverlapMapElementExtras<T, U>& operator=(OverlapMapElementExtras<T, U> const& other){
    OverlapMapElementExtras<T, U> tmp(other);
    swap(tmp);
    return *this;
  }
};

#define OVERLAPMAP_PARTICLES_JETS_UNLINKED_VARIABLES \
OVERLAPMAP_VARIABLE(cms3_listIndex_short_t, jet_match_index, 0) \
OVERLAPMAP_VARIABLE(cms3_listIndex_short_t, particle_match_index, 0)
#define OVERLAPMAP_PARTICLES_JETS_LINKED_VARIABLES \
OVERLAPMAP_VARIABLE(float, commonPFCandidates_sump4_pt, 0) \
OVERLAPMAP_VARIABLE(float, commonPFCandidates_sump4_eta, 0) \
OVERLAPMAP_VARIABLE(float, commonPFCandidates_sump4_phi, 0) \
OVERLAPMAP_VARIABLE(float, commonPFCandidates_sump4_mass, 0)
#define OVERLAPMAP_EGAMMAS_JETS_LINKED_VARIABLES \
OVERLAPMAP_VARIABLE(float, commonGoodMETPFMuons_sump4_px, 0) \
OVERLAPMAP_VARIABLE(float, commonGoodMETPFMuons_sump4_py, 0)

#define OVERLAPMAP_MUONS_JETS_UNLINKED_VARIABLES \
OVERLAPMAP_PARTICLES_JETS_UNLINKED_VARIABLES
#define OVERLAPMAP_MUONS_JETS_LINKED_VARIABLES \
OVERLAPMAP_PARTICLES_JETS_LINKED_VARIABLES
#define OVERLAPMAP_MUONS_JETS_VARIABLES \
OVERLAPMAP_MUONS_JETS_UNLINKED_VARIABLES \
OVERLAPMAP_MUONS_JETS_LINKED_VARIABLES

#define OVERLAPMAP_ELECTRONS_JETS_UNLINKED_VARIABLES \
OVERLAPMAP_PARTICLES_JETS_UNLINKED_VARIABLES
#define OVERLAPMAP_ELECTRONS_JETS_LINKED_VARIABLES \
OVERLAPMAP_PARTICLES_JETS_LINKED_VARIABLES \
OVERLAPMAP_EGAMMAS_JETS_LINKED_VARIABLES
#define OVERLAPMAP_ELECTRONS_JETS_VARIABLES \
OVERLAPMAP_ELECTRONS_JETS_UNLINKED_VARIABLES \
OVERLAPMAP_ELECTRONS_JETS_LINKED_VARIABLES

#define OVERLAPMAP_PHOTONS_JETS_UNLINKED_VARIABLES \
OVERLAPMAP_PARTICLES_JETS_UNLINKED_VARIABLES
#define OVERLAPMAP_PHOTONS_JETS_LINKED_VARIABLES \
OVERLAPMAP_PARTICLES_JETS_LINKED_VARIABLES \
OVERLAPMAP_EGAMMAS_JETS_LINKED_VARIABLES
#define OVERLAPMAP_PHOTONS_JETS_VARIABLES \
OVERLAPMAP_PHOTONS_JETS_UNLINKED_VARIABLES \
OVERLAPMAP_PHOTONS_JETS_LINKED_VARIABLES

template<> struct OverlapMapElementExtras<MuonObject, AK4JetObject>{
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) TYPE NAME;
  OVERLAPMAP_MUONS_JETS_LINKED_VARIABLES;
#undef OVERLAPMAP_VARIABLE

  OverlapMapElementExtras<MuonObject, AK4JetObject>(){
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) NAME = DEFVAL;
    OVERLAPMAP_MUONS_JETS_LINKED_VARIABLES;
#undef OVERLAPMAP_VARIABLE
  }
  OverlapMapElementExtras<MuonObject, AK4JetObject>(OverlapMapElementExtras<MuonObject, AK4JetObject> const& other){
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) NAME = other.NAME;
    OVERLAPMAP_MUONS_JETS_LINKED_VARIABLES;
#undef OVERLAPMAP_VARIABLE
  }
  void swap(OverlapMapElementExtras<MuonObject, AK4JetObject>& other){
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) std::swap(NAME, other.NAME);
    OVERLAPMAP_MUONS_JETS_LINKED_VARIABLES;
#undef OVERLAPMAP_VARIABLE
  }
  OverlapMapElementExtras<MuonObject, AK4JetObject>& operator=(const OverlapMapElementExtras<MuonObject, AK4JetObject>& other){
    OverlapMapElementExtras<MuonObject, AK4JetObject> tmp(other);
    swap(tmp);
    return *this;
  }
};
template<> struct OverlapMapElementExtras<MuonObject, AK8JetObject>{
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) TYPE NAME;
  OVERLAPMAP_MUONS_JETS_LINKED_VARIABLES;
#undef OVERLAPMAP_VARIABLE

  OverlapMapElementExtras<MuonObject, AK8JetObject>(){
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) NAME = DEFVAL;
    OVERLAPMAP_MUONS_JETS_LINKED_VARIABLES;
#undef OVERLAPMAP_VARIABLE
  }
  OverlapMapElementExtras<MuonObject, AK8JetObject>(OverlapMapElementExtras<MuonObject, AK8JetObject> const& other){
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) NAME = other.NAME;
    OVERLAPMAP_MUONS_JETS_LINKED_VARIABLES;
#undef OVERLAPMAP_VARIABLE
  }
  void swap(OverlapMapElementExtras<MuonObject, AK8JetObject>& other){
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) std::swap(NAME, other.NAME);
    OVERLAPMAP_MUONS_JETS_LINKED_VARIABLES;
#undef OVERLAPMAP_VARIABLE
  }
  OverlapMapElementExtras<MuonObject, AK8JetObject>& operator=(const OverlapMapElementExtras<MuonObject, AK8JetObject>& other){
    OverlapMapElementExtras<MuonObject, AK8JetObject> tmp(other);
    swap(tmp);
    return *this;
  }
};

template<> struct OverlapMapElementExtras<ElectronObject, AK4JetObject>{
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) TYPE NAME;
  OVERLAPMAP_ELECTRONS_JETS_LINKED_VARIABLES;
#undef OVERLAPMAP_VARIABLE

  OverlapMapElementExtras<ElectronObject, AK4JetObject>(){
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) NAME = DEFVAL;
    OVERLAPMAP_ELECTRONS_JETS_LINKED_VARIABLES;
#undef OVERLAPMAP_VARIABLE
  }
  OverlapMapElementExtras<ElectronObject, AK4JetObject>(OverlapMapElementExtras<ElectronObject, AK4JetObject> const& other){
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) NAME = other.NAME;
    OVERLAPMAP_ELECTRONS_JETS_LINKED_VARIABLES;
#undef OVERLAPMAP_VARIABLE
  }
  void swap(OverlapMapElementExtras<ElectronObject, AK4JetObject>& other){
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) std::swap(NAME, other.NAME);
    OVERLAPMAP_ELECTRONS_JETS_LINKED_VARIABLES;
#undef OVERLAPMAP_VARIABLE
  }
  OverlapMapElementExtras<ElectronObject, AK4JetObject>& operator=(const OverlapMapElementExtras<ElectronObject, AK4JetObject>& other){
    OverlapMapElementExtras<ElectronObject, AK4JetObject> tmp(other);
    swap(tmp);
    return *this;
  }
};
template<> struct OverlapMapElementExtras<ElectronObject, AK8JetObject>{
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) TYPE NAME;
  OVERLAPMAP_ELECTRONS_JETS_LINKED_VARIABLES;
#undef OVERLAPMAP_VARIABLE

  OverlapMapElementExtras<ElectronObject, AK8JetObject>(){
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) NAME = DEFVAL;
    OVERLAPMAP_ELECTRONS_JETS_LINKED_VARIABLES;
#undef OVERLAPMAP_VARIABLE
  }
  OverlapMapElementExtras<ElectronObject, AK8JetObject>(OverlapMapElementExtras<ElectronObject, AK8JetObject> const& other){
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) NAME = other.NAME;
    OVERLAPMAP_ELECTRONS_JETS_LINKED_VARIABLES;
#undef OVERLAPMAP_VARIABLE
  }
  void swap(OverlapMapElementExtras<ElectronObject, AK8JetObject>& other){
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) std::swap(NAME, other.NAME);
    OVERLAPMAP_ELECTRONS_JETS_LINKED_VARIABLES;
#undef OVERLAPMAP_VARIABLE
  }
  OverlapMapElementExtras<ElectronObject, AK8JetObject>& operator=(const OverlapMapElementExtras<ElectronObject, AK8JetObject>& other){
    OverlapMapElementExtras<ElectronObject, AK8JetObject> tmp(other);
    swap(tmp);
    return *this;
  }
};

template<> struct OverlapMapElementExtras<PhotonObject, AK4JetObject>{
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) TYPE NAME;
  OVERLAPMAP_PHOTONS_JETS_LINKED_VARIABLES;
#undef OVERLAPMAP_VARIABLE

  OverlapMapElementExtras<PhotonObject, AK4JetObject>(){
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) NAME = DEFVAL;
    OVERLAPMAP_PHOTONS_JETS_LINKED_VARIABLES;
#undef OVERLAPMAP_VARIABLE
  }
  OverlapMapElementExtras<PhotonObject, AK4JetObject>(OverlapMapElementExtras<PhotonObject, AK4JetObject> const& other){
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) NAME = other.NAME;
    OVERLAPMAP_PHOTONS_JETS_LINKED_VARIABLES;
#undef OVERLAPMAP_VARIABLE
  }
  void swap(OverlapMapElementExtras<PhotonObject, AK4JetObject>& other){
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) std::swap(NAME, other.NAME);
    OVERLAPMAP_PHOTONS_JETS_LINKED_VARIABLES;
#undef OVERLAPMAP_VARIABLE
  }
  OverlapMapElementExtras<PhotonObject, AK4JetObject>& operator=(const OverlapMapElementExtras<PhotonObject, AK4JetObject>& other){
    OverlapMapElementExtras<PhotonObject, AK4JetObject> tmp(other);
    swap(tmp);
    return *this;
  }
};
template<> struct OverlapMapElementExtras<PhotonObject, AK8JetObject>{
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) TYPE NAME;
  OVERLAPMAP_PHOTONS_JETS_LINKED_VARIABLES;
#undef OVERLAPMAP_VARIABLE

  OverlapMapElementExtras<PhotonObject, AK8JetObject>(){
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) NAME = DEFVAL;
    OVERLAPMAP_PHOTONS_JETS_LINKED_VARIABLES;
#undef OVERLAPMAP_VARIABLE
  }
  OverlapMapElementExtras<PhotonObject, AK8JetObject>(OverlapMapElementExtras<PhotonObject, AK8JetObject> const& other){
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) NAME = other.NAME;
    OVERLAPMAP_PHOTONS_JETS_LINKED_VARIABLES;
#undef OVERLAPMAP_VARIABLE
  }
  void swap(OverlapMapElementExtras<PhotonObject, AK8JetObject>& other){
#define OVERLAPMAP_VARIABLE(TYPE, NAME, DEFVAL) std::swap(NAME, other.NAME);
    OVERLAPMAP_PHOTONS_JETS_LINKED_VARIABLES;
#undef OVERLAPMAP_VARIABLE
  }
  OverlapMapElementExtras<PhotonObject, AK8JetObject>& operator=(const OverlapMapElementExtras<PhotonObject, AK8JetObject>& other){
    OverlapMapElementExtras<PhotonObject, AK8JetObject> tmp(other);
    swap(tmp);
    return *this;
  }
};


template<typename T, typename U> class OverlapMapElement : public OverlapMapElementBase, public OverlapMapElementExtras<T, U>{
protected:
  std::pair<T*, U*> linkedElementPair;

public:
  OverlapMapElement() : OverlapMapElementBase(), OverlapMapElementExtras<T, U>(), linkedElementPair(nullptr, nullptr){}
  OverlapMapElement(OverlapMapElement const& other) : OverlapMapElementBase(other), OverlapMapElementExtras<T, U>(other), linkedElementPair(other.linkedElementPair){}
  virtual ~OverlapMapElement(){}

  void swap(OverlapMapElement<T, U>& other){
    std::swap(index_pair, other.index_pair);
    OverlapMapElementExtras<T, U>::swap(other);
    std::swap(linkedElementPair, other.linkedElementPair);
  }
  OverlapMapElement<T, U>& operator=(OverlapMapElement<T, U> const& other){
    OverlapMapElement<T, U> tmp(other);
    swap(tmp);
    return *this;
  }

  bool isLinked() const{ return (linkedElementPair.first && linkedElementPair.second); }

  bool linkFirstElement(std::vector<T*> const& elist);
  bool linkSecondElement(std::vector<U*> const& elist);

  bool hasIdenticalElements(T* firstElement, U* secondElement) const;

  ParticleObject::LorentzVector_t p4_common() const{ return ParticleObject::LorentzVector_t(0, 0, 0, 0); }
  ParticleObject::LorentzVector_t p4_commonMuCands_goodMET() const{ return ParticleObject::LorentzVector_t(0, 0, 0, 0); }

};
template<typename T, typename U> bool OverlapMapElement<T, U>::hasIdenticalElements(T* firstElement, U* secondElement) const{
  if (this->isLinked()) return (firstElement == linkedElementPair.first && secondElement == linkedElementPair.second);
  else if (firstElement && secondElement && this->isValid()) return (
    (ParticleObject::UniqueId_t) index_pair.first == firstElement->getUniqueIdentifier()
    &&
    (ParticleObject::UniqueId_t) index_pair.second == secondElement->getUniqueIdentifier()
    );
  else return false;
}

template<typename T, typename U> bool OverlapMapElement<T, U>::linkFirstElement(std::vector<T*> const& elist){
  if (index_pair.first<0) return false;
  for (auto const& ee:elist){
    if ((ParticleObject::UniqueId_t) index_pair.first == ee->getUniqueIdentifier()){
      linkedElementPair.first = ee;
      return true;
    }
  }
  return false;
}
template<typename T, typename U> bool OverlapMapElement<T, U>::linkSecondElement(std::vector<U*> const& elist){
  if (index_pair.second<0) return false;
  for (auto const& ee:elist){
    if ((ParticleObject::UniqueId_t) index_pair.second == ee->getUniqueIdentifier()){
      linkedElementPair.second = ee;
      return true;
    }
  }
  return false;
}


#define OVERLAPMAP_SPECIALIZATIONS \
OVERLAPMAP_SPECIALIZATION(MuonObject, AK4JetObject) \
OVERLAPMAP_SPECIALIZATION(MuonObject, AK8JetObject) \
OVERLAPMAP_SPECIALIZATION(ElectronObject, AK4JetObject) \
OVERLAPMAP_SPECIALIZATION(ElectronObject, AK8JetObject) \
OVERLAPMAP_SPECIALIZATION(PhotonObject, AK4JetObject) \
OVERLAPMAP_SPECIALIZATION(PhotonObject, AK8JetObject)

#define OVERLAPMAP_SPECIALIZATION(T1, T2) \
template<> ParticleObject::LorentzVector_t OverlapMapElement<T1, T2>::p4_common() const; \
template<> ParticleObject::LorentzVector_t OverlapMapElement<T1, T2>::p4_commonMuCands_goodMET() const;

OVERLAPMAP_SPECIALIZATIONS;

#undef OVERLAPMAP_SPECIALIZATION
#undef OVERLAPMAP_SPECIALIZATIONS


#endif
