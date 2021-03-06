#ifndef TRIGGEROBJECT_H
#define TRIGGEROBJECT_H

#include <vector>
#include <utility>
#include <DataFormats/HLTReco/interface/TriggerTypeDefs.h>
#include <CMS3/Dictionaries/interface/EgammaFiduciality.h>
#include "ParticleObject.h"
#include "ParticleObjectHelpers.h"


#define TRIGGEROBJECT_EXTRA_VARIABLES \
TRIGGEROBJECT_VARIABLE(std::vector<cms3_triggerIndex_t>, associatedTriggers, std::vector<cms3_triggerIndex_t>()) \
TRIGGEROBJECT_VARIABLE(std::vector<cms3_triggerIndex_t>, passedTriggers, std::vector<cms3_triggerIndex_t>())
#define TRIGGEROBJECT_MOMENTUM_VARIABLES \
TRIGGEROBJECT_VARIABLE(cms3_triggertype_t, type, 0) \
TRIGGEROBJECT_VARIABLE(float, pt, 0) \
TRIGGEROBJECT_VARIABLE(float, eta, 0) \
TRIGGEROBJECT_VARIABLE(float, phi, 0) \
TRIGGEROBJECT_VARIABLE(float, mass, 0)

#define TRIGGEROBJECT_VARIABLES \
TRIGGEROBJECT_EXTRA_VARIABLES \
TRIGGEROBJECT_MOMENTUM_VARIABLES


class TriggerObjectVariables{
public:
#define TRIGGEROBJECT_VARIABLE(TYPE, NAME, DEFVAL) TYPE NAME;
  TRIGGEROBJECT_VARIABLES;
#undef TRIGGEROBJECT_VARIABLE

  TriggerObjectVariables();
  TriggerObjectVariables(TriggerObjectVariables const& other);
  TriggerObjectVariables& operator=(const TriggerObjectVariables& other);

  void swap(TriggerObjectVariables& other);

};

class TriggerObject : public ParticleObject{
public:
  TriggerObjectVariables extras;

  TriggerObject();
  TriggerObject(cms3_triggertype_t id_);
  TriggerObject(cms3_triggertype_t id_, LorentzVector_t const& mom_);
  TriggerObject(const TriggerObject& other);
  TriggerObject& operator=(const TriggerObject& other);
  ~TriggerObject(){}

  void swap(TriggerObject& other);

  trigger::TriggerObjectType getTriggerObjectType() const{ return (trigger::TriggerObjectType) id; }
  bool isTriggerObjectType(trigger::TriggerObjectType const& TOtype) const{ return (this->getTriggerObjectType() == TOtype); }

  template<typename T> static void getMatchedPhysicsObjects(
    std::vector<TriggerObject const*> const& allPassedTOs, std::vector<trigger::TriggerObjectType> const& TOreqs, double const& dRmatch,
    std::vector<T*> const& objlist, std::vector<T*>& res,
    std::vector< std::pair<T*, TriggerObject const*> >* matchingLegs=nullptr
  );

};

template<typename T> void TriggerObject::getMatchedPhysicsObjects(
  std::vector<TriggerObject const*> const& allPassedTOs, std::vector<trigger::TriggerObjectType> const& TOreqs, double const& dRmatch,
  std::vector<T*> const& objlist, std::vector<T*>& res,
  std::vector< std::pair<T*, TriggerObject const*> >* matchingLegs
){
  res.reserve(objlist.size());
  if (matchingLegs) matchingLegs->reserve(objlist.size());

  std::vector<TriggerObject const*> hltFlaggedObjs;
  for (auto const& to_obj:allPassedTOs){
    for (trigger::TriggerObjectType const& TOreq:TOreqs){
      if (to_obj->isTriggerObjectType(TOreq)){ hltFlaggedObjs.push_back(to_obj); break; }
    }
  }

  std::unordered_map<TriggerObject const*, T*> TO_physobj_map;
  ParticleObjectHelpers::matchParticles(
    ParticleObjectHelpers::kMatchBy_DeltaR, dRmatch,
    hltFlaggedObjs.cbegin(), hltFlaggedObjs.cend(),
    objlist.cbegin(), objlist.cend(),
    TO_physobj_map
  );
  for (auto const& it:TO_physobj_map){
    if (it.second){
      res.push_back(it.second);
      if (matchingLegs) matchingLegs->emplace_back(it.second, it.first);
    }
  }

  ParticleObjectHelpers::sortByGreaterPt(res);
}


#endif
