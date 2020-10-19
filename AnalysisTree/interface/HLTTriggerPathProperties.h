#ifndef CMS3_HLTTRIGGERPATHPROPERTIES_H
#define CMS3_HLTTRIGGERPATHPROPERTIES_H

#include <string>
#include <vector>
#include <unordered_map>
#include "HLTObjectProperties.h"
#include "MuonObject.h"
#include "ElectronObject.h"
#include "PhotonObject.h"
#include "AK4JetObject.h"
#include "AK8JetObject.h"
#include "TNumericUtil.hh"


class HLTTriggerPathProperties{
protected:
  std::string name;
  std::unordered_map< HLTObjectProperties::TriggerObjectType, std::vector<HLTObjectProperties> > triggerObjectProperties;
  bool hasCompositeFilters;

  std::vector< std::pair<unsigned int, unsigned int> > excluded_runRange_list;

  void setupName();

public:
  HLTTriggerPathProperties();
  HLTTriggerPathProperties(std::string const& name_);
  HLTTriggerPathProperties(std::string const& name_, std::vector<HLTObjectProperties> const& triggerObjectProperties_);
  HLTTriggerPathProperties(HLTTriggerPathProperties const& other);

  void setup();
  void addObjectProperties(HLTObjectProperties const& props);
  void setExcludedRunRanges(std::vector< std::pair<unsigned int, unsigned int> > const& rangelist){ excluded_runRange_list = rangelist; }

  void resetCuts();
  void resetExcludedRunRanges(){ excluded_runRange_list.clear(); }

  bool testCuts(
    std::vector<MuonObject const*> const& muons,
    std::vector<ElectronObject const*> const& electrons,
    std::vector<PhotonObject const*> const& photons,
    std::vector<AK4JetObject const*> const& ak4jets,
    std::vector<AK8JetObject const*> const& ak8jets,
    ParticleObject::LorentzVector_t const& pfmet_p4,
    ParticleObject::LorentzVector_t const& pfmet_nomu_p4,
    ParticleObject::LorentzVector_t const& ht_p4,
    ParticleObject::LorentzVector_t const& ht_nomu_p4
  ) const;
  template<typename T> bool testCutSet(HLTObjectProperties::TriggerObjectType const& type, std::vector<HLTObjectProperties> const& props, std::vector<T const*> const& objects, std::vector<std::vector<T const*>>* perms_objects_passing = nullptr) const;
  bool testRun(unsigned int const& run) const;

  bool isSameTrigger(std::string const& name_) const;
  bool hasExcludedRunRanges() const{ return !excluded_runRange_list.empty(); }

  std::string const& getName() const{ return name; }
  std::unordered_map< HLTObjectProperties::TriggerObjectType, std::vector<HLTObjectProperties> > const& getObjectProperties() const{ return triggerObjectProperties; }

};

template<typename T> bool HLTTriggerPathProperties::testCutSet(HLTObjectProperties::TriggerObjectType const& type, std::vector<HLTObjectProperties> const& props, std::vector<T const*> const& objects, std::vector<std::vector<T const*>>* perms_objects_passing) const{
  if (objects.size()<props.size()) return false;

  std::vector<std::vector<int>> perms;
  TNumericUtil::PermutationGenerator(objects.size(), props.size(), perms, 0, 1);

  bool passAnyPerm = false;
  for (auto const& perm:perms){
    bool res = true;
    unsigned int iprop=0;
    std::vector<T const*> perm_objects_passing;
    if (perms_objects_passing) perm_objects_passing.reserve(perm.size());

    for (auto const& ipos:perm){
      bool cutres = props.at(iprop).testCuts(objects.at(ipos)->p4(), type);
      if (perms_objects_passing && cutres) perm_objects_passing.push_back(objects.at(ipos));
      res = res && cutres;
      iprop++;
    }

    if (res){
      if (perms_objects_passing) perms_objects_passing->push_back(perm_objects_passing);
      else return res;
    }
    passAnyPerm |= res;
  }
  return passAnyPerm;
}
template bool HLTTriggerPathProperties::testCutSet<MuonObject>(HLTObjectProperties::TriggerObjectType const& type, std::vector<HLTObjectProperties> const& props, std::vector<MuonObject const*> const& objects, std::vector<std::vector<MuonObject const*>>* perms_objects_passing) const;
template bool HLTTriggerPathProperties::testCutSet<ElectronObject>(HLTObjectProperties::TriggerObjectType const& type, std::vector<HLTObjectProperties> const& props, std::vector<ElectronObject const*> const& objects, std::vector<std::vector<ElectronObject const*>>* perms_objects_passing) const;
template bool HLTTriggerPathProperties::testCutSet<PhotonObject>(HLTObjectProperties::TriggerObjectType const& type, std::vector<HLTObjectProperties> const& props, std::vector<PhotonObject const*> const& objects, std::vector<std::vector<PhotonObject const*>>* perms_objects_passing) const;
template bool HLTTriggerPathProperties::testCutSet<AK4JetObject>(HLTObjectProperties::TriggerObjectType const& type, std::vector<HLTObjectProperties> const& props, std::vector<AK4JetObject const*> const& objects, std::vector<std::vector<AK4JetObject const*>>* perms_objects_passing) const;
template bool HLTTriggerPathProperties::testCutSet<AK8JetObject>(HLTObjectProperties::TriggerObjectType const& type, std::vector<HLTObjectProperties> const& props, std::vector<AK8JetObject const*> const& objects, std::vector<std::vector<AK8JetObject const*>>* perms_objects_passing) const;
template bool HLTTriggerPathProperties::testCutSet<ParticleObject>(HLTObjectProperties::TriggerObjectType const& type, std::vector<HLTObjectProperties> const& props, std::vector<ParticleObject const*> const& objects, std::vector<std::vector<ParticleObject const*>>* perms_objects_passing) const;


#endif
