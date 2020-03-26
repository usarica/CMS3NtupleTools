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

  void setupName();

public:
  HLTTriggerPathProperties();
  HLTTriggerPathProperties(std::string const& name_);
  HLTTriggerPathProperties(std::string const& name_, std::vector<HLTObjectProperties> const& triggerObjectProperties_);
  HLTTriggerPathProperties(HLTTriggerPathProperties const& other);

  void setup();
  void addObjectProperties(HLTObjectProperties const& props);

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
  template<typename T> bool testCutSet(HLTObjectProperties::TriggerObjectType const& type, std::vector<HLTObjectProperties> const& props, std::vector<T const*> const& objects) const;

  bool isSameTrigger(std::string const& name_) const;

  std::string const& getName() const{ return name; }
  std::unordered_map< HLTObjectProperties::TriggerObjectType, std::vector<HLTObjectProperties> > const& getObjectProperties() const{ return triggerObjectProperties; }

};

template<typename T> bool HLTTriggerPathProperties::testCutSet(HLTObjectProperties::TriggerObjectType const& type, std::vector<HLTObjectProperties> const& props, std::vector<T const*> const& objects) const{
  if (objects.size()<props.size()) return false;

  std::vector<std::vector<int>> perms;
  TNumericUtil::PermutationGenerator(objects.size(), props.size(), perms, 0, 1);

  for (auto const& perm:perms){
    bool res=true;
    unsigned int iprop=0;
    for (auto const& ipos:perm){
      res = res && props.at(iprop).testCuts(objects.at(ipos)->p4(), type);
      iprop++;
    }
    if (res) return res;
  }

  return false;
}
template bool HLTTriggerPathProperties::testCutSet<MuonObject>(HLTObjectProperties::TriggerObjectType const& type, std::vector<HLTObjectProperties> const& props, std::vector<MuonObject const*> const& objects) const;
template bool HLTTriggerPathProperties::testCutSet<ElectronObject>(HLTObjectProperties::TriggerObjectType const& type, std::vector<HLTObjectProperties> const& props, std::vector<ElectronObject const*> const& objects) const;
template bool HLTTriggerPathProperties::testCutSet<PhotonObject>(HLTObjectProperties::TriggerObjectType const& type, std::vector<HLTObjectProperties> const& props, std::vector<PhotonObject const*> const& objects) const;
template bool HLTTriggerPathProperties::testCutSet<AK4JetObject>(HLTObjectProperties::TriggerObjectType const& type, std::vector<HLTObjectProperties> const& props, std::vector<AK4JetObject const*> const& objects) const;
template bool HLTTriggerPathProperties::testCutSet<AK8JetObject>(HLTObjectProperties::TriggerObjectType const& type, std::vector<HLTObjectProperties> const& props, std::vector<AK8JetObject const*> const& objects) const;


#endif
