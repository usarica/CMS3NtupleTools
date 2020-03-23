#ifndef CMS3_HLTTRIGGERPATHPROPERTIES_H
#define CMS3_HLTTRIGGERPATHPROPERTIES_H

#include <string>
#include <vector>
#include "HLTObjectProperties.h"


class HLTTriggerPathProperties{
protected:
  std::string name;
  std::vector<HLTObjectProperties> triggerObjectProperties;

  void setupName();

public:
  HLTTriggerPathProperties();
  HLTTriggerPathProperties(std::string const& name_);
  HLTTriggerPathProperties(std::string const& name_, std::vector<HLTObjectProperties> const& triggerObjectProperties_);
  HLTTriggerPathProperties(HLTTriggerPathProperties const& other);

  void setup();
  void addObjectProperties(HLTObjectProperties const& props);

  bool testCuts(ParticleObject::LorentzVector_t const& p4, HLTObjectProperties::TriggerObjectType const& type_) const;

  bool isSameTrigger(std::string const& name_) const;

  std::string const& getName() const{ return name; }
  std::vector<HLTObjectProperties> const& getObjectProperties() const{ return triggerObjectProperties; }

};


#endif
