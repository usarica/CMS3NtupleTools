#ifndef HLTTRIGGERPATHOBJECT_H
#define HLTTRIGGERPATHOBJECT_H

#include <string>


#define HLTTRIGGERPATH_VARIABLES \
HLTTRIGGERPATH_VARIABLE(std::string, name, "") \
HLTTRIGGERPATH_VARIABLE(bool, passTrigger, false) \
HLTTRIGGERPATH_VARIABLE(int, L1prescale, 1) \
HLTTRIGGERPATH_VARIABLE(int, HLTprescale, 1)


struct HLTTriggerPathObject{
#define HLTTRIGGERPATH_VARIABLE(TYPE, NAME, DEFVAL) TYPE NAME;
  HLTTRIGGERPATH_VARIABLES;
#undef HLTTRIGGERPATH_VARIABLE

  HLTTriggerPathObject();
  HLTTriggerPathObject(HLTTriggerPathObject const& other);
  HLTTriggerPathObject& operator=(const HLTTriggerPathObject& other);

  void swap(HLTTriggerPathObject& other);

};


#endif
