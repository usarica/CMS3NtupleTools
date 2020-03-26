#ifndef CMS3_OFFSHELL_TRIGGERHELPERS_H
#define CMS3_OFFSHELL_TRIGGERHELPERS_H

#include "TriggerHelpersCore.h"


namespace TriggerHelpers{
  void configureHLTmap();

  std::vector<std::string> getHLTMenus(TriggerHelpers::TriggerType type);
  std::vector<std::string> getHLTMenus(std::vector<TriggerHelpers::TriggerType> const& types);

  std::vector< std::unordered_map< TriggerHelpers::TriggerType, std::vector<HLTTriggerPathProperties> >::const_iterator > getHLTMenuProperties(TriggerHelpers::TriggerType type);
  std::vector< std::unordered_map< TriggerHelpers::TriggerType, std::vector<HLTTriggerPathProperties> >::const_iterator > getHLTMenuProperties(std::vector<TriggerHelpers::TriggerType> const& types);
}

#endif
