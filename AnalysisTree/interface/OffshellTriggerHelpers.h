#ifndef CMS3_OFFSHELL_TRIGGERHELPERS_H
#define CMS3_OFFSHELL_TRIGGERHELPERS_H

#include <utility>
#include "TriggerHelpersCore.h"


namespace TriggerHelpers{
  void configureHLTmap();

  std::vector<std::string> getHLTMenus(TriggerHelpers::TriggerType type);
  std::vector<std::string> getHLTMenus(std::vector<TriggerHelpers::TriggerType> const& types);

  std::vector< std::pair<TriggerHelpers::TriggerType, HLTTriggerPathProperties const*> > getHLTMenuProperties(TriggerHelpers::TriggerType type);
  std::vector< std::pair<TriggerHelpers::TriggerType, HLTTriggerPathProperties const*> > getHLTMenuProperties(std::vector<TriggerHelpers::TriggerType> const& types);

  void dropSelectionCuts(TriggerHelpers::TriggerType type); // Drops only the pt, eta etc. cuts, not the type matching requirements
}

#endif
