#ifndef CMS3_OFFSHELL_TRIGGERHELPERS_H
#define CMS3_OFFSHELL_TRIGGERHELPERS_H

#include "OffshellSampleHelpers.h"


namespace OffshellTriggerHelpers{
  enum TriggerType{
    kDoubleMu=0,
    kDoubleEle,
    kMuEle,
    kSingleMu,
    kSingleEle,
    kSinglePho,
    kTripleLep,
    nTriggerTypes
  };

  void configureHLTmap();
  std::vector<std::string> getHLTMenus(OffshellTriggerHelpers::TriggerType type);
  std::vector<std::string> getHLTMenus(std::vector<OffshellTriggerHelpers::TriggerType> const& types);

}

#endif
