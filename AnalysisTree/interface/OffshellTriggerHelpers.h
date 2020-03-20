#ifndef CMS3_OFFSHELL_TRIGGERHELPERS_H
#define CMS3_OFFSHELL_TRIGGERHELPERS_H

#include "OffshellSampleHelpers.h"


namespace OffshellTriggerHelpers{
  enum TriggerType{
    kTripleLep=0,

    kDoubleMu,
    kDoubleEle,
    kMuEle,

    kSingleMu,
    kSingleEle,
    kSinglePho,

    kSingleMu_Control,
    kSingleMu_Control_HighPt,
    kSingleEle_Control,

    kAK8PFJet_Control,
    kPFHT_Control,
    kPFMET_MHT_Control,

    nTriggerTypes
  };

  void configureHLTmap();
  std::vector<std::string> getHLTMenus(OffshellTriggerHelpers::TriggerType type);
  std::vector<std::string> getHLTMenus(std::vector<OffshellTriggerHelpers::TriggerType> const& types);

}

#endif
