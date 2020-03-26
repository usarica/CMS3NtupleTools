#ifndef CMS3_TRIGGERHELPERS_H
#define CMS3_TRIGGERHELPERS_H

#include <unordered_map>
#include "HLTTriggerPathProperties.h"


namespace TriggerHelpers{
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
}

#endif
