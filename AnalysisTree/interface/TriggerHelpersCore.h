#ifndef CMS3_TRIGGERHELPERS_H
#define CMS3_TRIGGERHELPERS_H

#include <unordered_map>
#include "HLTTriggerPathProperties.h"


namespace TriggerHelpers{
  enum TriggerType{
    kTripleLep=0,

    kDoubleMu,
    kDoubleMu_Prescaled,

    kDoubleEle,
    kDoubleEle_HighPt,

    kMuEle,

    kSingleMu,
    kSingleMu_Prescaled,
    kSingleMu_HighPt,
    kSingleMu_Control,
    // Subdivisions of kSingleMu_Control
    kSingleMu_Control_NoIso,
    kSingleMu_Control_Iso,

    kSingleEle,
    kSingleEle_Prescaled,
    kSingleEle_HighPt,
    kSingleEle_Control,
    // Subdivisions of kSingleEle_Control
    kSingleEle_Control_NoIso,
    kSingleEle_Control_Iso,

    kSinglePho,

    kAK8PFJet_Control,
    kPFHT_Control,
    kPFMET_MHT_Control,

    nTriggerTypes
  };
}

#endif
