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
    kSingleEle_L1EG, // We should use these for double-ele or egamma, not for single electron

    kMuEle,
    kMuEle_Extra,

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
    kPFMET_Control,
    kPFMET_MHT_Control,

    nTriggerTypes
  };

  bool hasRunRangeExclusions(std::string const& name, HLTTriggerPathProperties const** out_hltprop = nullptr); // This is to allow a string-based recognition of run range exclusions.
}

#endif
