#include <cassert>
#include "OffshellSampleHelpers.h"
#include "OffshellTriggerHelpers.h"
#include "MELAStreamHelpers.hh"


namespace TriggerHelpers{
  std::unordered_map< TriggerHelpers::TriggerType, std::vector<std::string> > HLT_type_list_map;
  std::unordered_map< TriggerHelpers::TriggerType, std::vector<HLTTriggerPathProperties> > HLT_type_proplist_map;
}


using namespace std;
using namespace MELAStreamHelpers;


std::vector<std::string> TriggerHelpers::getHLTMenus(TriggerHelpers::TriggerType type){ return getHLTMenus(std::vector<TriggerHelpers::TriggerType>{ type }); }
std::vector<std::string> TriggerHelpers::getHLTMenus(std::vector<TriggerHelpers::TriggerType> const& types){
  std::vector<std::string> res;
  for (auto const& type:types){
    auto it = HLT_type_list_map.find(type);
    if (it != HLT_type_list_map.end()) HelperFunctions::appendVector(res, it->second);
  }
  return res;
}
std::vector< std::pair<TriggerHelpers::TriggerType, HLTTriggerPathProperties const*> > TriggerHelpers::getHLTMenuProperties(TriggerHelpers::TriggerType type){ return getHLTMenuProperties(std::vector<TriggerHelpers::TriggerType>{ type }); }
std::vector< std::pair<TriggerHelpers::TriggerType, HLTTriggerPathProperties const*> > TriggerHelpers::getHLTMenuProperties(std::vector<TriggerHelpers::TriggerType> const& types){
  unsigned int isize=0;
  for (auto const& type:types){
    std::unordered_map< TriggerHelpers::TriggerType, std::vector<HLTTriggerPathProperties> >::const_iterator it = HLT_type_proplist_map.find(type);
    isize += it->second.size();
  }

  std::vector< std::pair<TriggerHelpers::TriggerType, HLTTriggerPathProperties const*> > res; res.reserve(isize);
  for (auto const& type:types){
    std::unordered_map< TriggerHelpers::TriggerType, std::vector<HLTTriggerPathProperties> >::const_iterator it = HLT_type_proplist_map.find(type);
    for (auto const& hltprop:it->second) res.emplace_back(type, &hltprop);
  }

  return res;
}

void TriggerHelpers::dropSelectionCuts(TriggerHelpers::TriggerType type){
  std::unordered_map< TriggerHelpers::TriggerType, std::vector<HLTTriggerPathProperties> >::iterator it = HLT_type_proplist_map.find(type);
  if (it != HLT_type_proplist_map.end()){ for (auto& props:it->second) props.resetCuts(); }
  else{
    MELAerr << "TriggerHelpers::dropSelectionCuts: Trigger type " << type << " is not defined." << endl;
    assert(0);
  }
}

void TriggerHelpers::configureHLTmap(){
  if (!SampleHelpers::runConfigure){
    MELAerr << "TriggerHelpers::configureHLTmap: Need to call SampleHelpers::configure(period, tag) first!" << endl;
    assert(0);
  }
  // Notice that the triggers that require cuts are ORDERED!
  // Ordering within each list is important.
  // FIXME: Triggers with prescales need to include higher-pT thresholds as well. SinglePhoton is done, but the rest needs revision either.
  switch (SampleHelpers::theDataYear){
  case 2016:
    HLT_type_proplist_map[kDoubleMu] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kMuon } } },
      { "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kMuon } } }
    };
    HLT_type_proplist_map[kDoubleMu_Prescaled] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*", { { HLTObjectProperties::kMuon, { { HLTObjectProperties::kPt, 17.f*1.1f } } }, { HLTObjectProperties::kMuon, { { HLTObjectProperties::kPt, 8.f*1.1f } } } } },
      { "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*", { { HLTObjectProperties::kMuon, { { HLTObjectProperties::kPt, 17.f*1.1f } } }, { HLTObjectProperties::kMuon, { { HLTObjectProperties::kPt, 8.f*1.1f } } } } }
    };
    HLT_type_proplist_map[kDoubleEle] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*", { { HLTObjectProperties::kElectron }, { HLTObjectProperties::kElectron } } },
      { "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*", { { HLTObjectProperties::kElectron }, { HLTObjectProperties::kElectron } } },
      { "HLT_DoubleEle33_CaloIdL_MW_v*", { { HLTObjectProperties::kElectron }, { HLTObjectProperties::kElectron } } },
      { "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*", { { HLTObjectProperties::kElectron }, { HLTObjectProperties::kElectron } } }
    };
    HLT_type_proplist_map[kDoubleEle_HighPt] = std::vector<HLTTriggerPathProperties>{
      { "HLT_DoublePhoton60_v*", { { HLTObjectProperties::kElectron }, { HLTObjectProperties::kElectron } } }
    };
    // These MuEle triggers are tricky to deal with because they are disabled for part of the run...
    HLT_type_proplist_map[kMuEle] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kElectron } } },
      { "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kElectron } } },
      { "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kElectron } } },
      { "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kElectron } } },
      { "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kElectron } } },
      { "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kElectron } } },
      { "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kElectron } } }
    };
    HLT_type_proplist_map[kSingleMu] = std::vector<HLTTriggerPathProperties>{
      { "HLT_IsoMu24_v*", { { HLTObjectProperties::kMuon } } },
      { "HLT_IsoTkMu24_v*", { { HLTObjectProperties::kMuon } } }
    };
    // Prescales are 0 or 1, so ignore pT thresholds
    HLT_type_proplist_map[kSingleMu_Prescaled] = std::vector<HLTTriggerPathProperties>{
      { "HLT_IsoMu22_v*", { { HLTObjectProperties::kMuon } } },
      { "HLT_IsoTkMu22_v*", { { HLTObjectProperties::kMuon } } },
      { "HLT_IsoMu20_v*", { { HLTObjectProperties::kMuon } } },
      { "HLT_IsoTkMu20_v*", { { HLTObjectProperties::kMuon } } },
    };
    HLT_type_proplist_map[kSingleMu_HighPt] = std::vector<HLTTriggerPathProperties>{
      { "HLT_TkMu50_v*", { { HLTObjectProperties::kMuon} } },
      { "HLT_Mu50_v*", { { HLTObjectProperties::kMuon } } }
    };
    HLT_type_proplist_map[kSingleEle] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Ele25_eta2p1_WPTight_Gsf_v*", { { HLTObjectProperties::kElectron } } },
      { "HLT_Ele27_WPTight_Gsf_v*", { { HLTObjectProperties::kElectron } } },
      { "HLT_Ele27_eta2p1_WPLoose_Gsf_v*", { { HLTObjectProperties::kElectron } } },
      { "HLT_Ele32_eta2p1_WPTight_Gsf_v*", { { HLTObjectProperties::kElectron } } }
    };
    HLT_type_proplist_map[kSingleEle_HighPt] = std::vector<HLTTriggerPathProperties>{ { "HLT_Photon175_v*", { { HLTObjectProperties::kElectron } } } };
    HLT_type_proplist_map[kSinglePho] = std::vector<HLTTriggerPathProperties>{
      //{ "HLT_Photon300_NoHE_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 300.f+40.f } } } } }, // Unprescaled
      //{ "HLT_Photon175_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 200.f }, { HLTObjectProperties::kPtHigh, 300.f+40.f } } } } }, // Unprescaled
      { "HLT_Photon175_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 175.f+45.f } } } } }, // Unprescaled
      { "HLT_Photon165_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 165.f+35.f }, { HLTObjectProperties::kPtHigh, 175.f+45.f } } } } },
      //{ "HLT_Photon165_HE10_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 165.f*1.1f }, { HLTObjectProperties::kPtHigh, 175.f*1.1f } } } } },
      //{ "HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_VBF_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 120.f*1.1f }, { HLTObjectProperties::kPtHigh, 165.f*1.1f } } } } },
      { "HLT_Photon120_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 120.f+15.f }, { HLTObjectProperties::kPtHigh, 165.f+35.f } } } } },
      //{ "HLT_Photon120_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 120.f*1.1f }, { HLTObjectProperties::kPtHigh, 165.f*1.1f } } } } },
      //{ "HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_VBF_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 90.f*1.1f }, { HLTObjectProperties::kPtHigh, 120.f*1.1f } } } } },
      { "HLT_Photon90_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 90.f*1.1f }, { HLTObjectProperties::kPtHigh, 120.f+15.f } } } } },
      //{ "HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_VBF_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 75.f*1.1f }, { HLTObjectProperties::kPtHigh, 90.f*1.1f } } } } },
      { "HLT_Photon75_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 75.f*1.1f }, { HLTObjectProperties::kPtHigh, 90.f*1.1f } } } } },
      //{ "HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_VBF_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 50.f*1.1f }, { HLTObjectProperties::kPtHigh, 75.f*1.1f } } } } },
      { "HLT_Photon50_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 50.f*1.1f }, { HLTObjectProperties::kPtHigh, 75.f*1.1f } } } } },
      //{ "HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_VBF_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 36.f*1.1f }, { HLTObjectProperties::kPtHigh, 50.f*1.1f } } } } },
      { "HLT_Photon36_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 36.f*1.1f }, { HLTObjectProperties::kPtHigh, 50.f*1.1f } } } } },
      { "HLT_Photon30_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 30.f*1.1f }, { HLTObjectProperties::kPtHigh, 36.f*1.1f } } } } },
      //{ "HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_VBF_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 22.f*1.1f }, { HLTObjectProperties::kPtHigh, 30.f*1.1f } } } } },
      { "HLT_Photon22_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 22.f*1.1f }, { HLTObjectProperties::kPtHigh, 30.f*1.1f } } } } }
    };
    HLT_type_proplist_map[kTripleLep] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kElectron }, { HLTObjectProperties::kElectron } } },
      { "HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kMuon }, { HLTObjectProperties::kElectron } } },
      { "HLT_TripleMu_12_10_5_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kMuon }, { HLTObjectProperties::kMuon } } },
      { "HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v*", { { HLTObjectProperties::kElectron }, { HLTObjectProperties::kElectron }, { HLTObjectProperties::kElectron } } }
    };
    HLT_type_proplist_map[kSingleMu_Control] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Mu17_TrkIsoVVL_v*", { { HLTObjectProperties::kMuon, { { HLTObjectProperties::kPt, 17.f } } } } }, // Has L1 and HLT prescales
      { "HLT_Mu17_v*", { { HLTObjectProperties::kMuon, { { HLTObjectProperties::kPt, 17.f } } } } }, // Has L1 and HLT prescales
      { "HLT_Mu8_TrkIsoVVL_v*", { { HLTObjectProperties::kMuon, { { HLTObjectProperties::kPt, 8.f } } } } }, // Has L1 and HLT prescales
      { "HLT_Mu8_v*", { { HLTObjectProperties::kMuon, { { HLTObjectProperties::kPt, 8.f } } } } } // Has L1 and HLT prescales
    };
    HLT_type_proplist_map[kSingleMu_Control_NoIso] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Mu17_v*", { { HLTObjectProperties::kMuon, { { HLTObjectProperties::kPt, 17.f } } } } }, // Has L1 and HLT prescales
      { "HLT_Mu8_v*", { { HLTObjectProperties::kMuon, { { HLTObjectProperties::kPt, 8.f } } } } } // Has L1 and HLT prescales
    };
    HLT_type_proplist_map[kSingleMu_Control_Iso] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Mu17_TrkIsoVVL_v*", { { HLTObjectProperties::kMuon, { { HLTObjectProperties::kPt, 17.f } } } } }, // Has L1 and HLT prescales
      { "HLT_Mu8_TrkIsoVVL_v*", { { HLTObjectProperties::kMuon, { { HLTObjectProperties::kPt, 8.f } } } } } // Has L1 and HLT prescales
    };
    HLT_type_proplist_map[kSingleEle_Control] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v*", { { HLTObjectProperties::kElectron, { { HLTObjectProperties::kPt, 17.f } } }, { HLTObjectProperties::kAK4Jet, { { HLTObjectProperties::kPt, 30.f } } } } }, // Has L1 and HLT prescales
      { "HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v*", { { HLTObjectProperties::kElectron, { { HLTObjectProperties::kPt, 17.f } } }, { HLTObjectProperties::kAK4Jet, { { HLTObjectProperties::kPt, 30.f } } } } }, // Has L1 and HLT prescales
      { "HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v*", { { HLTObjectProperties::kElectron, { { HLTObjectProperties::kPt, 8.f } } }, { HLTObjectProperties::kAK4Jet, { { HLTObjectProperties::kPt, 30.f } } } } }, // Has L1 prescales
      { "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v*", { { HLTObjectProperties::kElectron, { { HLTObjectProperties::kPt, 8.f } } }, { HLTObjectProperties::kAK4Jet, { { HLTObjectProperties::kPt, 30.f } } } } } // Has L1 prescales
    };
    HLT_type_proplist_map[kSingleEle_Control_NoIso] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v*", { { HLTObjectProperties::kElectron, { { HLTObjectProperties::kPt, 17.f } } }, { HLTObjectProperties::kAK4Jet, { { HLTObjectProperties::kPt, 30.f } } } } }, // Has L1 and HLT prescales
      { "HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v*", { { HLTObjectProperties::kElectron, { { HLTObjectProperties::kPt, 8.f } } }, { HLTObjectProperties::kAK4Jet, { { HLTObjectProperties::kPt, 30.f } } } } } // Has L1 prescales
    };
    HLT_type_proplist_map[kSingleEle_Control_Iso] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v*", { { HLTObjectProperties::kElectron, { { HLTObjectProperties::kPt, 17.f } } }, { HLTObjectProperties::kAK4Jet, { { HLTObjectProperties::kPt, 30.f } } } } }, // Has L1 and HLT prescales
      { "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v*", { { HLTObjectProperties::kElectron, { { HLTObjectProperties::kPt, 8.f } } }, { HLTObjectProperties::kAK4Jet, { { HLTObjectProperties::kPt, 30.f } } } } } // Has L1 prescales
    };
    HLT_type_proplist_map[kAK8PFJet_Control] = std::vector<HLTTriggerPathProperties>{
      { "HLT_AK8PFJet360_TrimMass30_v*", { { HLTObjectProperties::kAK8Jet, { { HLTObjectProperties::kPt, 360.f }, { HLTObjectProperties::kMass, 30.f } } } } }
    };
    HLT_type_proplist_map[kPFHT_Control] = std::vector<HLTTriggerPathProperties>{
      { "HLT_PFHT900_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 900.f } } } } }, // Prescaled
      { "HLT_PFHT800_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 800.f } } } } }, // Prescaled, might not exist?
      { "HLT_PFHT650_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 650.f } } } } }, // Prescaled
      { "HLT_PFHT600_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 600.f } } } } }, // Prescaled
      { "HLT_PFHT475_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 475.f } } } } }, // Prescaled
      { "HLT_PFHT400_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 400.f } } } } }, // Prescaled
      { "HLT_PFHT350_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 350.f } } } } }, // Prescaled
      { "HLT_PFHT300_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 300.f } } } } }, // Prescaled
      { "HLT_PFHT250_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 250.f } } } } }, // Prescaled
      { "HLT_PFHT200_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 200.f } } } } }, // Prescaled
      { "HLT_PFHT125_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 125.f } } } } } // Prescaled
    };
    HLT_type_proplist_map[kPFMET_MHT_Control] = std::vector<HLTTriggerPathProperties>{
      { "HLT_PFMET170_HBHECleaned_v*", { { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 170.f } } } } },
      { "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v*", { { HLTObjectProperties::kMET_NoMu, { { HLTObjectProperties::kPt, 120.f } } }, { HLTObjectProperties::kHT_NoMu, { { HLTObjectProperties::kMass, 120.f } } } } },
      { "HLT_PFMET120_PFMHT120_IDTight_v*", { { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 120.f } } }, { HLTObjectProperties::kHT, { { HLTObjectProperties::kMass, 120.f } } } } }
    };
    break;
  case 2017:
    HLT_type_proplist_map[kDoubleMu] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kMuon } } },
      { "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kMuon } } }
    };
    HLT_type_proplist_map[kDoubleMu_Prescaled] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*", { { HLTObjectProperties::kMuon, { { HLTObjectProperties::kPt, 17.f*1.1f } } }, { HLTObjectProperties::kMuon, { { HLTObjectProperties::kPt, 8.f*1.1f } } } } }
    };
    HLT_type_proplist_map[kDoubleEle] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*", { { HLTObjectProperties::kElectron }, { HLTObjectProperties::kElectron } } },
      { "HLT_DoubleEle33_CaloIdL_MW_v*", { { HLTObjectProperties::kElectron }, { HLTObjectProperties::kElectron } } }
    };
    HLT_type_proplist_map[kDoubleEle_HighPt] = std::vector<HLTTriggerPathProperties>{
      { "HLT_DoublePhoton70_v*", { { HLTObjectProperties::kElectron }, { HLTObjectProperties::kElectron } } }
    };
    HLT_type_proplist_map[kMuEle] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kElectron } } },
      { "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kElectron } } },
      { "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kElectron } } },
      { "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kElectron } } }
    };
    HLT_type_proplist_map[kSingleMu] = std::vector<HLTTriggerPathProperties>{ { "HLT_IsoMu27_v*", { { HLTObjectProperties::kMuon } } } };
    // HLT_IsoMu20_v* have L1 and HLT prescales
    HLT_type_proplist_map[kSingleMu_Prescaled] = std::vector<HLTTriggerPathProperties>{
      { "HLT_IsoMu24_v*", { { HLTObjectProperties::kMuon, { { HLTObjectProperties::kPt, 24.f*1.1f } } } } },
      { "HLT_IsoMu20_v*", { { HLTObjectProperties::kMuon, { { HLTObjectProperties::kPt, 20.f*1.1f }, { HLTObjectProperties::kPtHigh, 24.f*1.1f } } } } }
    };
    HLT_type_proplist_map[kSingleMu_HighPt] = std::vector<HLTTriggerPathProperties>{
      { "HLT_TkMu100_v*", { { HLTObjectProperties::kMuon } } },
      { "HLT_OldMu100_v*", { { HLTObjectProperties::kMuon } } },
      { "HLT_Mu50_v*", { { HLTObjectProperties::kMuon } } }
    };
    HLT_type_proplist_map[kSingleEle] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Ele35_WPTight_Gsf_v*", { { HLTObjectProperties::kElectron } } },
      { "HLT_Ele38_WPTight_Gsf_v*", { { HLTObjectProperties::kElectron } } },
      { "HLT_Ele40_WPTight_Gsf_v*", { { HLTObjectProperties::kElectron } } }
    };
    HLT_type_proplist_map[kSingleEle_HighPt] = std::vector<HLTTriggerPathProperties>{ { "HLT_Photon200_v*", { { HLTObjectProperties::kElectron } } } };
    HLT_type_proplist_map[kSinglePho] = std::vector<HLTTriggerPathProperties>{
      //{ "HLT_Photon300_NoHE_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 300.f+40.f } } } } }, // Unprescaled
      //{ "HLT_Photon200_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 200.f*1.15f }, { HLTObjectProperties::kPtHigh, 300.f+40.f } } } } }, // Unprescaled
      //{ "HLT_Photon175_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 175.f*1.1f }, { HLTObjectProperties::kPtHigh, 200.f*1.15f } } } } },
      //{ "HLT_Photon165_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 165.f*1.1f }, { HLTObjectProperties::kPtHigh, 175.f*1.15f } } } } },
      { "HLT_Photon200_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 200.f*1.15f } } } } }, // Unprescaled
      { "HLT_Photon165_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 165.f*1.1f }, { HLTObjectProperties::kPtHigh, 200.f*1.15f } } } } },
      { "HLT_Photon120_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 120.f*1.1f }, { HLTObjectProperties::kPtHigh, 165.f*1.1f } } } } },
      { "HLT_Photon90_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 90.f*1.1f }, { HLTObjectProperties::kPtHigh, 120.f*1.1f } } } } },
      { "HLT_Photon75_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 75.f*1.1f }, { HLTObjectProperties::kPtHigh, 90.f*1.1f } } } } },
      { "HLT_Photon50_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 50.f*1.1f }, { HLTObjectProperties::kPtHigh, 75.f*1.1f } } } } },
      { "HLT_Photon33_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 33.f*1.1f }, { HLTObjectProperties::kPtHigh, 50.f*1.1f } } } } }//,
      //{ "HLT_Photon30_HoverELoose_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 30.f*1.1f }, { HLTObjectProperties::kPtHigh, 33.f*1.1f } } } } }, // Only runs 299368 - 306460 (36.75 fb-1)
      //{ "HLT_Photon25_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 25.f*1.1f }, { HLTObjectProperties::kPtHigh, 30.f*1.1f } } } } }, // Only runs 302026 - 306460 (27.13 fb-1)
      //{ "HLT_Photon20_HoverELoose_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 20.f*1.1f }, { HLTObjectProperties::kPtHigh, 25.f*1.1f } } } } } // Only runs 299368 - 306460 (36.75 fb-1)
    };
    HLT_type_proplist_map[kTripleLep] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kElectron }, { HLTObjectProperties::kElectron } } },
      { "HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kMuon }, { HLTObjectProperties::kElectron } } },
      { "HLT_TripleMu_10_5_5_DZ_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kMuon }, { HLTObjectProperties::kMuon } } },
      { "HLT_TripleMu_12_10_5_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kMuon }, { HLTObjectProperties::kMuon } } },
      { "HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v*", { { HLTObjectProperties::kElectron }, { HLTObjectProperties::kElectron }, { HLTObjectProperties::kElectron } } }
    };
    HLT_type_proplist_map[kSingleMu_Control] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Mu17_TrkIsoVVL_v*", { { HLTObjectProperties::kMuon, { { HLTObjectProperties::kPt, 17.f } } } } }, // Has L1 and HLT prescales
      { "HLT_Mu17_v*", { { HLTObjectProperties::kMuon, { { HLTObjectProperties::kPt, 17.f } } } } }, // Has L1 and HLT prescales
      { "HLT_Mu8_TrkIsoVVL_v*", { { HLTObjectProperties::kMuon, { { HLTObjectProperties::kPt, 8.f } } } } }, // Has L1 and HLT prescales
      { "HLT_Mu8_v*", { { HLTObjectProperties::kMuon, { { HLTObjectProperties::kPt, 8.f } } } } } // Has L1 and HLT prescales
    };
    HLT_type_proplist_map[kSingleMu_Control_NoIso] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Mu17_v*", { { HLTObjectProperties::kMuon, { { HLTObjectProperties::kPt, 17.f } } } } }, // Has L1 and HLT prescales
      { "HLT_Mu8_v*", { { HLTObjectProperties::kMuon, { { HLTObjectProperties::kPt, 8.f } } } } } // Has L1 and HLT prescales
    };
    HLT_type_proplist_map[kSingleMu_Control_Iso] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Mu17_TrkIsoVVL_v*", { { HLTObjectProperties::kMuon, { { HLTObjectProperties::kPt, 17.f } } } } }, // Has L1 and HLT prescales
      { "HLT_Mu8_TrkIsoVVL_v*", { { HLTObjectProperties::kMuon, { { HLTObjectProperties::kPt, 8.f } } } } } // Has L1 and HLT prescales
    };
    HLT_type_proplist_map[kSingleEle_Control] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v*", { { HLTObjectProperties::kElectron, { { HLTObjectProperties::kPt, 17.f } } }, { HLTObjectProperties::kAK4Jet, { { HLTObjectProperties::kPt, 30.f } } } } }, // Has L1 and HLT prescales
      { "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v*", { { HLTObjectProperties::kElectron, { { HLTObjectProperties::kPt, 12.f } } }, { HLTObjectProperties::kAK4Jet, { { HLTObjectProperties::kPt, 30.f } } } } }, // Has L1 and HLT prescales
      { "HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v*", { { HLTObjectProperties::kElectron, { { HLTObjectProperties::kPt, 8.f } } }, { HLTObjectProperties::kAK4Jet, { { HLTObjectProperties::kPt, 30.f } } } } }, // Has L1 prescales
      { "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v*", { { HLTObjectProperties::kElectron, { { HLTObjectProperties::kPt, 8.f } } }, { HLTObjectProperties::kAK4Jet, { { HLTObjectProperties::kPt, 30.f } } } } } // Has L1 prescales
    };
    HLT_type_proplist_map[kSingleEle_Control_NoIso] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v*", { { HLTObjectProperties::kElectron, { { HLTObjectProperties::kPt, 17.f } } }, { HLTObjectProperties::kAK4Jet, { { HLTObjectProperties::kPt, 30.f } } } } }, // Has L1 and HLT prescales
      { "HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v*", { { HLTObjectProperties::kElectron, { { HLTObjectProperties::kPt, 8.f } } }, { HLTObjectProperties::kAK4Jet, { { HLTObjectProperties::kPt, 30.f } } } } } // Has L1 prescales
    };
    HLT_type_proplist_map[kSingleEle_Control_Iso] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v*", { { HLTObjectProperties::kElectron, { { HLTObjectProperties::kPt, 12.f } } }, { HLTObjectProperties::kAK4Jet, { { HLTObjectProperties::kPt, 30.f } } } } }, // Has L1 and HLT prescales
      { "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v*", { { HLTObjectProperties::kElectron, { { HLTObjectProperties::kPt, 8.f } } }, { HLTObjectProperties::kAK4Jet, { { HLTObjectProperties::kPt, 30.f } } } } } // Has L1 prescales
    };
    HLT_type_proplist_map[kAK8PFJet_Control] = std::vector<HLTTriggerPathProperties>{
      { "HLT_AK8PFJet360_TrimMass30_v*", { { HLTObjectProperties::kAK8Jet, { { HLTObjectProperties::kPt, 360.f }, { HLTObjectProperties::kMass, 30.f } } } } }
    };
    HLT_type_proplist_map[kPFHT_Control] = std::vector<HLTTriggerPathProperties>{
      { "HLT_PFHT1050_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 1050.f } } } } }, // Prescaled
      { "HLT_PFHT890_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 890.f } } } } }, // Prescaled
      { "HLT_PFHT780_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 780.f } } } } }, // Prescaled
      { "HLT_PFHT680_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 680.f } } } } }, // Prescaled
      { "HLT_PFHT590_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 590.f } } } } }, // Prescaled
      { "HLT_PFHT510_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 510.f } } } } }, // Prescaled
      { "HLT_PFHT430_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 430.f } } } } }, // Prescaled
      { "HLT_PFHT370_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 370.f } } } } }, // Prescaled
      { "HLT_PFHT350_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 350.f } } } } }, // Prescaled
      { "HLT_PFHT250_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 250.f } } } } }, // Prescaled
      { "HLT_PFHT180_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 180.f } } } } } // Prescaled
    };
    HLT_type_proplist_map[kPFMET_MHT_Control] = std::vector<HLTTriggerPathProperties>{
      { "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v*", { { HLTObjectProperties::kMET_NoMu, { { HLTObjectProperties::kPt, 120.f } } }, { HLTObjectProperties::kHT_NoMu, { { HLTObjectProperties::kMass, 120.f } } }, { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 60.f } } } } },
      { "HLT_PFMET120_PFMHT120_IDTight_PFHT60_v*", { { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 120.f } } }, { HLTObjectProperties::kHT, { { HLTObjectProperties::kMass, 120.f }, { HLTObjectProperties::kPt, 60.f } } } } },
      { "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v*", { { HLTObjectProperties::kMET_NoMu, { { HLTObjectProperties::kPt, 120.f } } }, { HLTObjectProperties::kHT_NoMu, { { HLTObjectProperties::kMass, 120.f } } } } },
      { "HLT_PFMET120_PFMHT120_IDTight_v*", { { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 120.f } } }, { HLTObjectProperties::kHT, { { HLTObjectProperties::kMass, 120.f } } } } },
    };
    break;
  case 2018:
    HLT_type_proplist_map[kDoubleMu] = std::vector<HLTTriggerPathProperties>{
      // No gain from adding the Mass8 version of the one below
      { "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kMuon } } }
    };
    HLT_type_proplist_map[kDoubleMu_Prescaled] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*", { { HLTObjectProperties::kMuon, { { HLTObjectProperties::kPt, 17.f*1.1f } } }, { HLTObjectProperties::kMuon, { { HLTObjectProperties::kPt, 8.f*1.1f } } } } }
    };
    HLT_type_proplist_map[kDoubleEle] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*", { { HLTObjectProperties::kElectron }, { HLTObjectProperties::kElectron } } },
      //{ "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*", { { HLTObjectProperties::kElectron }, { HLTObjectProperties::kElectron } } },
      { "HLT_DoubleEle25_CaloIdL_MW_v*", { { HLTObjectProperties::kElectron }, { HLTObjectProperties::kElectron } } }
    };
    HLT_type_proplist_map[kDoubleEle_HighPt] = std::vector<HLTTriggerPathProperties>{
      { "HLT_DoublePhoton70_v*", { { HLTObjectProperties::kElectron }, { HLTObjectProperties::kElectron } } }
    };
    HLT_type_proplist_map[kMuEle] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kElectron } } },
      { "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kElectron } } },
      { "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kElectron } } },
      { "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kElectron } } }
    };
    HLT_type_proplist_map[kSingleMu] = std::vector<HLTTriggerPathProperties>{ { "HLT_IsoMu24_v*", { { HLTObjectProperties::kMuon } } } };
    HLT_type_proplist_map[kSingleMu_Prescaled] = std::vector<HLTTriggerPathProperties>{ { "HLT_IsoMu20_v*", { { HLTObjectProperties::kMuon } } } };
    HLT_type_proplist_map[kSingleMu_HighPt] = std::vector<HLTTriggerPathProperties>{
      { "HLT_TkMu100_v*", { { HLTObjectProperties::kMuon } } },
      { "HLT_OldMu100_v*", { { HLTObjectProperties::kMuon } } },
      { "HLT_Mu50_v*", { { HLTObjectProperties::kMuon } } }
    };
    HLT_type_proplist_map[kSingleEle] = std::vector<HLTTriggerPathProperties>{ { "HLT_Ele32_WPTight_Gsf_v*", { { HLTObjectProperties::kElectron } } } };
    HLT_type_proplist_map[kSingleEle_HighPt] = std::vector<HLTTriggerPathProperties>{ { "HLT_Photon200_v*", { { HLTObjectProperties::kElectron } } } };
    HLT_type_proplist_map[kSinglePho] = std::vector<HLTTriggerPathProperties>{
      //{ "HLT_Photon300_NoHE_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 300.f+40.f } } } } }, // Unprescaled
      //{ "HLT_Photon200_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 200.f*1.15f }, { HLTObjectProperties::kPtHigh, 300.f+40.f } } } } }, // Unprescaled
      //{ "HLT_Photon175_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 175.f*1.1f }, { HLTObjectProperties::kPtHigh, 200.f*1.15f } } } } },
      //{ "HLT_Photon165_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 165.f*1.1f }, { HLTObjectProperties::kPtHigh, 175.f*1.15f } } } } },
      { "HLT_Photon200_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 200.f*1.15f } } } } }, // Unprescaled
      { "HLT_Photon165_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 165.f*1.1f }, { HLTObjectProperties::kPtHigh, 200.f*1.15f } } } } },
      { "HLT_Photon120_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 120.f*1.1f }, { HLTObjectProperties::kPtHigh, 165.f*1.1f } } } } },
      { "HLT_Photon90_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 90.f*1.1f }, { HLTObjectProperties::kPtHigh, 120.f*1.1f } } } } },
      { "HLT_Photon75_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 90.f }, { HLTObjectProperties::kPtHigh, 90.f*1.1f } } } } },
      { "HLT_Photon50_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 50.f*1.1f }, { HLTObjectProperties::kPtHigh, 90.f } } } } },
      { "HLT_Photon33_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 33.f*1.1f }, { HLTObjectProperties::kPtHigh, 50.f*1.1f } } } } },
      { "HLT_Photon30_HoverELoose_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 30.f*1.1f }, { HLTObjectProperties::kPtHigh, 33.f*1.1f } } } } },
      { "HLT_Photon20_HoverELoose_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 20.f*1.1f }, { HLTObjectProperties::kPtHigh, 20.f*1.1f } } } } }
    };
    HLT_type_proplist_map[kTripleLep] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kElectron }, { HLTObjectProperties::kElectron } } },
      { "HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kMuon }, { HLTObjectProperties::kElectron } } },
      { "HLT_TripleMu_10_5_5_DZ_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kMuon }, { HLTObjectProperties::kMuon } } },
      { "HLT_TripleMu_12_10_5_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kMuon }, { HLTObjectProperties::kMuon } } }//,
      //{ "HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v*", { { HLTObjectProperties::kElectron }, { HLTObjectProperties::kElectron }, { HLTObjectProperties::kElectron } } } // Somehow its effective lumi is not the same as active lumi, and it does not show up in the spreadsheets
    };
    HLT_type_proplist_map[kSingleMu_Control] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Mu17_TrkIsoVVL_v*", { { HLTObjectProperties::kMuon, { { HLTObjectProperties::kPt, 17.f } } } } }, // Has L1 and HLT prescales
      { "HLT_Mu17_v*", { { HLTObjectProperties::kMuon, { { HLTObjectProperties::kPt, 17.f } } } } }, // Has L1 and HLT prescales
      { "HLT_Mu8_TrkIsoVVL_v*", { { HLTObjectProperties::kMuon, { { HLTObjectProperties::kPt, 8.f } } } } }, // Has L1 and HLT prescales
      { "HLT_Mu8_v*", { { HLTObjectProperties::kMuon, { { HLTObjectProperties::kPt, 8.f } } } } } // Has L1 and HLT prescales
    };
    HLT_type_proplist_map[kSingleMu_Control_NoIso] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Mu17_v*", { { HLTObjectProperties::kMuon, { { HLTObjectProperties::kPt, 17.f } } } } }, // Has L1 and HLT prescales
      { "HLT_Mu8_v*", { { HLTObjectProperties::kMuon, { { HLTObjectProperties::kPt, 8.f } } } } } // Has L1 and HLT prescales
    };
    HLT_type_proplist_map[kSingleMu_Control_Iso] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Mu17_TrkIsoVVL_v*", { { HLTObjectProperties::kMuon, { { HLTObjectProperties::kPt, 17.f } } } } }, // Has L1 and HLT prescales
      { "HLT_Mu8_TrkIsoVVL_v*", { { HLTObjectProperties::kMuon, { { HLTObjectProperties::kPt, 8.f } } } } } // Has L1 and HLT prescales
    };
    HLT_type_proplist_map[kSingleEle_Control] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v*", { { HLTObjectProperties::kElectron, { { HLTObjectProperties::kPt, 17.f } } }, { HLTObjectProperties::kAK4Jet, { { HLTObjectProperties::kPt, 30.f } } } } }, // Has L1 and HLT prescales
      { "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v*", { { HLTObjectProperties::kElectron, { { HLTObjectProperties::kPt, 12.f } } }, { HLTObjectProperties::kAK4Jet, { { HLTObjectProperties::kPt, 30.f } } } } }, // Has L1 and HLT prescales
      { "HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v*", { { HLTObjectProperties::kElectron, { { HLTObjectProperties::kPt, 8.f } } }, { HLTObjectProperties::kAK4Jet, { { HLTObjectProperties::kPt, 30.f } } } } }, // Has L1 prescales
      { "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v*", { { HLTObjectProperties::kElectron, { { HLTObjectProperties::kPt, 8.f } } }, { HLTObjectProperties::kAK4Jet, { { HLTObjectProperties::kPt, 30.f } } } } } // Has L1 prescales
    };
    HLT_type_proplist_map[kSingleEle_Control_NoIso] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v*", { { HLTObjectProperties::kElectron, { { HLTObjectProperties::kPt, 17.f } } }, { HLTObjectProperties::kAK4Jet, { { HLTObjectProperties::kPt, 30.f } } } } }, // Has L1 and HLT prescales
      { "HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v*", { { HLTObjectProperties::kElectron, { { HLTObjectProperties::kPt, 8.f } } }, { HLTObjectProperties::kAK4Jet, { { HLTObjectProperties::kPt, 30.f } } } } } // Has L1 prescales
    };
    HLT_type_proplist_map[kSingleEle_Control_Iso] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v*", { { HLTObjectProperties::kElectron, { { HLTObjectProperties::kPt, 12.f } } }, { HLTObjectProperties::kAK4Jet, { { HLTObjectProperties::kPt, 30.f } } } } }, // Has L1 and HLT prescales
      { "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v*", { { HLTObjectProperties::kElectron, { { HLTObjectProperties::kPt, 8.f } } }, { HLTObjectProperties::kAK4Jet, { { HLTObjectProperties::kPt, 30.f } } } } } // Has L1 prescales
    };
    HLT_type_proplist_map[kAK8PFJet_Control] = std::vector<HLTTriggerPathProperties>{
      { "HLT_AK8PFJet400_TrimMass30_v*", { { HLTObjectProperties::kAK8Jet, { { HLTObjectProperties::kPt, 400.f }, { HLTObjectProperties::kMass, 30.f } } } } }
    };
    // FIXME: 200313 production contains only HLT_PFHT1050_v*
    HLT_type_proplist_map[kPFHT_Control] = std::vector<HLTTriggerPathProperties>{
      { "HLT_PFHT1050_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 1050.f } } } } },
      { "HLT_PFHT890_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 890.f } } } } }, // Prescaled
      { "HLT_PFHT780_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 780.f } } } } }, // Prescaled
      { "HLT_PFHT680_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 680.f } } } } }, // Prescaled
      { "HLT_PFHT590_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 590.f } } } } }, // Prescaled
      { "HLT_PFHT510_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 510.f } } } } }, // Prescaled
      { "HLT_PFHT430_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 430.f } } } } }, // Prescaled
      { "HLT_PFHT370_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 370.f } } } } }, // Prescaled
      { "HLT_PFHT350_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 350.f } } } } }, // Prescaled
      { "HLT_PFHT250_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 250.f } } } } }, // Prescaled
      { "HLT_PFHT180_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 180.f } } } } } // Prescaled
    };
    HLT_type_proplist_map[kPFMET_MHT_Control] = std::vector<HLTTriggerPathProperties>{
      { "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v*", { { HLTObjectProperties::kMET_NoMu, { { HLTObjectProperties::kPt, 120.f } } }, { HLTObjectProperties::kHT_NoMu, { { HLTObjectProperties::kMass, 120.f } } }, { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 60.f } } } } },
      { "HLT_PFMET120_PFMHT120_IDTight_PFHT60_v*", { { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 120.f } } }, { HLTObjectProperties::kHT, { { HLTObjectProperties::kMass, 120.f }, { HLTObjectProperties::kPt, 60.f } } } } },
      { "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v*", { { HLTObjectProperties::kMET_NoMu, { { HLTObjectProperties::kPt, 120.f } } }, { HLTObjectProperties::kHT_NoMu, { { HLTObjectProperties::kMass, 120.f } } } } },
      { "HLT_PFMET120_PFMHT120_IDTight_v*", { { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 120.f } } }, { HLTObjectProperties::kHT, { { HLTObjectProperties::kMass, 120.f } } } } }
    };
    break;
  default:
    break;
  }

  // Check that all triggers are defined
  for (int itt=(int) kTripleLep; itt!=(int) nTriggerTypes; itt++){
    if (HLT_type_proplist_map.find((TriggerType) itt)==HLT_type_proplist_map.cend()){
      MELAerr << "TriggerHelpers::configureHLTmap: Triggers for type " << itt << " are not defined for year " << SampleHelpers::theDataYear << ". Please fix the vectors." << endl;
      assert(0);
    }
  }

  // Fill the name map as well for simpler checking functionality
  for (auto const& it:HLT_type_proplist_map){
    auto const& props = it.second;
    std::vector<std::string> tmplist; tmplist.reserve(props.size());
    for (auto const& prop:props) tmplist.push_back(prop.getName());
    HLT_type_list_map[it.first] = tmplist;
  }
}
