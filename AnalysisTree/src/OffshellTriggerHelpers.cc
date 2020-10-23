#include <cassert>
#include "OffshellSampleHelpers.h"
#include "OffshellTriggerHelpers.h"
#include "MELAStreamHelpers.hh"


namespace TriggerHelpers{
  std::unordered_map< TriggerHelpers::TriggerType, std::vector<std::string> > HLT_type_list_map;
  std::unordered_map< TriggerHelpers::TriggerType, std::vector<HLTTriggerPathProperties> > HLT_type_proplist_map;
  std::vector<HLTTriggerPathProperties const*> runRangeExcluded_HLTprop_list; // This list is only for ease

  void assignRunRangeExclusions(std::string const& name, std::vector< std::pair<int, int> > const& rangelist);
  void assignTriggerObjectCheckException(std::string const& name, HLTTriggerPathProperties::TriggerObjectExceptionType const& flag);
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
    if (it == HLT_type_proplist_map.cend()){
      MELAerr << "TriggerHelpers::getHLTMenuProperties: Trigger type " << type << " is not defined in the HLT type-properties map." << endl;
      assert(0);
    }
    isize += it->second.size();
  }

  std::vector< std::pair<TriggerHelpers::TriggerType, HLTTriggerPathProperties const*> > res; res.reserve(isize);
  for (auto const& type:types){
    std::unordered_map< TriggerHelpers::TriggerType, std::vector<HLTTriggerPathProperties> >::const_iterator it = HLT_type_proplist_map.find(type);
    if (it == HLT_type_proplist_map.cend()){
      MELAerr << "TriggerHelpers::getHLTMenuProperties: Trigger type " << type << " is not defined in the HLT type-properties map." << endl;
      assert(0);
    }
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

void TriggerHelpers::assignRunRangeExclusions(std::string const& name, std::vector< std::pair<int, int> > const& rangelist){
  if (rangelist.empty()) return;

  std::pair<unsigned int, unsigned int> runrange_dp = SampleHelpers::getRunRangeFromDataPeriod(SampleHelpers::getDataPeriod());
  std::pair<unsigned int, unsigned int> runrange_year = SampleHelpers::getRunRangeFromDataPeriod(Form("%i", SampleHelpers::getDataYear()));
  std::vector< std::pair<unsigned int, unsigned int> > rangelist_eff; rangelist_eff.reserve(rangelist.size());

  // Place run range exclusions only when needed!
  // Do not place them if the excluded run range falls outside the data period run range!
  // Also modify run ranges in case they contain -1.
  for (auto pp:rangelist){
    if (pp.first == -1) pp.first = runrange_year.first;
    if (pp.second == -1) pp.second = runrange_year.second;
    if ((unsigned int) pp.first>runrange_dp.second || (unsigned int) pp.second<runrange_dp.first) continue;
    rangelist_eff.emplace_back(pp.first, pp.second);
  }

  if (rangelist_eff.empty()) return;

  for (int itt=0; itt!=(int) nTriggerTypes; itt++){
    // No need to check the iterator since it was already checked.
    auto it = HLT_type_proplist_map.find((TriggerType) itt);
    for (auto& hltprop:it->second){
      if (!hltprop.isSameTrigger(name)) continue;

      MELAout << "TriggerHelpers::assignRunRangeExclusions: Adding run range exclusion to " << name << ". Excluded ranges = " << rangelist_eff << endl;

      hltprop.setExcludedRunRanges(rangelist_eff);
      runRangeExcluded_HLTprop_list.push_back(&hltprop);
    }
  }
}
bool TriggerHelpers::hasRunRangeExclusions(std::string const& name, HLTTriggerPathProperties const** out_hltprop){
  if (out_hltprop) *out_hltprop = nullptr;
  for (auto const& hltprop:runRangeExcluded_HLTprop_list){
    if (!hltprop->isSameTrigger(name)) continue;
    if (out_hltprop) *out_hltprop = hltprop;
    return true;
  }
  return false;
}

void TriggerHelpers::assignTriggerObjectCheckException(std::string const& name, HLTTriggerPathProperties::TriggerObjectExceptionType const& flag){
  for (int itt=0; itt!=(int) nTriggerTypes; itt++){
    // No need to check the iterator since it was already checked.
    auto it = HLT_type_proplist_map.find((TriggerType) itt);
    for (auto& hltprop:it->second){
      if (!hltprop.isSameTrigger(name)) continue;

      MELAout << "TriggerHelpers::assignTriggerObjectCheckException: Adding TO exception of type " << flag << " to " << name << "." << endl;

      hltprop.setTOException(flag);
    }
  }
}

void TriggerHelpers::configureHLTmap(){
  if (!SampleHelpers::runConfigure){
    MELAerr << "TriggerHelpers::configureHLTmap: Need to call SampleHelpers::configure(period, tag) first!" << endl;
    assert(0);
  }

  // Clear the maps
  HLT_type_list_map.clear();
  HLT_type_proplist_map.clear();

  // Clear the list of HLT properties with run range exclusions
  runRangeExcluded_HLTprop_list.clear();

  // Notice that the triggers that require cuts are ORDERED!
  // Ordering within each list is important.
  // FIXME: Triggers with prescales need to include higher-pT thresholds as well. SinglePhoton is done, but the rest needs revision either.
  switch (SampleHelpers::getDataYear()){
  case 2016:
    HLT_type_proplist_map[kDoubleMu] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kMuon } } },
      { "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kMuon } } }
    };
    HLT_type_proplist_map[kDoubleMu_Extra] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Mu30_TkMu11_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kMuon } } },
      { "HLT_Mu40_TkMu11_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kMuon } } } // Probably no need for this
    };
    HLT_type_proplist_map[kDoubleMu_Prescaled] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kMuon } } },
      { "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kMuon } } }
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
    HLT_type_proplist_map[kSingleEle_L1EG] = std::vector<HLTTriggerPathProperties>();
    // These MuEle triggers are tricky to deal with because they are disabled for part of the run...
    HLT_type_proplist_map[kMuEle] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kElectron } } },
      { "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kElectron } } },
      { "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kElectron } } },
      { "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kElectron } } }
    };
    HLT_type_proplist_map[kMuEle_Extra] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kElectron } } },
      { "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kElectron } } },
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
      { "HLT_Ele27_eta2p1_WPLoose_Gsf_v*", { { HLTObjectProperties::kElectron } } }, // Prescale>1 for 2016H, but we include here with a run range exclusion below
      { "HLT_Ele27_WPTight_Gsf_v*", { { HLTObjectProperties::kElectron } } },
      { "HLT_Ele32_eta2p1_WPTight_Gsf_v*", { { HLTObjectProperties::kElectron } } }
    };
    HLT_type_proplist_map[kSingleEle_Prescaled] = std::vector<HLTTriggerPathProperties>();
    HLT_type_proplist_map[kSingleEle_HighPt] = std::vector<HLTTriggerPathProperties>{ { "HLT_Photon175_v*", { { HLTObjectProperties::kElectron } } } };
    HLT_type_proplist_map[kSinglePho] = std::vector<HLTTriggerPathProperties>{
      //{ "HLT_Photon300_NoHE_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 300.f+40.f } } } } }, // Unprescaled
      //{ "HLT_Photon175_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 200.f }, { HLTObjectProperties::kPtHigh, 300.f+40.f } } } } }, // Unprescaled
      { "HLT_Photon175_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 175.f+45.f } } } } }, // Unprescaled
      { "HLT_Photon165_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 165.f+35.f }, { HLTObjectProperties::kPtHigh, 175.f+45.f } } } } },
      //{ "HLT_Photon165_HE10_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 165.f*1.1f }, { HLTObjectProperties::kPtHigh, 175.f*1.1f } } } } },
      //{ "HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_VBF_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 120.f*1.1f }, { HLTObjectProperties::kPtHigh, 165.f*1.1f } } } } },
      { "HLT_Photon120_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 120.f+15.f }, { HLTObjectProperties::kPtHigh, 165.f+35.f } } } } },
      //{ "HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_VBF_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 90.f*1.1f }, { HLTObjectProperties::kPtHigh, 120.f*1.1f } } } } },
      { "HLT_Photon90_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 90.f*1.1f }, { HLTObjectProperties::kPtHigh, 120.f+15.f } } } } },
      //{ "HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_VBF_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 75.f*1.1f }, { HLTObjectProperties::kPtHigh, 90.f*1.1f } } } } },
      { "HLT_Photon75_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 75.f*1.1f }, { HLTObjectProperties::kPtHigh, 90.f*1.1f } } } } },
      //{ "HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_VBF_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 50.f*1.1f }, { HLTObjectProperties::kPtHigh, 75.f*1.1f } } } } },
      { "HLT_Photon50_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 50.f*1.1f }, { HLTObjectProperties::kPtHigh, 75.f*1.1f } } } } },
      //{ "HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_VBF_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 36.f*1.1f }, { HLTObjectProperties::kPtHigh, 50.f*1.1f } } } } },
      { "HLT_Photon36_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 36.f*1.1f }, { HLTObjectProperties::kPtHigh, 50.f*1.1f } } } } }, // Has L1 prescales as well
      { "HLT_Photon30_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 30.f*1.1f }, { HLTObjectProperties::kPtHigh, 36.f*1.1f } } } } }, // Has L1 prescales as well
      //{ "HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_VBF_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 22.f*1.1f }, { HLTObjectProperties::kPtHigh, 30.f*1.1f } } } } }, // Has L1 prescales as well
      { "HLT_Photon22_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 22.f*1.1f }, { HLTObjectProperties::kPtHigh, 30.f*1.1f } } } } } // Has L1 prescales as well
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
    HLT_type_proplist_map[kVBFJets_Control] = std::vector<HLTTriggerPathProperties>{
      { "HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu140_v*", { { HLTObjectProperties::kAK4Jet, { { HLTObjectProperties::kPt, 40.f } } }, { HLTObjectProperties::kAK4Jet, { { HLTObjectProperties::kPt, 40.f } } }, { HLTObjectProperties::kAK4DiJetSumWithDEtaDPhi, { { HLTObjectProperties::kEtaLow, 3.5f }, { HLTObjectProperties::kMass, 600.f } } }, { HLTObjectProperties::kMET_NoMu, { { HLTObjectProperties::kPt, 140.f } } } } }
    };
    HLT_type_proplist_map[kPFHT_Control] = std::vector<HLTTriggerPathProperties>{
      { "HLT_PFHT900_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 1000.f } } } } }, // Prescaled
      //{ "HLT_PFHT800_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 800.f } } } } }, // Prescaled, but it doesn't exist in 2016H. Therefore, it is disabled.
      { "HLT_PFHT650_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 740.f }, { HLTObjectProperties::kPtHigh, 1000.f } } } } }, // Prescaled
      { "HLT_PFHT600_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 690.f }, { HLTObjectProperties::kPtHigh, 740.f } } } } }, // Prescaled
      { "HLT_PFHT475_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 600.f }, { HLTObjectProperties::kPtHigh, 690.f } } } } }, // Prescaled
      { "HLT_PFHT400_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 550.f }, { HLTObjectProperties::kPtHigh, 600.f } } } } }, // Prescaled
      { "HLT_PFHT350_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 520.f }, { HLTObjectProperties::kPtHigh, 550.f } } } } }, // Prescaled
      { "HLT_PFHT300_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 500.f }, { HLTObjectProperties::kPtHigh, 520.f } } } } }, // Prescaled
      { "HLT_PFHT250_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 480.f }, { HLTObjectProperties::kPtHigh, 500.f } } } } }, // Prescaled
      { "HLT_PFHT200_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 440.f }, { HLTObjectProperties::kPtHigh, 480.f } } } } }, // Prescaled
      { "HLT_PFHT125_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 350.f }, { HLTObjectProperties::kPtHigh, 440.f } } } } } // Prescaled
    };
    HLT_type_proplist_map[kMET_Control] = std::vector<HLTTriggerPathProperties>{
      { "HLT_MET600_v*", { { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 600.f } } } } },
      { "HLT_MET300_v*", { { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 300.f } } } } },
      { "HLT_MET250_v*", { { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 250.f } } } } },
      { "HLT_MET200_v*", { { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 200.f } } } } }
    };
    HLT_type_proplist_map[kPFMET_Control] = std::vector<HLTTriggerPathProperties>{
      { "HLT_PFMET600_v*", { { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 600.f } } } } },
      { "HLT_PFMET500_v*", { { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 500.f } } } } },
      { "HLT_PFMET400_v*", { { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 400.f } } } } },
      { "HLT_PFMET300_v*", { { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 300.f } } } } },
      { "HLT_PFMET170_HBHECleaned_v*", { { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 170.f } } } } }
    };
    HLT_type_proplist_map[kPFHT_PFMET_Control] = std::vector<HLTTriggerPathProperties>{
      { "HLT_PFHT300_PFMET110_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 300.f } } }, { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 110.f } } } } }
    };
    HLT_type_proplist_map[kPFMET_MHT_Control] = std::vector<HLTTriggerPathProperties>{
      { "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v*", { { HLTObjectProperties::kMET_NoMu, { { HLTObjectProperties::kPt, 120.f } } }, { HLTObjectProperties::kHT_NoMu, { { HLTObjectProperties::kMass, 120.f } } } } },
      { "HLT_PFMET120_PFMHT120_IDTight_v*", { { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 120.f } } }, { HLTObjectProperties::kHT, { { HLTObjectProperties::kMass, 120.f } } } } }
    };
    HLT_type_proplist_map[kPFHT_PFMET_MHT_Control] = std::vector<HLTTriggerPathProperties>();

    // Assign TO check exceptions to recover passing legs from the failed ones
    assignTriggerObjectCheckException("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*", HLTTriggerPathProperties::toRecoverObjectsFromFailing);
    assignTriggerObjectCheckException("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v*", HLTTriggerPathProperties::toRecoverObjectsFromFailing);
    assignTriggerObjectCheckException("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*", HLTTriggerPathProperties::toRecoverObjectsFromFailing);
    assignTriggerObjectCheckException("HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v*", HLTTriggerPathProperties::toRecoverObjectsFromFailing);
    assignTriggerObjectCheckException("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v*", HLTTriggerPathProperties::toRecoverObjectsFromFailing);
    assignTriggerObjectCheckException("HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v*", HLTTriggerPathProperties::toRecoverObjectsFromFailing);
    assignTriggerObjectCheckException("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v*", HLTTriggerPathProperties::toRecoverObjectsFromFailing);

    // Assign run range exclusions
    assignRunRangeExclusions(
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*", {
        { 280919/*281613*/, -1 }
      }
    );
    assignRunRangeExclusions(
      "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*", {
        { 280919/*281613*/, -1 }
      }
    );
    assignRunRangeExclusions(
      "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*", {
        { 274968, -1 } // Prescales in 2016D-H are crazy, they interchange between 0 and 1.
      }
    );
    assignRunRangeExclusions(
      "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*", {
        { 278873, -1 }
      }
    );
    /*
    Note on HLT_DoubleEle33_CaloIdL_MW_v*:
    - Egamma POG recommends 276453 as the lower boundary, but the golden JSON has 276454 currently as the first run after that.
    They also mention this trigger 'has a bug in that only one of the electrons is required to pass the medium window (MW) pixel requirement before 276453', but
    since this will create a mismatch between data and MC, we exclude all runs before 276453 as well.
    - HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v* and this trigger are able to cover all of Run 2016 for pT1,2>33 GeV since the first good run after 278822 is 278873.
    */
    assignRunRangeExclusions(
      "HLT_DoubleEle33_CaloIdL_MW_v*", {
        { -1/*276453*//*276454*/, 278822 }
      }
    );
    assignRunRangeExclusions(
      "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*", {
        { 274968, -1 }
      }
    );
    assignRunRangeExclusions(
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*", {
        { 280919/*281613*/, -1 }
      }
    );
    assignRunRangeExclusions(
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*", {
        { -1, 278240 }
      }
    );
    assignRunRangeExclusions(
      "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*", {
        { 274968, -1 }
      }
    );
    assignRunRangeExclusions(
      "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v*", {
        { 280919/*281613*/, -1 }
      }
    );
    assignRunRangeExclusions(
      "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v*", {
        { -1, 278240 }
      }
    );
    assignRunRangeExclusions(
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*", {
        { 280919/*281613*/, -1 }
      }
    );
    assignRunRangeExclusions(
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*", {
        { -1, 278240 }
      }
    );
    assignRunRangeExclusions(
      "HLT_TkMu50_v*", {
        { -1, 274442 }
      }
    );
    assignRunRangeExclusions(
      "HLT_Ele27_eta2p1_WPLoose_Gsf_v*", {
        { 280919/*281613*/, -1 } // Prescale>1 for 2016H
      }
    );
    // HLT_PFHT800_v* seems disabled in 2016H. However, this would create a gap in the HT spectrum, so it is excluded from the set of control triggers altogether.
    assignRunRangeExclusions(
      "HLT_PFHT800_v*", {
        { 280919/*281613*/, -1 }
      }
    );
    break;
  case 2017:
    HLT_type_proplist_map[kDoubleMu] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kMuon } } },
      { "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kMuon } } }
    };
    HLT_type_proplist_map[kDoubleMu_Extra] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Mu37_TkMu27_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kMuon } } }
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
    HLT_type_proplist_map[kSingleEle_L1EG] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Ele32_WPTight_Gsf_L1DoubleEG_v*", { { HLTObjectProperties::kElectron }, { HLTObjectProperties::kElectron } } }
    };
    HLT_type_proplist_map[kMuEle] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kElectron } } },
      { "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kElectron } } },
      { "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kElectron } } },
      { "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kElectron } } }
    };
    HLT_type_proplist_map[kMuEle_Extra] = std::vector<HLTTriggerPathProperties>();
    HLT_type_proplist_map[kSingleMu] = std::vector<HLTTriggerPathProperties>{ { "HLT_IsoMu27_v*", { { HLTObjectProperties::kMuon } } } };
    // HLT_IsoMu20_v* have L1 and HLT prescales
    // HLT_IsoMu24_v* is disabled for a large portion of 2017, so omit it.
    HLT_type_proplist_map[kSingleMu_Prescaled] = std::vector<HLTTriggerPathProperties>{
      //{ "HLT_IsoMu24_v*", { { HLTObjectProperties::kMuon, { { HLTObjectProperties::kPt, 24.f*1.1f } } } } },
      //{ "HLT_IsoMu20_v*", { { HLTObjectProperties::kMuon, { { HLTObjectProperties::kPt, 20.f*1.1f }, { HLTObjectProperties::kPtHigh, 24.f*1.1f } } } } }
      { "HLT_IsoMu20_v*", { { HLTObjectProperties::kMuon, { { HLTObjectProperties::kPt, 20.f*1.1f } } } } }
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
    HLT_type_proplist_map[kSingleEle_Prescaled] = std::vector<HLTTriggerPathProperties>();
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
    // HLT_AK8PFJet360_TrimMass30_v* is disabled for most of 2017, so omit it altogether.
    HLT_type_proplist_map[kAK8PFJet_Control] = std::vector<HLTTriggerPathProperties>{
      //{ "HLT_AK8PFJet360_TrimMass30_v*", { { HLTObjectProperties::kAK8Jet, { { HLTObjectProperties::kPt, 360.f }, { HLTObjectProperties::kMass, 30.f } } } } }
    };
    HLT_type_proplist_map[kVBFJets_Control] = std::vector<HLTTriggerPathProperties>{
      { "HLT_DiJet110_35_Mjj650_PFMET110_v*", { { HLTObjectProperties::kAK4Jet, { { HLTObjectProperties::kPt, 110.f } } }, { HLTObjectProperties::kAK4Jet, { { HLTObjectProperties::kPt, 35.f } } }, { HLTObjectProperties::kAK4DiJetSumWithDEtaDPhi, { { HLTObjectProperties::kMass, 650.f } } }, { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 110.f } } } } },
      { "HLT_DiJet110_35_Mjj650_PFMET120_v*", { { HLTObjectProperties::kAK4Jet, { { HLTObjectProperties::kPt, 110.f } } }, { HLTObjectProperties::kAK4Jet, { { HLTObjectProperties::kPt, 35.f } } }, { HLTObjectProperties::kAK4DiJetSumWithDEtaDPhi, { { HLTObjectProperties::kMass, 650.f } } }, { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 120.f } } } } },
      { "HLT_DiJet110_35_Mjj650_PFMET130_v*", { { HLTObjectProperties::kAK4Jet, { { HLTObjectProperties::kPt, 110.f } } }, { HLTObjectProperties::kAK4Jet, { { HLTObjectProperties::kPt, 35.f } } }, { HLTObjectProperties::kAK4DiJetSumWithDEtaDPhi, { { HLTObjectProperties::kMass, 650.f } } }, { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 130.f } } } } }
    };
    HLT_type_proplist_map[kPFHT_Control] = std::vector<HLTTriggerPathProperties>{
      { "HLT_PFHT1050_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 1150.f } } } } },
      { "HLT_PFHT890_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 1000.f }, { HLTObjectProperties::kPtHigh, 1150.f } } } } },
      { "HLT_PFHT780_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 900.f }, { HLTObjectProperties::kPtHigh, 1000.f } } } } }, // Prescaled
      { "HLT_PFHT680_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 800.f }, { HLTObjectProperties::kPtHigh, 900.f } } } } }, // Prescaled
      { "HLT_PFHT590_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 750.f }, { HLTObjectProperties::kPtHigh, 800.f } } } } }, // Prescaled
      { "HLT_PFHT510_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 700.f }, { HLTObjectProperties::kPtHigh, 750.f } } } } }, // Prescaled
      { "HLT_PFHT430_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 680.f }, { HLTObjectProperties::kPtHigh, 700.f } } } } }, // Prescaled
      { "HLT_PFHT370_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 600.f }, { HLTObjectProperties::kPtHigh, 680.f } } } } }, // Prescaled
      { "HLT_PFHT350_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 580.f }, { HLTObjectProperties::kPtHigh, 600.f } } } } }, // Prescaled
      { "HLT_PFHT250_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 550.f }, { HLTObjectProperties::kPtHigh, 580.f } } } } }, // Prescaled
      { "HLT_PFHT180_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 500.f }, { HLTObjectProperties::kPtHigh, 550.f } } } } } // Prescaled
    };
    HLT_type_proplist_map[kMET_Control] = std::vector<HLTTriggerPathProperties>();
    HLT_type_proplist_map[kPFMET_Control] = std::vector<HLTTriggerPathProperties>{
      { "HLT_PFMET300_HBHECleaned_v*", { { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 300.f } } } } },
      { "HLT_PFMET250_HBHECleaned_v*", { { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 250.f } } } } },
      { "HLT_PFMET200_HBHE_BeamHaloCleaned_v*", { { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 200.f } } } } }
    };
    HLT_type_proplist_map[kPFHT_PFMET_Control] = std::vector<HLTTriggerPathProperties>();
    HLT_type_proplist_map[kPFMET_MHT_Control] = std::vector<HLTTriggerPathProperties>{
      { "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v*", { { HLTObjectProperties::kMET_NoMu, { { HLTObjectProperties::kPt, 120.f } } }, { HLTObjectProperties::kHT_NoMu, { { HLTObjectProperties::kMass, 120.f } } }, { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 60.f } } } } },
      { "HLT_PFMET120_PFMHT120_IDTight_PFHT60_v*", { { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 120.f } } }, { HLTObjectProperties::kHT, { { HLTObjectProperties::kMass, 120.f }, { HLTObjectProperties::kPt, 60.f } } } } },
      { "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v*", { { HLTObjectProperties::kMET_NoMu, { { HLTObjectProperties::kPt, 120.f } } }, { HLTObjectProperties::kHT_NoMu, { { HLTObjectProperties::kMass, 120.f } } } } },
      { "HLT_PFMET120_PFMHT120_IDTight_v*", { { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 120.f } } }, { HLTObjectProperties::kHT, { { HLTObjectProperties::kMass, 120.f } } } } }
    };
    HLT_type_proplist_map[kPFHT_PFMET_MHT_Control] = std::vector<HLTTriggerPathProperties>{
      { "HLT_PFHT500_PFMET100_PFMHT100_IDTight_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 500.f }, { HLTObjectProperties::kMass, 100.f } } }, { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 100.f } } } } },
      { "HLT_PFHT500_PFMET110_PFMHT110_IDTight_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 500.f }, { HLTObjectProperties::kMass, 110.f } } }, { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 110.f } } } } },
      { "HLT_PFHT700_PFMET85_PFMHT85_IDTight_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 700.f }, { HLTObjectProperties::kMass, 85.f } } }, { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 85.f } } } } },
      { "HLT_PFHT700_PFMET95_PFMHT95_IDTight_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 700.f }, { HLTObjectProperties::kMass, 95.f } } }, { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 95.f } } } } },
      { "HLT_PFHT800_PFMET75_PFMHT75_IDTight_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 800.f }, { HLTObjectProperties::kMass, 75.f } } }, { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 75.f } } } } },
      { "HLT_PFHT800_PFMET85_PFMHT85_IDTight_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 800.f }, { HLTObjectProperties::kMass, 85.f } } }, { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 85.f } } } } }
    };

    assignTriggerObjectCheckException("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*", HLTTriggerPathProperties::toRecoverObjectsFromFailing);
    assignTriggerObjectCheckException("HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v*", HLTTriggerPathProperties::toRecoverObjectsFromFailing);
    assignTriggerObjectCheckException("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v*", HLTTriggerPathProperties::toRecoverObjectsFromFailing);
    assignTriggerObjectCheckException("HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v*", HLTTriggerPathProperties::toRecoverObjectsFromFailing);
    assignTriggerObjectCheckException("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v*", HLTTriggerPathProperties::toRecoverObjectsFromFailing);

    assignRunRangeExclusions(
      "HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v*", {
        { 303824/*303832*/, -1/*306462*//*306456*/ } // The actual boundary read off is 303832-306456, which is a large subset of 2017E+F. We exclude entire 2017E+F.
      }
    );
    // These triggers are replaced by the HT60 versions after 2017B.
    assignRunRangeExclusions(
      "HLT_PFMET120_PFMHT120_IDTight_v*", {
        { 305586, -1 }
      }
    );
    assignRunRangeExclusions(
      "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v*", {
        { 305586, -1 }
      }
    );
    // These triggers need to exclude 2017B because prescale=0.
    assignRunRangeExclusions(
      "HLT_PFMET250_HBHECleaned_v*", {
        { -1, 299329 }
      }
    );
    assignRunRangeExclusions(
      "HLT_PFMET300_HBHECleaned_v*", {
        { -1, 299329 }
      }
    );
    assignRunRangeExclusions(
      "HLT_PFMET120_PFMHT120_IDTight_PFHT60_v*", {
        { -1, 299329 }
      }
    );
    assignRunRangeExclusions(
      "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v*", {
        { -1, 299329 }
      }
    );
    assignRunRangeExclusions(
      "HLT_DiJet110_35_Mjj650_PFMET110_v*", {
        { -1, 299329 }
      }
    );
    assignRunRangeExclusions(
      "HLT_DiJet110_35_Mjj650_PFMET120_v*", {
        { -1, 299329 }
      }
    );
    assignRunRangeExclusions(
      "HLT_DiJet110_35_Mjj650_PFMET130_v*", {
        { -1, 299329 }
      }
    );
    // Mu50 is unprescaled for the entire run, but TkMu100 is not.
    assignRunRangeExclusions(
      "HLT_TkMu100_v*", {
        { -1, 299329 }
      }
    );
    assignRunRangeExclusions(
      "HLT_OldMu100_v*", {
        { -1, 299329 }
      }
    );
    // HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v* is unprescaled in 2017B.
    // We don't need to use it explicitly as a replacement for HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*
    // because we also have HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v*, which is fully unprescaled.
    //assignRunRangeExclusions(
    //  "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*", {
    //    { 299368, -1 }
    //  }
    //);
    assignRunRangeExclusions(
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*", {
        { -1, 299329 }
      }
    );
    assignRunRangeExclusions(
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*", {
        { -1, 299329 }
      }
    );
    assignRunRangeExclusions(
      "HLT_Mu37_TkMu27_v*", {
        { -1, 302029/*302019*/ } // No prescale=0 in the last few valid runs of 2017C, but we will ignore that and exclude all of 2017B+C.
      }
    );
    assignRunRangeExclusions(
      "HLT_PFMET200_HBHE_BeamHaloCleaned_v*", {
        { -1, 304797 }
      }
    );
    break;
  case 2018:
    HLT_type_proplist_map[kDoubleMu] = std::vector<HLTTriggerPathProperties>{
      // No gain from adding the Mass8 version of the one below
      { "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kMuon } } }
    };
    HLT_type_proplist_map[kDoubleMu_Extra] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Mu37_TkMu27_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kMuon } } }
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
    HLT_type_proplist_map[kSingleEle_L1EG] = std::vector<HLTTriggerPathProperties>();
    HLT_type_proplist_map[kMuEle] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kElectron } } },
      { "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kElectron } } },
      { "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kElectron } } },
      { "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kElectron } } }
    };
    HLT_type_proplist_map[kMuEle_Extra] = std::vector<HLTTriggerPathProperties>();
    HLT_type_proplist_map[kSingleMu] = std::vector<HLTTriggerPathProperties>{ { "HLT_IsoMu24_v*", { { HLTObjectProperties::kMuon } } } };
    HLT_type_proplist_map[kSingleMu_Prescaled] = std::vector<HLTTriggerPathProperties>{ { "HLT_IsoMu20_v*", { { HLTObjectProperties::kMuon } } } };
    HLT_type_proplist_map[kSingleMu_HighPt] = std::vector<HLTTriggerPathProperties>{
      { "HLT_TkMu100_v*", { { HLTObjectProperties::kMuon } } },
      { "HLT_OldMu100_v*", { { HLTObjectProperties::kMuon } } },
      { "HLT_Mu50_v*", { { HLTObjectProperties::kMuon } } }
    };
    HLT_type_proplist_map[kSingleEle] = std::vector<HLTTriggerPathProperties>{ { "HLT_Ele32_WPTight_Gsf_v*", { { HLTObjectProperties::kElectron } } } };
    HLT_type_proplist_map[kSingleEle_Prescaled] = std::vector<HLTTriggerPathProperties>();
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
    HLT_type_proplist_map[kVBFJets_Control] = std::vector<HLTTriggerPathProperties>{
      { "HLT_DiJet110_35_Mjj650_PFMET110_v*", { { HLTObjectProperties::kAK4Jet, { { HLTObjectProperties::kPt, 110.f } } }, { HLTObjectProperties::kAK4Jet, { { HLTObjectProperties::kPt, 35.f } } }, { HLTObjectProperties::kAK4DiJetSumWithDEtaDPhi, { { HLTObjectProperties::kMass, 650.f } } }, { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 110.f } } } } },
      { "HLT_DiJet110_35_Mjj650_PFMET120_v*", { { HLTObjectProperties::kAK4Jet, { { HLTObjectProperties::kPt, 110.f } } }, { HLTObjectProperties::kAK4Jet, { { HLTObjectProperties::kPt, 35.f } } }, { HLTObjectProperties::kAK4DiJetSumWithDEtaDPhi, { { HLTObjectProperties::kMass, 650.f } } }, { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 120.f } } } } },
      { "HLT_DiJet110_35_Mjj650_PFMET130_v*", { { HLTObjectProperties::kAK4Jet, { { HLTObjectProperties::kPt, 110.f } } }, { HLTObjectProperties::kAK4Jet, { { HLTObjectProperties::kPt, 35.f } } }, { HLTObjectProperties::kAK4DiJetSumWithDEtaDPhi, { { HLTObjectProperties::kMass, 650.f } } }, { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 130.f } } } } }
    };
    HLT_type_proplist_map[kPFHT_Control] = std::vector<HLTTriggerPathProperties>{
      { "HLT_PFHT1050_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 1150.f } } } } },
      { "HLT_PFHT890_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 1000.f }, { HLTObjectProperties::kPtHigh, 1150.f } } } } }, 
      { "HLT_PFHT780_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 900.f }, { HLTObjectProperties::kPtHigh, 1000.f } } } } }, // Prescaled
      { "HLT_PFHT680_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 800.f }, { HLTObjectProperties::kPtHigh, 900.f } } } } }, // Prescaled
      { "HLT_PFHT590_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 750.f }, { HLTObjectProperties::kPtHigh, 800.f } } } } }, // Prescaled
      { "HLT_PFHT510_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 700.f }, { HLTObjectProperties::kPtHigh, 750.f } } } } }, // Prescaled
      { "HLT_PFHT430_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 680.f }, { HLTObjectProperties::kPtHigh, 700.f } } } } }, // Prescaled
      { "HLT_PFHT370_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 600.f }, { HLTObjectProperties::kPtHigh, 680.f } } } } }, // Prescaled
      { "HLT_PFHT350_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 580.f }, { HLTObjectProperties::kPtHigh, 600.f } } } } }, // Prescaled
      { "HLT_PFHT250_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 550.f }, { HLTObjectProperties::kPtHigh, 580.f } } } } }, // Prescaled
      { "HLT_PFHT180_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 500.f }, { HLTObjectProperties::kPtHigh, 550.f } } } } } // Prescaled
    };
    HLT_type_proplist_map[kMET_Control] = std::vector<HLTTriggerPathProperties>();
    HLT_type_proplist_map[kPFMET_Control] = std::vector<HLTTriggerPathProperties>{
      { "HLT_PFMET300_HBHECleaned_v*", { { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 300.f } } } } },
      { "HLT_PFMET250_HBHECleaned_v*", { { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 250.f } } } } },
      { "HLT_PFMET200_HBHE_BeamHaloCleaned_v*", { { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 200.f } } } } }
    };
    HLT_type_proplist_map[kPFHT_PFMET_Control] = std::vector<HLTTriggerPathProperties>();
    HLT_type_proplist_map[kPFMET_MHT_Control] = std::vector<HLTTriggerPathProperties>{
      { "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v*", { { HLTObjectProperties::kMET_NoMu, { { HLTObjectProperties::kPt, 120.f } } }, { HLTObjectProperties::kHT_NoMu, { { HLTObjectProperties::kMass, 120.f } } }, { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 60.f } } } } },
      { "HLT_PFMET120_PFMHT120_IDTight_PFHT60_v*", { { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 120.f } } }, { HLTObjectProperties::kHT, { { HLTObjectProperties::kMass, 120.f }, { HLTObjectProperties::kPt, 60.f } } } } },
      { "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v*", { { HLTObjectProperties::kMET_NoMu, { { HLTObjectProperties::kPt, 120.f } } }, { HLTObjectProperties::kHT_NoMu, { { HLTObjectProperties::kMass, 120.f } } } } },
      { "HLT_PFMET120_PFMHT120_IDTight_v*", { { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 120.f } } }, { HLTObjectProperties::kHT, { { HLTObjectProperties::kMass, 120.f } } } } }
    };
    HLT_type_proplist_map[kPFHT_PFMET_MHT_Control] = std::vector<HLTTriggerPathProperties>{
      { "HLT_PFHT500_PFMET100_PFMHT100_IDTight_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 500.f }, { HLTObjectProperties::kMass, 100.f } } }, { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 100.f } } } } },
      { "HLT_PFHT500_PFMET110_PFMHT110_IDTight_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 500.f }, { HLTObjectProperties::kMass, 110.f } } }, { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 110.f } } } } },
      { "HLT_PFHT700_PFMET85_PFMHT85_IDTight_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 700.f }, { HLTObjectProperties::kMass, 85.f } } }, { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 85.f } } } } },
      { "HLT_PFHT700_PFMET95_PFMHT95_IDTight_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 700.f }, { HLTObjectProperties::kMass, 95.f } } }, { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 95.f } } } } },
      { "HLT_PFHT800_PFMET75_PFMHT75_IDTight_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 800.f }, { HLTObjectProperties::kMass, 75.f } } }, { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 75.f } } } } },
      { "HLT_PFHT800_PFMET85_PFMHT85_IDTight_v*", { { HLTObjectProperties::kHT, { { HLTObjectProperties::kPt, 800.f }, { HLTObjectProperties::kMass, 85.f } } }, { HLTObjectProperties::kMET, { { HLTObjectProperties::kPt, 85.f } } } } }
    };

    assignTriggerObjectCheckException("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*", HLTTriggerPathProperties::toRecoverObjectsFromFailing);
    assignTriggerObjectCheckException("HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v*", HLTTriggerPathProperties::toRecoverObjectsFromFailing);
    assignTriggerObjectCheckException("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v*", HLTTriggerPathProperties::toRecoverObjectsFromFailing);
    assignTriggerObjectCheckException("HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v*", HLTTriggerPathProperties::toRecoverObjectsFromFailing);
    assignTriggerObjectCheckException("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v*", HLTTriggerPathProperties::toRecoverObjectsFromFailing);

    // Prescale=0 for a small portion of 2018A
    assignRunRangeExclusions(
      "HLT_DiJet110_35_Mjj650_PFMET110_v*", {
        { 316153, 316239 }
      }
    );
    assignRunRangeExclusions(
      "HLT_DiJet110_35_Mjj650_PFMET120_v*", {
        { 316153, 316239 }
      }
    );
    assignRunRangeExclusions(
      "HLT_DiJet110_35_Mjj650_PFMET130_v*", {
        { 316153, 316239 }
      }
    );
    break;
  default:
    break;
  }

  // Check that all triggers are defined
  for (int itt=0; itt!=(int) nTriggerTypes; itt++){
    if (HLT_type_proplist_map.find((TriggerType) itt)==HLT_type_proplist_map.cend()){
      MELAerr << "TriggerHelpers::configureHLTmap: Triggers for type " << itt << " are not defined for year " << SampleHelpers::getDataYear() << ". Please fix the map HLT_type_proplist_map." << endl;
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
