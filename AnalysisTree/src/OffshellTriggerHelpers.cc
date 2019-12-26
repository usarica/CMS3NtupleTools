#include <cassert>
#include "OffshellTriggerHelpers.h"
#include "MELAStreamHelpers.hh"


namespace OffshellTriggerHelpers{
  std::unordered_map<OffshellTriggerHelpers::TriggerType, std::vector<std::string>> HLT_type_list_map;
}


using namespace std;
using namespace MELAStreamHelpers;


std::vector<std::string> OffshellTriggerHelpers::getHLTMenus(OffshellTriggerHelpers::TriggerType type){ return getHLTMenus(std::vector<OffshellTriggerHelpers::TriggerType>{ type }); }
std::vector<std::string> OffshellTriggerHelpers::getHLTMenus(std::vector<OffshellTriggerHelpers::TriggerType> const& types){
  std::vector<std::string> res;
  for (auto const& type:types){
    auto it = HLT_type_list_map.find(type);
    if (it != HLT_type_list_map.end()) HelperFunctions::appendVector(res, it->second);
  }
  return res;
}
void OffshellTriggerHelpers::configureHLTmap(){
  if (!SampleHelpers::runConfigure){
    MELAerr << "OffshellTriggerHelpers::configureHLTmap: Need to call SampleHelpers::configure(period, tag) first!" << endl;
    assert(0);
  }
  // All triggers except the single photon ones below are unprescaled.
  switch (SampleHelpers::theDataYear){
  case 2016:
    HLT_type_list_map[kDoubleMu] = std::vector<std::string>{
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*",
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*"
    };
    HLT_type_list_map[kDoubleEle] = std::vector<std::string>{
      "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
      "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*"
    };
    HLT_type_list_map[kMuEle] = std::vector<std::string>{
      "HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v*","HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v*",
      "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*", "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*", "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*",
      "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*", "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*", "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v*"
    };
    HLT_type_list_map[kSingleMu] = std::vector<std::string>{ "HLT_IsoMu20_v*", "HLT_IsoTkMu20_v*", "HLT_IsoMu22_v*", "HLT_IsoTkMu22_v*", "HLT_IsoMu24_v*", "HLT_IsoTkMu24_v*" };
    HLT_type_list_map[kSingleEle] = std::vector<std::string>{ "HLT_Ele25_eta2p1_WPTight_Gsf_v*","HLT_Ele27_WPTight_Gsf_v*","HLT_Ele27_eta2p1_WPLoose_Gsf_v*", "HLT_Ele32_eta2p1_WPTight_Gsf_v*" };
    HLT_type_list_map[kSinglePho] = std::vector<std::string>{
      "HLT_Photon22_R9Id90_HE10_IsoM_v*", "HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_VBF_v*",
      "HLT_Photon30_R9Id90_HE10_IsoM_v*",
      "HLT_Photon36_R9Id90_HE10_IsoM_v*", "HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_VBF_v*",
      "HLT_Photon50_R9Id90_HE10_IsoM_v*", "HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_VBF_v*",
      "HLT_Photon75_R9Id90_HE10_IsoM_v*", "HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_VBF_v*",
      "HLT_Photon90_R9Id90_HE10_IsoM_v*", "HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_VBF_v*",
      /*"HLT_Photon120_v*", */"HLT_Photon120_R9Id90_HE10_IsoM_v*", "HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_VBF_v*",
      "HLT_Photon165_HE10_v*", "HLT_Photon165_R9Id90_HE10_IsoM_v*",
      "HLT_Photon175_v*",
      "HLT_Photon300_NoHE_v*" // Unprescaled
    };
    HLT_type_list_map[kTripleLep] = std::vector<std::string>{
      "HLT_TripleMu_12_10_5_v*", "HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v*"
    };
    break;
  case 2017:
    HLT_type_list_map[kDoubleMu] = std::vector<std::string>{ "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v*" };
    HLT_type_list_map[kDoubleEle] = std::vector<std::string>{ "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*","HLT_DoubleEle33_CaloIdL_MW_v*" };
    HLT_type_list_map[kMuEle] = std::vector<std::string>{
      "HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v*","HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v*","HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_v*",
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*", "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*", "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*"
    };
    HLT_type_list_map[kSingleMu] = std::vector<std::string>{ "HLT_IsoMu27_v*" };
    HLT_type_list_map[kSingleEle] = std::vector<std::string>{ "HLT_Ele35_WPTight_Gsf_v*","HLT_Ele38_WPTight_Gsf_v*","HLT_Ele40_WPTight_Gsf_v*" };
    HLT_type_list_map[kSinglePho] = std::vector<std::string>{
      //"HLT_Photon20_HoverELoose_v*", // Only runs 299368 - 306460 (36.75 fb-1) // FIXME: We should get this trigger?
      //"HLT_Photon25_v*", // Only runs 302026 - 306460 (27.13 fb-1) // FIXME: We should get this trigger?
      //"HLT_Photon30_HoverELoose_v*", // Only runs 299368 - 306460 (36.75 fb-1) // FIXME: We should get this trigger?
      "HLT_Photon33_v*",
      "HLT_Photon50_R9Id90_HE10_IsoM_v*",
      "HLT_Photon75_R9Id90_HE10_IsoM_v*",
      "HLT_Photon90_R9Id90_HE10_IsoM_v*",
      "HLT_Photon120_v*", "HLT_Photon120_R9Id90_HE10_IsoM_v*",
      "HLT_Photon165_R9Id90_HE10_IsoM_v*",
      "HLT_Photon175_v*",
      "HLT_Photon200_v*", "HLT_Photon300_NoHE_v*" // Unprescaled
    };
    HLT_type_list_map[kTripleLep] = std::vector<std::string>{
      "HLT_TripleMu_10_5_5_DZ_v*","HLT_TripleMu_12_10_5_v*",
      "HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v*"
    };
    break;
  case 2018:
    HLT_type_list_map[kDoubleMu] = std::vector<std::string>{ "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*" };
    HLT_type_list_map[kDoubleEle] = std::vector<std::string>{ "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*", "HLT_DoubleEle25_CaloIdL_MW_v*" };
    HLT_type_list_map[kMuEle] = std::vector<std::string>{
      "HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v*", "HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_v*",
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*", "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*","HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*"
    };
    HLT_type_list_map[kSingleMu] = std::vector<std::string>{ "HLT_IsoMu24_v*" };
    HLT_type_list_map[kSingleEle] = std::vector<std::string>{ "HLT_Ele32_WPTight_Gsf_v*" };
    HLT_type_list_map[kSinglePho] = std::vector<std::string>{
      "HLT_Photon20_HoverELoose_v*",
      "HLT_Photon30_HoverELoose_v*",
      "HLT_Photon33_v*",
      "HLT_Photon50_R9Id90_HE10_IsoM_v*",
      "HLT_Photon75_R9Id90_HE10_IsoM_v*",
      "HLT_Photon90_R9Id90_HE10_IsoM_v*",
      "HLT_Photon120_v*","HLT_Photon120_R9Id90_HE10_IsoM_v*",
      "HLT_Photon165_R9Id90_HE10_IsoM_v*",
      "HLT_Photon175_v*",
      "HLT_Photon200_v*", "HLT_Photon300_NoHE_v*" // Unprescaled
    };
    HLT_type_list_map[kTripleLep] = std::vector<std::string>{
      "HLT_TripleMu_10_5_5_DZ_v*", "HLT_TripleMu_12_10_5_v*"
    };
    break;
  default:
    break;
  }
}
