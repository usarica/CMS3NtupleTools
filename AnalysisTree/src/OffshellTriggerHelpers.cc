#include <cassert>
#include "HLTTriggerPathProperties.h"
#include "OffshellTriggerHelpers.h"
#include "MELAStreamHelpers.hh"


namespace OffshellTriggerHelpers{
  std::unordered_map< OffshellTriggerHelpers::TriggerType, std::vector<std::string> > HLT_type_list_map;
  std::unordered_map< OffshellTriggerHelpers::TriggerType, std::vector<HLTTriggerPathProperties> > HLT_type_proplist_map;
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
    HLT_type_proplist_map[kDoubleMu] = std::vector<HLTTriggerPathProperties>{
      //{ "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kMuon } } },
      //{ "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kMuon } } },
      { "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*", { { HLTObjectProperties::kMuon },{ HLTObjectProperties::kMuon } } },
      { "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*", { { HLTObjectProperties::kMuon },{ HLTObjectProperties::kMuon } } }
    };
    HLT_type_proplist_map[kDoubleEle] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*", { { HLTObjectProperties::kElectron },{ HLTObjectProperties::kElectron } } },
      { "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*", { { HLTObjectProperties::kElectron },{ HLTObjectProperties::kElectron } } },
      { "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*", { { HLTObjectProperties::kElectron },{ HLTObjectProperties::kElectron } } }
    };
    HLT_type_proplist_map[kMuEle] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*",{ { HLTObjectProperties::kMuon },{ HLTObjectProperties::kElectron } } },
      { "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*",{ { HLTObjectProperties::kMuon },{ HLTObjectProperties::kElectron } } },
      { "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*",{ { HLTObjectProperties::kMuon },{ HLTObjectProperties::kElectron } } },
      { "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*",{ { HLTObjectProperties::kMuon },{ HLTObjectProperties::kElectron } } },
      { "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v*",{ { HLTObjectProperties::kMuon },{ HLTObjectProperties::kElectron } } },
      { "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*",{ { HLTObjectProperties::kMuon },{ HLTObjectProperties::kElectron } } },
      { "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",{ { HLTObjectProperties::kMuon },{ HLTObjectProperties::kElectron } } }
    };
    HLT_type_proplist_map[kSingleMu] = std::vector<HLTTriggerPathProperties>{
      { "HLT_IsoMu20_v*",{ { HLTObjectProperties::kMuon } } },
      { "HLT_IsoTkMu20_v*",{ { HLTObjectProperties::kMuon } } },
      { "HLT_IsoMu22_v*",{ { HLTObjectProperties::kMuon } } },
      { "HLT_IsoTkMu22_v*",{ { HLTObjectProperties::kMuon } } },
      { "HLT_IsoMu24_v*",{ { HLTObjectProperties::kMuon } } },
      { "HLT_IsoTkMu24_v*",{ { HLTObjectProperties::kMuon } } }
    };
    HLT_type_proplist_map[kSingleEle] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Ele25_eta2p1_WPTight_Gsf_v*",{ { HLTObjectProperties::kElectron } } },
      { "HLT_Ele27_WPTight_Gsf_v*",{ { HLTObjectProperties::kElectron } } },
      { "HLT_Ele27_eta2p1_WPLoose_Gsf_v*",{ { HLTObjectProperties::kElectron } } },
      { "HLT_Ele32_eta2p1_WPTight_Gsf_v*",{ { HLTObjectProperties::kElectron } } }
    };
    HLT_type_proplist_map[kSinglePho] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Photon22_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 22.f } } } } },
      //{ "HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_VBF_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 22.f } } } } },
      { "HLT_Photon30_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 30.f } } } } },
      { "HLT_Photon36_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 36.f } } } } },
      //{ "HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_VBF_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 36.f } } } } },
      { "HLT_Photon50_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 50.f } } } } },
      //{ "HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_VBF_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 50.f } } } } },
      { "HLT_Photon75_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 75.f } } } } },
      //{ "HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_VBF_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 75.f } } } } },
      { "HLT_Photon90_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 90.f } } } } },
      //{ "HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_VBF_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 90.f } } } } },
      { "HLT_Photon120_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 120.f } } } } },
      { "HLT_Photon120_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 120.f } } } } },
      //{ "HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_VBF_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 120.f } } } } },
      { "HLT_Photon165_HE10_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 165.f } } } } },
      { "HLT_Photon165_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 165.f } } } } },
      { "HLT_Photon175_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 175.f } } } } },
      { "HLT_Photon300_NoHE_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 300.f } } } } } // Unprescaled
    };
    HLT_type_proplist_map[kTripleLep] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v*" },
      { "HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v*" },
      { "HLT_TripleMu_12_10_5_v*" },
      { "HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v*" }
    };
    break;
  case 2017:
    HLT_type_proplist_map[kDoubleMu] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*",{ { HLTObjectProperties::kMuon },{ HLTObjectProperties::kMuon } } },
      { "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v*",{ { HLTObjectProperties::kMuon },{ HLTObjectProperties::kMuon } } }
    };
    HLT_type_proplist_map[kDoubleEle] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*",{ { HLTObjectProperties::kElectron },{ HLTObjectProperties::kElectron } } },
      { "HLT_DoubleEle33_CaloIdL_MW_v*",{ { HLTObjectProperties::kElectron },{ HLTObjectProperties::kElectron } } }
    };
    HLT_type_proplist_map[kMuEle] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*",{ { HLTObjectProperties::kMuon },{ HLTObjectProperties::kElectron } } },
      { "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",{ { HLTObjectProperties::kMuon },{ HLTObjectProperties::kElectron } } },
      { "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*",{ { HLTObjectProperties::kMuon },{ HLTObjectProperties::kElectron } } },
      { "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*",{ { HLTObjectProperties::kMuon },{ HLTObjectProperties::kElectron } } }
    };
    HLT_type_proplist_map[kSingleMu] = std::vector<HLTTriggerPathProperties>{ { "HLT_IsoMu27_v*",{ { HLTObjectProperties::kMuon } } } };
    HLT_type_proplist_map[kSingleEle] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Ele35_WPTight_Gsf_v*",{ { HLTObjectProperties::kElectron } } },
      { "HLT_Ele38_WPTight_Gsf_v*",{ { HLTObjectProperties::kElectron } } },
      { "HLT_Ele40_WPTight_Gsf_v*",{ { HLTObjectProperties::kElectron } } }
    };
    HLT_type_proplist_map[kSinglePho] = std::vector<HLTTriggerPathProperties>{
      //{ "HLT_Photon20_HoverELoose_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 20.f } } } } }, // Only runs 299368 - 306460 (36.75 fb-1) // FIXME: We should get this trigger?
      //{ "HLT_Photon25_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 25.f } } } } }, // Only runs 302026 - 306460 (27.13 fb-1) // FIXME: We should get this trigger?
      //{ "HLT_Photon30_HoverELoose_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 30.f } } } } }, // Only runs 299368 - 306460 (36.75 fb-1) // FIXME: We should get this trigger?
      { "HLT_Photon33_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 33.f } } } } },
      { "HLT_Photon50_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 50.f } } } } },
      { "HLT_Photon75_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 75.f } } } } },
      { "HLT_Photon90_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 90.f } } } } },
      { "HLT_Photon120_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 120.f } } } } },
      { "HLT_Photon120_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 120.f } } } } },
      { "HLT_Photon165_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 165.f } } } } },
      { "HLT_Photon175_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 175.f } } } } },
      { "HLT_Photon200_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 200.f } } } } },
      { "HLT_Photon300_NoHE_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 300.f } } } } } // Unprescaled
    };
    HLT_type_proplist_map[kTripleLep] = std::vector<HLTTriggerPathProperties>{
      { "HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v*" },
      { "HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v*" }, { "HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_v*" },
      { "HLT_TripleMu_10_5_5_DZ_v*" }, { "HLT_TripleMu_12_10_5_v*" },
      { "HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v*" }
    };
    break;
  case 2018:
    HLT_type_proplist_map[kDoubleMu] = std::vector<HLTTriggerPathProperties>{
      // No gain from adding the Mass8 version of the one below
      { "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*", { { HLTObjectProperties::kMuon }, { HLTObjectProperties::kMuon } } }
    };
    HLT_type_proplist_map[kDoubleEle] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*", { { HLTObjectProperties::kElectron },{ HLTObjectProperties::kElectron } } },
      //{ "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*", { { HLTObjectProperties::kElectron },{ HLTObjectProperties::kElectron } } },
      { "HLT_DoubleEle25_CaloIdL_MW_v*",{ { HLTObjectProperties::kElectron },{ HLTObjectProperties::kElectron } } }
    };
    HLT_type_proplist_map[kMuEle] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*",{ { HLTObjectProperties::kMuon },{ HLTObjectProperties::kElectron } } },
      { "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",{ { HLTObjectProperties::kMuon },{ HLTObjectProperties::kElectron } } },
      { "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*",{ { HLTObjectProperties::kMuon },{ HLTObjectProperties::kElectron } } },
      { "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*",{ { HLTObjectProperties::kMuon },{ HLTObjectProperties::kElectron } } }
    };
    HLT_type_proplist_map[kSingleMu] = std::vector<HLTTriggerPathProperties>{ { "HLT_IsoMu24_v*",{ { HLTObjectProperties::kMuon } } } };
    HLT_type_proplist_map[kSingleEle] = std::vector<HLTTriggerPathProperties>{ { "HLT_Ele32_WPTight_Gsf_v*",{ { HLTObjectProperties::kElectron } } } };
    HLT_type_proplist_map[kSinglePho] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Photon20_HoverELoose_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 20.f } } } } },
      { "HLT_Photon30_HoverELoose_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 30.f } } } } },
      { "HLT_Photon33_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 33.f } } } } },
      { "HLT_Photon50_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 50.f } } } } },
      { "HLT_Photon75_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 75.f } } } } },
      { "HLT_Photon90_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 90.f } } } } },
      { "HLT_Photon120_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 120.f } } } } },
      { "HLT_Photon120_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 120.f } } } } },
      { "HLT_Photon165_R9Id90_HE10_IsoM_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 165.f } } } } },
      { "HLT_Photon175_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 175.f } } } } },
      { "HLT_Photon200_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 200.f } } } } },
      { "HLT_Photon300_NoHE_v*", { { HLTObjectProperties::kPhoton, { { HLTObjectProperties::kPt, 300.f } } } } } // Unprescaled
    };
    HLT_type_proplist_map[kTripleLep] = std::vector<HLTTriggerPathProperties>{
      { "HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v*" },
      { "HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_v*" },
      //{ "HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v*" }, // Somehow its effective lumi is not the same as active lumi, and it does not show up in the spreadsheets
      { "HLT_TripleMu_10_5_5_DZ_v*" },
      { "HLT_TripleMu_12_10_5_v*" }
    };
    HLT_type_proplist_map[kSingleMu_Control] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Mu8_v*",{ { HLTObjectProperties::kMuon, { { HLTObjectProperties::kPt, 8.f } } } } }, // Has L1 and HLT prescales
      { "HLT_Mu8_TrkIsoVVL_v*",{ { HLTObjectProperties::kMuon,{ { HLTObjectProperties::kPt, 8.f } } } } }, // Has L1 and HLT prescales
      { "HLT_Mu17_v*",{ { HLTObjectProperties::kMuon,{ { HLTObjectProperties::kPt, 17.f } } } } }, // Has L1 and HLT prescales
      { "HLT_Mu17_TrkIsoVVL_v*",{ { HLTObjectProperties::kMuon,{ { HLTObjectProperties::kPt, 17.f } } } } } // Has L1 and HLT prescales
    };
    // FIXME: New production needed for kSingleMu_Control_HighPt triggers in 2018
    HLT_type_proplist_map[kSingleMu_Control_HighPt] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Mu50_v*",{ { HLTObjectProperties::kMuon,{ { HLTObjectProperties::kPt, 50.f } } } } }, // For high-pT single muon eff.
      { "HLT_OldMu100_v*",{ { HLTObjectProperties::kMuon,{ { HLTObjectProperties::kPt, 100.f } } } } }, // For high-pT single muon eff.
      { "HLT_TkMu100_v*",{ { HLTObjectProperties::kMuon,{ { HLTObjectProperties::kPt, 100.f } } } } } // For high-pT single muon eff.
    };
    HLT_type_proplist_map[kSingleEle_Control] = std::vector<HLTTriggerPathProperties>{
      { "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v*" }, // Has L1 prescales
      { "HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v*" }, // Has L1 prescales
      { "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v*" }, // Has L1 and HLT prescales
      { "HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v*" } // Has L1 and HLT prescales
    };
    HLT_type_proplist_map[kAK8PFJet_Control] = std::vector<HLTTriggerPathProperties>{
      { "HLT_AK8PFJet400_TrimMass30_v*" }
    };
    // FIXME: 200313 production contains only HLT_PFHT1050_v*
    HLT_type_proplist_map[kPFHT_Control] = std::vector<HLTTriggerPathProperties>{
      { "HLT_PFHT250_v*" }, // Prescaled
      { "HLT_PFHT350_v*" }, // Prescaled
      { "HLT_PFHT370_v*" }, // Prescaled
      { "HLT_PFHT430_v*" }, // Prescaled
      { "HLT_PFHT510_v*" }, // Prescaled
      { "HLT_PFHT590_v*" }, // Prescaled
      { "HLT_PFHT680_v*" }, // Prescaled
      { "HLT_PFHT780_v*" }, // Prescaled
      { "HLT_PFHT890_v*" }, // Prescaled
      { "HLT_PFHT1050_v*" }
    };
    HLT_type_proplist_map[kPFMET_MHT_Control] = std::vector<HLTTriggerPathProperties>{
      { "HLT_PFMET120_PFMHT120_IDTight_v*" },
      { "HLT_PFMET120_PFMHT120_IDTight_PFHT60_v*" },
      { "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v*" },
      { "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v*" }
    };
    break;
  default:
    break;
  }

  for (auto const& it:HLT_type_proplist_map){
    auto const& props = it.second;
    std::vector<std::string> tmplist; tmplist.reserve(props.size());
    for (auto const& prop:props) tmplist.push_back(prop.getName());
    HLT_type_list_map[it.first] = tmplist;
  }
}
