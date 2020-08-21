#include <cassert>
#include <CMSDataTools/AnalysisTree/interface/HostHelpersCore.h>
#include "BtagHelpers.h"
#include "SamplesCore.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace SampleHelpers;
using namespace MELAStreamHelpers;


namespace BtagHelpers{
  BtagWPType btagWPType = nBtagWPTypes;
}


void BtagHelpers::setBtagWPType(BtagWPType type){ btagWPType=type; }

float BtagHelpers::getBtagWP(bool is80X){ return getBtagWP(btagWPType, is80X); }
float BtagHelpers::getBtagWP(BtagHelpers::BtagWPType type, bool is80X){
  // From https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation
  float WP_DEEPFLAV_TIGHT(0.), WP_DEEPFLAV_MEDIUM(0.), WP_DEEPFLAV_LOOSE(0.);
  float WP_DEEPCSV_TIGHT(0.), WP_DEEPCSV_MEDIUM(0.), WP_DEEPCSV_LOOSE(0.);
  if (theDataYear == 2016 && is80X){
    // From https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco
    WP_DEEPCSV_TIGHT = 0.8958;
    WP_DEEPCSV_MEDIUM = 0.6324;
    WP_DEEPCSV_LOOSE = 0.2219;
  }
  else if (theDataYear == 2016){
    // From https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation2016Legacy
    WP_DEEPFLAV_TIGHT  = 0.7221;
    WP_DEEPFLAV_MEDIUM = 0.3093;
    WP_DEEPFLAV_LOOSE  = 0.0614;

    WP_DEEPCSV_TIGHT = 0.8953;
    WP_DEEPCSV_MEDIUM = 0.6321;
    WP_DEEPCSV_LOOSE = 0.2217;
  }
  else if (theDataYear == 2017){
    // From https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
    WP_DEEPFLAV_TIGHT  = 0.7489;
    WP_DEEPFLAV_MEDIUM = 0.3033;
    WP_DEEPFLAV_LOOSE  = 0.0521;

    WP_DEEPCSV_TIGHT = 0.8001;
    WP_DEEPCSV_MEDIUM = 0.4941;
    WP_DEEPCSV_LOOSE = 0.1522;
  }
  else if (theDataYear == 2018){
    // From https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
    WP_DEEPFLAV_TIGHT  = 0.7264;
    WP_DEEPFLAV_MEDIUM = 0.2770;
    WP_DEEPFLAV_LOOSE  = 0.0494;

    WP_DEEPCSV_TIGHT = 0.7527;
    WP_DEEPCSV_MEDIUM = 0.4184;
    WP_DEEPCSV_LOOSE = 0.1241;
  }
  else{
    MELAerr << "BtagHelpers::getBtagWP: Cannot determine the b tag SF file name for year " << theDataYear << ". Aborting..." << endl;
    assert(0);
  }
  switch (type){
  case kDeepCSV_Loose:
    return WP_DEEPCSV_LOOSE;
  case kDeepCSV_Medium:
    return WP_DEEPCSV_MEDIUM;
  case kDeepCSV_Tight:
    return WP_DEEPCSV_TIGHT;
  case kDeepFlav_Loose:
    return WP_DEEPFLAV_LOOSE;
  case kDeepFlav_Medium:
    return WP_DEEPFLAV_MEDIUM;
  case kDeepFlav_Tight:
    return WP_DEEPFLAV_TIGHT;
  default:
    MELAerr << "BtagHelpers::getBtagWP: No implementation for the b tagging WP type " << type << ". Aborting..." << endl;
    assert(0);
    return 0; // Just return to avoid warnings
  }
}

std::vector<float> BtagHelpers::getBtagWPs(bool is80X){
  std::vector<float> res; res.reserve(3);
  switch (btagWPType){
  case kDeepCSV_Loose:
  case kDeepCSV_Medium:
  case kDeepCSV_Tight:
    res.push_back(getBtagWP(kDeepCSV_Loose, is80X));
    res.push_back(getBtagWP(kDeepCSV_Medium, is80X));
    res.push_back(getBtagWP(kDeepCSV_Tight, is80X));
    break;
  case kDeepFlav_Loose:
  case kDeepFlav_Medium:
  case kDeepFlav_Tight:
    res.push_back(getBtagWP(kDeepFlav_Loose, is80X));
    res.push_back(getBtagWP(kDeepFlav_Medium, is80X));
    res.push_back(getBtagWP(kDeepFlav_Tight, is80X));
    break;
  default:
    MELAerr << "BtagHelpers::getBtagWPs: No implementation for the b tagging WP type " << btagWPType << ". Aborting..." << endl;
    assert(0);
  }
  return res;
}

TString BtagHelpers::getBtagSFFileName(BtagWPType type){
  TString res;

  if (theDataYear == 2016){
    switch (type){
    case kDeepCSV_Loose:
    case kDeepCSV_Medium:
    case kDeepCSV_Tight:
      res = "DeepCSV_2016LegacySF_WP_V1.csv";
      break;
    case kDeepFlav_Loose:
    case kDeepFlav_Medium:
    case kDeepFlav_Tight:
      res = "DeepJet_2016LegacySF_WP_V1.csv";
      break;
    default:
      break;
    }
  }
  else if (theDataYear == 2017){
    switch (type){
    case kDeepCSV_Loose:
    case kDeepCSV_Medium:
    case kDeepCSV_Tight:
      res = "DeepCSV_94XSF_WP_V4_B_F.csv";
      break;
    case kDeepFlav_Loose:
    case kDeepFlav_Medium:
    case kDeepFlav_Tight:
      res = "DeepFlavour_94XSF_WP_V3_B_F.csv";
      break;
    default:
      break;
    }
  }
  else if (theDataYear == 2018){
    switch (type){
    case kDeepCSV_Loose:
    case kDeepCSV_Medium:
    case kDeepCSV_Tight:
      res = "DeepCSV_102XSF_WP_V1.csv";
      break;
    case kDeepFlav_Loose:
    case kDeepFlav_Medium:
    case kDeepFlav_Tight:
      res = "DeepJet_102XSF_WP_V1.csv";
      break;
    default:
      break;
    }
  }
  if (res == ""){
    MELAerr << "BtagHelpers::getBtagSFFileName: WP " << type << " is not implemented for year " << theDataYear << "." << endl;
    assert(0);
  }

  res = ANALYSISTREEPKGDATAPATH + Form("ScaleFactors/bTagging/%i/", theDataYear) + res;
  HostHelpers::ExpandEnvironmentVariables(res);
  if (!HostHelpers::FileReadable(res.Data())){
    MELAerr << "BtagHelpers::getBtagSFFileName: File " << res << " is not readable." << endl;
    assert(0);
  }

  return res;
}
TString BtagHelpers::getBtagEffFileName(){
  TString res = ANALYSISTREEPKGDATAPATH + Form("ScaleFactors/bTagging/%i/Final_bTag_Efficiencies_AllMC.root", theDataYear);
  HostHelpers::ExpandEnvironmentVariables(res);
  if (!HostHelpers::FileReadable(res.Data())){
    MELAerr << "BtagHelpers::getBtagEffFileName: File " << res << " is not readable." << endl;
    assert(0);
  }
  return res;
}
TString BtagHelpers::getBtagEffHistName(BtagWPType type, const char* jet_type){
  switch (type){
  case kDeepCSV_Loose:
    return Form("DeepCSV_LooseJets_%s", jet_type);
  case kDeepCSV_Medium:
    return Form("DeepCSV_MediumJets_%s", jet_type);
  case kDeepCSV_Tight:
    return Form("DeepCSV_TightJets_%s", jet_type);
  case kDeepFlav_Loose:
    return Form("DeepFlavor_LooseJets_%s", jet_type);
  case kDeepFlav_Medium:
    return Form("DeepFlavor_MediumJets_%s", jet_type);
  case kDeepFlav_Tight:
    return Form("DeepFlavor_TightJets_%s", jet_type);
  default:
    MELAerr << "BtagHelpers::getBtagEffHistName: WP " << type << " is not implemented." << endl;
    assert(0);
    break;
  }
  return "";
}
