#include <cassert>
#include "BtagHelpers.h"
#include "SamplesCore.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace SampleHelpers;
using namespace MELAStreamHelpers;


namespace BtagHelpers{
  BtagWPType btagWPType = kDeepCSV_Medium;
}


void BtagHelpers::setBtagWPType(BtagWPType type){ btagWPType=type; }

float BtagHelpers::getBtagWP(bool is80X){
  // From https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation
  float WP_DEEPCSV_TIGHT, WP_DEEPCSV_MEDIUM, WP_DEEPCSV_LOOSE;
  if (theDataYear == 2016 && is80X){
    // From https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco
    WP_DEEPCSV_TIGHT = 0.8958;
    WP_DEEPCSV_MEDIUM = 0.6324;
    WP_DEEPCSV_LOOSE = 0.2219;
  }
  else if (theDataYear == 2016){
    // From https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation2016Legacy
    WP_DEEPCSV_TIGHT = 0.8953;
    WP_DEEPCSV_MEDIUM = 0.6321;
    WP_DEEPCSV_LOOSE = 0.2217;
  }
  else if (theDataYear == 2017){
    // From https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
    WP_DEEPCSV_TIGHT = 0.8001;
    WP_DEEPCSV_MEDIUM = 0.4941;
    WP_DEEPCSV_LOOSE = 0.1522;
  }
  else if (theDataYear == 2018){
    // From https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
    WP_DEEPCSV_TIGHT = 0.7527;
    WP_DEEPCSV_MEDIUM = 0.4184;
    WP_DEEPCSV_LOOSE = 0.1241;
  }
  else{
    MELAerr << "BtagHelpers::getBtagWP: Cannot determine the b tag SF file name for year " << theDataYear << ". Aborting..." << endl;
    assert(0);
  }
  switch (btagWPType){
  case kDeepCSV_Loose:
    return WP_DEEPCSV_LOOSE;
  case kDeepCSV_Medium:
    return WP_DEEPCSV_MEDIUM;
  case kDeepCSV_Tight:
    return WP_DEEPCSV_TIGHT;
  default:
    MELAerr << "BtagHelpers::getBtagWP: No implementation for the b tagging WP type " << btagWPType << ". Aborting..." << endl;
    assert(0);
    return 0; // Just return to avoid warnings
  }
}
