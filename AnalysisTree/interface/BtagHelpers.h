#ifndef BTAGHELPERS_H
#define BTAGHELPERS_H

#include <vector>


namespace BtagHelpers{
  enum BtagWPType{
    kDeepCSV_Loose=0,
    kDeepCSV_Medium,
    kDeepCSV_Tight,
    kDeepFlav_Loose,
    kDeepFlav_Medium,
    kDeepFlav_Tight,
    nBtagWPTypes
  };
  extern BtagWPType btagWPType;

  void setBtagWPType(BtagWPType type);
  float getBtagWP(bool is80X);
  float getBtagWP(BtagHelpers::BtagWPType type, bool is80X);
  std::vector<float> getBtagWPs(bool is80X);

  TString getBtagSFFileName(BtagWPType type);
  TString getBtagEffFileName();
  TString getBtagEffHistName(BtagWPType type, const char* jet_type);



}


#endif
