#ifndef BTAGHELPERS_H
#define BTAGHELPERS_H


namespace BtagHelpers{
  enum BtagWPType{
    kDeepCSV_Loose,
    kDeepCSV_Medium,
    kDeepCSV_Tight
  };
  extern BtagWPType btagWPType;

  void setBtagWPType(BtagWPType type);
  float getBtagWP(bool is80X);

}


#endif
