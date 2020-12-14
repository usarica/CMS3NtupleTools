#ifndef DISCRIMINANTCLASSES_H
#define DISCRIMINANTCLASSES_H

#include "HelperFunctionsCore.h"
#include "SimpleDiscriminant.h"
#include "PA2PB1Discriminant.h"
#include "PA2PB2Discriminant.h"
#include "PA1PB1PBp1Discriminant.h"
#include "PA1PB2PBp2Discriminant.h"
#include "VHProdDecACDiscriminant.h"
#include "VHProdDecACGsDiscriminant.h"
#include "JJEWQCDBkgDiscriminant.h"
#include "JJEWQCDBkgM4LDiscriminant.h"
#include "SimpleInterferenceDiscriminant.h"
#include "VHProdIntACDiscriminant.h"
#include "JJEWQCDInterferenceDiscriminant.h"
#include "SimpleInterferenceTrigPhase.h"
#include "SimpleAverageInterferenceTrigPhase.h"


namespace DiscriminantClasses{
  enum Type{
    kDbkgkin,
    kDbkgdec,
    kDbkgm4l,

    kDggbkgkin,
    kDggint,
    kCggint,

    kDjVBF,
    kDjjVBF,
    kDjjZH,
    kDjjWH,

    kDjjVBFL1,
    kDjjZHL1,
    kDjjWHL1,

    kDjjVBFa2,
    kDjjZHa2,
    kDjjWHa2,

    kDjjVBFa3,
    kDjjZHa3,
    kDjjWHa3,

    kDjjVBFL1ZGs,
    kDjjZHL1ZGs,

    kDbkgjjEWQCD,
    kDbkgm4ljjEWQCD,
    kDintjjEWQCD,
    kCjjVBFint,
    kCjjVHint,

    kDL1dec,
    kDL1decint,
    kCL1decint,
    kDa2dec,
    kDa2decint,
    kCa2decint,
    kDa3dec,
    kDa3decint,
    kCa3decint,
    kDL1ZGsdec,
    kDL1ZGsdecint,
    kCL1ZGsdecint,

    kDL1jjVBFdec,
    kDL1jjVBFint,
    kCL1jjVBFint,
    kDa2jjVBFdec,
    kDa2jjVBFint,
    kCa2jjVBFint,
    kDa3jjVBFdec,
    kDa3jjVBFint,
    kCa3jjVBFint,
    kDL1ZGsjjVBFdec,
    kDL1ZGsjjVBFint,
    kCL1ZGsjjVBFint,

    kDL1jjVHdec,
    kDL1jjVHint,
    kCL1jjVHint,
    kDa2jjVHdec,
    kDa2jjVHint,
    kCa2jjVHint,
    kDa3jjVHdec,
    kDa3jjVHint,
    kCa3jjVHint,
    kDL1ZGsjjVHdec,
    kDL1ZGsjjVHint,
    kCL1ZGsjjVHint,

    kNTypes
  };

  struct KDspecs{
    Type KDtype;
    TString KDname;
    TString KDlabel;
    std::vector<TString> KDvars;

    Discriminant* KD;

    KDspecs();
    KDspecs(DiscriminantClasses::Type type);
    KDspecs(TString strname);
    bool isValid() const;
    void resetKD();
  };

  extern const std::unordered_map<TString, DiscriminantClasses::Type> mapKDNameType;
  std::unordered_map<TString, DiscriminantClasses::Type> getKDNameTypeMap();

  DiscriminantClasses::Type getKDType(const TString name);
  TString getKDName(DiscriminantClasses::Type type);
  TString getKDLabel(DiscriminantClasses::Type type);
  TString getKDLabel(const TString name);
  float getKDWP(DiscriminantClasses::Type type);
  float getKDWP(const TString name);

  Discriminant* constructKDFromType(
    const Type type,
    const TString cfilename="", const TString splinename="",
    const TString gfilename="", const TString gsplinename="", const float gscale=1
  );
  Discriminant* constructKDFromType(
    const TString name,
    const TString cfilename="", const TString splinename="",
    const TString gfilename="", const TString gsplinename="", const float gscale=1
  );

  std::vector<TString> getKDVars(const Type type);
  std::vector<TString> getKDVars(const TString name);

  bool isCPSensitive(const Type type);
  bool isCPSensitive(const TString name);

  bool usesDecInfo(const Type type);
  bool usesDecInfo(const TString name);

  bool usesProdInfo(const Type type);
  bool usesProdInfo(const TString name);

  // Channel index should be 0 unless constructing discriminants specifically for a particular 4l final state.
  // In that case, index should be the product of lepton ids.
  // strCategory can be "", "JJVBFTagged" or "HadVHTagged"
  void constructDiscriminants(std::vector<DiscriminantClasses::KDspecs>& KDlist, unsigned int idx_channel=0, TString strCategory="");

}


#endif
