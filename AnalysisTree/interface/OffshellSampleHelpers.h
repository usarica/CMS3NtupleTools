#ifndef CMS3_OFFSHELL_SAMPLEHELPERS_H
#define CMS3_OFFSHELL_SAMPLEHELPERS_H

#include "SampleHelpersCore.h"
#include "SamplesCore.h"
#include "SystematicVariations.h"


// Still quite general stuff
namespace SampleHelpers{
  extern bool runConfigure;
  extern TString theSamplesTag;

  void configure(TString period, TString stag, HostHelpers::Hosts input_host = HostHelpers::kUCSDT2); // Run this before doing anything else!

  std::string getDatasetDirectoryCoreName(std::string sname);
  TString getDatasetDirectoryName(std::string sname, bool ignoreDNE=false);
  TString getDatasetDirectoryName(TString snam, bool ignoreDNE=false);

  TString getDatasetFileName(std::string sname, bool ignoreDNE=false);
  TString getDatasetFileName(TString sname, bool ignoreDNE=false);

  std::vector<TString> getDatasetFileNames(std::string sname, bool ignoreDNE=false);
  std::vector<TString> getDatasetFileNames(TString sname, bool ignoreDNE=false);

  TString getDatasetFirstFileName(std::string sname, bool ignoreDNE=false);
  TString getDatasetFirstFileName(TString sname, bool ignoreDNE=false);

  // Check whether the file exists locally or in the Worker directory.
  // If it exists in the Worker directory, modify the file name to point to that.
  bool checkFileOnWorker(TString& fname, HostHelpers::Hosts input_host = HostHelpers::kUCSDT2);

}
// Here begins the more specialized functions for off-shell analysis
namespace SampleHelpers{
  void constructSamplesList(TString const& sname, SystematicsHelpers::SystematicVariationTypes syst, std::vector<TString>& samples);
  void getSamplesList(std::vector<TString> const& s, std::vector<TString>& vs, SystematicsHelpers::SystematicVariationTypes syst, std::vector<size_t>* ns=nullptr);

  enum HiggsSampleDecayMode{
    // 4l samples
    kZZTo4L,
    kZZTo2L2X,

    // ZZ 2l2nu samples
    kZZTo2L2Nu,
    kZZTo2Nu2X,
    kZZTo2L2Q,
    kZZTo4Q,

    // WW 2l2nu samples
    kWWTo2L2Nu,
    kWWToLNuQQ,
    kWWToLNuXX,

    nHiggsSampleDecayModes
  };
  HiggsSampleDecayMode getHiggsSampleDecayMode(TString const& sname);
  bool isHiggsToWWDecay(SampleHelpers::HiggsSampleDecayMode const& dkmode);

  double calculateAdjustedHiggsBREff(TString const& sname, double const& sum_wgts_defaultMemberZero, double const& sum_wgts_defaultLHEEventWeight, bool hasTaus);
}

#endif
