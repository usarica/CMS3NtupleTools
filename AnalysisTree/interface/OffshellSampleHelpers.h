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

}
// Here begins the more specialized functions for off-shell analysis
namespace SampleHelpers{
  void constructSamplesList(TString const& sname, SystematicsHelpers::SystematicVariationTypes syst, std::vector<TString>& samples);
  void getSamplesList(std::vector<TString> const& s, std::vector<TString>& vs, SystematicsHelpers::SystematicVariationTypes syst, std::vector<size_t>* ns=nullptr);

}

#endif
