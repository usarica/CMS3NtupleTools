#ifndef CMS3_OFFSHELL_SAMPLEHELPERS_H
#define CMS3_OFFSHELL_SAMPLEHELPERS_H

#include "SampleHelpersCore.h"
#include "SamplesCore.h"


namespace SampleHelpers{
  extern TString theSamplesTag;

  void configure(TString period, TString stag); // Run this before doing anything else!

  TString getDatasetDirectoryName(std::string sname);
  TString getDatasetDirectoryName(TString sname);
}


#endif
