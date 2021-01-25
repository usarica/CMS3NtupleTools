#ifndef CMS3_SAMPLEEXCEPTIONS_H
#define CMS3_SAMPLEEXCEPTIONS_H

#include <vector>
#include <utility>
#include "TString.h"

namespace SampleHelpers{
  enum GenWeightExceptionType{
    kLargeDefaultGenWeight,

    nGenWeightExceptionType
  };

  std::vector<TString> getPUExceptions(int const& year);
  bool hasPUException(TString const& sampleIdentifier, int const& year);

  std::vector< std::pair<TString, SampleHelpers::GenWeightExceptionType> > getGenWeightExceptions(int const& year);
  bool hasGenWeightException(TString const& sampleIdentifier, int const& year, SampleHelpers::GenWeightExceptionType& type);

  std::vector< std::pair<TString, float> > getXSecExceptions(int const& year);
  bool hasXSecException(TString const& sampleIdentifier, int const& year, float* xsec_mult=nullptr);

}

#endif
