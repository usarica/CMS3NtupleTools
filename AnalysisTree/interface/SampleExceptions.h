#ifndef CMS3_SAMPLEEXCEPTIONS_H
#define CMS3_SAMPLEEXCEPTIONS_H

#include <vector>
#include "TString.h"

namespace SampleHelpers{
  std::vector<TString> getPUExceptions(int const& year);
  bool hasPUException(TString const& sampleIdentifier, int const& year);

}

#endif
