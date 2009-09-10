// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      CommonUtils
//
/**\class CommonUtils CommonUtils.cc CMS2/CommonUtils/src/CommonUtils.cc

 Description: make associations between muons and tracks

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Thu Sep  10 11:07:38 CDT 2008
//
//


#include <Math/VectorUtil.h>
#include <iostream>
#include <limits>
#include "TMath.h"

namespace CommonUtils {
  template<typename T> inline bool isinf(T value)
    {
      value = TMath::Abs(value);
      return std::numeric_limits<T>::has_infinity && value == std::numeric_limits<T>::infinity();
    }
}
