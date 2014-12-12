// -*- C++ -*-
//
// Package:    JetUtilities
// Class:      JetUtilities
// 
/**\class JetUtilities JetUtilities.h CMS2/NtupleMaker/interface/JetUtilities.h

Description: Jet utilities

*/
//
// Original Author:  Oliver Gutsche
// Tue Jul 22 22:41:55 UTC 2008
// $Id: JetUtilities.h,v 1.2 2009/09/10 10:51:37 fgolf Exp $
//
//
#ifndef CMS2_JETUTILITIES_H
#define CMS2_JETUTILITIES_H

#include <utility>
#include <vector>
#include <algorithm>
#include <iostream>

#include "DataFormats/Math/interface/LorentzVector.h"

class JetUtilities {
public:
  JetUtilities();
  ~JetUtilities();

  static bool testJetForElectrons(const math::XYZTLorentzVectorF& jetP4, const math::XYZTLorentzVectorF& elP4);

private:

};

#endif
