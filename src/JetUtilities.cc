// -*- C++ -*-
//
// Package:    JetUtilities
// Class:      JetUtilities
// 
/**\class JetUtilities JetUtilities.cc CMS2/NtupleMaker/src/JetUtilities.cc

Description: Jet utilities

*/
//
// Original Author:  Oliver Gutsche
// Tue Jul 22 22:41:55 UTC 2008
// $Id: JetUtilities.cc,v 1.2 2009/09/10 10:51:43 fgolf Exp $
//
//
#include "CMS2/NtupleMaker/interface/JetUtilities.h"
#include "TMath.h"

typedef math::XYZTLorentzVectorF LorentzVector;

JetUtilities::JetUtilities() {
}

JetUtilities::~JetUtilities() {
}

bool JetUtilities::testJetForElectrons(const LorentzVector& jetP4, const LorentzVector& elP4) {
  
  
  bool matched = false;
  float elphi  = elP4.Phi();
  float jetphi = jetP4.Phi();
   
  float eleta  = elP4.Eta();
  float jeteta = jetP4.Eta();
   
  float dphi = elphi - jetphi;
  float deta = eleta - jeteta;
  if(fabs(dphi) > TMath::Pi() ) dphi = 2*TMath::Pi() - fabs(dphi);
   
  double dR = sqrt(dphi*dphi + deta*deta);
  if (dR < 0.4) 
    matched = true;
  
  return !matched;
}
