// -*- C++ -*-
//
// Package:    ElUtilities
// Class:      ElUtilities
// 
/**\class ElUtilities ElUtilities.h CMS2/NtupleMaker/interface/ElUtilities.h

Description: MC utilities

*/
//
// Original Author:  Puneeth Kalavase
// Wed Oct 15 12:55:46 UTC 2008
// $Id: ElUtilities.h,v 1.1 2008/10/21 16:58:55 kalavase Exp $
//
//
#ifndef CMS2_ELUTILITIES_H
#define CMS2_ELUTILITIES_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

class ElUtilities {
public:
  ElUtilities();
  ~ElUtilities();

  static std::vector<const reco::GsfElectron* > getElectrons(const edm::Event&, const edm::InputTag);
  static void removeElectrons(const std::vector<const reco::GsfElectron*>* );
private:

};

#endif
