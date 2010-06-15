// -*- C++ -*-
//
// Package:    PATElectronMaker
// Class:      PATElectronMaker
// 
/**\class PATElectronMaker PATElectronMaker.h CMS2/NtupleMaker/interface/PATElectronMaker.h

Description: copy additional PAT electron variables in simple data structures into the EDM event tree

Implementation:
- take PAT electrons
- extract and fill variables

*/
//
// Original Author:  Puneeth Kalavase
// Thu Jun 12 22:55:46 UTC 2008
// $Id: PATElectronMaker.h,v 1.5 2010/06/15 10:08:36 fgolf Exp $
//
//
#ifndef CMS2_PATElECTRONMAKER_H
#define CMS2_PATELECTRONMAKER_H

// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

//
// class decleration
//

class PATElectronMaker : public edm::EDProducer {
public:
  explicit PATElectronMaker (const edm::ParameterSet&);
  ~PATElectronMaker();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  // ----------member data ---------------------------
  edm::InputTag patElectronsInputTag_;
  edm::InputTag recoElectronsInputTag_;
	std::string aliasprefix_;
};

#endif
