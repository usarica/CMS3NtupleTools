// -*- C++ -*-
//
// Package:    TQElectronMaker
// Class:      TQElectronMaker
// 
/**\class TQElectronMaker TQElectronMaker.h CMS2/NtupleMaker/interface/TQElectronMaker.h

Description: copy additional TQAF electron variables in simple data structures into the EDM event tree

Implementation:
- take TQAF electrons
- extract and fill variables

*/
//
// Original Author:  Puneeth Kalavase
// Thu Jun 12 22:55:46 UTC 2008
// $Id: TQElectronMaker.h,v 1.1 2008/06/14 16:27:26 kalavase Exp $
//
//
#ifndef CMS2_TQElECTRONMAKER_H
#define CMS2_TQELECTRONMAKER_H

// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

//
// class decleration
//

class TQElectronMaker : public edm::EDProducer {
public:
  explicit TQElectronMaker (const edm::ParameterSet&);
  ~TQElectronMaker();

private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  // ----------member data ---------------------------
  edm::InputTag tqElectronsInputTag;

};


#endif
