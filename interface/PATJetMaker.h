// -*- C++ -*-
//
// Package:    PATJetMaker
// Class:      PATJetMaker
// 
/**\class PATJetMaker PATJetMaker.h CMS2/NtupleMaker/interface/PATJetMaker.h

Description: copy additional PAT jet variables in simple data structures into the EDM event tree

Implementation:
- take PAT jets
- extract and fill variables

*/
//
// Original Author:  Oliver Gutsche
// Thu Jun 12 22:55:46 UTC 2008
// $Id: PATJetMaker.h,v 1.1 2008/10/21 07:23:46 kalavase Exp $
//
//
#ifndef CMS2_PATJETMAKER_H
#define CMS2_PATJETMAKER_H

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

class PATJetMaker : public edm::EDProducer {
public:
  explicit PATJetMaker (const edm::ParameterSet&);
  ~PATJetMaker();

private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  // ----------member data ---------------------------
  edm::InputTag patJetsInputTag;

};


#endif
