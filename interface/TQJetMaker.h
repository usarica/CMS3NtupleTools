// -*- C++ -*-
//
// Package:    TQJetMaker
// Class:      TQJetMaker
// 
/**\class TQJetMaker TQJetMaker.h CMS2/NtupleMaker/interface/TQJetMaker.h

Description: copy additional TQAF jet variables in simple data structures into the EDM event tree

Implementation:
- take TQAF jets
- extract and fill variables

*/
//
// Original Author:  Oliver Gutsche
// Thu Jun 12 22:55:46 UTC 2008
// $Id: TQJetMaker.h,v 1.1 2008/06/12 23:40:09 gutsche Exp $
//
//
#ifndef CMS2_TQJETMAKER_H
#define CMS2_TQJETMAKER_H

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

class TQJetMaker : public edm::EDProducer {
public:
  explicit TQJetMaker (const edm::ParameterSet&);
  ~TQJetMaker();

private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  // ----------member data ---------------------------
  edm::InputTag tqJetsInputTag;

};


#endif
