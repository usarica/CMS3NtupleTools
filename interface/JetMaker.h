// -*- C++ -*-
//
// Package:    JetMaker
// Class:      JetMaker
// 
/**\class JetMaker JetMaker.h CMS2/JetMaker/interface/JetMaker.h

Description: copy some jet variables in simble data structures into the EDM event tree

Implementation:
- take TQAF jets
- extract and fill variables
*/
//
// Original Author:  Oliver Gutsche
//         Created:  Tue Jun  9 11:07:38 CDT 2008
// $Id: JetMaker.h,v 1.1 2008/06/12 00:11:51 gutsche Exp $
//
//
#ifndef CMS2_JETMAKER_H
#define CMS2_JETMAKER_H

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

class JetMaker : public edm::EDProducer {
public:
  explicit JetMaker (const edm::ParameterSet&);
  ~JetMaker();

private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  // ----------member data ---------------------------
  edm::InputTag jetsInputTag;
  edm::InputTag genJetsInputTag;
  edm::InputTag genParticlesInputTag;

};


#endif
