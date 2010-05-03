// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      JPTtoCaloJetAssMaker
// 
/**\class JPTtoCaloJetAssMaker JPTtoCaloJetAssMaker.h CMS2/NtupleMaker/interface/JPTtoCaloJetAssMaker.h

Description: make associations between jets and muons

*/
//
// Original Author:  Oliver Gutsche
//         Created:  Tue Jun 17 20:40:42 UTC 2008
// $Id: JPTtoCaloJetAssMaker.h,v 1.1 2010/05/03 23:20:45 kalavase Exp $
//
//
#ifndef CMS2_JETTOELASSMAKER_H
#define CMS2_JETTOELASSMAKER_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class declaration
//

class JPTtoCaloJetAssMaker : public edm::EDProducer {
public:
  explicit JPTtoCaloJetAssMaker (const edm::ParameterSet&);

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  // ----------member data ---------------------------
  edm::InputTag jptInputTag_;
  edm::InputTag cms2CaloJetInputTag_;

};

#endif
