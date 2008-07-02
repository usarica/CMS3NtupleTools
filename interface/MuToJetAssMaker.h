// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      MuToJetAssMaker
// 
/**\class MuToJetAssMaker MuToJetAssMaker.h CMS2/NtupleMaker/interface/MuToJetAssMaker.h

Description: make associations between muons and jets

*/
//
// Original Author:  Frank Golf
//         Created:  Wed Jun 25 18:32:24 UTC 2008
// $Id: MuToJetAssMaker.h,v 1.1 2008/07/02 02:28:31 jmuelmen Exp $
//
//
#ifndef CMS2_MUTOJETASSMAKER_H
#define CMS2_MUTOJETASSMAKER_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class declaration
//

class MuToJetAssMaker : public edm::EDProducer {
public:
  explicit MuToJetAssMaker (const edm::ParameterSet&);

private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  // ----------member data ---------------------------
  double minDR;
  
};

#endif
