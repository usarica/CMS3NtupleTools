// -*- C++ -*-
//
// Package:    AllTrkJetMaker
// Class:      AllTrkJetMaker
//
/**\class AllTrkJetMaker AllTrkJetMaker.h CMS2/NtupleMaker/interface/AllTrkJetMaker.h

Description: Produces TrkJet Collection

Implementation:
*/
//
//
// Original Author:  Sanjay Padhi
//         Created:  Mon Jun 23 03:57:47 CEST 2008
// $Id: AllTrkJetMaker.h,v 1.1 2008/09/05 00:22:16 fkw Exp $
//
//

#ifndef CMS2_ALLTRKJETMAKER_H
#define CMS2_ALLTRKJETMAKER_H

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

class AllTrkJetMaker : public edm::EDProducer {
public:
  explicit AllTrkJetMaker (const edm::ParameterSet&);
  ~AllTrkJetMaker();

private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  edm::InputTag trkJetsInputTag;

};


#endif
