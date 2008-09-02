#ifndef NTUPLEMAKER_TrackerMuonMaker_H
#define NTUPLEMAKER_TrackerMuonMaker_H
// -*- C++ -*-
// $Id: TrackerMuonMaker.h,v 1.1 2008/09/02 00:55:30 dmytro Exp $


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

class TrackerMuonMaker : public edm::EDProducer {
 public:
   explicit TrackerMuonMaker(const edm::ParameterSet&);
   ~TrackerMuonMaker();

 private:
   void produce(edm::Event&, const edm::EventSetup&);
   edm::InputTag m_inputTag;
   edm::InputTag m_referenceTag;
};

#endif
